module simdata3D
#include "macro.h"
  use unit
  implicit none

  integer :: ld,lu,kd,ku,jd,ju,lv_min,lv_max
  

  real(4),allocatable :: x(:,:),y(:,:),z(:,:),&
       vol3D(:,:,:,:)

  real(4) :: tms(1)
  
#ifdef FUGAKU
  real(8),parameter :: time_unit_h5= 1d0 ! in second
  real(8),parameter :: vel_unit_h5 = 2.99792458d10
#else
  real(8),parameter :: time_unit_h5= 1d-3 ! in millisecond
  real(8),parameter :: vel_unit_h5 = 1d0
#endif
  
  real(4),allocatable :: &
       qrho(:,:,:,:),&
       ye  (:,:,:,:),&
       tem (:,:,:,:),&
       ut  (:,:,:,:),&
       qb  (:,:,:,:),&
       sen (:,:,:,:),&
       vlx (:,:,:,:),&
       vly (:,:,:,:),&
       vlz (:,:,:,:),&
       ! hhh (:,:,:,:),&
       rne (:,:,:,:),&
       rae (:,:,:,:),&
       deptn(:,:,:,:),&
       depta(:,:,:,:)
       
  ! integer,allocatable :: &
  !      flag_active(:,:,:,:)
  
  real(4),allocatable :: &
       vlx_b (:,:,:,:),&
       vly_b (:,:,:,:),&
       vlz_b (:,:,:,:)

!!! secondary
  real(4),allocatable :: &
       qe  (:,:,:,:),&
       pres(:,:,:,:),&
       eps (:,:,:,:),&
       hhh (:,:,:,:)


  integer,allocatable :: &
       ip_ejecta_vol(:,:,:,:)

contains
  subroutine get_ngrid_info(file_id)
    use hdf5
    use h5lt
    
    INTEGER(HID_T),intent(in) :: file_id
    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: dataset_id    ! Dataset identifier
    INTEGER(HID_T) :: dataspace_id  ! Data space identifier
    INTEGER(HSIZE_T) :: npoints
    logical :: link_exists
    character(10) :: str1
    integer :: lv

!!! open dataset
    call h5dopen_f(file_id,"/level1/x",dataset_id,error)
!!! get dataspace
    call h5dget_space_f(dataset_id, dataspace_id, error)
!!! get number of points
    call h5sget_simple_extent_npoints_f(dataspace_id, npoints, error) 
!!! close dataspace
    call h5sclose_f(dataspace_id, error)
!!! close dataset
    call h5dclose_f(dataset_id, error)

    if(npoints<10)then
       write(*,*) "Npoints does not set properly!"
       stop
    endif

#ifdef STAGGERED
    jd = -npoints/2
    ju = +npoints/2-1
    kd = -npoints/2
    ku = +npoints/2-1
#ifdef FULL
    ld = -npoints/2
#else
    ! ld = 0
    ld = -1
#endif
    lu = +npoints/2-1

#else
    jd = -(npoints-1)/2
    ju = +(npoints-1)/2
    kd = -(npoints-1)/2
    ku = +(npoints-1)/2
#ifdef FULL
    ld = -(npoints-1)/2
#else
    ld = 0
#endif
    lu = +(npoints-1)/2
#endif

    lv_min = 1
    lv = 0
    setlevel: do
       lv = lv + 1
       write(str1,'(i10)') lv
       call h5lexists_f(file_id,"/level"//trim(adjustl(str1)),link_exists,error)
       if(.not. link_exists)then
          lv_max = lv - 1
          exit setlevel
       endif

    enddo setlevel

    lv_max = 13
    write(6,*)
    write(6,*) "level min, max = ",lv_min,lv_max
    write(6,'("# of grid points = ",99i5)') npoints, jd, ju, kd, ku, ld, lu
    
  end subroutine get_ngrid_info

  subroutine get_coor(file_id)
    use hdf5
    use h5lt
    
    INTEGER(HID_T),intent(in) :: file_id
    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: dataset_id    ! Dataset identifier
    INTEGER(HID_T) :: dataspace_id  ! Data space identifier
    INTEGER(HSIZE_T) :: npoints
    integer(HSIZE_T) :: dims1(1)

    character(10) :: str1
    integer :: j,k,l,lv, jdat,kdat,ldat, ld_read

    real(8) :: dx,dy,dz,fac1,fac,vol,vol_tot,vol_tot_analytic,vol_tot_analytic2
    
    ldat = (lu-ld+1)
    kdat = (ku-kd+1)
    jdat = (ju-jd+1)

    ld_read=ld
#ifdef STAGGERED
#ifndef FULL
    ld_read=ld+1
#endif
#endif

    do lv = lv_min,lv_max
       
       write(str1,'(i10)') lv
       dims1(1) = jdat
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/x",x(:,lv),dims1,error)
     
       dims1(1) = kdat
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/y",y(:,lv),dims1,error)

       dims1(1) = ldat
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/z",z(ld_read:lu,lv),dims1,error)

    enddo

#ifdef STAGGERED
#ifndef FULL
    z(ld,lv_min:lv_max) = -z(ld+1,lv_min:lv_max)
#endif
#endif

    write(6,*)
    write(6,'("setting 3D volume element...")')
    vol3D(:,:,:,:) = 0.d0
    do lv = lv_min,lv_max

       ! dx = dlx0*2.d0**( lv_max - lv )
       ! dy = dlx0*2.d0**( lv_max - lv )
       ! dz = dlx0*2.d0**( lv_max - lv )
       dx = x(1,lv) - x(0,lv)
       dy = y(1,lv) - y(0,lv)
       dz = z(1,lv) - z(0,lv)

       !$omp parallel &
       !$omp default(none) &
       !$omp shared(ld,lu,kd,ku,jd,ju,vol3d,dx,dy,dz,lv_max,lv) &
       !$omp private(vol,fac,fac1)
       !$omp do 
       do l = ld,lu
          do k = kd,ku
             do j = jd,ju

                vol = dx*dy*dz
#ifdef STAGGERED
                if(lv/=lv_max.and.jd/2<=j.and.j<=(ju-1)/2.and.kd/2<=k.and.k<=(ku-1)/2.and.ld/2<=l.and.l<=(lu-1)/2)then
                   vol=0.d0
                endif
#ifdef FULL
                vol3D(j,k,l,lv) = vol
#else
                if(l==ld) vol=0.d0
                vol3D(j,k,l,lv) = vol*2.d0
#endif

#else

                fac1= 0.d0
                if(lv.ne.lv_max.and.jd/2.le.j.and.j.le.ju/2.and.kd/2.le.k.and.k.le.ku/2.and.ld/2.le.l.and.l.le.lu/2)then
                   fac1= 1.d0

                   if(j.eq.ju/2.or.j.eq.jd/2)then
                      fac1 = fac1*0.5d0
                   endif
                   if(k.eq.ku/2.or.k.eq.kd/2)then
                      fac1 = fac1*0.5d0
                   endif
#ifdef FULL
                   if(l.eq.lu/2.or.l.eq.ld/2)then
                      fac1 = fac1*0.5d0
                   endif
#else
                   if(l.eq.lu/2)then
                      fac1 = fac1*0.5d0
                   endif
#endif
                endif

                fac = 1.d0
                if (j.eq.ju.or.j.eq.jd)then
                   fac = fac * 0.5d0
                endif
                if (k.eq.ku.or.k.eq.kd)then
                   fac = fac * 0.5d0
                endif
                if (l.eq.lu.or.l.eq.ld)then
                   fac = fac * 0.5d0
                endif

#ifdef FULL
                vol3D(j,k,l,lv) = vol*fac*(1.d0-fac1)
#else
                vol3D(j,k,l,lv) = vol*fac*(1.d0-fac1)*2.d0
#endif
#endif

             end do
          end do
       end do
       !$omp end do
       !$omp end parallel
    end do


    vol_tot = 0.d0
    do lv = lv_min,lv_max
       do l = ld,lu
          do k = kd,ku
             do j = jd,ju

                vol_tot = vol_tot + vol3D(j,k,l,lv)

             end do
          end do
       end do
    enddo
    ! vol_tot = sum(vol3D(:,:,:,lv_max))
#ifdef STAGGERED

    lv=lv_min
    dx = x(1,lv) - x(0,lv)
    dy = y(1,lv) - y(0,lv)
    dz = z(1,lv) - z(0,lv)

#ifdef FULL
    vol_tot_analytic =  (x(ju,lv)-x(jd,lv)+1d0*dx)*(y(ku,lv)-y(kd,lv)+1d0*dy)*(z(lu,lv)+0.75d0*dz)
    vol_tot_analytic2= dx*dy*dz*dble(ju-jd+1)*dble(ku-kd+1)*dble(lu-ld+1)
#else
    vol_tot_analytic =  (x(ju,lv)-x(jd,lv)+1d0*dx)*(y(ku,lv)-y(kd,lv)+1d0*dy)*(z(lu,lv)+0.75d0*dz)*2d0 !* 7d0/8d0
    vol_tot_analytic2= dx*dy*dz*dble(ju-jd+1)*dble(ku-kd+1)*dble(lu-(ld+1)+1)*2d0
#endif
    
#else

    lv=lv_min
    dx = x(1,lv) - x(0,lv)
    dy = y(1,lv) - y(0,lv)
    dz = z(1,lv) - z(0,lv)

#ifdef FULL
    vol_tot_analytic =  (x(ju,lv)-x(jd,lv))*(y(ku,lv)-y(kd,lv))*(z(lu,lv)-z(ld,lv))
    vol_tot_analytic2= dx*dy*dz*dble(ju-jd+1)*dble(ku-kd+1)*dble(lu-ld+1)
#else
    vol_tot_analytic =  (x(ju,lv)-x(jd,lv))*(y(ku,lv)-y(kd,lv))*(z(lu,lv)-z(ld,lv))*2.d0
    vol_tot_analytic2= dx*dy*dz*dble(ju-jd+1)*dble(ku-kd+1)*dble(lu-ld+1)*2d0
#endif

#endif

    write(6,'("x,y,z = ",99es13.5)') x(ju,lv_min),x(jd,lv_min),y(ku,lv_min),y(kd,lv_min),z(lu,lv_min),z(ld,lv_min)
    write(6,'("total volume            = ",99es13.5)') vol_tot,vol_tot_analytic,vol_tot_analytic2

  end subroutine get_coor

  subroutine allocate_simdata
    
    allocate(x(jd:ju,lv_min:lv_max),y(kd:ku,lv_min:lv_max),z(ld:lu,lv_min:lv_max),&
         qrho(jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         ye  (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         tem (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         ut  (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         qb  (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         sen (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vlx (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vly (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vlz (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         rne (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         rae (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         deptn(jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         depta(jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vol3D(jd:ju,kd:ku,ld:lu,lv_min:lv_max), &
         ! flag_active(jd:ju,kd:ku,ld:lu,lv_min:lv_max), &
         vlx_b (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vly_b (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vlz_b (jd:ju,kd:ku,ld:lu,lv_min:lv_max)   )

    allocate( &
         qe  (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         pres(jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         eps (jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         hhh (jd:ju,kd:ku,ld:lu,lv_min:lv_max) )

    allocate( &
         ip_ejecta_vol(jd:ju,kd:ku,ld:lu,lv_min:lv_max) )

    
  end subroutine allocate_simdata

  subroutine read_simdata(file_id,it,t)
    use hdf5
    use h5lt
    use unit

    INTEGER(HID_T),intent(in) :: file_id
    integer,intent(in)  :: it
    real(8),intent(out) :: t

    integer :: error, sum_err
    integer(HSIZE_T) :: dims1(1),dims3(3)
    real(4),allocatable :: buf3d_real4_1(:,:,:), buf3d_real4_2(:,:,:), buf3d_real4_3(:,:,:), buf3d_real4_4(:,:,:), buf3d_real4_5(:,:,:), buf3d_real4_6(:,:,:), buf3d_real4_7(:,:,:), buf3d_real4_8(:,:,:)

    integer :: lv,jdat,kdat,ldat,ld_read
    character(10) :: str1,str2

    integer :: lv_max_present
    logical :: link_exists


    ld_read=ld

#ifdef STAGGERED
#ifndef FULL
    ld_read=ld+1
#endif
#endif

    write(str2,'(i10)') it

    call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str2))//"/time",tms,dims1,error)
    t = tms(1)*time_unit_h5
    
    ldat = (lu-ld+1)
    kdat = (ku-kd+1)
    jdat = (ju-jd+1)

    dims3(1) = jdat
    dims3(2) = kdat
    dims3(3) = ldat

    allocate(buf3d_real4_1(jdat,kdat,ldat), buf3d_real4_2(jdat,kdat,ldat), buf3d_real4_3(jdat,kdat,ldat), buf3d_real4_4(jdat,kdat,ldat), buf3d_real4_5(jdat,kdat,ldat), buf3d_real4_6(jdat,kdat,ldat), buf3d_real4_7(jdat,kdat,ldat), buf3d_real4_8(jdat,kdat,ldat))
    
    ! !$omp parallel default(none) &
    ! !$omp num_threads(1) &
    ! !$omp shared(lv_max,lv_min,ld_read,lu,file_id,str2,dims3,qrho,ut,ye,sen,tem,qb,vlx,vly,vlz) &
    ! !$omp private(rbuf3,error,sum_err,str1)
    ! !$omp do
    loop_read:do lv=lv_min,lv_max

       write(str1,'(i10)') lv

       call h5lexists_f(file_id,"/level"//trim(adjustl(str1)),link_exists,error)
       if(link_exists)then
          lv_max_present=lv
       else
          write(6,*) "lv=",lv,"does not exists. skip."
          exit loop_read
       endif
        
       !write(6,*) "/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/density"
       sum_err = 0
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/density",buf3d_real4_1,dims3,error); sum_err = sum_err + error
       qrho(:,:,ld_read:lu,lv) = buf3d_real4_1(:,:,:)
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/u_t"    ,buf3d_real4_2,dims3,error); sum_err = sum_err + error
       ut  (:,:,ld_read:lu,lv) = buf3d_real4_2(:,:,:)
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/ye"     ,buf3d_real4_3,dims3,error); sum_err = sum_err + error
       ye  (:,:,ld_read:lu,lv) = buf3d_real4_3(:,:,:)
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/entropy",buf3d_real4_4,dims3,error); sum_err = sum_err + error
       sen (:,:,ld_read:lu,lv) = buf3d_real4_4(:,:,:)
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/temperature",buf3d_real4_5,dims3,error); sum_err = sum_err + error
       tem (:,:,ld_read:lu,lv) = buf3d_real4_5(:,:,:)
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/vx"     ,buf3d_real4_6,dims3,error); sum_err = sum_err + error
       vlx (:,:,ld_read:lu,lv) = buf3d_real4_6(:,:,:)/vel_unit_h5
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/vy"     ,buf3d_real4_7,dims3,error); sum_err = sum_err + error
       vly (:,:,ld_read:lu,lv) = buf3d_real4_7(:,:,:)/vel_unit_h5
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/vz"     ,buf3d_real4_8,dims3,error); sum_err = sum_err + error
       vlz (:,:,ld_read:lu,lv) = buf3d_real4_8(:,:,:)/vel_unit_h5

       call h5lexists_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/rho_star",link_exists,error)
       if(link_exists)then
          call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))//"/rho_star",buf3d_real4_1,dims3,error); sum_err = sum_err + error
          qb(:,:,ld_read:lu,lv) = buf3d_real4_1(:,:,:)*rho_uni
          ! qb(:,:,:,lv) = qb(:,:,:,lv) *rho_uni
       else
          !$omp parallel
          !$omp do
          do l=ld,lu
             do k=kd,ku
                do j=jd,ju
                   qb(j,k,l,lv) = qrho(j,k,l,lv)/sqrt(1d0 - ( vlx(j,k,l,lv)**2 + vly(j,k,l,lv)**2 + vlz(j,k,l,lv)**2 ) )
                enddo
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       
       if(sum_err>0)then
          write(6,*) "Error in reading hdf5 file. STOP.",lv,sum_err
          stop
       endif
       
       ! write(6,*) lv

    enddo loop_read
    ! !$omp end do
    ! !$omp end parallel

    deallocate(buf3d_real4_1, buf3d_real4_2, buf3d_real4_3, buf3d_real4_4, buf3d_real4_5, buf3d_real4_6, buf3d_real4_7, buf3d_real4_8)
    
#ifdef STAGGERED
#ifndef FULL
    call zboundary
#endif
#endif

    write(6,*) lv_max, lv_max_present
    call interp_finer(lv_max_present)
    ! stop
    ! block
    !   integer :: j,k,l
    !   lv=lv_max
    !   l=ld+1
    !   do k=kd,ku,3
    !      do j=jd,ju,3
    !         write(99,'(99es12.4)') x(j,lv), y(k,lv), qb(j,k,l,lv), qrho(j,k,l,lv), ut(j,k,l,lv)
    !      enddo
    !   enddo
    ! end block

  end subroutine read_simdata

  subroutine interp_finer(lv_max_present)
    ! data exists for lv <= lv_max_present
    integer,intent(in) :: lv_max_present
    
    integer :: j,k,l,lv,lv0, j1,k1,l1, ld_read
    real(8) :: xx, yy, zz

    ld_read=ld
#ifdef STAGGERED
#ifndef FULL
    ld_read=ld+1
#endif
#endif
    
    if(lv_max>lv_max_present)then
       do lv=lv_max_present+1,lv_max
          do l=ld_read,lu
             do k=kd,ku
                do j=jd,ju
                   xx=x(j,lv)
                   yy=y(k,lv)
                   zz=z(l,lv)
                   call coorindex3D(xx,yy,zz,j1,k1,l1,lv0,lv_max_present)
                   write(6,*) j1,k1,l1,lv0
                   write(6,*) x(j1-1,lv0), x(j,lv), x(j1,lv0)
                   write(6,*) y(k1-1,lv0), y(k,lv), y(k1,lv0)
                   write(6,*) z(l1-1,lv0), z(l,lv), z(l1,lv0)
                   stop
                   
                enddo
             enddo
          enddo
          
       enddo
    endif


  end subroutine interp_finer

  subroutine zboundary
    integer :: j,k,l,lv

    do lv=lv_min,lv_max
       !$omp parallel default(none) &
       !$omp shared(kd,ku,jd,ju,ld,lv,qrho,qb,ut,ye,sen,tem,vlx,vly,vlz)
       !$omp do
       do k=kd,ku
          do j=jd,ju
             qrho(j,k,ld,lv) = qrho(j,k,ld+1,lv)
             qb  (j,k,ld,lv) = qb  (j,k,ld+1,lv)
             ut  (j,k,ld,lv) = ut  (j,k,ld+1,lv)
             ye  (j,k,ld,lv) = ye  (j,k,ld+1,lv)
             tem (j,k,ld,lv) = tem (j,k,ld+1,lv)
             sen (j,k,ld,lv) = sen (j,k,ld+1,lv)
             vlx (j,k,ld,lv) = vlx (j,k,ld+1,lv)
             vly (j,k,ld,lv) = vly (j,k,ld+1,lv)
             vlz (j,k,ld,lv) =-vlz (j,k,ld+1,lv)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    
  end subroutine zboundary
  
  subroutine set_secondary
    use module_eos
    use unit

    integer :: j,k,l,lv, irho,irho1,iye,iye1,itemp,itemp1
    real(8) :: rhot,fyet,temt, ss,ssp,uu,uup,tt,ttp

    do lv=lv_min,lv_max
       !$omp parallel default(none) &
       !$omp shared(lv,ld,lu,kd,ku,jd,ju,qrho,ye,tem,rho_e_min,rho_e,nrho,drhoi,ye_e_min,ye_e,nye,dyei,tem_e_min,tem_e,ntemp,dtemi, &
       !$omp     pres_e, eps_e,pres,eps,hhh) &
       !$omp private(j,k,l,rhot,fyet,temt,irho,irho1,uu,uup,iye,iye1,tt,ttp,itemp,itemp1,ss,ssp)
       !$omp do
       do l=ld,lu
          do k=kd,ku
             do j=jd,ju

                rhot = qrho(j,k,l,lv)
                fyet = ye  (j,k,l,lv)
                temt = tem (j,k,l,lv)

                irho = max(1   , min(nrho-1, int((log10(rhot)-rho_e_min)*drhoi)+1))
                irho1= irho+1
                uu   = max(0.d0, min(1.d0, (log10(rhot)-rho_e(irho))*drhoi))
                uup  = 1.d0-uu

                iye  = max(1 , min(nye-1, int((fyet-ye_e_min  )*dyei) ))
                iye1 = iye +1
                tt   = max(0.d0, min(1.d0 ,     (fyet-ye_e(iye))*dyei))
                ttp  = 1.d0-tt

                itemp = max(1 , min(ntemp-1, int((log10(temt)-tem_e_min  )*dtemi)+1))
                itemp1=itemp+1
                ss    = max(0.d0, min(1.d0 ,     (log10(temt)-tem_e(itemp))*dtemi))
                ssp   = 1.d0-ss

                if(irho<1.or.nrho-1<irho.or. &
                     iye<1.or.nye-1<iye.or. &
                     itemp<1.or.ntemp-1<itemp)then
                   write(6,'(4i5)') j,k,l,lv
                   write(6,'(3i5)') irho,itemp,iye
                   write(6,'(99es12.4)') rhot, fyet, temt
                endif
                
                pres(j,k,l,lv) = ssp *ttp *uup *pres_e(itemp ,iye ,irho )   &
                               + ss  *ttp *uup *pres_e(itemp1,iye ,irho )   &
                               + ssp *tt  *uup *pres_e(itemp ,iye1,irho )   &
                               + ssp *ttp *uu  *pres_e(itemp ,iye ,irho1)   &
                               + ss  *tt  *uup *pres_e(itemp1,iye1,irho )   &
                               + ss  *ttp *uu  *pres_e(itemp1,iye ,irho1)   &
                               + ssp *tt  *uu  *pres_e(itemp ,iye1,irho1)   &
                               + ss  *tt  *uu  *pres_e(itemp1,iye1,irho1)
                eps (j,k,l,lv) = ssp *ttp *uup * eps_e(itemp ,iye ,irho )   &
                               + ss  *ttp *uup * eps_e(itemp1,iye ,irho )   &
                               + ssp *tt  *uup * eps_e(itemp ,iye1,irho )   &
                               + ssp *ttp *uu  * eps_e(itemp ,iye ,irho1)   &
                               + ss  *tt  *uup * eps_e(itemp1,iye1,irho )   &
                               + ss  *ttp *uu  * eps_e(itemp1,iye ,irho1)   &
                               + ssp *tt  *uu  * eps_e(itemp ,iye1,irho1)   &
                               + ss  *tt  *uu  * eps_e(itemp1,iye1,irho1)

                eps (j,k,l,lv) = 1.d1**eps (j,k,l,lv) - 1.d0
                pres(j,k,l,lv) = 1.d1**pres(j,k,l,lv)
                hhh (j,k,l,lv) = 1.d0 + eps(j,k,l,lv) + pres(j,k,l,lv)/rhot/v_uni**2
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo


    return
  end subroutine set_secondary

  subroutine all_proc(fn)
    use hdf5
    use h5lt
    character(*),intent(in) :: fn
    INTEGER        :: hdf_err     ! Error flag
    INTEGER(HID_T) :: file_id     ! File identifier
    integer(HSIZE_T) :: dims1(1),dims2(2),dims3(3)
    real(4),allocatable :: buf1(:), buf2(:,:), buf3(:,:,:)
    
    real(8) :: t
    integer :: it

    it = 18
    write(6,*) fn
    call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
!!! obtain information on the grid
    call get_ngrid_info(file_id)
!!! allocate simulation data variables
    call allocate_simdata
!!! set coordinate (x,z)
    call get_coor(file_id)
stop    
    call read_simdata(file_id,it,t)
    call h5fclose_f(file_id, hdf_err)
  end subroutine all_proc

  subroutine coorindex3D(xx,yy,zz,j1,k1,l1,lv,lv_max_search)
!    use simdata3D
!    implicit none
    real(8),intent(in) :: xx,yy,zz
    integer,intent(out) :: lv,j1,k1,l1
    integer,intent(in),optional :: lv_max_search

    lv = lv_max
    if(present(lv_max_search))then
       if(lv_max_search > lv_max)then
          write(6,*) "lv_max_search > lv_max", lv_max_search, lv_max
          stop
       else
          lv = lv_max_search
       endif
    endif

    do
       ! write(6,'(99es11.3)') abs(xx),x(ju,lv), abs(yy),y(ku,lv), abs(zz),z(lu,lv)
       if((abs(xx)<x(ju,lv).and.abs(yy)<y(ku,lv).and.abs(zz)<z(lu,lv)).or.lv==lv_min)then
          exit
       else
          lv = lv - 1
       endif
    enddo

    j1 = jd + 1
    do while( ( xx - x(j1-1,lv) )*( xx - x(j1,lv) ) > 0.d0 .and.j1<ju)
       j1 = j1 + 1
    end do

    k1 = kd + 1
    do while( ( yy - y(k1-1,lv) )*( yy - y(k1,lv) ) > 0.d0 .and.k1<ku)
       k1 = k1 + 1
    end do

    l1 = ld + 1
    do while( ( zz - z(l1-1,lv) )*( zz - z(l1,lv) ) > 0.d0 .and.l1<lu)
       l1 = l1 + 1
    end do

    if(lv<lv_min .or. lv>lv_max .or. &
         j1 < jd + 1 .or. ju < j1 .or. &
         k1 < kd + 1 .or. ku < k1 .or. &
         l1 < ld + 1 .or. lu < l1 )then
       write(6,'(a,99i10)') "index out of range", j1,k1,l1,lv
    endif

  end subroutine coorindex3D

end module simdata3D


