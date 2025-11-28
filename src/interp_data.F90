module interp_data3D
#include "macro.h"
  use unit
  implicit none

  private
  
  public :: interp_get_grid_info, interp_get_coor, interp_allocate_simdata, interp_data, interp_output_hdf
  
  integer :: ld,lu,kd,ku,jd,ju,lv_min,lv_max
  
  real(4),allocatable :: x_ip(:,:),y_ip(:,:),z_ip(:,:)

  real(4) :: tms(1)
  
#ifdef FUGAKU
  real(8),parameter :: time_unit_h5= 1d0 ! in second
  real(8),parameter :: vel_unit_h5 = 2.99792458d10
#else
  real(8),parameter :: time_unit_h5= 1d-3 ! in millisecond
  real(8),parameter :: vel_unit_h5 = 1d0
#endif
  
  real(4),allocatable :: &
       vlx_ip (:,:,:,:),&
       vly_ip (:,:,:,:),&
       vlz_ip (:,:,:,:) 
contains

  subroutine interp_get_grid_info(file_id)
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

    write(6,*)
    write(6,*) "level min, max = ",lv_min,lv_max
    write(6,'("# of grid points = ",99i5)') npoints, jd, ju, kd, ku, ld, lu
    
  end subroutine interp_get_grid_info

  subroutine interp_get_coor(file_id)
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
    
    ld_read=ld
#ifdef STAGGERED
#ifndef FULL
    ld_read=ld+1
#endif
#endif

    ldat = (lu-ld_read+1)
    kdat = (ku-kd+1)
    jdat = (ju-jd+1)


    do lv = lv_min,lv_max
       
       write(str1,'(i10)') lv
       dims1(1) = jdat
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/x",x_ip(:,lv),dims1,error)
     
       dims1(1) = kdat
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/y",y_ip(:,lv),dims1,error)

       dims1(1) = ldat
       call H5LTread_dataset_float_f(file_id,"/level"//trim(adjustl(str1))//"/z",z_ip(ld_read:lu,lv),dims1,error)
       
    enddo

#ifdef STAGGERED
#ifndef FULL
    z_ip(ld,lv_min:lv_max) = -z_ip(ld+1,lv_min:lv_max)
#endif
#endif

    write(6,'("x,y,z = ",99es13.5)') x_ip(ju,lv_min),x_ip(jd,lv_min),y_ip(ku,lv_min),y_ip(kd,lv_min),z_ip(lu,lv_min),z_ip(ld,lv_min)


  end subroutine interp_get_coor

  subroutine interp_allocate_simdata
    
    allocate(x_ip(jd:ju,lv_min:lv_max),y_ip(kd:ku,lv_min:lv_max),z_ip(ld:lu,lv_min:lv_max),&
         vlx_ip(jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vly_ip(jd:ju,kd:ku,ld:lu,lv_min:lv_max),&
         vlz_ip(jd:ju,kd:ku,ld:lu,lv_min:lv_max))

  end subroutine interp_allocate_simdata


  subroutine interp_data()
    use simdata3D,only : x,y,z,vlx,vly,vlz,coorindex3D

    integer :: j,k,l,lv,lv0, j1,k1,l1, ld_read, j0,k0,l0
    real(8) :: xx, yy, zz
    real(8) :: x0,x1,y0,y1,z0,z1
    real(8) :: vlx_i, vly_i, vlz_i

    ld_read=ld
#ifdef STAGGERED
#ifndef FULL
    ld_read=ld+1
#endif
#endif
    
    do lv=lv_min,lv_max
       !$omp parallel default(none)&
!        !$omp num_threads(1) &
       !$omp shared(lv,ld_read,lu,kd,ku,jd,ju,x,y,z,x_ip,y_ip,z_ip,vlx,vly,vlz,vlx_ip,vly_ip,vlz_ip) &
       !$omp private(xx,yy,zz,lv0,j0,k0,l0,j1,k1,l1,vlx_i,vly_i,vlz_i,x0,y0,z0,x1,y1,z1)
       !$omp do
       do l=ld_read,lu
          do k=kd,ku
             do j=jd,ju

                xx=x_ip(j,lv)
                yy=y_ip(k,lv)
                zz=z_ip(l,lv)

                call coorindex3D(xx,yy,zz,j1,k1,l1,lv0)

                j0=j1-1
                k0=k1-1
                l0=l1-1

                x1 = (xx-x(j0,lv0))/(x(j1,lv0)-x(j0,lv0))
                x0 = 1.d0-x1
                y1 = (yy-y(k0,lv0))/(y(k1,lv0)-y(k0,lv0))
                y0 = 1.d0-y1
                z1 = (zz-z(l0,lv0))/(z(l1,lv0)-z(l0,lv0))
                z0 = 1.d0-z1

                ! if(mod(l-ld_read,10)==0.and.mod(k-kd,10)==0.and.mod(j-jd,10)==0)then
                !    write(6,*) x1,y1,z1
                ! endif
                vlx_i = x1*y1*z1* vlx(j1,k1,l1,lv0) &
                     + x0*y1*z1* vlx(j0,k1,l1,lv0) &
                     + x1*y0*z1* vlx(j1,k0,l1,lv0) &
                     + x0*y0*z1* vlx(j0,k0,l1,lv0) &
                     + x1*y1*z0* vlx(j1,k1,l0,lv0) &
                     + x0*y1*z0* vlx(j0,k1,l0,lv0) &
                     + x1*y0*z0* vlx(j1,k0,l0,lv0) &
                     + x0*y0*z0* vlx(j0,k0,l0,lv0)
                vly_i = x1*y1*z1* vly(j1,k1,l1,lv0) &
                     + x0*y1*z1* vly(j0,k1,l1,lv0) &
                     + x1*y0*z1* vly(j1,k0,l1,lv0) &
                     + x0*y0*z1* vly(j0,k0,l1,lv0) &
                     + x1*y1*z0* vly(j1,k1,l0,lv0) &
                     + x0*y1*z0* vly(j0,k1,l0,lv0) &
                     + x1*y0*z0* vly(j1,k0,l0,lv0) &
                     + x0*y0*z0* vly(j0,k0,l0,lv0)
                vlz_i = x1*y1*z1* vlz(j1,k1,l1,lv0) &
                     + x0*y1*z1* vlz(j0,k1,l1,lv0) &
                     + x1*y0*z1* vlz(j1,k0,l1,lv0) &
                     + x0*y0*z1* vlz(j0,k0,l1,lv0) &
                     + x1*y1*z0* vlz(j1,k1,l0,lv0) &
                     + x0*y1*z0* vlz(j0,k1,l0,lv0) &
                     + x1*y0*z0* vlz(j1,k0,l0,lv0) &
                     + x0*y0*z0* vlz(j0,k0,l0,lv0)

                vlx_ip(j,k,l,lv) = vlx_i
                vly_ip(j,k,l,lv) = vly_i
                vlz_ip(j,k,l,lv) = vlz_i

             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       write(6,*) "lv=",lv,"done"
    enddo



  end subroutine interp_data


  subroutine interp_output_hdf(fn_interp, t, it)
#include "macro.h"
    use hdf5
    use h5lt
    implicit none

    character(*),intent(in) :: fn_interp
    real(8),intent(in)  :: t
    integer,intent(in)  :: it

    character(256) :: fn, str1,str2,str3

    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: group_id, group2_id      ! Group identifier
    integer :: hdf_err
    integer(HSIZE_T) :: dims1(1),dims2(2),dims3(3)
    integer :: jdat,kdat,ldat, ld_write

    integer :: j,k,l,lv

    real(4),allocatable :: buf3d_real4_1(:,:,:)

    fn = trim(fn_interp)
    call h5fcreate_f(fn, H5F_ACC_TRUNC_F, file_id, hdf_err)

    ld_write=ld
#ifdef STAGGERED
#ifndef FULL
    ld_write=ld+1
#endif
#endif

    jdat = ju-jd+1
    kdat = ku-kd+1
    ldat = lu-ld_write+1

    allocate(buf3d_real4_1(jdat,kdat,ldat))

    do lv=lv_min,lv_max

       write(str1,'(i10)') lv
       call h5gcreate_f(file_id, "level"//trim(adjustl(str1)) , group_id, hdf_err)

       dims1(1) = jdat
       call h5ltmake_dataset_float_f(group_id, "x", 1, dims1, x_ip(:,lv), hdf_err)
       dims1(1) = kdat
       call h5ltmake_dataset_float_f(group_id, "y", 1, dims1, y_ip(:,lv), hdf_err)
       dims1(1) = ldat
       call h5ltmake_dataset_float_f(group_id, "z", 1, dims1, z_ip(ld_write:lu,lv), hdf_err)

       ! t-dependent data

       write(str2,'(i10)') it
       call h5gcreate_f(group_id, "data"//trim(adjustl(str2)) , group2_id, hdf_err)

       dims1(1) = 1
       tms(1) = t/time_unit_h5
       call h5ltmake_dataset_float_f(group2_id,"time", 1, dims1, tms, hdf_err)

       dims3(1) = jdat
       dims3(2) = kdat
       dims3(3) = ldat

       buf3d_real4_1(:,:,:) = vlx_ip(:,:,ld_write:lu,lv)*vel_unit_h5
       call h5ltmake_dataset_float_f(group2_id, "vx", 3, dims3, buf3d_real4_1, hdf_err)
       buf3d_real4_1(:,:,:) = vly_ip(:,:,ld_write:lu,lv)*vel_unit_h5
       call h5ltmake_dataset_float_f(group2_id, "vy", 3, dims3, buf3d_real4_1, hdf_err)
       buf3d_real4_1(:,:,:) = vlz_ip(:,:,ld_write:lu,lv)*vel_unit_h5
       call h5ltmake_dataset_float_f(group2_id, "vz", 3, dims3, buf3d_real4_1, hdf_err)

       call h5gclose_f(group_id, hdf_err)
       call h5gclose_f(group2_id, hdf_err)

    enddo

    call h5fclose_f(file_id, hdf_err)

    deallocate(buf3d_real4_1)
  end subroutine interp_output_hdf

end module interp_data3D


