module module_flux_sample

  implicit none

  private
  public :: init_angle_sample, sample_compare

  integer :: n_points
  real(8),allocatable :: th(:),ph(:), dom(:)
  real(8) :: dth


  real(8),allocatable :: &
       ut_sample(:), &
       hhh_sample(:), &
       vlr_sample(:), &
       qrho_sample(:), &
       ye_sample(:), &
       sen_sample(:), &
       dm_sample(:)
       

contains


  subroutine init_angle_sample(nside)
#include "macro.h"
    use const
    use io
    integer,intent(in) :: nside
    
    integer :: ip
    integer :: p,i,j,pp
    real(8) :: z,phi,s,phh
    integer :: npix, nnorth, neq, nring

    npix = 12*nside**2
    neq = 4*nside
    nnorth = int((npix - neq)/2) + neq
    nring = 4*nside - 1
    write(6,*) "Npix, Nnorth = ", npix, nnorth

    dth = 2d0*pi/dble(neq)
    write(6,*) "dth/pi = ", dth/pi
    
#ifdef MIRROR
    n_points = nnorth
#else
    n_points = npix
#endif

    allocate(th(n_points),ph(n_points),dom(n_points))
    
    do ip=1,nnorth
       
       p = ip-1
       
       phh = dble(p+1)/2d0
       i = int(sqrt(phh - sqrt(dble(int(phh))) )) + 1

       if(i<nside)then

          j = p + 1 - 2*i*(i-1)

          z = 1d0 -dble(i)**2/(3d0*dble(nside)**2)
          s = 1d0
          phi = pi/(2d0*dble(i)) * (dble(j) - s/2d0)
          ! write(6,'(3i5,2es12.4," ",a)') p,i,j,z,phi, "north-pole"
       else

          pp = p - 2*nside*(nside-1)
          i = int( dble(pp)/dble(4*nside) ) + nside

          if(nside <= i .and. i <= 2*nside)then

             j = mod(pp, (4*nside)) + 1
             
             z=4d0/3d0 - 2d0*dble(i)/(3d0*dble(nside))
             s = dble(mod((i-nside+1), 2))
             phi = pi/(2d0*dble(nside)) * (dble(j) - s/2d0)

             !write(6,'(3i5,2es12.4," ",a)') p,i,j,z,phi,"equatorial"
          else
        
             write(6,*) "Something is wrong, stop."
             stop
             
          endif
       endif
       
       th(ip) = acos(z)
       ph(ip) = phi
    enddo

    do ip = 1,n_points
       dom(ip) = 4d0*pi/dble(npix)
    enddo

#ifdef MIRROR
    do ip = 1,nnorth-neq
       dom(ip) = dom(ip) * 2d0
    enddo
#endif


    allocate(ut_sample(n_points), &
       hhh_sample(n_points), &
       vlr_sample(n_points), &
       qrho_sample(n_points), &
       ye_sample(n_points), &
       sen_sample(n_points), &
       dm_sample(n_points) )
    
    block
      integer :: unit
      open(newunit=unit,file=trim(dir_out)//"/angle_info_sampling.dat",status="replace")
      write(unit,*) "#",sum(dom(:))/4.d0/pi,n_points
      do ip=1,n_points
         write(unit,'(i5,99es12.4)') ip,th(ip),ph(ip),sin(th(ip))*cos(ph(ip)),sin(th(ip))*sin(ph(ip)),cos(th(ip)),dom(ip)
      enddo
      close(unit)
    end block
    
  end subroutine init_angle_sample


  subroutine sample_compare(rfl,ipu,np_set,job,it,time)
    use particle_data, only : var_p, index_dm, index_ye, index_sen, index_hhh, index_ut, index_rho, index_x, index_y, index_z
    use module_eos, only : hhh_min
    use io, only : dir_out
    !$use omp_lib

    real(8),intent(in) :: rfl,time
    integer,intent(in) :: ipu, np_set, job,it
    
    integer :: iye,ientr
    integer,parameter :: nye=61,nentr=101
    
    real(8) :: ye_val (nye),logentr_val (nentr)
    real(8) :: ye_hist_particle(nye),entr_hist_particle(nentr)
    real(8) :: ye_hist_sample(nye),entr_hist_sample(nentr)
    real(8),parameter :: ye_min = 0.005d0, ye_max = 0.605d0, dye = (ye_max-ye_min)/dble(nye-1)
    real(8),parameter :: logentr_min = -1d0, logentr_max = 3d0, dlogentr = (logentr_max-logentr_min)/dble(nentr-1)

    integer :: ip

    integer :: nunit
    character(200) :: str1,str2

    do iye = 1,nye
       ye_val(iye) = ye_min + (ye_max-ye_min)*dble(iye-1)/dble(nye-1)
       !write(6,*) iye,ye_val(iye)
    enddo
    do ientr = 1,nentr
       logentr_val(ientr) = logentr_min + (logentr_max-logentr_min)*dble(ientr-1)/dble(nentr-1)
       !write(6,*) ientr, logentr_val(ientr)
    enddo

    call sample_flux(rfl)

    ye_hist_particle(:) = 0d0
    entr_hist_particle(:) = 0d0
    do ip = ipu-np_set+1,ipu
       iye = min(nye,max(1,int( (var_p(index_ye,ip)-ye_min)/dye ) + 1))
       ientr = min(nentr,max(1,int( (log10(var_p(index_sen,ip))-logentr_min)/dlogentr ) + 1))
       
       ye_hist_particle(iye) = ye_hist_particle(iye) + var_p(index_dm,ip)
       entr_hist_particle(ientr) = entr_hist_particle(ientr) + var_p(index_dm,ip)
    enddo

    ye_hist_sample(:) = 0d0
    entr_hist_sample(:) = 0d0
    !$omp parallel default(none) &
    !$omp shared(n_points,hhh_sample,ut_sample,hhh_min,ye_sample, &
    !$omp   sen_sample,dm_sample) &
    !$omp private(iye,ientr) &
    !$omp reduction(+:ye_hist_sample, entr_hist_sample)
    !$omp do
    do ip = 1,n_points
       
       if( hhh_sample(ip)*ut_sample(ip) + hhh_min < 0d0 )then
          iye = min(nye,max(1,int( (ye_sample(ip)-ye_min)/dye ) + 1))
          ientr = min(nentr,max(1,int( (log10(sen_sample(ip))-logentr_min)/dlogentr ) + 1))
          
          ye_hist_sample(iye) = ye_hist_sample(iye) + dm_sample(ip)
          entr_hist_sample(ientr) = entr_hist_sample(ientr) + dm_sample(ip)
       endif
    enddo
    !$omp end do
    !$omp end parallel

    write(str1,'(i3.3)') job
    write(str2,'(i6.6)') it

    open(newunit=nunit, file=trim(dir_out)//"/ye_hist_flux_"//trim(str1)//"_"//trim(str2)//".dat",status="replace",action="write")
    write(nunit,'("#",99es15.7)') time, sum(ye_hist_particle(:)), sum(ye_hist_sample(:))
    do iye=1,nye
       write(nunit,'(99es15.7)') ye_val(iye),ye_hist_particle(iye)/sum(ye_hist_particle(:)), ye_hist_sample(iye)/sum(ye_hist_sample(:))
    enddo
    close(nunit)

    open(newunit=nunit, file=trim(dir_out)//"/sen_hist_flux_"//trim(str1)//"_"//trim(str2)//".dat",status="replace",action="write")
    write(nunit,'("#",99es15.7)') time, sum(entr_hist_particle(:)), sum(entr_hist_sample(:))
    do ientr=1,nentr
       write(nunit,'(99es15.7)') 10d0**logentr_val(ientr),entr_hist_particle(ientr)/sum(entr_hist_particle(:)), entr_hist_sample(ientr)/sum(entr_hist_sample(:))
    enddo
    close(nunit)
    
    open(newunit=nunit, file=trim(dir_out)//"/flux_sample_"//trim(str1)//"_"//trim(str2)//".dat",status="replace",action="write")
    write(nunit,'("#",99es15.7)') time
    do ip=1,n_points
       write(nunit,'(i8,99es15.7)') ip,cos(ph(ip)),cos(th(ip)),&
            ut_sample(ip), &
            hhh_sample(ip), &
            qrho_sample(ip), &
            ye_sample(ip), &
            sen_sample(ip), &
            dm_sample(ip)
    enddo
    close(nunit)

    open(newunit=nunit, file=trim(dir_out)//"/flux_particle_"//trim(str1)//"_"//trim(str2)//".dat",status="replace",action="write")
    write(nunit,'("#",99es15.7)') time
    do ip=ipu-np_set+1,ipu
       write(nunit,'(i8,99es15.7)') ip,var_p(index_x,ip)/sqrt(var_p(index_x,ip)**2+var_p(index_y,ip)**2), var_p(index_z,ip)/sqrt(var_p(index_x,ip)**2 + var_p(index_y,ip)**2 + var_p(index_z,ip)**2 ), &
            var_p(index_ut,ip), &
            var_p(index_hhh,ip), &
            var_p(index_rho,ip), &
            var_p(index_ye,ip), &
            var_p(index_sen,ip), &
            var_p(index_dm,ip)
    enddo
    close(nunit)
    

  end subroutine sample_compare



  subroutine sample_flux(rfl)
    
    use const, only : pi
    use simdata3D
    use unit
    use module_eos
    !$use omp_lib

    real(8),intent(in) :: rfl
    
    integer :: j,k,l,lv,j1,k1,l1,ip
    real(8) :: xi,yi,zi, x0,x1,y0,y1,z0,z1
    real(8) :: ut_p,hhh_p,vlx_p,vly_p,vlz_p,vlr_p,qb_p, qrho_p,ye_p,sen_p, dm_p

    !$omp parallel default(none) &
    !$omp shared(n_points,rfl,th,ph,x,y,z,ut,hhh,vlx,vly,vlz,qb,qrho,ye,sen,dom, &
    !$omp   ut_sample,hhh_sample,vlr_sample,ye_sample,sen_sample,dm_sample,qrho_sample) &
    !$omp private(xi,yi,zi,j,k,l,lv,j1,k1,l1,x1,x0,y1,y0,z1,z0, &
    !$omp   ut_p,hhh_p,vlx_p,vly_p,vlz_p,vlr_p,qb_p,qrho_p,ye_p,sen_p,dm_p)
    !$omp do
    do ip = 1, n_points

       xi = rfl*sin(th(ip))*cos(ph(ip))
       yi = rfl*sin(th(ip))*sin(ph(ip))
       zi = rfl*cos(th(ip))

!!! interpolate variables
       call coorindex3D(xi,yi,zi,j1,k1,l1,lv)
       j=j1-1
       k=k1-1
       l=l1-1
       
       x1 = (xi-x(j,lv))/(x(j1,lv)-x(j,lv))
       x0 = 1.d0-x1
       y1 = (yi-y(k,lv))/(y(k1,lv)-y(k,lv))
       y0 = 1.d0-y1
       z1 = (zi-z(l,lv))/(z(l1,lv)-z(l,lv))
       z0 = 1.d0-z1

       ut_p = x1*y1*z1* ut   (j1,k1,l1,lv) &
            + x0*y1*z1* ut   (j ,k1,l1,lv) &
            + x1*y0*z1* ut   (j1,k ,l1,lv) &
            + x0*y0*z1* ut   (j ,k ,l1,lv) &
            + x1*y1*z0* ut   (j1,k1,l ,lv) &
            + x0*y1*z0* ut   (j ,k1,l ,lv) &
            + x1*y0*z0* ut   (j1,k ,l ,lv) &
            + x0*y0*z0* ut   (j ,k ,l ,lv)

       hhh_p = x1*y1*z1* hhh   (j1,k1,l1,lv) &
             + x0*y1*z1* hhh   (j ,k1,l1,lv) &
             + x1*y0*z1* hhh   (j1,k ,l1,lv) &
             + x0*y0*z1* hhh   (j ,k ,l1,lv) &
             + x1*y1*z0* hhh   (j1,k1,l ,lv) &
             + x0*y1*z0* hhh   (j ,k1,l ,lv) &
             + x1*y0*z0* hhh   (j1,k ,l ,lv) &
             + x0*y0*z0* hhh   (j ,k ,l ,lv)

       vlx_p = x1*y1*z1* vlx   (j1,k1,l1,lv) &
             + x0*y1*z1* vlx   (j ,k1,l1,lv) &
             + x1*y0*z1* vlx   (j1,k ,l1,lv) &
             + x0*y0*z1* vlx   (j ,k ,l1,lv) &
             + x1*y1*z0* vlx   (j1,k1,l ,lv) &
             + x0*y1*z0* vlx   (j ,k1,l ,lv) &
             + x1*y0*z0* vlx   (j1,k ,l ,lv) &
             + x0*y0*z0* vlx   (j ,k ,l ,lv)
       vly_p = x1*y1*z1* vly   (j1,k1,l1,lv) &
             + x0*y1*z1* vly   (j ,k1,l1,lv) &
             + x1*y0*z1* vly   (j1,k ,l1,lv) &
             + x0*y0*z1* vly   (j ,k ,l1,lv) &
             + x1*y1*z0* vly   (j1,k1,l ,lv) &
             + x0*y1*z0* vly   (j ,k1,l ,lv) &
             + x1*y0*z0* vly   (j1,k ,l ,lv) &
             + x0*y0*z0* vly   (j ,k ,l ,lv)
       vlz_p = x1*y1*z1* vlz   (j1,k1,l1,lv) &
             + x0*y1*z1* vlz   (j ,k1,l1,lv) &
             + x1*y0*z1* vlz   (j1,k ,l1,lv) &
             + x0*y0*z1* vlz   (j ,k ,l1,lv) &
             + x1*y1*z0* vlz   (j1,k1,l ,lv) &
             + x0*y1*z0* vlz   (j ,k1,l ,lv) &
             + x1*y0*z0* vlz   (j1,k ,l ,lv) &
             + x0*y0*z0* vlz   (j ,k ,l ,lv)
       vlr_p = vlx_p*sin(th(ip))*cos(ph(ip)) &
             + vly_p*sin(th(ip))*sin(ph(ip)) &
             + vlz_p*cos(th(ip))
       
       qb_p = x1*y1*z1* qb (j1,k1,l1,lv) &
            + x0*y1*z1* qb (j ,k1,l1,lv) &
            + x1*y0*z1* qb (j1,k ,l1,lv) &
            + x0*y0*z1* qb (j ,k ,l1,lv) &
            + x1*y1*z0* qb (j1,k1,l ,lv) &
            + x0*y1*z0* qb (j ,k1,l ,lv) &
            + x1*y0*z0* qb (j1,k ,l ,lv) &
            + x0*y0*z0* qb (j ,k ,l ,lv)

       qrho_p= x1*y1*z1* qrho  (j1,k1,l1,lv) &
             + x0*y1*z1* qrho  (j ,k1,l1,lv) &
             + x1*y0*z1* qrho  (j1,k ,l1,lv) &
             + x0*y0*z1* qrho  (j ,k ,l1,lv) &
             + x1*y1*z0* qrho  (j1,k1,l ,lv) &
             + x0*y1*z0* qrho  (j ,k1,l ,lv) &
             + x1*y0*z0* qrho  (j1,k ,l ,lv) &
             + x0*y0*z0* qrho  (j ,k ,l ,lv)
       ye_p  = x1*y1*z1* ye    (j1,k1,l1,lv) &
             + x0*y1*z1* ye    (j ,k1,l1,lv) &
             + x1*y0*z1* ye    (j1,k ,l1,lv) &
             + x0*y0*z1* ye    (j ,k ,l1,lv) &
             + x1*y1*z0* ye    (j1,k1,l ,lv) &
             + x0*y1*z0* ye    (j ,k1,l ,lv) &
             + x1*y0*z0* ye    (j1,k ,l ,lv) &
             + x0*y0*z0* ye    (j ,k ,l ,lv)
       sen_p = x1*y1*z1* sen   (j1,k1,l1,lv) &
             + x0*y1*z1* sen   (j ,k1,l1,lv) &
             + x1*y0*z1* sen   (j1,k ,l1,lv) &
             + x0*y0*z1* sen   (j ,k ,l1,lv) &
             + x1*y1*z0* sen   (j1,k1,l ,lv) &
             + x0*y1*z0* sen   (j ,k1,l ,lv) &
             + x1*y0*z0* sen   (j1,k ,l ,lv) &
             + x0*y0*z0* sen   (j ,k ,l ,lv)
       
       dm_p  = rfl**2*dom(ip) * qb_p*abs(vlr_p)*v_uni

       
       ut_sample(ip)  = ut_p
       hhh_sample(ip) = hhh_p
       vlr_sample(ip) = vlr_p
       qrho_sample(ip)= qrho_p
       ye_sample(ip)  = ye_p
       sen_sample(ip) = sen_p
       dm_sample(ip)  = dm_p
       
    enddo
    !$omp end do
    !$omp end parallel
    
  end subroutine sample_flux
  
end module module_flux_sample
