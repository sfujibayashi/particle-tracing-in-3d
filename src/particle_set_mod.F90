module particle_set

  implicit none
  
  integer :: n_points
  real(8),allocatable :: th(:),ph(:), dom(:)
  real(8) :: dth
contains
  
  subroutine init_angle(nside,n_pset)
#include "macro.h"
    use const
    use io
    integer,intent(in) :: nside
    integer,intent(out) :: n_pset
    
    integer :: ip
    integer :: p,i,j,pp
    real(8) :: z,phi,s,phh
    integer :: npix, nnorth, neq, nring

    write(6,*) 
    write(6,*) "Initializing for flux-based particle setting..."

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
    n_pset = n_points

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

             ! write(6,'(3i5,2es12.4," ",a)') p,i,j,z,phi,"equatorial"
          else
        
             write(6,*) "Something is wrong, stop."
             stop
             
          endif
       endif
       
       th(ip) = acos(z)
       ph(ip) = phi
    enddo

    
#ifdef FULL
    ! set southern points
    do ip=nnorth+1,n_points
       p = ip-1
       p = npix-1-p

       phh = dble(p+1)/2d0
       i = int(sqrt(phh - sqrt(dble(int(phh))) )) + 1
       
       if(i<nside)then

          j = p + 1 - 2*i*(i-1)

          z =-1d0 +dble(i)**2/(3d0*dble(nside)**2)
          s = 1d0
          phi = - pi/(2d0*dble(i)) * (dble(j) - s/2d0) + 2d0*pi
          write(99,'(3i5,2es12.4," ",a)') p,i,j,z,phi, "south-pole"
       else

          pp = p - 2*nside*(nside-1)
          i = int( dble(pp)/dble(4*nside) ) + nside

          ! check whether it is actually on the equatorial belt
          if(nside <= i .and. i <= 2*nside)then

             j = mod(pp, (4*nside)) + 1
             
             z=-4d0/3d0 + 2d0*dble(i)/(3d0*dble(nside))
             s = dble(mod((i-nside+1), 2))
             phi =-pi/(2d0*dble(nside)) * (dble(j) - s/2d0)  + 2d0*pi

             write(99,'(3i5,2es12.4," ",a)') p,i,j,z,phi,"south-equatorial"
          else
             
             write(6,*) "Something is wrong, stop."
             stop
             
          endif
       endif
       
       th(ip) = acos(z)
       ph(ip) = phi
       
    enddo
#endif


    do ip = 1,n_points
       dom(ip) = 4d0*pi/dble(npix)
    enddo

#ifdef MIRROR
    do ip = 1,nnorth-neq
       dom(ip) = dom(ip) * 2d0
    enddo
#endif
    
    block
      integer :: unit
      open(newunit=unit,file=trim(dir_out)//"/angle_info.dat",status="replace")
      write(unit,*) "#",sum(dom(:))/4.d0/pi,n_points
      do ip=1,n_points
         write(unit,'(i5,99es12.4)') ip,th(ip),ph(ip),sin(th(ip))*cos(ph(ip)),sin(th(ip))*sin(ph(ip)),cos(th(ip)),dom(ip)
      enddo
      close(unit)
    end block


  end subroutine init_angle


  !subroutine set_new_particle(rfl,ipu,n_th,dt,it_skip,it_skip_out,it_skip_pset)
  subroutine set_new_particle(ips,rfl,dt,v_average,it_skip,it_skip_out,it_skip_pset,m_max,m_min,m_average,v_max,v_min,np_set)
    
    use const, only : pi
    use simdata3D
    use particle_data
    use unit
    use module_eos

    implicit none

    integer,intent(inout) :: ips
    real(8),intent(in) :: rfl,dt
    
    real(8),intent(out) :: v_average,m_average,m_max,m_min,v_max,v_min
    integer,intent(in) :: it_skip_out,it_skip
    integer,intent(out) :: it_skip_pset,np_set
    
    integer :: ipnt,j,k,l,lv,j1,k1,l1,ip
    real(8) :: x0,x1,y0,y1,z0,z1,ut_i,hhh_i,vlx_i,vly_i,vlz_i,vlr_i,qb_i,xi,yi,zi,dx,dy,dz

    real(8) :: dt_pset
    !integer :: np_set

    np_set = 0
    v_average=0d0
    v_max=0d0
    v_min=1d99

    do ipnt = 1, n_points

       xi = rfl*sin(th(ipnt))*cos(ph(ipnt))
       yi = rfl*sin(th(ipnt))*sin(ph(ipnt))
       zi = rfl*cos(th(ipnt))

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

       ut_i = x1*y1*z1* ut   (j1,k1,l1,lv) &
            + x0*y1*z1* ut   (j ,k1,l1,lv) &
            + x1*y0*z1* ut   (j1,k ,l1,lv) &
            + x0*y0*z1* ut   (j ,k ,l1,lv) &
            + x1*y1*z0* ut   (j1,k1,l ,lv) &
            + x0*y1*z0* ut   (j ,k1,l ,lv) &
            + x1*y0*z0* ut   (j1,k ,l ,lv) &
            + x0*y0*z0* ut   (j ,k ,l ,lv)

       hhh_i = x1*y1*z1* hhh   (j1,k1,l1,lv) &
             + x0*y1*z1* hhh   (j ,k1,l1,lv) &
             + x1*y0*z1* hhh   (j1,k ,l1,lv) &
             + x0*y0*z1* hhh   (j ,k ,l1,lv) &
             + x1*y1*z0* hhh   (j1,k1,l ,lv) &
             + x0*y1*z0* hhh   (j ,k1,l ,lv) &
             + x1*y0*z0* hhh   (j1,k ,l ,lv) &
             + x0*y0*z0* hhh   (j ,k ,l ,lv)

       vlx_i = x1*y1*z1* vlx   (j1,k1,l1,lv) &
             + x0*y1*z1* vlx   (j ,k1,l1,lv) &
             + x1*y0*z1* vlx   (j1,k ,l1,lv) &
             + x0*y0*z1* vlx   (j ,k ,l1,lv) &
             + x1*y1*z0* vlx   (j1,k1,l ,lv) &
             + x0*y1*z0* vlx   (j ,k1,l ,lv) &
             + x1*y0*z0* vlx   (j1,k ,l ,lv) &
             + x0*y0*z0* vlx   (j ,k ,l ,lv)
       vly_i = x1*y1*z1* vly   (j1,k1,l1,lv) &
             + x0*y1*z1* vly   (j ,k1,l1,lv) &
             + x1*y0*z1* vly   (j1,k ,l1,lv) &
             + x0*y0*z1* vly   (j ,k ,l1,lv) &
             + x1*y1*z0* vly   (j1,k1,l ,lv) &
             + x0*y1*z0* vly   (j ,k1,l ,lv) &
             + x1*y0*z0* vly   (j1,k ,l ,lv) &
             + x0*y0*z0* vly   (j ,k ,l ,lv)
       vlz_i = x1*y1*z1* vlz   (j1,k1,l1,lv) &
             + x0*y1*z1* vlz   (j ,k1,l1,lv) &
             + x1*y0*z1* vlz   (j1,k ,l1,lv) &
             + x0*y0*z1* vlz   (j ,k ,l1,lv) &
             + x1*y1*z0* vlz   (j1,k1,l ,lv) &
             + x0*y1*z0* vlz   (j ,k1,l ,lv) &
             + x1*y0*z0* vlz   (j1,k ,l ,lv) &
             + x0*y0*z0* vlz   (j ,k ,l ,lv)
       vlr_i = vlx_i*sin(th(ipnt))*cos(ph(ipnt)) &
             + vly_i*sin(th(ipnt))*sin(ph(ipnt)) &
             + vlz_i*cos(th(ipnt))

       ! if( -ut_i-1.d0 > 0.d0 .and. vlr_i > 0.d0)then
       ! if( -ut_i*hhh_i - hhh_min > 0.d0 .and. vlr_i > 0.d0)then
       if( vlr_i > 0.d0)then

          ips = ips + 1
          ip = ips
          np_set = np_set + 1

          flag_evol(ip) = 1
          
          qb_i = x1*y1*z1* qb (j1,k1,l1,lv) &
               + x0*y1*z1* qb (j ,k1,l1,lv) &
               + x1*y0*z1* qb (j1,k ,l1,lv) &
               + x0*y0*z1* qb (j ,k ,l1,lv) &
               + x1*y1*z0* qb (j1,k1,l ,lv) &
               + x0*y1*z0* qb (j ,k1,l ,lv) &
               + x1*y0*z0* qb (j1,k ,l ,lv) &
               + x0*y0*z0* qb (j ,k ,l ,lv)

          var_p(index_dm,ip) = abs(dt) *rfl**2*dom(ipnt) * qb_i*abs(vlr_i)*v_uni
          
          var_p(index_ut1,ip) = ut_i + 1.d0
          var_p(index_hut,ip) = ut_i*hhh_i + hhh_min
          
          var_p(index_x,ip) = xi
          var_p(index_y,ip) = yi
          var_p(index_z,ip) = zi
          

          v_average = v_average + vlr_i
          v_max=max(v_max,vlr_i)
          v_min=min(v_min,vlr_i)

       endif

    enddo
    

    if(np_set > 0)then
       v_average = v_average / dble(np_set)
       
       dt_pset = rfl*dth / abs(max(0.05,v_average)*v_uni)
       it_skip_pset = int(dt_pset/abs(dt))


       ! do while(mod(it_skip_pset,it_skip_out)>0)
       !    it_skip_pset=it_skip_pset+1
       ! enddo
       
       !write(6,*) it_skip_pset

       m_max=0d0
       m_min=1d99
       m_average=0d0
       do ip=ips-np_set+1,ips
          !dm_p(ip) = dm_p(ip) * dble(it_skip_pset)/dble(it_skip)
          var_p(index_dm,ip) = var_p(index_dm,ip) * dble(it_skip_pset)
          m_average=m_average+var_p(index_dm,ip)
          m_max=max(m_max,var_p(index_dm,ip))
          m_min=min(m_min,var_p(index_dm,ip))
       enddo
       m_average=m_average/dble(np_set)

    else

       write(6,*) "No particles set."
       v_average = 0.05d0
       dt_pset = rfl*dth / (max(0.05,v_average)*v_uni)
       it_skip_pset = int(dt_pset/abs(dt))

       ! do while(mod(it_skip_pset,it_skip_out)>0)
       !    it_skip_pset=it_skip_pset+1
       ! enddo
       m_max=0d0
       m_min=0d0
       m_average=0d0
       v_max=0d0
       v_min=0d0
       v_average=0d0
    endif

  end subroutine set_new_particle
  
  subroutine set_particle_data(ipu)
    use simdata3D
    use particle_data
    integer,intent(in) :: ipu
    integer :: ip
    real(8) :: xi,yi,zi

    integer :: j,k,l,lv, j1,k1,l1
    real(8) :: x0,y0,z0,x1,y1,z1

    do ip=1,ipu
       if(flag_evol(ip)==1)then

          xi=var_p(index_x,ip)
          yi=var_p(index_y,ip)
          zi=var_p(index_z,ip)
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
          
          var_p(index_vlx,ip) = &
                 x1*y1*z1* vlx   (j1,k1,l1,lv) &
               + x0*y1*z1* vlx   (j ,k1,l1,lv) &
               + x1*y0*z1* vlx   (j1,k ,l1,lv) &
               + x0*y0*z1* vlx   (j ,k ,l1,lv) &
               + x1*y1*z0* vlx   (j1,k1,l ,lv) &
               + x0*y1*z0* vlx   (j ,k1,l ,lv) &
               + x1*y0*z0* vlx   (j1,k ,l ,lv) &
               + x0*y0*z0* vlx   (j ,k ,l ,lv)
          var_p(index_vly,ip) = &
                 x1*y1*z1* vly   (j1,k1,l1,lv) &
               + x0*y1*z1* vly   (j ,k1,l1,lv) &
               + x1*y0*z1* vly   (j1,k ,l1,lv) &
               + x0*y0*z1* vly   (j ,k ,l1,lv) &
               + x1*y1*z0* vly   (j1,k1,l ,lv) &
               + x0*y1*z0* vly   (j ,k1,l ,lv) &
               + x1*y0*z0* vly   (j1,k ,l ,lv) &
               + x0*y0*z0* vly   (j ,k ,l ,lv)
          var_p(index_vlz,ip) = &
                 x1*y1*z1* vlz   (j1,k1,l1,lv) &
               + x0*y1*z1* vlz   (j ,k1,l1,lv) &
               + x1*y0*z1* vlz   (j1,k ,l1,lv) &
               + x0*y0*z1* vlz   (j ,k ,l1,lv) &
               + x1*y1*z0* vlz   (j1,k1,l ,lv) &
               + x0*y1*z0* vlz   (j ,k1,l ,lv) &
               + x1*y0*z0* vlz   (j1,k ,l ,lv) &
               + x0*y0*z0* vlz   (j ,k ,l ,lv)
          var_p(index_rho,ip)= &
                 x1*y1*z1* qrho  (j1,k1,l1,lv) &
               + x0*y1*z1* qrho  (j ,k1,l1,lv) &
               + x1*y0*z1* qrho  (j1,k ,l1,lv) &
               + x0*y0*z1* qrho  (j ,k ,l1,lv) &
               + x1*y1*z0* qrho  (j1,k1,l ,lv) &
               + x0*y1*z0* qrho  (j ,k1,l ,lv) &
               + x1*y0*z0* qrho  (j1,k ,l ,lv) &
               + x0*y0*z0* qrho  (j ,k ,l ,lv)
          var_p (index_ye,ip) = &
                 x1*y1*z1* ye    (j1,k1,l1,lv) &
               + x0*y1*z1* ye    (j ,k1,l1,lv) &
               + x1*y0*z1* ye    (j1,k ,l1,lv) &
               + x0*y0*z1* ye    (j ,k ,l1,lv) &
               + x1*y1*z0* ye    (j1,k1,l ,lv) &
               + x0*y1*z0* ye    (j ,k1,l ,lv) &
               + x1*y0*z0* ye    (j1,k ,l ,lv) &
               + x0*y0*z0* ye    (j ,k ,l ,lv)
          var_p(index_tem,ip) = &
                 x1*y1*z1* tem   (j1,k1,l1,lv) &
               + x0*y1*z1* tem   (j ,k1,l1,lv) &
               + x1*y0*z1* tem   (j1,k ,l1,lv) &
               + x0*y0*z1* tem   (j ,k ,l1,lv) &
               + x1*y1*z0* tem   (j1,k1,l ,lv) &
               + x0*y1*z0* tem   (j ,k1,l ,lv) &
               + x1*y0*z0* tem   (j1,k ,l ,lv) &
               + x0*y0*z0* tem   (j ,k ,l ,lv)
          var_p(index_sen,ip) = &
                 x1*y1*z1* sen   (j1,k1,l1,lv) &
               + x0*y1*z1* sen   (j ,k1,l1,lv) &
               + x1*y0*z1* sen   (j1,k ,l1,lv) &
               + x0*y0*z1* sen   (j ,k ,l1,lv) &
               + x1*y1*z0* sen   (j1,k1,l ,lv) &
               + x0*y1*z0* sen   (j ,k1,l ,lv) &
               + x1*y0*z0* sen   (j1,k ,l ,lv) &
               + x0*y0*z0* sen   (j ,k ,l ,lv)
          
          var_p(index_hhh,ip) = &
                 x1*y1*z1* hhh   (j1,k1,l1,lv) &
               + x0*y1*z1* hhh   (j ,k1,l1,lv) &
               + x1*y0*z1* hhh   (j1,k ,l1,lv) &
               + x0*y0*z1* hhh   (j ,k ,l1,lv) &
               + x1*y1*z0* hhh   (j1,k1,l ,lv) &
               + x0*y1*z0* hhh   (j ,k1,l ,lv) &
               + x1*y0*z0* hhh   (j1,k ,l ,lv) &
               + x0*y0*z0* hhh   (j ,k ,l ,lv)
          var_p (index_ut,ip) = &
                 x1*y1*z1* ut   (j1,k1,l1,lv) &
               + x0*y1*z1* ut   (j ,k1,l1,lv) &
               + x1*y0*z1* ut   (j1,k ,l1,lv) &
               + x0*y0*z1* ut   (j ,k ,l1,lv) &
               + x1*y1*z0* ut   (j1,k1,l ,lv) &
               + x0*y1*z0* ut   (j ,k1,l ,lv) &
               + x1*y0*z0* ut   (j1,k ,l ,lv) &
               + x0*y0*z0* ut   (j ,k ,l ,lv)

          var_p (index_qb,ip) = &
                 x1*y1*z1* qb   (j1,k1,l1,lv) &
               + x0*y1*z1* qb   (j ,k1,l1,lv) &
               + x1*y0*z1* qb   (j1,k ,l1,lv) &
               + x0*y0*z1* qb   (j ,k ,l1,lv) &
               + x1*y1*z0* qb   (j1,k1,l ,lv) &
               + x0*y1*z0* qb   (j ,k1,l ,lv) &
               + x1*y0*z0* qb   (j1,k ,l ,lv) &
               + x0*y0*z0* qb   (j ,k ,l ,lv)
          
       endif
    enddo
    
    
    
  end subroutine set_particle_data
  
end module particle_set
