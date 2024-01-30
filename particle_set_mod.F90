module particle_set

  implicit none
  
  integer :: n_points
  real(8),allocatable :: th(:),ph(:), dom(:)
  real(8) :: dth
contains
  
  subroutine init_angle(nth_pset,n_pset)
#include "macro.h"
    use const
    use io
    integer,intent(in) :: nth_pset
    integer,intent(out) :: n_pset
    integer :: ith,iph,ip,n_th, n_ph, ip_prev_th
    real(8) :: th1,th2
    real(8) :: dph,dom_min,dph_prev
    real(8) :: dcos_prev,dcos
    integer :: unit

    n_th = nth_pset
#ifdef FULL
    dth = pi/dble(n_th)
#else
    dth = pi/2.d0/dble(n_th)
#endif

    dom_min = 4.d0*pi/dble(n_th**2)
    n_ph=2
    ! dcos_prev=0.d0
    n_points = 0
    do ith=1,n_th
       th1 = dth*dble(ith-1)
       th2 = dth*dble(ith  )
       
       dcos= cos(th1)-cos(th2)
       ! if ( dcos_prev/dble(n_ph) < dcos/dble(n_ph+1) )then
       !n_ph = max(n_ph+1,int(2.d0*pi*(cos(th1)-cos(th2))/dom_min)+1)
       ! else
       !    n_ph = n_ph
       ! endif
       n_ph = int(2.d0*pi*(cos(th1)-cos(th2))/dom_min)+1
       
       dph = 2.d0*pi/dble(n_ph)
       ! write(6,'(2i5,99es12.4)') ith,n_ph,dph,dcos*dph
       n_points = n_points + n_ph

       !dcos_prev=dcos
    enddo

    n_pset = n_points
    allocate(th(n_points),ph(n_points), dom(n_points))
    
    n_ph=2
    dph_prev=0.d0
    ip = 1
    do ith=1,n_th
      
       
       th1 = dth*dble(ith-1)
       th2 = dth*dble(ith  )
       
       !n_ph = max(n_ph+1,int(2.d0*pi*(cos(th1)-cos(th2))/dom_min)+1)
       n_ph = int(2.d0*pi*(cos(th1)-cos(th2))/dom_min)+1
       dph = 2.d0*pi/dble(n_ph)
       do iph=1,n_ph
          th(ip) = dth*(dble(ith-1)+0.5d0)
          ph(ip) = dph*dble(iph-1)
#ifdef FULL
          dom(ip)=dph*(cos(th1)-cos(th2))
#else
          dom(ip)=dph*(cos(th1)-cos(th2))*2.d0
#endif
          ip = ip + 1
       enddo
       dph_prev = dph
    enddo

    open(newunit=unit,file=trim(dir_out)//"/angle_info.dat",status="replace")
    write(unit,*) "#",sum(dom(:))/4.d0/pi,n_points
    do ip=1,n_points
       write(unit,'(i5,99es12.4)') ip,th(ip),ph(ip),sin(th(ip))*cos(ph(ip)),sin(th(ip))*sin(ph(ip)),cos(th(ip)),dom(ip)
    enddo
    close(unit)

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

          dm_p(ip) = abs(dt) *rfl**2*dom(ipnt) * qb_i*abs(vlr_i)*v_uni
          
          ut1_p(ip) = ut_i + 1.d0
          hut_p(ip) = ut_i*hhh_i + hhh_min
          
          x_p(ip) = xi
          y_p(ip) = yi
          z_p(ip) = zi
          

          v_average = v_average + vlr_i
          v_max=max(v_max,vlr_i)
          v_min=min(v_min,vlr_i)

       endif

    enddo
    

    if(np_set > 0)then
       v_average = v_average / dble(np_set)
       
       dt_pset = rfl*dth / abs(v_average*v_uni)
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
          dm_p(ip) = dm_p(ip) * dble(it_skip_pset)
          m_average=m_average+dm_p(ip)
          m_max=max(m_max,dm_p(ip))
          m_min=min(m_min,dm_p(ip))
       enddo
       m_average=m_average/dble(np_set)

    else

       write(6,*) "No particles set."
       v_average = 0.05d0
       dt_pset = rfl*dth / (v_average*v_uni)
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
    real(8) :: www

    do ip=1,ipu
       if(flag_evol(ip)==1)then

          xi=x_p(ip)
          yi=y_p(ip)
          zi=z_p(ip)
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
          
          vlx_p(ip) = x1*y1*z1* vlx   (j1,k1,l1,lv) &
                    + x0*y1*z1* vlx   (j ,k1,l1,lv) &
                    + x1*y0*z1* vlx   (j1,k ,l1,lv) &
                    + x0*y0*z1* vlx   (j ,k ,l1,lv) &
                    + x1*y1*z0* vlx   (j1,k1,l ,lv) &
                    + x0*y1*z0* vlx   (j ,k1,l ,lv) &
                    + x1*y0*z0* vlx   (j1,k ,l ,lv) &
                    + x0*y0*z0* vlx   (j ,k ,l ,lv)
          vly_p(ip) = x1*y1*z1* vly   (j1,k1,l1,lv) &
                    + x0*y1*z1* vly   (j ,k1,l1,lv) &
                    + x1*y0*z1* vly   (j1,k ,l1,lv) &
                    + x0*y0*z1* vly   (j ,k ,l1,lv) &
                    + x1*y1*z0* vly   (j1,k1,l ,lv) &
                    + x0*y1*z0* vly   (j ,k1,l ,lv) &
                    + x1*y0*z0* vly   (j1,k ,l ,lv) &
                    + x0*y0*z0* vly   (j ,k ,l ,lv)
          vlz_p(ip) = x1*y1*z1* vlz   (j1,k1,l1,lv) &
                    + x0*y1*z1* vlz   (j ,k1,l1,lv) &
                    + x1*y0*z1* vlz   (j1,k ,l1,lv) &
                    + x0*y0*z1* vlz   (j ,k ,l1,lv) &
                    + x1*y1*z0* vlz   (j1,k1,l ,lv) &
                    + x0*y1*z0* vlz   (j ,k1,l ,lv) &
                    + x1*y0*z0* vlz   (j1,k ,l ,lv) &
                    + x0*y0*z0* vlz   (j ,k ,l ,lv)
          qrho_p(ip)= x1*y1*z1* qrho  (j1,k1,l1,lv) &
                    + x0*y1*z1* qrho  (j ,k1,l1,lv) &
                    + x1*y0*z1* qrho  (j1,k ,l1,lv) &
                    + x0*y0*z1* qrho  (j ,k ,l1,lv) &
                    + x1*y1*z0* qrho  (j1,k1,l ,lv) &
                    + x0*y1*z0* qrho  (j ,k1,l ,lv) &
                    + x1*y0*z0* qrho  (j1,k ,l ,lv) &
                    + x0*y0*z0* qrho  (j ,k ,l ,lv)
          ye_p (ip) = x1*y1*z1* ye    (j1,k1,l1,lv) &
                    + x0*y1*z1* ye    (j ,k1,l1,lv) &
                    + x1*y0*z1* ye    (j1,k ,l1,lv) &
                    + x0*y0*z1* ye    (j ,k ,l1,lv) &
                    + x1*y1*z0* ye    (j1,k1,l ,lv) &
                    + x0*y1*z0* ye    (j ,k1,l ,lv) &
                    + x1*y0*z0* ye    (j1,k ,l ,lv) &
                    + x0*y0*z0* ye    (j ,k ,l ,lv)
          tem_p(ip) = x1*y1*z1* tem   (j1,k1,l1,lv) &
                    + x0*y1*z1* tem   (j ,k1,l1,lv) &
                    + x1*y0*z1* tem   (j1,k ,l1,lv) &
                    + x0*y0*z1* tem   (j ,k ,l1,lv) &
                    + x1*y1*z0* tem   (j1,k1,l ,lv) &
                    + x0*y1*z0* tem   (j ,k1,l ,lv) &
                    + x1*y0*z0* tem   (j1,k ,l ,lv) &
                    + x0*y0*z0* tem   (j ,k ,l ,lv)
          sen_p(ip) = x1*y1*z1* sen   (j1,k1,l1,lv) &
                    + x0*y1*z1* sen   (j ,k1,l1,lv) &
                    + x1*y0*z1* sen   (j1,k ,l1,lv) &
                    + x0*y0*z1* sen   (j ,k ,l1,lv) &
                    + x1*y1*z0* sen   (j1,k1,l ,lv) &
                    + x0*y1*z0* sen   (j ,k1,l ,lv) &
                    + x1*y0*z0* sen   (j1,k ,l ,lv) &
                    + x0*y0*z0* sen   (j ,k ,l ,lv)
          tem_p(ip) = x1*y1*z1* tem   (j1,k1,l1,lv) &
                    + x0*y1*z1* tem   (j ,k1,l1,lv) &
                    + x1*y0*z1* tem   (j1,k ,l1,lv) &
                    + x0*y0*z1* tem   (j ,k ,l1,lv) &
                    + x1*y1*z0* tem   (j1,k1,l ,lv) &
                    + x0*y1*z0* tem   (j ,k1,l ,lv) &
                    + x1*y0*z0* tem   (j1,k ,l ,lv) &
                    + x0*y0*z0* tem   (j ,k ,l ,lv)
          
          hhh_p(ip) = x1*y1*z1* hhh   (j1,k1,l1,lv) &
                    + x0*y1*z1* hhh   (j ,k1,l1,lv) &
                    + x1*y0*z1* hhh   (j1,k ,l1,lv) &
                    + x0*y0*z1* hhh   (j ,k ,l1,lv) &
                    + x1*y1*z0* hhh   (j1,k1,l ,lv) &
                    + x0*y1*z0* hhh   (j ,k1,l ,lv) &
                    + x1*y0*z0* hhh   (j1,k ,l ,lv) &
                    + x0*y0*z0* hhh   (j ,k ,l ,lv)
          ut_p (ip) = x1*y1*z1* ut   (j1,k1,l1,lv) &
                    + x0*y1*z1* ut   (j ,k1,l1,lv) &
                    + x1*y0*z1* ut   (j1,k ,l1,lv) &
                    + x0*y0*z1* ut   (j ,k ,l1,lv) &
                    + x1*y1*z0* ut   (j1,k1,l ,lv) &
                    + x0*y1*z0* ut   (j ,k1,l ,lv) &
                    + x1*y0*z0* ut   (j1,k ,l ,lv) &
                    + x0*y0*z0* ut   (j ,k ,l ,lv)

          www = 1d0/sqrt(1d0-vlx_p(ip)*vlx_p(ip)-vly_p(ip)*vly_p(ip)-vlz_p(ip)*vlz_p(ip))
          qb_p(ip) = www*qrho_p(ip)
       endif
    enddo
    
    
    
  end subroutine set_particle_data
  
end module particle_set
