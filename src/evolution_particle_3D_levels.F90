subroutine evolution_particle_3D_levels(ipu,substep_max,lv_evolved_min,lv_evolved_max,np_evolved,evolution_finished)
#include "macro.h"
  use simdata3D
  use particle_data
  use unit
  !$ use omp_lib

  implicit none
  ! real(8),intent(in) :: time,time_prv
  integer,intent(in) :: ipu
  integer,intent(out) :: substep_max
  integer,intent(in) :: lv_evolved_min,lv_evolved_max
  logical,intent(inout) :: evolution_finished(ipu)
  integer,intent(out) :: np_evolved

  integer :: ip
  real(8) :: dt,t,xtmp,ytmp,ztmp,vxtmp,vytmp,vztmp
  integer :: j,k,l,lv,j1,k1,l1,lv_p
  real(8) :: x0,y0,z0,x1,y1,z1,t0,t1
  real(8) :: dx,dy,dz,dt_p,dtt
  integer :: nt_sub,itr,its
  real(8) :: xm,ym,zm,vxm,vym,vzm,xp,yp,zp,xc,yc,zc,vxc,vyc,vzc,vxdx,vxdy,vxdz,vydx,vydy,vydz,vzdx,vzdy,vzdz,dxc,dyc,dzc
  real(8) :: fx,fy,fz,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,det,fixx,fixy,fixz,fiyx,fiyy,fiyz,fizx,fizy,fizz

  real(8) :: dxdt, dydt, dzdt, dtb, dt_xb, dt_yb, dt_zb

  logical,allocatable :: evolved_in_this_level(:)
  ! integer :: np_evolved

  integer :: ip_test
  
  ip_test = 50

  ! dt = time-time_prv

  allocate(evolved_in_this_level(ipu))

  evolved_in_this_level(:) = .false.
  do ip=1,ipu

     if(flag_evol(ip)==1)then
        ! evolve the particle
        t    =  t_p(ip)
        xtmp =  x_p(ip)
        ytmp =  y_p(ip)
        ztmp =  z_p(ip)

        call coorindex3D(xtmp,ytmp,ztmp,j1,k1,l1,lv_p)
        if(lv_evolved_min<=lv_p .and. lv_p<=lv_evolved_max)then
           evolved_in_this_level(ip) = .true.
        elseif(remain(ip)==1)then
           evolved_in_this_level(ip) = .true.
        endif
        
     endif

     ! write(6,*) ip, evolved_in_this_level(ip)
     
  enddo
  
  evolution_finished(:) = .false.
  remain(:) = 0
  
  substep_max = 0
  !$omp parallel default(none) &
  !$omp shared(ipu,flag_evol,x_p,y_p,z_p,x,y,z,vlx_b,vly_b,vlz_b,vlx,vly,vlz,ju,jd,ku,kd,lu,ld,lv_min,&
  !$omp   lv_evolved_min,lv_evolved_max,t_p,time_level,remain,evolved_in_this_level,evolution_finished,ip_test) &
  !$omp private(t,xtmp,ytmp,ztmp,j1,k1,l1,lv,j,k,l,x1,x0,y1,y0,z1,z0,t1,t0, &
  !$omp   vxtmp,vytmp,vztmp,dx,dy,dz,dt_p,nt_sub,dtt,xm,ym,zm,vxm,vym,vzm,xp,yp,zp, &
  !$omp   xc,yc,zc,vxc,vyc,vzc,vxdx,vxdy,vxdz,vydx,vydy,vydz,vzdx,vzdy,vzdz, &
  !$omp   fx,fy,fz,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,det,fixx,fixy,fixz,fiyx,fiyy,fiyz,fizx,fizy,fizz, &
  !$omp   dxc,dyc,dzc,&
  !$omp   lv_p,dxdt, dydt, dzdt, dtb, dt_xb, dt_yb, dt_zb, dt) &
  !$omp reduction(max:substep_max)
  !$omp do
  do ip=1,ipu
     if(evolved_in_this_level(ip))then
        ! evolve the particle
        t    =  t_p(ip)
        xtmp =  x_p(ip)
        ytmp =  y_p(ip)
        ztmp =  z_p(ip)

        call coorindex3D(xtmp,ytmp,ztmp,j1,k1,l1,lv_p)
        ! if(lv_p < lv_evolved_min .or. lv_evolved_max < lv_p) goto 11
        
        lv = lv_p
        dt = time_level(lv)-t
        
        j=j1-1
        k=k1-1
        l=l1-1

        x1 = (xtmp-x(j,lv))/(x(j1,lv)-x(j,lv))
        x0 = 1.d0-x1
        y1 = (ytmp-y(k,lv))/(y(k1,lv)-y(k,lv))
        y0 = 1.d0-y1
        z1 = (ztmp-z(l,lv))/(z(l1,lv)-z(l,lv))
        z0 = 1.d0-z1
        
        vxtmp = x1*y1*z1* vlx_b(j1,k1,l1,lv) &
              + x0*y1*z1* vlx_b(j ,k1,l1,lv) &
              + x1*y0*z1* vlx_b(j1,k ,l1,lv) &
              + x0*y0*z1* vlx_b(j ,k ,l1,lv) &
              + x1*y1*z0* vlx_b(j1,k1,l ,lv) &
              + x0*y1*z0* vlx_b(j ,k1,l ,lv) &
              + x1*y0*z0* vlx_b(j1,k ,l ,lv) &
              + x0*y0*z0* vlx_b(j ,k ,l ,lv)
        vytmp = x1*y1*z1* vly_b(j1,k1,l1,lv) &
              + x0*y1*z1* vly_b(j ,k1,l1,lv) &
              + x1*y0*z1* vly_b(j1,k ,l1,lv) &
              + x0*y0*z1* vly_b(j ,k ,l1,lv) &
              + x1*y1*z0* vly_b(j1,k1,l ,lv) &
              + x0*y1*z0* vly_b(j ,k1,l ,lv) &
              + x1*y0*z0* vly_b(j1,k ,l ,lv) &
              + x0*y0*z0* vly_b(j ,k ,l ,lv)
        vztmp = x1*y1*z1* vlz_b(j1,k1,l1,lv) &
              + x0*y1*z1* vlz_b(j ,k1,l1,lv) &
              + x1*y0*z1* vlz_b(j1,k ,l1,lv) &
              + x0*y0*z1* vlz_b(j ,k ,l1,lv) &
              + x1*y1*z0* vlz_b(j1,k1,l ,lv) &
              + x0*y1*z0* vlz_b(j ,k1,l ,lv) &
              + x1*y0*z0* vlz_b(j1,k ,l ,lv) &
              + x0*y0*z0* vlz_b(j ,k ,l ,lv)

        dx = x(j1,lv)-x(j,lv)
        dy = y(k1,lv)-y(k,lv)
        dz = z(l1,lv)-z(l,lv)
        dt_p = min(dx/(abs(vxtmp)+1.d-20)/v_uni, &
                   dy/(abs(vytmp)+1.d-20)/v_uni, &
                   dz/(abs(vztmp)+1.d-20)/v_uni  )

        nt_sub = max(1,int(2.d0*abs(dt)/dt_p))

        ! write(1000+ip,*)
        ! write(1000+ip,'(3i4,99es12.4)') ip,lv_p, nt_sub, t, xtmp,ytmp,ztmp, time_level(lv),dt
        ! write(1000+ip,*)
        
        ! if(ip==ip_test)then
        !    write(6,*)
        !    write(6,'(3i4,99es12.4)') ip,lv_p, nt_sub, t, xtmp,ytmp,ztmp, time_level(lv),dt
        !    write(6,*)
        ! endif
        
        !write(6,*) ip, dt_p, dt, nt_sub

        substep_max = max( substep_max, nt_sub)
        dtt  = dt/dble(nt_sub)
        
        do its=1,nt_sub

           t  = t + dtt
!!! predictor-corrector method (2nd order)
           xm = xtmp
           ym = ytmp
           zm = ztmp
           vxm= vxtmp
           vym= vytmp
           vzm= vztmp

!!! predictor
           xp = xtmp + dtt*vxtmp*v_uni
           yp = ytmp + dtt*vytmp*v_uni
           zp = ztmp + dtt*vztmp*v_uni

!!! derive corrector
           xc = xp
           yc = yp
           zc = zp

!!! the method is reduced to Euler method if this iteration is commented-out.
           ! goto 10
           do itr = 1,10
              !interpolate variables
              call coorindex3D(xc,yc,zc,j1,k1,l1,lv)
              j=j1-1
              k=k1-1
              l=l1-1

              x1 = (xc-x(j,lv))/(x(j1,lv)-x(j,lv))
              x0 = 1.d0-x1
              y1 = (yc-y(k,lv))/(y(k1,lv)-y(k,lv))
              y0 = 1.d0-y1
              z1 = (zc-z(l,lv))/(z(l1,lv)-z(l,lv))
              z0 = 1.d0-z1
              t1 = (time_level(lv)-t)/dt
              t0 = 1.d0-t1

              vxc = t1*x1*y1*z1* vlx_b(j1,k1,l1,lv) &
                  + t1*x0*y1*z1* vlx_b(j ,k1,l1,lv) &
                  + t1*x1*y0*z1* vlx_b(j1,k ,l1,lv) &
                  + t1*x0*y0*z1* vlx_b(j ,k ,l1,lv) &
                  + t1*x1*y1*z0* vlx_b(j1,k1,l ,lv) &
                  + t1*x0*y1*z0* vlx_b(j ,k1,l ,lv) &
                  + t1*x1*y0*z0* vlx_b(j1,k ,l ,lv) &
                  + t1*x0*y0*z0* vlx_b(j ,k ,l ,lv) &
                  + t0*x1*y1*z1* vlx  (j1,k1,l1,lv) &
                  + t0*x0*y1*z1* vlx  (j ,k1,l1,lv) &
                  + t0*x1*y0*z1* vlx  (j1,k ,l1,lv) &
                  + t0*x0*y0*z1* vlx  (j ,k ,l1,lv) &
                  + t0*x1*y1*z0* vlx  (j1,k1,l ,lv) &
                  + t0*x0*y1*z0* vlx  (j ,k1,l ,lv) &
                  + t0*x1*y0*z0* vlx  (j1,k ,l ,lv) &
                  + t0*x0*y0*z0* vlx  (j ,k ,l ,lv)

              vyc = t1*x1*y1*z1* vly_b(j1,k1,l1,lv) &
                  + t1*x0*y1*z1* vly_b(j ,k1,l1,lv) &
                  + t1*x1*y0*z1* vly_b(j1,k ,l1,lv) &
                  + t1*x0*y0*z1* vly_b(j ,k ,l1,lv) &
                  + t1*x1*y1*z0* vly_b(j1,k1,l ,lv) &
                  + t1*x0*y1*z0* vly_b(j ,k1,l ,lv) &
                  + t1*x1*y0*z0* vly_b(j1,k ,l ,lv) &
                  + t1*x0*y0*z0* vly_b(j ,k ,l ,lv) &
                  + t0*x1*y1*z1* vly  (j1,k1,l1,lv) &
                  + t0*x0*y1*z1* vly  (j ,k1,l1,lv) &
                  + t0*x1*y0*z1* vly  (j1,k ,l1,lv) &
                  + t0*x0*y0*z1* vly  (j ,k ,l1,lv) &
                  + t0*x1*y1*z0* vly  (j1,k1,l ,lv) &
                  + t0*x0*y1*z0* vly  (j ,k1,l ,lv) &
                  + t0*x1*y0*z0* vly  (j1,k ,l ,lv) &
                  + t0*x0*y0*z0* vly  (j ,k ,l ,lv)

              vzc = t1*x1*y1*z1* vlz_b(j1,k1,l1,lv) &
                  + t1*x0*y1*z1* vlz_b(j ,k1,l1,lv) &
                  + t1*x1*y0*z1* vlz_b(j1,k ,l1,lv) &
                  + t1*x0*y0*z1* vlz_b(j ,k ,l1,lv) &
                  + t1*x1*y1*z0* vlz_b(j1,k1,l ,lv) &
                  + t1*x0*y1*z0* vlz_b(j ,k1,l ,lv) &
                  + t1*x1*y0*z0* vlz_b(j1,k ,l ,lv) &
                  + t1*x0*y0*z0* vlz_b(j ,k ,l ,lv) &
                  + t0*x1*y1*z1* vlz  (j1,k1,l1,lv) &
                  + t0*x0*y1*z1* vlz  (j ,k1,l1,lv) &
                  + t0*x1*y0*z1* vlz  (j1,k ,l1,lv) &
                  + t0*x0*y0*z1* vlz  (j ,k ,l1,lv) &
                  + t0*x1*y1*z0* vlz  (j1,k1,l ,lv) &
                  + t0*x0*y1*z0* vlz  (j ,k1,l ,lv) &
                  + t0*x1*y0*z0* vlz  (j1,k ,l ,lv) &
                  + t0*x0*y0*z0* vlz  (j ,k ,l ,lv)

!!! derivative of velocity wrt x,y,z
              vxdx=(t1*( 1.d0)*y1*z1* vlx_b(j1,k1,l1,lv) &
                  + t1*(-1.d0)*y1*z1* vlx_b(j ,k1,l1,lv) &
                  + t1*( 1.d0)*y0*z1* vlx_b(j1,k ,l1,lv) &
                  + t1*(-1.d0)*y0*z1* vlx_b(j ,k ,l1,lv) &
                  + t1*( 1.d0)*y1*z0* vlx_b(j1,k1,l ,lv) &
                  + t1*(-1.d0)*y1*z0* vlx_b(j ,k1,l ,lv) &
                  + t1*( 1.d0)*y0*z0* vlx_b(j1,k ,l ,lv) &
                  + t1*(-1.d0)*y0*z0* vlx_b(j ,k ,l ,lv) &
                  + t0*( 1.d0)*y1*z1* vlx  (j1,k1,l1,lv) &
                  + t0*(-1.d0)*y1*z1* vlx  (j ,k1,l1,lv) &
                  + t0*( 1.d0)*y0*z1* vlx  (j1,k ,l1,lv) &
                  + t0*(-1.d0)*y0*z1* vlx  (j ,k ,l1,lv) &
                  + t0*( 1.d0)*y1*z0* vlx  (j1,k1,l ,lv) &
                  + t0*(-1.d0)*y1*z0* vlx  (j ,k1,l ,lv) &
                  + t0*( 1.d0)*y0*z0* vlx  (j1,k ,l ,lv) &
                  + t0*(-1.d0)*y0*z0* vlx  (j ,k ,l ,lv) )/(x(j1,lv)-x(j,lv))
              vxdy=(t1*x1*( 1.d0)*z1* vlx_b(j1,k1,l1,lv) &
                  + t1*x0*( 1.d0)*z1* vlx_b(j ,k1,l1,lv) &
                  + t1*x1*(-1.d0)*z1* vlx_b(j1,k ,l1,lv) &
                  + t1*x0*(-1.d0)*z1* vlx_b(j ,k ,l1,lv) &
                  + t1*x1*( 1.d0)*z0* vlx_b(j1,k1,l ,lv) &
                  + t1*x0*( 1.d0)*z0* vlx_b(j ,k1,l ,lv) &
                  + t1*x1*(-1.d0)*z0* vlx_b(j1,k ,l ,lv) &
                  + t1*x0*(-1.d0)*z0* vlx_b(j ,k ,l ,lv) &
                  + t0*x1*( 1.d0)*z1* vlx  (j1,k1,l1,lv) &
                  + t0*x0*( 1.d0)*z1* vlx  (j ,k1,l1,lv) &
                  + t0*x1*(-1.d0)*z1* vlx  (j1,k ,l1,lv) &
                  + t0*x0*(-1.d0)*z1* vlx  (j ,k ,l1,lv) &
                  + t0*x1*( 1.d0)*z0* vlx  (j1,k1,l ,lv) &
                  + t0*x0*( 1.d0)*z0* vlx  (j ,k1,l ,lv) &
                  + t0*x1*(-1.d0)*z0* vlx  (j1,k ,l ,lv) &
                  + t0*x0*(-1.d0)*z0* vlx  (j ,k ,l ,lv) )/(y(k1,lv)-y(k,lv))
              vxdz=(t1*x1*y1*( 1.d0)* vlx_b(j1,k1,l1,lv) &
                  + t1*x0*y1*( 1.d0)* vlx_b(j ,k1,l1,lv) &
                  + t1*x1*y0*( 1.d0)* vlx_b(j1,k ,l1,lv) &
                  + t1*x0*y0*( 1.d0)* vlx_b(j ,k ,l1,lv) &
                  + t1*x1*y1*(-1.d0)* vlx_b(j1,k1,l ,lv) &
                  + t1*x0*y1*(-1.d0)* vlx_b(j ,k1,l ,lv) &
                  + t1*x1*y0*(-1.d0)* vlx_b(j1,k ,l ,lv) &
                  + t1*x0*y0*(-1.d0)* vlx_b(j ,k ,l ,lv) &
                  + t0*x1*y1*( 1.d0)* vlx  (j1,k1,l1,lv) &
                  + t0*x0*y1*( 1.d0)* vlx  (j ,k1,l1,lv) &
                  + t0*x1*y0*( 1.d0)* vlx  (j1,k ,l1,lv) &
                  + t0*x0*y0*( 1.d0)* vlx  (j ,k ,l1,lv) &
                  + t0*x1*y1*(-1.d0)* vlx  (j1,k1,l ,lv) &
                  + t0*x0*y1*(-1.d0)* vlx  (j ,k1,l ,lv) &
                  + t0*x1*y0*(-1.d0)* vlx  (j1,k ,l ,lv) &
                  + t0*x0*y0*(-1.d0)* vlx  (j ,k ,l ,lv) )/(z(l1,lv)-z(l,lv))

              vydx=(t1*( 1.d0)*y1*z1* vly_b(j1,k1,l1,lv) &
                  + t1*(-1.d0)*y1*z1* vly_b(j ,k1,l1,lv) &
                  + t1*( 1.d0)*y0*z1* vly_b(j1,k ,l1,lv) &
                  + t1*(-1.d0)*y0*z1* vly_b(j ,k ,l1,lv) &
                  + t1*( 1.d0)*y1*z0* vly_b(j1,k1,l ,lv) &
                  + t1*(-1.d0)*y1*z0* vly_b(j ,k1,l ,lv) &
                  + t1*( 1.d0)*y0*z0* vly_b(j1,k ,l ,lv) &
                  + t1*(-1.d0)*y0*z0* vly_b(j ,k ,l ,lv) &
                  + t0*( 1.d0)*y1*z1* vly  (j1,k1,l1,lv) &
                  + t0*(-1.d0)*y1*z1* vly  (j ,k1,l1,lv) &
                  + t0*( 1.d0)*y0*z1* vly  (j1,k ,l1,lv) &
                  + t0*(-1.d0)*y0*z1* vly  (j ,k ,l1,lv) &
                  + t0*( 1.d0)*y1*z0* vly  (j1,k1,l ,lv) &
                  + t0*(-1.d0)*y1*z0* vly  (j ,k1,l ,lv) &
                  + t0*( 1.d0)*y0*z0* vly  (j1,k ,l ,lv) &
                  + t0*(-1.d0)*y0*z0* vly  (j ,k ,l ,lv) )/(x(j1,lv)-x(j,lv))
              vydy=(t1*x1*( 1.d0)*z1* vly_b(j1,k1,l1,lv) &
                  + t1*x0*( 1.d0)*z1* vly_b(j ,k1,l1,lv) &
                  + t1*x1*(-1.d0)*z1* vly_b(j1,k ,l1,lv) &
                  + t1*x0*(-1.d0)*z1* vly_b(j ,k ,l1,lv) &
                  + t1*x1*( 1.d0)*z0* vly_b(j1,k1,l ,lv) &
                  + t1*x0*( 1.d0)*z0* vly_b(j ,k1,l ,lv) &
                  + t1*x1*(-1.d0)*z0* vly_b(j1,k ,l ,lv) &
                  + t1*x0*(-1.d0)*z0* vly_b(j ,k ,l ,lv) &
                  + t0*x1*( 1.d0)*z1* vly  (j1,k1,l1,lv) &
                  + t0*x0*( 1.d0)*z1* vly  (j ,k1,l1,lv) &
                  + t0*x1*(-1.d0)*z1* vly  (j1,k ,l1,lv) &
                  + t0*x0*(-1.d0)*z1* vly  (j ,k ,l1,lv) &
                  + t0*x1*( 1.d0)*z0* vly  (j1,k1,l ,lv) &
                  + t0*x0*( 1.d0)*z0* vly  (j ,k1,l ,lv) &
                  + t0*x1*(-1.d0)*z0* vly  (j1,k ,l ,lv) &
                  + t0*x0*(-1.d0)*z0* vly  (j ,k ,l ,lv) )/(y(k1,lv)-y(k,lv))
              vydz=(t1*x1*y1*( 1.d0)* vly_b(j1,k1,l1,lv) &
                  + t1*x0*y1*( 1.d0)* vly_b(j ,k1,l1,lv) &
                  + t1*x1*y0*( 1.d0)* vly_b(j1,k ,l1,lv) &
                  + t1*x0*y0*( 1.d0)* vly_b(j ,k ,l1,lv) &
                  + t1*x1*y1*(-1.d0)* vly_b(j1,k1,l ,lv) &
                  + t1*x0*y1*(-1.d0)* vly_b(j ,k1,l ,lv) &
                  + t1*x1*y0*(-1.d0)* vly_b(j1,k ,l ,lv) &
                  + t1*x0*y0*(-1.d0)* vly_b(j ,k ,l ,lv) &
                  + t0*x1*y1*( 1.d0)* vly  (j1,k1,l1,lv) &
                  + t0*x0*y1*( 1.d0)* vly  (j ,k1,l1,lv) &
                  + t0*x1*y0*( 1.d0)* vly  (j1,k ,l1,lv) &
                  + t0*x0*y0*( 1.d0)* vly  (j ,k ,l1,lv) &
                  + t0*x1*y1*(-1.d0)* vly  (j1,k1,l ,lv) &
                  + t0*x0*y1*(-1.d0)* vly  (j ,k1,l ,lv) &
                  + t0*x1*y0*(-1.d0)* vly  (j1,k ,l ,lv) &
                  + t0*x0*y0*(-1.d0)* vly  (j ,k ,l ,lv) )/(z(l1,lv)-z(l,lv))

              vzdx=(t1*( 1.d0)*y1*z1* vlz_b(j1,k1,l1,lv) &
                  + t1*(-1.d0)*y1*z1* vlz_b(j ,k1,l1,lv) &
                  + t1*( 1.d0)*y0*z1* vlz_b(j1,k ,l1,lv) &
                  + t1*(-1.d0)*y0*z1* vlz_b(j ,k ,l1,lv) &
                  + t1*( 1.d0)*y1*z0* vlz_b(j1,k1,l ,lv) &
                  + t1*(-1.d0)*y1*z0* vlz_b(j ,k1,l ,lv) &
                  + t1*( 1.d0)*y0*z0* vlz_b(j1,k ,l ,lv) &
                  + t1*(-1.d0)*y0*z0* vlz_b(j ,k ,l ,lv) &
                  + t0*( 1.d0)*y1*z1* vlz  (j1,k1,l1,lv) &
                  + t0*(-1.d0)*y1*z1* vlz  (j ,k1,l1,lv) &
                  + t0*( 1.d0)*y0*z1* vlz  (j1,k ,l1,lv) &
                  + t0*(-1.d0)*y0*z1* vlz  (j ,k ,l1,lv) &
                  + t0*( 1.d0)*y1*z0* vlz  (j1,k1,l ,lv) &
                  + t0*(-1.d0)*y1*z0* vlz  (j ,k1,l ,lv) &
                  + t0*( 1.d0)*y0*z0* vlz  (j1,k ,l ,lv) &
                  + t0*(-1.d0)*y0*z0* vlz  (j ,k ,l ,lv) )/(x(j1,lv)-x(j,lv))
              vzdy=(t1*x1*( 1.d0)*z1* vlz_b(j1,k1,l1,lv) &
                  + t1*x0*( 1.d0)*z1* vlz_b(j ,k1,l1,lv) &
                  + t1*x1*(-1.d0)*z1* vlz_b(j1,k ,l1,lv) &
                  + t1*x0*(-1.d0)*z1* vlz_b(j ,k ,l1,lv) &
                  + t1*x1*( 1.d0)*z0* vlz_b(j1,k1,l ,lv) &
                  + t1*x0*( 1.d0)*z0* vlz_b(j ,k1,l ,lv) &
                  + t1*x1*(-1.d0)*z0* vlz_b(j1,k ,l ,lv) &
                  + t1*x0*(-1.d0)*z0* vlz_b(j ,k ,l ,lv) &
                  + t0*x1*( 1.d0)*z1* vlz  (j1,k1,l1,lv) &
                  + t0*x0*( 1.d0)*z1* vlz  (j ,k1,l1,lv) &
                  + t0*x1*(-1.d0)*z1* vlz  (j1,k ,l1,lv) &
                  + t0*x0*(-1.d0)*z1* vlz  (j ,k ,l1,lv) &
                  + t0*x1*( 1.d0)*z0* vlz  (j1,k1,l ,lv) &
                  + t0*x0*( 1.d0)*z0* vlz  (j ,k1,l ,lv) &
                  + t0*x1*(-1.d0)*z0* vlz  (j1,k ,l ,lv) &
                  + t0*x0*(-1.d0)*z0* vlz  (j ,k ,l ,lv) )/(y(k1,lv)-y(k,lv))
              vzdz=(t1*x1*y1*( 1.d0)* vlz_b(j1,k1,l1,lv) &
                  + t1*x0*y1*( 1.d0)* vlz_b(j ,k1,l1,lv) &
                  + t1*x1*y0*( 1.d0)* vlz_b(j1,k ,l1,lv) &
                  + t1*x0*y0*( 1.d0)* vlz_b(j ,k ,l1,lv) &
                  + t1*x1*y1*(-1.d0)* vlz_b(j1,k1,l ,lv) &
                  + t1*x0*y1*(-1.d0)* vlz_b(j ,k1,l ,lv) &
                  + t1*x1*y0*(-1.d0)* vlz_b(j1,k ,l ,lv) &
                  + t1*x0*y0*(-1.d0)* vlz_b(j ,k ,l ,lv) &
                  + t0*x1*y1*( 1.d0)* vlz  (j1,k1,l1,lv) &
                  + t0*x0*y1*( 1.d0)* vlz  (j ,k1,l1,lv) &
                  + t0*x1*y0*( 1.d0)* vlz  (j1,k ,l1,lv) &
                  + t0*x0*y0*( 1.d0)* vlz  (j ,k ,l1,lv) &
                  + t0*x1*y1*(-1.d0)* vlz  (j1,k1,l ,lv) &
                  + t0*x0*y1*(-1.d0)* vlz  (j ,k1,l ,lv) &
                  + t0*x1*y0*(-1.d0)* vlz  (j1,k ,l ,lv) &
                  + t0*x0*y0*(-1.d0)* vlz  (j ,k ,l ,lv) )/(z(l1,lv)-z(l,lv))

              fx = xc - xm - 0.5d0*dtt*(vxm+vxc)*v_uni
              fy = yc - ym - 0.5d0*dtt*(vym+vyc)*v_uni
              fz = zc - zm - 0.5d0*dtt*(vzm+vzc)*v_uni

              fxx = 1.d0 - 0.5d0*dtt*vxdx*v_uni
              fxy =      - 0.5d0*dtt*vxdy*v_uni
              fxz =      - 0.5d0*dtt*vxdz*v_uni

              fyx =      - 0.5d0*dtt*vydx*v_uni
              fyy = 1.d0 - 0.5d0*dtt*vydy*v_uni
              fyz =      - 0.5d0*dtt*vydz*v_uni

              fzx =      - 0.5d0*dtt*vzdx*v_uni
              fzy =      - 0.5d0*dtt*vzdy*v_uni
              fzz = 1.d0 - 0.5d0*dtt*vzdz*v_uni

              det = fxx*(fyy*fzz-fyz*fzy) - fxy*(fyz*fzx-fyx*fzz) + fxz*(fyx*fzy-fyy*fzx)

              fixx =  fyy*fzz - fyz*fzy
              fixy =-(fxy*fzz - fxz*fzy)
              fixz =  fxy*fyz - fxz*fyy
              fiyx =-(fyx*fzz - fyz*fzx)
              fiyy =  fxx*fzz - fxz*fzx
              fiyz =-(fxx*fyz - fxz*fyx)
              fizx =  fyx*fzy - fyy*fzx
              fizy =-(fxx*fzy - fxy*fzx)
              fizz =  fxx*fyy - fxy*fyx

              dxc = -(fixx*fx +fixy*fy +fixz*fz)/det
              dyc = -(fiyx*fx +fiyy*fy +fiyz*fz)/det
              dzc = -(fizx*fx +fizy*fy +fizz*fz)/det
              
              ! if(ip==ip_test.or.abs(dxc/xc)>1d0.or.abs(dyc/yc)>1d0.or.abs(dzc/zc)>1d0)then
              !    write(6,'(2i4,99es12.4)') ip, itr, x1,y1,z1,t1,dxc,dyc,dzc, xc,yc,zc
              ! endif
              ! write(1000+ip,'(2i4,99es12.4)') ip, itr, x1,y1,z1,t1,dxc,dyc,dzc, xc,yc,zc

              if(abs(dxc) < 1.d-8*abs(xc) .and. &
                 abs(dyc) < 1.d-8*abs(yc) .and. &
                 abs(dzc) < 1.d-8*abs(zc)       ) goto 10
              
              xc = xc + dxc
              yc = yc + dyc
              zc = zc + dzc

              xc = max(min(xc,x(ju,lv_min)),x(jd,lv_min))
              yc = max(min(yc,y(ku,lv_min)),y(kd,lv_min))
#ifdef FULL
              zc = max(min(zc,z(lu,lv_min)),z(ld,lv_min))
#else
              zc = max(min(zc,z(lu,lv_min)),0.d0)
#endif
              
           enddo
10         continue
           
           ! condition to cross the FMR boundary
           dxdt = (xc-xtmp)/dtt
           dydt = (yc-ytmp)/dtt
           dzdt = (zc-ztmp)/dtt

           if    (dxdt>0d0)then
              dt_xb = (x(ju,lv_p)-xtmp)/dxdt
           elseif(dxdt<0d0)then
              dt_xb = (x(jd,lv_p)-xtmp)/dxdt
           else
              dt_xb = 1d99
           endif

           if    (dydt>0d0)then
              dt_yb = (y(ku,lv_p)-ytmp)/dydt
           elseif(dydt<0d0)then
              dt_yb = (y(kd,lv_p)-ytmp)/dydt
           else
              dt_yb = 1d99
           endif

           if    (dzdt>0d0)then
              dt_zb = (z(lu,lv_p)-ztmp)/dzdt
           elseif(dzdt<0d0)then
              dt_zb = (z(ld,lv_p)-ztmp)/dzdt
           else
              dt_zb = 1d99
           endif

           dtb = min(dt_xb,dt_yb,dt_zb)
           

           if(dtb < dtt)then

              write(6,*) ip,"reached a boundary"
              t_p(ip) = t - dtt + dtb
              x_p(ip) = xtmp + dtb*dxdt
              y_p(ip) = ytmp + dtb*dydt
              z_p(ip) = ztmp + dtb*dzdt
              remain(ip) = 1
              
              goto 11
           endif
           
           ! update current position
           xtmp = xc
           ytmp = yc
           ztmp = zc
           

!!! interpolate variables
           call coorindex3D(xtmp,ytmp,ztmp,j1,k1,l1,lv)
           j=j1-1
           k=k1-1
           l=l1-1

           x1 = (xtmp-x(j,lv))/(x(j1,lv)-x(j,lv))
           x0 = 1.d0-x1
           y1 = (ytmp-y(k,lv))/(y(k1,lv)-y(k,lv))
           y0 = 1.d0-y1
           z1 = (ztmp-z(l,lv))/(z(l1,lv)-z(l,lv))
           z0 = 1.d0-z1

           vxtmp = x1*y1*z1* vlx  (j1,k1,l1,lv) &
                 + x0*y1*z1* vlx  (j ,k1,l1,lv) &
                 + x1*y0*z1* vlx  (j1,k ,l1,lv) &
                 + x0*y0*z1* vlx  (j ,k ,l1,lv) &
                 + x1*y1*z0* vlx  (j1,k1,l ,lv) &
                 + x0*y1*z0* vlx  (j ,k1,l ,lv) &
                 + x1*y0*z0* vlx  (j1,k ,l ,lv) &
                 + x0*y0*z0* vlx  (j ,k ,l ,lv)

           vytmp = x1*y1*z1* vly  (j1,k1,l1,lv) &
                 + x0*y1*z1* vly  (j ,k1,l1,lv) &
                 + x1*y0*z1* vly  (j1,k ,l1,lv) &
                 + x0*y0*z1* vly  (j ,k ,l1,lv) &
                 + x1*y1*z0* vly  (j1,k1,l ,lv) &
                 + x0*y1*z0* vly  (j ,k1,l ,lv) &
                 + x1*y0*z0* vly  (j1,k ,l ,lv) &
                 + x0*y0*z0* vly  (j ,k ,l ,lv)

           vztmp = x1*y1*z1* vlz  (j1,k1,l1,lv) &
                 + x0*y1*z1* vlz  (j ,k1,l1,lv) &
                 + x1*y0*z1* vlz  (j1,k ,l1,lv) &
                 + x0*y0*z1* vlz  (j ,k ,l1,lv) &
                 + x1*y1*z0* vlz  (j1,k1,l ,lv) &
                 + x0*y1*z0* vlz  (j ,k1,l ,lv) &
                 + x1*y0*z0* vlz  (j1,k ,l ,lv) &
                 + x0*y0*z0* vlz  (j ,k ,l ,lv)

        enddo
        
        t_p(ip) = t
        x_p(ip) = xtmp
        y_p(ip) = ytmp
        z_p(ip) = ztmp
        evolution_finished(ip) = .true.
     endif

11   continue
     
  enddo ! particle evolution end
  !$omp end do
  !$omp end parallel

  np_evolved = 0
  do ip=1,ipu
     if(evolution_finished(ip))then
        np_evolved = np_evolved + 1
     endif
  enddo
  
end subroutine evolution_particle_3D_levels
