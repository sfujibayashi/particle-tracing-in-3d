subroutine set_ejecta_uniform(rfl,rin,mass_crit,mass_min,npv)
  
  use simdata3D
  use particle_data
  use module_eos
  use divide
  use io
  use const,only:pi

  implicit none

  real(8),intent(in) :: rfl,rin,mass_crit,mass_min
  integer,intent(out) :: npv
  
  integer :: n_r, n_theta, n_phi
  integer :: i_r, i_theta, i_phi
  real(8) :: ratio, dr, dtheta, dphi

  real(8),allocatable :: r(:), theta(:), phi(:)

  integer :: ip
  
  !n_r = 15; n_theta = 8
  !n_r = 28; n_theta = 16
  n_r = 41; n_theta = 16
  n_phi = n_theta*4
  
  npv = n_r * n_theta * n_phi

  dtheta = 0.5d0*pi/dble(n_theta)
  dphi = 2d0*pi/dble(n_phi)
  dr = dtheta*rin

  !dr = (rfl-rin)/dble(n_r-1)
  call find_ratio(rin,rfl,dr,n_r-1,ratio)

  allocate(r(n_r), theta(n_theta), phi(n_phi))
  
  do i_r = 1,n_r
     r(i_r) = rin + (ratio**(i_r-1) - 1d0) / (ratio - 1d0)*dr
     !r(i_r) = rin + dr*dble(i_r-1)
     
     !write(6,*) i_r, r(i_r)
  enddo

  do i_theta = 1,n_theta
     theta(i_theta) = dtheta * (dble(i_theta-1)+0.5d0)

     !write(6,*) i_theta, theta(i_theta)
  enddo


  do i_phi = 1,n_phi
     phi(i_phi) = dphi * (dble(i_phi-1)+0.5d0)
     
     !write(6,*) i_phi, phi(i_phi)
  enddo
  

  call allocate_particle_data(npv)

  open(101,file=trim(dir_out)//'/setting.dat',status='replace')
  !open(101,file="setting.dat",status='replace')
  write(101,'(99i5)') n_r, n_theta, n_phi
  ip = 1
  do i_phi = 1,n_phi
     do i_r = 1,n_r
        do i_theta = 1,n_theta
           
           var_p(index_x,ip) = r(i_r)*sin(theta(i_theta))*cos(phi(i_phi))
           var_p(index_y,ip) = r(i_r)*sin(theta(i_theta))*sin(phi(i_phi))
           var_p(index_z,ip) = r(i_r)*cos(theta(i_theta))

           write(101,'(4i5,99es15.7)') ip, i_r, i_theta, i_phi, &
                r(i_r), theta(i_theta), phi(i_phi), &
                var_p(index_x,ip),var_p(index_y,ip),var_p(index_z,ip)
           

           ip = ip + 1
        enddo
     enddo
  enddo

  close(101)


  block
    
    real(8) :: xx,yy,zz, x1,y1,z1,x0,y0,z0
    integer :: j1,k1,l1,lv0, j0,k0,l0
    real(8) :: hhh_i, ut_i, ye_i, vx_i, vy_i, vz_i, floss, gam_inf_r
    real(8) :: hhh_r

    hhh_r = 1d0
    open(101,file=trim(dir_out)//'/p_inside.dat',status='replace')
    write(101,*) "#", npv
    write(101,'(a1,99a14)') "#", &
         "ip","m","ut-1","hut","x","y","z","vx","vy","vz","ye","Gamma_inf(r)"
    do ip=1,npv
       
       flag_evol(ip) = 1
       xx=var_p(index_x,ip)
       yy=var_p(index_y,ip)
       zz=var_p(index_z,ip)
  
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
  
       hhh_i = x1*y1*z1* hhh   (j1,k1,l1,lv0) &
             + x0*y1*z1* hhh   (j0,k1,l1,lv0) &
             + x1*y0*z1* hhh   (j1,k0,l1,lv0) &
             + x0*y0*z1* hhh   (j0,k0,l1,lv0) &
             + x1*y1*z0* hhh   (j1,k1,l0,lv0) &
             + x0*y1*z0* hhh   (j0,k1,l0,lv0) &
             + x1*y0*z0* hhh   (j1,k0,l0,lv0) &
             + x0*y0*z0* hhh   (j0,k0,l0,lv0)
       ut_i  = x1*y1*z1* ut    (j1,k1,l1,lv0) &
             + x0*y1*z1* ut    (j0,k1,l1,lv0) &
             + x1*y0*z1* ut    (j1,k0,l1,lv0) &
             + x0*y0*z1* ut    (j0,k0,l1,lv0) &
             + x1*y1*z0* ut    (j1,k1,l0,lv0) &
             + x0*y1*z0* ut    (j0,k1,l0,lv0) &
             + x1*y0*z0* ut    (j1,k0,l0,lv0) &
             + x0*y0*z0* ut    (j0,k0,l0,lv0)
       ye_i  = x1*y1*z1* ye    (j1,k1,l1,lv0) &
             + x0*y1*z1* ye    (j0,k1,l1,lv0) &
             + x1*y0*z1* ye    (j1,k0,l1,lv0) &
             + x0*y0*z1* ye    (j0,k0,l1,lv0) &
             + x1*y1*z0* ye    (j1,k1,l0,lv0) &
             + x0*y1*z0* ye    (j0,k1,l0,lv0) &
             + x1*y0*z0* ye    (j1,k0,l0,lv0) &
             + x0*y0*z0* ye    (j0,k0,l0,lv0)

       vx_i  = x1*y1*z1* vlx    (j1,k1,l1,lv0) &
             + x0*y1*z1* vlx    (j0,k1,l1,lv0) &
             + x1*y0*z1* vlx    (j1,k0,l1,lv0) &
             + x0*y0*z1* vlx    (j0,k0,l1,lv0) &
             + x1*y1*z0* vlx    (j1,k1,l0,lv0) &
             + x0*y1*z0* vlx    (j0,k1,l0,lv0) &
             + x1*y0*z0* vlx    (j1,k0,l0,lv0) &
             + x0*y0*z0* vlx    (j0,k0,l0,lv0)
       vy_i  = x1*y1*z1* vly    (j1,k1,l1,lv0) &
             + x0*y1*z1* vly    (j0,k1,l1,lv0) &
             + x1*y0*z1* vly    (j1,k0,l1,lv0) &
             + x0*y0*z1* vly    (j0,k0,l1,lv0) &
             + x1*y1*z0* vly    (j1,k1,l0,lv0) &
             + x0*y1*z0* vly    (j0,k1,l0,lv0) &
             + x1*y0*z0* vly    (j1,k0,l0,lv0) &
             + x0*y0*z0* vly    (j0,k0,l0,lv0)
       vz_i  = x1*y1*z1* vlz    (j1,k1,l1,lv0) &
             + x0*y1*z1* vlz    (j0,k1,l1,lv0) &
             + x1*y0*z1* vlz    (j1,k0,l1,lv0) &
             + x0*y0*z1* vlz    (j0,k0,l1,lv0) &
             + x1*y1*z0* vlz    (j1,k1,l0,lv0) &
             + x0*y1*z0* vlz    (j0,k1,l0,lv0) &
             + x1*y0*z0* vlz    (j1,k0,l0,lv0) &
             + x0*y0*z0* vlz    (j0,k0,l0,lv0)

       floss = max(0d0,0.0032d0 - 0.0085d0*ye_i)
       
       gam_inf_r = - hhh_i*ut_i/hhh_r * (1d0-floss)
       
       var_p(index_ut1,ip) = ut_i + 1.d0
       var_p(index_hut,ip) = ut_i*hhh_i! + hhh_at
       
       write(101,'(a1,i14,99es14.6)') " ",ip, var_p(index_dm,ip), var_p(index_ut1,ip), var_p(index_hut,ip), var_p(index_x,ip), var_p(index_y,ip), var_p(index_z,ip), vx_i, vy_i, vz_i, ye_i, gam_inf_r
    enddo
    
  end block
  
end subroutine set_ejecta_uniform

subroutine find_ratio(r_in,r_out,dr,n_r,ratio)
  implicit none
  
  integer,intent(in) :: n_r
  real(8),intent(in) :: r_in, r_out, dr
  real(8),intent(out) :: ratio

  integer,parameter :: itrlim = 30
  integer :: itr
  real(8),parameter :: ratio_min = 1d0 + 1d-5
  real(8) :: dratio, f, df

  ratio = 1.01d0

  !write(6,*) r_in,r_out,dr,n_r
  
  do itr=1,itrlim
     f  = r_out - r_in - dr*(ratio**n_r-1d0)/(ratio-1d0)
     df = -dr*( dble(n_r)*ratio**(n_r-1)/(ratio-1d0) - (ratio**n_r-1d0)/(ratio-1d0)**2 )

     dratio = -f/df
     !write(6,'(i5,99es12.4)') itr, ratio, dratio, f, df, dr*(ratio**n_r-1d0)/(ratio-1d0)
     dratio = max(min(dratio,1d-2),-1d-2)
     
     ratio = ratio + dratio
     ratio = max(ratio, ratio_min)

     if(abs(dratio/ratio) < 1d-10)exit
     
  enddo
end subroutine find_ratio
