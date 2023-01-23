subroutine set_particle_uniform(ib,ittot,n_r,n_phi,n_theta,r_max,npv)
  implicit none
  integer,intent(in) :: ib,ittot,n_r,n_phi,n_theta
  real(8),intent(in) :: r_max
  integer,intent(out) :: npv
  
  integer :: i_r, i_phi, i_theta, ip
  real(8) :: xx,yy,zz,r,phi,theta,rin,rout

  rin=5d7
  rout=5d8

  if(ib==0)then
     npv=n_r*n_phi*n_theta
  endif

  if(ib==1)then
     ip=1
     do i_r=1,n_r
        do i_theta=1,n_theta
           do i_phi=1,n_phi
              r = rin + (rout-rin)*dble(i_r-1)/dble(n_r-1)
              theta = pi/2d0*dble(i_theta)/dble(n_theta)
              phi = 2d0*pi*dble(i_phi)/dble(n_phi)

              xx=r*sin(theta)*cos(phi)
              yy=r*sin(theta)*sin(phi)
              zz=r*cos(theta)
              
              write(6,*) ip,xx,yy,zz
              
              ip=ip+1
           enddo
        enddo
     enddo

  endif
  
end subroutine set_particle_uniform
