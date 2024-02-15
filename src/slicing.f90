subroutine slicing_phi(nphi,lv)
  use const,only : pi
  use simdata3D
  implicit none
  integer,intent(in) :: nphi

  real(8) :: dphi,phi
  integer :: iphi
  
  real(8),allocatable :: x_interp(:), y_interp(:), z_interp(:)
  real(8),allocatable :: sen_interp(:)
  real(8) :: nr,nz

  nr=100
  nz=100
  allocate(x_interp(nr), y_interp(nr), z_interp(nz))
  allocate(sen_interp(nr,nz))
  rmin=-min(abs(x(ju,lv)),abs(x(jd,lv)))
  rmax=+min(abs(x(ju,lv)),abs(x(jd,lv)))
  
  dphi = pi/dble(nphi)
  do iphi=1,nphi
     phi = dphi*dble(iphi-1)

     do ir=1,nr
        x_interp(ir) = (rmin + (rmax-rmin)*dble(ir-1)/dble(nr-1))*cos(phi)
        y_interp(ir) = (rmin + (rmax-rmin)*dble(ir-1)/dble(nr-1))*sin(phi)
     enddo
     do iz=1,nz
        z_interp(iz) = (rmin + (rmax-rmin)*dble(iz-1)/dble(nz-1))
     enddo
     
     do iz=1,nz
        do ir=1,nr
           
        enddo
     enddo
     
  enddo
  
end subroutine slicing_phi
