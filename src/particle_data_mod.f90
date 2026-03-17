module particle_data

  implicit none

  integer,parameter :: nvar = 24, nvar_out = 20

  integer,parameter :: index_x = 1
  integer,parameter :: index_y = 2
  integer,parameter :: index_z = 3
  integer,parameter :: index_vlx = 4
  integer,parameter :: index_vly = 5
  integer,parameter :: index_vlz = 6
  integer,parameter :: index_rho = 7
  integer,parameter :: index_tem = 8
  integer,parameter :: index_ye = 9
  integer,parameter :: index_sen = 10
  integer,parameter :: index_rne = 11
  integer,parameter :: index_rae = 12
  integer,parameter :: index_deptn = 13
  integer,parameter :: index_depta = 14
  integer,parameter :: index_ut = 15
  integer,parameter :: index_hhh = 16
  integer,parameter :: index_ch_nuf = 17
  integer,parameter :: index_ch_naf = 18
  integer,parameter :: index_b2 = 19
  integer,parameter :: index_pres = 20

  integer,parameter :: index_dm = 21
  integer,parameter :: index_ut1 = 22
  integer,parameter :: index_hut = 23
  integer,parameter :: index_qb = 24

  real(8),parameter :: mev_to_kelvin  = 1.160445d10, clight = 2.99792458d10
  character(len=14),parameter :: label_var(nvar_out) = [character(len=14) :: "x [cm]", "y [cm]", "z [cm]", "Vx [cm/s]", "Vy [cm/s]", "Vz [cm/s]", "rho [g/cm^3]", "T [K]", "Ye", "S [k_b/nuc]", "Ee [erg/cm^3]", "Ea [erg/cm^3]", "tau_e", "tau_a", "u_t", "h/c^2", "eta(nu_e)", "eta(nu_ebar)", "b^2 [G^2]", "P [dyn/cm^2]"]

  real(8),parameter :: unit_var(nvar_out) = [1d0, 1d0, 1d0, clight, clight, clight, 1d0, mev_to_kelvin, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0]

  ! Their size is only np
  integer,allocatable :: flag_evol(:)
  real(8),allocatable :: var_p(:,:)
  
contains
  
  subroutine allocate_particle_data(np)
    integer,intent(in) :: np
    integer :: ip

    if(.not.allocated(flag_evol))then
       allocate( &
         flag_evol(np), &
         var_p(nvar,np) )
    endif
    
    ! initialize
    flag_evol(:) = 0
    
  end subroutine allocate_particle_data


  subroutine reallocate_particle_data(np_old,np_new)
    integer,intent(in) :: np_old,np_new
    integer :: ip
    
    integer,allocatable :: flag_evol_buf(:)
    real(8),allocatable :: var_buf(:,:)
    
    allocate( &
         flag_evol_buf(np_old), &
         var_buf   (nvar,np_old) )
    
    do ip=1,np_old
       flag_evol_buf(ip) = flag_evol(ip)
       var_buf(:,ip) = var_p(:,ip)
    enddo

    deallocate( &
         flag_evol, &
         var_p)

    call allocate_particle_data(np_new)

    do ip=1,np_old
       flag_evol(ip) = flag_evol_buf(ip)
       var_p(:,ip) = var_buf(:,ip)
    enddo
    
    deallocate( &
         flag_evol_buf, &
         var_buf)
    
  end subroutine reallocate_particle_data

  subroutine deallocate_particle_data
    
    if(allocated(flag_evol))then
       deallocate( &
         flag_evol, &
         var_p)
    else
       write(6,*) "particle data not allocated yet."
    endif
    
  end subroutine deallocate_particle_data
  
  subroutine output_position(fn,ipu,time)
    integer,intent(in) :: ipu
    character(*),intent(in) :: fn
    real(8),intent(in) :: time
    
    character(10) :: str1
    integer :: unit,ip
    
    open (newunit=unit,file=fn,status="replace",action="write")
    write(unit,'("# ",99es12.4)') time
    do ip = 1,ipu
       if(flag_evol(ip)==1)then
          ! write(unit,'(i10,99es10.2)') ip, x_p(ip),z_p(ip),qrho_p(ip),ye_p(ip),tem_p(ip),sen_p(ip),dm_p(ip)
          write(unit,'(i10,99es25.17)') ip, var_p(index_x,ip), var_p(index_z,ip), var_p(index_rho,ip), var_p(index_ye,ip), var_p(index_tem,ip), var_p(index_sen,ip), var_p(index_dm,ip)
       endif
    enddo
    close(unit)
    
  end subroutine output_position

end module particle_data
