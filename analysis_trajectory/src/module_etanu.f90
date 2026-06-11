module module_etanu
  implicit none
  
  private
  public :: etanu_init, etanu_eta_ene
  
  real(8) :: dlogf2,logf2_min,logf2_max
  integer :: nf2

  real(8),allocatable :: logf2(:),eta_from_f2(:),f3_from_f2(:)!,diff_fermi2(:)
  
  real(8),parameter :: pi = 4d0*atan(1d0) !3.14159265358979323846264338327950288d0
  real(8),parameter :: amu_cgs = 1.66053873d-24
  real(8),parameter :: na = 1d0/amu_cgs
  real(8),parameter :: mev_to_erg = 1.60217733d-6
  real(8),parameter :: kb_erg = 1.380658d-16
      
  real(8),parameter :: h       = 6.6260689633d-27, &
                       clight  = 2.99792458d10, &
                       hbar    = 0.5d0 * h/pi, &
                       hbc     = hbar*clight, &
                       hbc_mevfm = hbar*clight/mev_to_erg*1d39, &
                       hbc_mevcm = hbar*clight/mev_to_erg

contains
  subroutine etanu_init(fn)
    character(*),intent(in) :: fn

    integer :: access, iac

    iac = access(fn," ")
    if(iac==0)then
       call etanu_readtable(fn)
    else
       call etanu_maketable(fn)
    endif
    
  end subroutine etanu_init

  subroutine etanu_readtable(fn)
    character(*),intent(in) :: fn
    integer :: n_unit
    integer :: if2

    open(newunit=n_unit,file=fn,status="old")
    read(n_unit,*) nf2!,logf2_min,logf2_max
    allocate(logf2(nf2),eta_from_f2(nf2),f3_from_f2(nf2))!,diff_fermi2(neta))
    do if2=1,nf2
       read(n_unit,*) logf2(if2),eta_from_f2(if2),f3_from_f2(if2)
    enddo
    close(n_unit)
    
    logf2_min = logf2(1)
    logf2_max = logf2(nf2)
    dlogf2 = (logf2_max-logf2_min)/dble(nf2-1)
    ! write(6,*) logf2_min, logf2_max, dlogf2
  end subroutine etanu_readtable
  
  subroutine etanu_maketable(fn)
    character(*),intent(in) :: fn
    integer :: n_unit

    integer :: i
    real(8) :: dx
    real(8) :: h1,h2,h3,h4,f1,f2,f3,f4,fint,x,f_ref,fint_red

    integer :: ieta,if2

    real(8) :: deta,eta_min,eta_max
    integer :: neta
    real(8),allocatable :: eta(:),fermi2(:),fermi3(:)

    real(8) :: xx, xxp


    nf2=15001
    logf2_min=log10(1d-40)
    logf2_max=log10(1d3)
    dlogf2=(logf2_max-logf2_min)/dble(nf2-1)
    allocate(logf2(nf2),eta_from_f2(nf2),f3_from_f2(nf2))!,diff_fermi2(neta))
    
    neta = 20001
    eta_min = -1d2
    eta_max = 1d2
    deta = (eta_max-eta_min)/dble(neta-1)
    allocate(eta(neta),fermi2(neta),fermi3(neta))!,diff_fermi2(neta))

    write(6,*) "eta_min,eta_max : ", eta_min, eta_max
    write(6,*) "Neta, deta : ", neta, deta

    dx=1d-2
    do i=1,neta
       eta(i) = eta_min + deta*dble(i-1)
       
       fint = 0.d0
       x = 0.d0
       fint_red=0.d0
       if(eta(i)>3d2)then
          fint = eta(i)**3/3.d0 + 4.d0*eta(i) + 2.d0*exp(min(3d2,max(-3d2,-eta(i))))
       elseif(eta(i)<-3d2)then
          fint = 2.d0*exp(min(3d2,max(-3d2,eta(i))))
       else
          do while(x<100d0.or.x-eta(i)<100d0)
             f1 = dx* (x         )**2/(exp((x         )-eta(i))+1d0)
             f2 = dx* (x+0.5d0*dx)**2/(exp((x+0.5d0*dx)-eta(i))+1d0)
             f3 = dx* (x+0.5d0*dx)**2/(exp((x+0.5d0*dx)-eta(i))+1d0)
             f4 = dx* (x+      dx)**2/(exp((x+      dx)-eta(i))+1d0)
             x = x + dx
             fint = fint + (f1 +2d0*f2 +2d0*f3 +f4)/6.d0
          enddo
          fint_red = exp(eta(i)-x)*(x**2 + 2.d0*x + 2.d0)
       endif
       
       fermi2(i) = fint


       fint = 0.d0
       x = 0.d0
       fint_red=0.d0
       if(eta(i)>3d2)then
          fint = eta(i)**4/4.d0 + pi**2/2d0*eta(i)**2 + 12d0 - 6.d0*exp(min(3d2,max(-3d2,-eta(i))))
       elseif(eta(i)<-3d2)then
          fint = 6.d0*exp(min(3d2,max(-3d2,eta(i))))
       else
          do while(x<100d0.or.x-eta(i)<100d0)
             f1 = dx* (x         )**3/(exp((x         )-eta(i))+1d0)
             f2 = dx* (x+0.5d0*dx)**3/(exp((x+0.5d0*dx)-eta(i))+1d0)
             f3 = dx* (x+0.5d0*dx)**3/(exp((x+0.5d0*dx)-eta(i))+1d0)
             f4 = dx* (x+      dx)**3/(exp((x+      dx)-eta(i))+1d0)
             x = x + dx
             fint = fint + (f1 +2d0*f2 +2d0*f3 +f4)/6.d0
          enddo
          fint_red = exp(eta(i)-x)*(x**3 + 3.d0*x**2 + 6.d0*x + 6.d0)
       endif
       
       fermi3(i) = fint

       ! if(eta(i)>0d0)then
       !    f_ref = eta(i)**3/3.d0 + 4.d0*eta(i)+2.d0*exp(min(2d2,max(-2d2,-eta(i))))
       ! elseif(eta(i)<0d0)then
       !    f_ref = 2.d0*exp(min(2d2,max(-2d2,eta(i))))
       ! elseif(eta(i)==0d0)then
       !    f_ref = 1.803085354739391428099607d0
       ! endif
       ! diff_fermi2(i) = pi**2/3.d0*eta(i) + 1.d0/3.d0*eta(i)**3
       ! write(6,'("eta, F2, F3 = ",99es12.4)') eta(i),fermi2(i),fermi3(i)
    enddo
    !stop

    fermi2(:) = log10(fermi2(:))
    fermi3(:) = log10(fermi3(:))

    write(6,'("eta, F2 = ",99es12.4)') eta(1),fermi2(1),fermi2(1),logf2_min
    write(6,'("eta, F2 = ",99es12.4)') eta(neta),fermi2(neta),fermi2(nf2),logf2_max

    if(logf2_min<fermi2(1).or.fermi2(neta)<logf2_max)then
       write(6,*) "range of F2 is not enough!"
       stop
    endif


    do if2=1,nf2
       logf2(if2) = logf2_min + dlogf2*dble(if2-1)
    enddo

    ieta=1
    do if2=1,nf2
       do while(logf2(if2)>fermi2(ieta+1))
          ieta=ieta+1
          if(ieta>neta)stop "ieta exceeded limit"
       enddo
       ! write(6,*) fermi2(ieta),logf2(if2),fermi2(ieta+1)
       xxp= (logf2(if2)-fermi2(ieta))/(fermi2(ieta+1)-fermi2(ieta))
       xx = 1.d0-xxp

       eta_from_f2(if2) = xx*eta(ieta) + xxp*eta(ieta+1)
       f3_from_f2 (if2) = xx*fermi3(ieta) + xxp*fermi3(ieta+1)
       
       ! write(6,'(99es15.7)') logf2(if2),eta_from_f2(if2),f3_from_f2(if2), eta(ieta),eta(ieta+1)
    enddo

    open(newunit=n_unit,file=fn,status="new")
    write(n_unit,*) nf2,logf2_min,logf2_max
    do if2=1,nf2
       write(n_unit,'(99es15.7)') logf2(if2),eta_from_f2(if2),f3_from_f2(if2)
    enddo
    close(n_unit)
    
    ! do ieta=1,neta
    !    write(13,'(99es15.7)') fermi2(ieta), eta(ieta), fermi3(ieta)
    ! enddo
    
    deallocate(eta,fermi2,fermi3)!,diff_fermi2(neta))
    
  end subroutine etanu_maketable
  
  subroutine etanu_eta_f3(f2,eta,f3)
    real(8),intent(in) :: f2
    real(8),intent(out) :: eta,f3
    
    integer :: if2
    real(8) :: xx, xxp, f2_tmp

    !write(6,*) f2, logf2_max
    f2_tmp = max(f2,1d1**logf2_min)
    if2 = int((log10(f2_tmp)-logf2_min)/dlogf2)+1
    xxp= (log10(f2_tmp)-logf2(if2))/dlogf2
    xx = 1d0-xxp

    eta = xx*eta_from_f2(if2) + xxp*eta_from_f2(if2+1)
    f3  = xx*f3_from_f2(if2) + xxp*f3_from_f2(if2+1)
    
    f3 = 1d1**f3
  end subroutine etanu_eta_f3

  subroutine etanu_eta_ene(rhoynu,tem,eta,ene)
    real(8),intent(in) :: rhoynu,tem
    real(8),intent(out) :: eta,ene

    !real(8) :: const = na*hbc_mevfm**3*2d0*pi
    real(8) :: f2,f3
    integer :: if2
    real(8) :: xx, xxp

    f2 = rhoynu/tem**3*na*hbc_mevcm**3*2d0*pi

    !write(6,*) f2, rhoynu, tem

    call etanu_eta_f3(f2,eta,f3)

    ene = f3/(2d0*pi*pi)*(tem/hbc_mevcm)**3 *tem*mev_to_erg
    
  end subroutine etanu_eta_ene

end module module_etanu
