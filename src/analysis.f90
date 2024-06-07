subroutine analysis(dir_out,time)
  implicit none
  character(*),intent(in) :: dir_out
  real(8),intent(in) :: time

  write(6,*) "enter analysis routine"

  ! call angle_dependence(dir_out,time)
  call histograms(dir_out,time)

end subroutine analysis

subroutine histograms(dir_out,time)
  use simdata3D
  use module_eos, only: hhh_min
  implicit none
  character(*),intent(in) :: dir_out
  real(8),intent(in) :: time
  
  integer,parameter :: nrho = 131
  real(8) :: hist_v_lrho(nrho), hist_lrho(nrho)
  real(8),parameter :: lrho_max = log10(1d15), lrho_min = log10(1d2)
  real(8),parameter :: dlrho=(lrho_max-lrho_min)/dble(nrho)

  integer,parameter :: nbgam = 91
  real(8) :: hist_v_lbgam(nbgam), hist_lbgam(nbgam), hist_lbgam_hut(nbgam)
  real(8),parameter :: lbgam_max = log10(1d1), lbgam_min = log10(1d-2)
  real(8),parameter :: dlbgam=(lbgam_max-lbgam_min)/dble(nbgam)

  integer :: irho, ibgam
  integer :: j,k,l,lv
  real(8) :: dm, mass_tot, bgam, gam, bgam_hut, gam_hut
  integer :: ncell_bgam1

  integer :: nunit
  character(256) :: fn

  do irho=1,nrho
     hist_v_lrho(irho) = lrho_min + dlrho*dble(irho-1)
  enddo

  do ibgam=1,nbgam
     hist_v_lbgam(ibgam) = lbgam_min + dlbgam*dble(ibgam-1)
  enddo

  hist_lrho(:) = 0d0
  hist_lbgam(:) = 0d0
  hist_lbgam_hut(:) = 0d0
  ncell_bgam1 = 0
  
  do lv=lv_min,lv_max
     
     !$omp parallel default(none) &
     !$omp shared(ld,lu,kd,ku,jd,ju,lv,vol3d,qb,qrho,ut,hhh,hhh_min) &
     !$omp private(dm, irho, ibgam, gam, bgam, gam_hut, bgam_hut) &
     !$omp reduction(+: hist_lrho, hist_lbgam, mass_tot, ncell_bgam1, hist_lbgam_hut)
     !$omp do
     do l=ld,lu
        do k=kd,ku
           do j=jd,ju
              
              dm = qb(j,k,l,lv)*vol3D(j,k,l,lv)
              mass_tot = mass_tot + dm
              
              irho = max(1,min(nrho,int( (log10(qrho(j,k,l,lv))-lrho_min)/dlrho ) + 1 ))
              
              hist_lrho(irho) = hist_lrho(irho) + dm
              
              gam = -ut(j,k,l,lv)
              if(gam>1d0)then
                 bgam = sqrt(gam**2-1d0)
                 
                 ibgam = max(1,min(nbgam,int( (log10(bgam)-lbgam_min)/dlbgam ) + 1 ))

                 hist_lbgam(ibgam) = hist_lbgam(ibgam) + dm

                 if(bgam>=1d0.and.dm>0d0)then
                    ncell_bgam1 = ncell_bgam1 + 1
                 endif
              endif


              gam_hut = -hhh(j,k,l,lv)*ut(j,k,l,lv)/hhh_min
              if(gam_hut>1d0)then
                 bgam_hut = sqrt(gam_hut**2-1d0)
                 
                 ibgam = max(1,min(nbgam,int( (log10(bgam_hut)-lbgam_min)/dlbgam ) + 1 ))

                 hist_lbgam_hut(ibgam) = hist_lbgam_hut(ibgam) + dm
                 
              endif

           enddo
        enddo
     enddo
     !$omp end do
     !$omp end parallel

  enddo
  
  fn = trim(dir_out) // "/hist_rho.dat"
  open(newunit=nunit,file=fn, status="replace", action="write")
  write(nunit,'("#","Mtot=",99es12.4)') mass_tot
  write(nunit,'("#",99a15)') "rho(bottom)", "mass (g)", "cumulative"
  do irho=1,nrho
     write(nunit,'(" ",99es15.7)') 10d0**hist_v_lrho(irho), hist_lrho(irho), sum(hist_lrho(:irho))
  enddo
  close(nunit)

  fn = trim(dir_out) // "/hist_bgam.dat"
  open(newunit=nunit,file=fn, status="replace", action="write")
  write(nunit,'("#","Mtot=",es12.4, "Ncell(betaGamma>1)=",i8)') mass_tot, ncell_bgam1
  write(nunit,'("#",99a15)') "bGam(bottom)", "mass (g)", "cumulative"
  do ibgam=1,nbgam
     write(nunit,'(" ",99es15.7)') 10d0**hist_v_lbgam(ibgam), hist_lbgam(ibgam), sum(hist_lbgam(ibgam:))
  enddo
  close(nunit)

  fn = trim(dir_out) // "/hist_bgam_hut.dat"
  open(newunit=nunit,file=fn, status="replace", action="write")
  write(nunit,'("#","Mtot=",es12.4, "Ncell(betaGamma>1)=",i8)') mass_tot, ncell_bgam1
  write(nunit,'("#",99a15)') "bGam(bottom)", "mass (g)", "cumulative"
  do ibgam=1,nbgam
     write(nunit,'(" ",99es15.7)') 10d0**hist_v_lbgam(ibgam), hist_lbgam_hut(ibgam), sum(hist_lbgam_hut(ibgam:))
  enddo
  close(nunit)

end subroutine histograms

subroutine angle_dependence(dir_out,time)
  use simdata3D
  implicit none
  character(*),intent(in) :: dir_out
  real(8),intent(in) :: time

  real(8),parameter :: pi = atan(1d0)*4d0

  integer :: j,k,l,lv

  integer :: unum

  integer :: ir,nr
  real(8),allocatable :: r(:)
  real(8) :: r_min, r_max

  integer :: idir,ndir
  real(8),allocatable :: theta(:),phi(:)

  integer :: j1,k1,l1
  real(8) :: x_interp,y_interp,z_interp
  real(8) :: x0,x1,y0,y1,z0,z1

  real(8) :: vx_interp,vy_interp,vz_interp,qrho_interp,ye_interp,tem_interp,sen_interp,hhh_interp,ut_interp
  real(8) :: vr

  character(256) :: fn
  integer :: i


  nr = 100
  r_min = 1d7; r_max = 1d10
  allocate(r(nr))
  do ir=1,nr
     r(ir) = 10d0**(log10(r_min) + (log10(r_max)-log10(r_min))*dble(ir-1)/dble(nr-1))
  enddo
  
  ndir= 17
  allocate(theta(ndir), phi(ndir))
  theta( 1) = pi/2d0; phi( 1) = pi*0d0/4d0
  theta( 2) = pi/2d0; phi( 2) = pi*1d0/4d0
  theta( 3) = pi/2d0; phi( 3) = pi*2d0/4d0
  theta( 4) = pi/2d0; phi( 4) = pi*3d0/4d0
  theta( 5) = pi/2d0; phi( 5) = pi*4d0/4d0
  theta( 6) = pi/2d0; phi( 6) = pi*5d0/4d0
  theta( 7) = pi/2d0; phi( 7) = pi*6d0/4d0
  theta( 8) = pi/2d0; phi( 8) = pi*7d0/4d0
  theta( 9) = pi/4d0; phi( 9) = pi*0d0/4d0
  theta(10) = pi/4d0; phi(10) = pi*1d0/4d0
  theta(11) = pi/4d0; phi(11) = pi*2d0/4d0
  theta(12) = pi/4d0; phi(12) = pi*3d0/4d0
  theta(13) = pi/4d0; phi(13) = pi*4d0/4d0
  theta(14) = pi/4d0; phi(14) = pi*5d0/4d0
  theta(15) = pi/4d0; phi(15) = pi*6d0/4d0
  theta(16) = pi/4d0; phi(16) = pi*7d0/4d0
  theta(17) = 0d0   ; phi(17) = 0d0
  
  do idir=1,ndir
     write(fn,'(a,"/profile_along_dir_",i3.3,".dat")') trim(dir_out), idir
     !write(6,'(a)') trim(fn)

     open(newunit=unum,file=fn,status="replace",action="write")
     
     write(unum,'("# time  = ",99es15.7)') time
     write(unum,'("# theta = ",99es15.7)') theta(idir)
     write(unum,'("# phi   = ",99es15.7)') phi(idir)
     write(unum,'("# n_r   = ",99i15)') nr
     
     write(unum,'("#",99i15)') (i,i=1,6)
     write(unum,'("#",99a15)') "r (cm)", "v^r/c", "rho (g/cm^3)", "Ye", "s (k/nuc)", "u_t+1"

     do ir=1,nr
        x_interp = r(ir)*sin(theta(idir))*cos(phi(idir))
        y_interp = r(ir)*sin(theta(idir))*sin(phi(idir))
        z_interp = r(ir)*cos(theta(idir))
        call coorindex3D(x_interp,y_interp,z_interp,j1,k1,l1,lv)
        j = j1-1
        k = k1-1
        l = l1-1

        x1 = (x_interp-x(j,lv))/(x(j1,lv)-x(j,lv))
        x0 = 1.d0-x1
        y1 = (y_interp-y(k,lv))/(y(k1,lv)-y(k,lv))
        y0 = 1.d0-y1
        z1 = (z_interp-z(l,lv))/(z(l1,lv)-z(l,lv))
        z0 = 1.d0-z1
        
        vx_interp &
             = x1*y1*z1* vlx   (j1,k1,l1,lv) &
             + x0*y1*z1* vlx   (j ,k1,l1,lv) &
             + x1*y0*z1* vlx   (j1,k ,l1,lv) &
             + x0*y0*z1* vlx   (j ,k ,l1,lv) &
             + x1*y1*z0* vlx   (j1,k1,l ,lv) &
             + x0*y1*z0* vlx   (j ,k1,l ,lv) &
             + x1*y0*z0* vlx   (j1,k ,l ,lv) &
             + x0*y0*z0* vlx   (j ,k ,l ,lv)
        vy_interp &
             = x1*y1*z1* vly   (j1,k1,l1,lv) &
             + x0*y1*z1* vly   (j ,k1,l1,lv) &
             + x1*y0*z1* vly   (j1,k ,l1,lv) &
             + x0*y0*z1* vly   (j ,k ,l1,lv) &
             + x1*y1*z0* vly   (j1,k1,l ,lv) &
             + x0*y1*z0* vly   (j ,k1,l ,lv) &
             + x1*y0*z0* vly   (j1,k ,l ,lv) &
             + x0*y0*z0* vly   (j ,k ,l ,lv)
        vz_interp &
             = x1*y1*z1* vlz   (j1,k1,l1,lv) &
             + x0*y1*z1* vlz   (j ,k1,l1,lv) &
             + x1*y0*z1* vlz   (j1,k ,l1,lv) &
             + x0*y0*z1* vlz   (j ,k ,l1,lv) &
             + x1*y1*z0* vlz   (j1,k1,l ,lv) &
             + x0*y1*z0* vlz   (j ,k1,l ,lv) &
             + x1*y0*z0* vlz   (j1,k ,l ,lv) &
             + x0*y0*z0* vlz   (j ,k ,l ,lv)
        qrho_interp&
             = x1*y1*z1* qrho  (j1,k1,l1,lv) &
             + x0*y1*z1* qrho  (j ,k1,l1,lv) &
             + x1*y0*z1* qrho  (j1,k ,l1,lv) &
             + x0*y0*z1* qrho  (j ,k ,l1,lv) &
             + x1*y1*z0* qrho  (j1,k1,l ,lv) &
             + x0*y1*z0* qrho  (j ,k1,l ,lv) &
             + x1*y0*z0* qrho  (j1,k ,l ,lv) &
             + x0*y0*z0* qrho  (j ,k ,l ,lv)
        ye_interp &
             = x1*y1*z1* ye    (j1,k1,l1,lv) &
             + x0*y1*z1* ye    (j ,k1,l1,lv) &
             + x1*y0*z1* ye    (j1,k ,l1,lv) &
             + x0*y0*z1* ye    (j ,k ,l1,lv) &
             + x1*y1*z0* ye    (j1,k1,l ,lv) &
             + x0*y1*z0* ye    (j ,k1,l ,lv) &
             + x1*y0*z0* ye    (j1,k ,l ,lv) &
             + x0*y0*z0* ye    (j ,k ,l ,lv)
        tem_interp &
             = x1*y1*z1* tem   (j1,k1,l1,lv) &
             + x0*y1*z1* tem   (j ,k1,l1,lv) &
             + x1*y0*z1* tem   (j1,k ,l1,lv) &
             + x0*y0*z1* tem   (j ,k ,l1,lv) &
             + x1*y1*z0* tem   (j1,k1,l ,lv) &
             + x0*y1*z0* tem   (j ,k1,l ,lv) &
             + x1*y0*z0* tem   (j1,k ,l ,lv) &
             + x0*y0*z0* tem   (j ,k ,l ,lv)
        sen_interp &
             = x1*y1*z1* sen   (j1,k1,l1,lv) &
             + x0*y1*z1* sen   (j ,k1,l1,lv) &
             + x1*y0*z1* sen   (j1,k ,l1,lv) &
             + x0*y0*z1* sen   (j ,k ,l1,lv) &
             + x1*y1*z0* sen   (j1,k1,l ,lv) &
             + x0*y1*z0* sen   (j ,k1,l ,lv) &
             + x1*y0*z0* sen   (j1,k ,l ,lv) &
             + x0*y0*z0* sen   (j ,k ,l ,lv)
        hhh_interp &
             = x1*y1*z1* hhh   (j1,k1,l1,lv) &
             + x0*y1*z1* hhh   (j ,k1,l1,lv) &
             + x1*y0*z1* hhh   (j1,k ,l1,lv) &
             + x0*y0*z1* hhh   (j ,k ,l1,lv) &
             + x1*y1*z0* hhh   (j1,k1,l ,lv) &
             + x0*y1*z0* hhh   (j ,k1,l ,lv) &
             + x1*y0*z0* hhh   (j1,k ,l ,lv) &
             + x0*y0*z0* hhh   (j ,k ,l ,lv)
        ut_interp &
             = x1*y1*z1* ut   (j1,k1,l1,lv) &
             + x0*y1*z1* ut   (j ,k1,l1,lv) &
             + x1*y0*z1* ut   (j1,k ,l1,lv) &
             + x0*y0*z1* ut   (j ,k ,l1,lv) &
             + x1*y1*z0* ut   (j1,k1,l ,lv) &
             + x0*y1*z0* ut   (j ,k1,l ,lv) &
             + x1*y0*z0* ut   (j1,k ,l ,lv) &
             + x0*y0*z0* ut   (j ,k ,l ,lv)

        vr = (vx_interp*x_interp + vy_interp*y_interp + vz_interp*z_interp)/r(ir)

        write(unum,'(" ",99es15.7)') r(ir), vr, qrho_interp, ye_interp, sen_interp, ut_interp+1d0 
     enddo
     close(unum)
  enddo
  
  
end subroutine angle_dependence
