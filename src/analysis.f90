subroutine output_profile(dir_out,time)
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

  write(6,*) "enter analysis routine"

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

  write(6,*) "analysis finished"
  
end subroutine output_profile



subroutine analysis_3d_data(time,unum,rin,rfl)
  use simdata3D
  use module_eos, only: hhh_min
  use ejecta_mod
  implicit none
  real(8),intent(in) :: time
  integer,intent(in) :: unum
  real(8),intent(in) :: rin,rfl
  real(8) :: mass_ejecta_geo, mass_ejecta_bernoulli, mass_ejecta_hut1
  real(8) :: mass_total

  real(8) :: hhh_crit
  integer :: j,k,l,lv
  real(8) :: dm

  ! logical :: condition_ejecta_geo, condition_ejecta_bernoulli, condition_ejecta_hut1
  
  hhh_crit = hhh_min
  
  mass_ejecta_geo=0d0
  mass_ejecta_bernoulli=0d0
  mass_ejecta_hut1=0d0
  mass_total = 0d0

  do lv=lv_min,lv_max
     
     !$omp parallel default(none) &
     !$omp shared(lv,ld,lu,kd,ku,jd,ju,rfl,rin,hhh_crit,x,y,z,vol3D,qb) &
     !$omp private(dm) &
     !$omp reduction(+:mass_total,mass_ejecta_geo, mass_ejecta_bernoulli, mass_ejecta_hut1)
     !$omp do
     do l=ld,lu
        do k=kd,ku
           do j=jd,ju
              dm = vol3D(j,k,l,lv)*qb(j,k,l,lv)
              if(        sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
                   .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>=rin)then
                 mass_total = mass_total + dm
              endif

              if(condition_ejecta_geo(j,k,l,lv,rfl,rin))then
                 mass_ejecta_geo = mass_ejecta_geo + dm
              endif

              if(condition_ejecta_bernoulli(j,k,l,lv,rfl,rin,hhh_crit))then
                 mass_ejecta_bernoulli = mass_ejecta_bernoulli + dm
              endif

              if(condition_ejecta_hut1(j,k,l,lv,rfl,rin,hhh_crit))then
                 mass_ejecta_hut1 = mass_ejecta_hut1 + dm
              endif

           enddo
        enddo
     enddo
     !$omp end do
     !$omp end parallel
  enddo
  
  write(unum,'(" ",99es15.7)') time, mass_total, mass_ejecta_geo, mass_ejecta_bernoulli, mass_ejecta_hut1
  flush(unum)
  
end subroutine analysis_3d_data
