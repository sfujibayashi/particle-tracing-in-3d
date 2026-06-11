module module_weak_interaction

  implicit none
  
  ! mathematical constants
  real(8),parameter :: ln2 = log(2d0)
  real(8),parameter :: pi  = 4d0*atan(1d0)
  
  ! Basic physical constants
  real(8),parameter :: mev_to_erg = 1.60218d-6
  real(8),parameter :: erg_to_mev = 1.d0/mev_to_erg
  
  real(8),parameter :: h       = 6.6260689633d-27
  real(8),parameter :: clight  = 2.99792458d10
  real(8),parameter :: gfmev=1.1663787d-11
  real(8),parameter :: na = 6.02214076d23
  real(8),parameter :: alpha_fine = 7.2973525693d-3
  
  ! some composits
  real(8),parameter :: hbar    = h/2d0/pi
  real(8),parameter :: hbc = hbar*clight
  real(8),parameter :: gferg=gfmev*erg_to_mev**2
  real(8),parameter :: dmnp = 1.29333d0
  real(8),parameter :: memev = 0.51099996d0
  real(8),parameter :: mpmev = 938.272088d0
  real(8),parameter :: mumev = 931.49410242d0
  real(8),parameter :: mpi0mev = 134.9768d0
  real(8),parameter :: mpipmmev = 139.57039d0
  real(8),parameter :: mnmev = mpmev + dmnp
  real(8),parameter :: sigma0 = gfmev**2 * memev**2 * 4.d0/pi *(hbc*erg_to_mev)**2

  real(8),parameter :: kweak = 6146d0
  real(8),parameter :: ga = 1.26d0
  real(8),parameter :: wsin2 = 0.2319d0
  real(8),parameter :: wsin4 = wsin2*wsin2
  real(8),parameter :: c_v = 0.5d0 + 2d0*wsin2
  real(8),parameter :: c_a = 0.5d0
  
  real(8),parameter :: ga2= ga*ga
  real(8),parameter :: cap_nrate_base = ln2/kweak*(1d0+3d0*ga2)/memev**5
  real(8),parameter :: cap_erate_base = ln2/kweak*(1d0+3d0*ga2)/memev**5 * mev_to_erg

  ! pair
  real(8),parameter :: c_pair_e = 2d0*(c_v**2 + c_a**2)
  real(8),parameter :: c_pair_x = 2d0*((c_v-1d0)**2 + (c_a-1d0)**2)
  real(8),parameter :: pair_nrate_base_e = c_pair_e/36.d0/pi**4 *sigma0*clight/(memev*mev_to_erg)**2 /hbc**6 * mev_to_erg**8
  real(8),parameter :: pair_nrate_base_x = c_pair_x/36.d0/pi**4 *sigma0*clight/(memev*mev_to_erg)**2 /hbc**6 * mev_to_erg**8
  real(8),parameter :: pair_erate_base_e = c_pair_e/36.d0/pi**4 *sigma0*clight/(memev*mev_to_erg)**2 /hbc**6 * mev_to_erg**9
  real(8),parameter :: pair_erate_base_x = c_pair_x/36.d0/pi**4 *sigma0*clight/(memev*mev_to_erg)**2 /hbc**6 * mev_to_erg**9

  ! brems
  real(8),parameter :: f_pion = 1d0
  real(8),parameter :: zeta_brems = 0.5d0
  
  real(8),parameter :: brems_erate_base = 1d0/15d0/pi**4.5d0 * mnmev**1.5d0/memev**2 *(f_pion/mpi0mev)**4 *ga2 *zeta_brems *na**2 *sigma0*clight *1.66234d0 * mev_to_erg
  
  ! plasmon
  real(8),parameter :: plasm_nrate_base_e = c_v**2      /(192d0*pi**3*alpha_fine) * sigma0*clight /memev**2 *mev_to_erg**6 /hbc**6
  real(8),parameter :: plasm_nrate_base_x = (c_v-1d0)**2/(192d0*pi**3*alpha_fine) * sigma0*clight /memev**2 *mev_to_erg**6 /hbc**6
  ! some fermi integrals for eta=0
  real(8),parameter :: fermi3_0 = 7d0*pi**4/120d0
  real(8),parameter :: fermi4_0 = 45d0/2d0 * 1.0369277551433699d0
    
  ! 
  real(8),allocatable :: caprate_logeta(:)
  real(8),allocatable :: caprate_logtemp(:)
  
  real(8),allocatable :: caprate_num_e(:,:)
  real(8),allocatable :: caprate_num_p(:,:)
  real(8),allocatable :: caprate_ene_e(:,:)
  real(8),allocatable :: caprate_ene_p(:,:)

  real(8),allocatable :: caprate_neg_num_e(:,:)
  real(8),allocatable :: caprate_neg_num_p(:,:)
  real(8),allocatable :: caprate_neg_ene_e(:,:)
  real(8),allocatable :: caprate_neg_ene_p(:,:)

  integer :: n_eta, n_temp
  real(8) :: caprate_logeta_min,caprate_logtemp_min,caprate_dlogeta,caprate_dlogtemp
  

  
contains
  
  !subroutine capture_rate()
  !end subroutine capture_rate

  subroutine capture_integral_trap(temt,eta_e,eta_nu,eta_na,ih, &
       finte,fintp)
    ! use const
    implicit none
    
    real(8),intent(in) :: temt,eta_e,eta_nu,eta_na
    integer,intent(in) :: ih
    real(8),intent(out) :: finte,fintp
    
    !integer :: krho,krho1,itemt,itemt1,jye,jye1
    !real(8) :: uu,uup,ss,ssp,tt,ttp, ch_n,ch_p,ch_e
    real(8) :: temi,etaem,etapm,zeta_n,zeta_np
    
    real(8) :: etan,etaa,facec,facpc,fintn,finta
    
    temi = 1.d0/temt
    
    etaem = eta_e
    etapm = -etaem
    
    zeta_n = -dmnp*temi
    zeta_np=  dmnp*temi
      
    call fermicap(zeta_n ,etaem,temt,ih,finte)
    call fermicap(zeta_np,etapm,temt,ih,fintp)
    
    etan = eta_nu
    etaa = eta_na
    
    call fermicap(zeta_n ,etan-zeta_n  ,temt,ih,fintn)
    call fermicap(zeta_np,etaa-zeta_np ,temt,ih,finta)
    
    facec = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,etan-zeta_n -etaem)) ))
    facpc = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,etaa-zeta_np-etapm)) ))
    finte = facec * (finte-fintn)
    fintp = facpc * (fintp-finta)
    
    return
  end subroutine capture_integral_trap
  
  ! subroutine epcaprate_ntrap(rhot,temt,fyet,yn,ya, &
  !      finte,fintp)
  !   use eostab
  !   use const
  !   implicit none

  !   real(8),intent(in) :: rhot,temt,fyet,yn,ya
  !   real(8),intent(out) :: finte,fintp
    
  !   integer :: krho,krho1,itemt,itemt1,jye,jye1
  !   real(8) :: uu,uup,ss,ssp,tt,ttp, ch_n,ch_p,ch_e
  !   real(8) :: temi,etae,etaem,etapm,zeta_n,zeta_np
    
  !   integer :: inu,ina
  !   real(8) :: rynu,vv,vvp,ryna,vva,vvap,fynt,fyat
  !   real(8) :: ch_nu,ch_na
  !   real(8) :: etan,etaa,facec,facpc,fintn,finta
  !   real(8) :: ft_ecf,ft_pcf

  !   ft_ecf = ft_f
  !   ft_pcf = ft_f

  !   krho  = max(ked , min(keu-1, int((log10(rhot)-rho_e_min  )*drhoi)+1))
  !   krho1 = krho+1
  !   uu    = max(0.d0, min(1.d0 ,     (log10(rhot)-rho_e(krho))*drhoi)   )
  !   uup   = 1.d0-uu

  !   itemt = max(ied , min(ieu-1, int((log10(temt)-tem_e_min  )*dtemi)+1))
  !   itemt1=itemt+1
  !   ss    = max(0.d0, min(1.d0 ,     (log10(temt)-tem_e(itemt))*dtemi))
  !   ssp   = 1.d0-ss

  !   jye  = max(jed , min(jeu-1, int((fyet-ye_e_min )*dyei )+1))
  !   jye1 = jye +1
  !   tt   = max(0.d0, min(1.d0 ,     (fyet-ye_e(jye))*dyei))
  !   ttp  = 1.d0-tt

  !   ch_n = ssp *ttp *uup * chn_e(itemt ,jye ,krho )   &
  !        + ss  *ttp *uup * chn_e(itemt1,jye ,krho )   &
  !        + ssp *tt  *uup * chn_e(itemt ,jye1,krho )   &
  !        + ssp *ttp *uu  * chn_e(itemt ,jye ,krho1)   &
  !        + ss  *tt  *uup * chn_e(itemt1,jye1,krho )   &
  !        + ss  *ttp *uu  * chn_e(itemt1,jye ,krho1)   &
  !        + ssp *tt  *uu  * chn_e(itemt ,jye1,krho1)   &
  !        + ss  *tt  *uu  * chn_e(itemt1,jye1,krho1)
  !   ch_p = ssp *ttp *uup * chp_e(itemt ,jye ,krho )   &
  !        + ss  *ttp *uup * chp_e(itemt1,jye ,krho )   &
  !        + ssp *tt  *uup * chp_e(itemt ,jye1,krho )   &
  !        + ssp *ttp *uu  * chp_e(itemt ,jye ,krho1)   &
  !        + ss  *tt  *uup * chp_e(itemt1,jye1,krho )   &
  !        + ss  *ttp *uu  * chp_e(itemt1,jye ,krho1)   &
  !        + ssp *tt  *uu  * chp_e(itemt ,jye1,krho1)   &
  !        + ss  *tt  *uu  * chp_e(itemt1,jye1,krho1)
  !   ch_e = ssp *ttp *uup * che_e(itemt ,jye ,krho )   &
  !        + ss  *ttp *uup * che_e(itemt1,jye ,krho )   &
  !        + ssp *tt  *uup * che_e(itemt ,jye1,krho )   &
  !        + ssp *ttp *uu  * che_e(itemt ,jye ,krho1)   &
  !        + ss  *tt  *uup * che_e(itemt1,jye1,krho )   &
  !        + ss  *ttp *uu  * che_e(itemt1,jye ,krho1)   &
  !        + ssp *tt  *uu  * che_e(itemt ,jye1,krho1)   &
  !        + ss  *tt  *uu  * che_e(itemt1,jye1,krho1)
    
  !   temi = 1.d0/temt
  !   etae= ch_e -memev*temi
    
  !   etaem = etae + memev*temi
  !   etapm = -etaem
    
  !   zeta_n = -dmnp*temi
  !   zeta_np=  dmnp*temi
      
  !   call fermicap(zeta_n ,etaem,temt,0,finte)
  !   call fermicap(zeta_np,etapm,temt,0,fintp)
    
  !   ! additional part for neutrino blocking
  !   fynt = max(1.d-99,yn)
  !   rynu = log10(rhot*fynt/rho_uni/(max(tem01,temt)**3))
  !   inu  = max( min( int((rynu-escn0    )*dscni)+1, inut-1), 1)
  !   vv   = max( min(     (rynu-escn(inu))*dscni   , 1.d0  ), 0.d0)
  !   vvp  = 1.d0-vv
  !   ch_nu=vvp*ech_nn(inu) + vv*ech_nn(inu+1)
    
  !   fyat = max(1.d-99,ya)
  !   ryna = log10(rhot*fyat/rho_uni/(max(tem01,temt)**3))
  !   ina  = max( min( int((ryna-escn0)*dscni)+1, inut-1), 1)
  !   vva  = max( min( (ryna-escn(ina))*dscni, 1.d0), 0.d0)
  !   vvap = 1.d0-vva
  !   ch_na= vvap*ech_nn(ina) + vva*ech_nn(ina+1)
  
  !   etan = ch_nu
  !   etaa = ch_na
    
  !   call fermicap(zeta_n ,etan-zeta_n  ,temt,0,fintn)
  !   call fermicap(zeta_np,etaa-zeta_np ,temt,0,finta)
    
  !   facec = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,etan-zeta_n -etaem)) ))
  !   facpc = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,etaa-zeta_np-etapm)) ))
  !   finte = facec * (finte-fintn)
  !   fintp = facpc * (fintp-finta)
    
  !   return
  ! end subroutine epcaprate_ntrap
  
  subroutine fermicap(zeta,eta,tem,ih,fint)
    implicit none

    real(8),intent(in) :: zeta,eta,tem
    integer,intent(in) :: ih
    real(8),intent(out) :: fint
    
    real(8) :: dt,et,xl, xfin,dx
    real(8) :: x, x_1,x_2,x_3,x_4, ex_1,ex_2,ex_3,ex_4, f_1,f_2,f_3,f_4, residual
    
    integer :: i

    ! --- *
    ! A + e^- -> B + nu_e
    ! zeta = (m_A - m_B)/T
    ! eta is the chemical potential of electron (including mass!!) divided by T
    ! ih = 1 => energy, 0 => number
    ! --- *
    
    dt = dmnp/tem
    et = memev/tem

    xl = max(-zeta,et)

    xfin = max(1000.d0,10.d0*(eta-xl))
    !xfin = 10000d0
    dx = 0.01d0

    fint = 0.d0
    do i=0,nint(xfin/dx)
       x = dx*dble(i)
       
       ! x_1 = x;          ex_1 = max(-3d2,min(3d2,x_1+xl-eta))
       ! f_1 = (x_1 + xl + zeta)**(2+ih) * ( (x_1 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_1 + xl)**2) *1.d0/(exp(ex_1)+1d0)
       ! x_2 = x+0.5d0*dx; ex_2 = max(-3d2,min(3d2,x_2+xl-eta))
       ! f_2 = (x_2 + xl + zeta)**(2+ih) * ( (x_2 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_2 + xl)**2) *1.d0/(exp(ex_2)+1d0)
       ! x_3 = x+0.5d0*dx; ex_3 = max(-3d2,min(3d2,x_3+xl-eta))
       ! f_3 = (x_3 + xl + zeta)**(2+ih) * ( (x_3 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_3 + xl)**2) *1.d0/(exp(ex_3)+1d0)
       ! x_4 = x+dx;       ex_4 = max(-3d2,min(3d2,x_4+xl-eta))
       ! f_4 = (x_4 + xl + zeta)**(2+ih) * ( (x_4 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_4 + xl)**2) *1.d0/(exp(ex_4)+1d0)
       
       x_1 = x
       x_2 = x+0.5d0*dx
       x_3 = x+0.5d0*dx
       x_4 = x+dx
       f_1 = capture_rate(zeta,tem,eta,x_1,ih)
       f_2 = capture_rate(zeta,tem,eta,x_2,ih)
       f_3 = f_2!capture_rate(zeta,tem,eta,x_3,ih)
       f_4 = capture_rate(zeta,tem,eta,x_4,ih)
       
       fint = fint + (f_1 + 2d0*f_2 + 2d0*f_3 + f_4)/6d0*dx
       !fint = fint + f_1*dx
       
       ex_4 = max(-3d2,min(3d2,x_4+xl-eta))
       
       residual = exp(-ex_4)* ( (x_4**4 + 4.d0*x_4**3 + 12.d0*x_4**2 + 24.d0*x_4 + 24.d0)  * (5d0*x_4)**ih + 120d0*dble(ih) )
       
       !if(x>10.d0*abs(eta-xl).and.x>10.d0*abs(zeta).and.x>10.d0*et.and.x>10.d0*xl.and.fint > 1.d15*residual) goto 10
       if(x > 10.d0*max(abs(eta-xl),xl) .and. fint > 1.d15*residual) goto 10
    enddo

10  continue

    return
  end subroutine fermicap

  real(8) function capture_rate(zeta,tem,eta,x,ih)
    real(8),intent(in) :: zeta,tem,eta,x
    integer,intent(in) :: ih

    real(8) :: et,xl
    real(8) :: ex

    et = memev/tem
    xl = max(-zeta,et)
    
    ex = max(-3d2,min(3d2,x+xl-eta))
    ! capture_rate = (x + xl + zeta)**(2+ih) * ( (x + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x + xl)**2) *1.d0/(exp(ex)+1d0)

    capture_rate = (x + xl + zeta)**(2+ih) * (x + xl)**2 * sqrt(1.d0 - et**2/(x + xl)**2) *1.d0/(exp(ex)+1d0)
    
    ! test function
    !ex = max(-3d2,min(3d2,x))
    !capture_rate = x**(4+ih)*exp(-ex)

    return
  end function capture_rate

  subroutine fermik(eta,k,fint)
    implicit none

    real(8),intent(in) :: eta
    integer,intent(in) :: k
    real(8),intent(out) :: fint
    
    real(8) :: dt,et,xl, xfin,dx
    real(8) :: x, x_1,x_2,x_3,x_4, ex_1,ex_2,ex_3,ex_4, f_1,f_2,f_3,f_4, residual
    
    integer :: i
    
    xfin = max(1000.d0,10.d0*eta)
    dx = 0.01d0
    fint = 0.d0
    do i=0,nint(xfin/dx)
       x = dx*dble(i)
       
       ! x_1 = x;          ex_1 = max(-3d2,min(3d2,x_1+xl-eta))
       ! f_1 = (x_1 + xl + zeta)**(2+ih) * ( (x_1 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_1 + xl)**2) *1.d0/(exp(ex_1)+1d0)
       ! x_2 = x+0.5d0*dx; ex_2 = max(-3d2,min(3d2,x_2+xl-eta))
       ! f_2 = (x_2 + xl + zeta)**(2+ih) * ( (x_2 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_2 + xl)**2) *1.d0/(exp(ex_2)+1d0)
       ! x_3 = x+0.5d0*dx; ex_3 = max(-3d2,min(3d2,x_3+xl-eta))
       ! f_3 = (x_3 + xl + zeta)**(2+ih) * ( (x_3 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_3 + xl)**2) *1.d0/(exp(ex_3)+1d0)
       ! x_4 = x+dx;       ex_4 = max(-3d2,min(3d2,x_4+xl-eta))
       ! f_4 = (x_4 + xl + zeta)**(2+ih) * ( (x_4 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_4 + xl)**2) *1.d0/(exp(ex_4)+1d0)
       
       x_1 = x
       x_2 = x+0.5d0*dx
       x_3 = x+0.5d0*dx
       x_4 = x+dx
       f_1 = x_1**k/(exp(max(-3d2,min(3d2,x_1-eta)))+1d0)
       f_2 = x_2**k/(exp(max(-3d2,min(3d2,x_2-eta)))+1d0)
       f_3 = f_2
       f_4 = x_4**k/(exp(max(-3d2,min(3d2,x_4-eta)))+1d0)
       
       fint = fint + (f_1 + 2d0*f_2 + 2d0*f_3 + f_4)/6d0*dx
       
       ex_4 = max(-3d2,min(3d2,x_4))
       
       !residual = exp(-ex_4)* ( (x_4**4 + 4.d0*x_4**3 + 12.d0*x_4**2 + 24.d0*x_4 + 24.d0)  * (5d0*x_4)**ih + 120d0*dble(ih) )
       residual = exp(-ex_4)*x_4**k
       
       !if(x>10.d0*abs(eta-xl).and.x>10.d0*abs(zeta).and.x>10.d0*et.and.x>10.d0*xl.and.fint > 1.d15*residual) goto 10
       if(x > 10.d0*abs(eta) .and. fint > 1.d15*residual) goto 10
    enddo

10  continue

    return
  end subroutine fermik
  
  subroutine rate_pairann(eta_e,tem,nemis_n,nemis_o,eemis_n,eemis_o,ene_n)
    real(8),intent(in)  :: eta_e,tem
    real(8),intent(out) :: nemis_n,nemis_o,eemis_n,eemis_o,ene_n

    real(8) :: fmaske,deta_e,fermi4e,fermi3e,fermi4p,fermi3p,ferm43e,ferm43p,teme

    real(8) :: fermi4e_int,fermi3e_int,fermi4p_int,fermi3p_int
    real(8),parameter :: f23 = 2d0/3d0
    real(8),parameter :: ppi = (4d0*atan(1d0))**2
    
    fmaske= (1.d0+sign(1.d0,eta_e))*0.5d0
    deta_e = exp(-abs(eta_e))
    
    fermi4e = ( eta_e**5 *0.2d0          &
               +eta_e**3 *f23 *ppi       &
               +eta_e    *48.d0          &
               +24.d0*deta_e      ) * fmaske &
               +24.d0*deta_e        *(1.d0-fmaske)
 
    fermi3e = ( eta_e**4 *0.25d0         &
               +eta_e**2 *0.5d0 *ppi     &
               +12.d0 -6.d0*deta_e )     &
                               * fmaske &
            + 6.d0 *deta_e *(1.d0-fmaske)
    ferm43e = fermi4e/(fermi3e +1.d-6) *fmaske  &
                          + 4.d0 *(1.d0-fmaske)

    fermi4p = (-eta_e**5 *0.2d0          &
               -eta_e**3 *f23 *ppi       &
               -eta_e    *48.d0          &
               +24.d0*deta_e      ) *(1.d0-fmaske) &
               +24.d0*deta_e        * fmaske
    fermi3p = ( eta_e**4 *0.25d0         &
               +eta_e**2 *0.5d0 *ppi     &
               +12.d0 -6.d0*deta_e )*(1.d0-fmaske) &
               + 6.d0 *deta_e       * fmaske
    ferm43p = fermi4p/(fermi3p +1.d-6) *(1.d0-fmaske) &
                                      + 4.d0 *fmaske

    call fermik( eta_e,3,fermi3e_int)
    call fermik( eta_e,4,fermi4e_int)
    call fermik(-eta_e,3,fermi3p_int)
    call fermik(-eta_e,4,fermi4p_int)
    ! write(6,'(99es15.7)') fermi3e_int, fermi4e_int
    ! write(6,'(99es15.7)') fermi3e, fermi4e
    ! write(6,'(99es15.7)') fermi3_0, fermi4_0

    fermi4e = fermi4e_int
    fermi3e = fermi3e_int
    fermi4p = fermi4p_int
    fermi3p = fermi3p_int
    
    ferm43e = fermi4e/fermi3e
    ferm43p = fermi4p/fermi3p

    teme = tem/memev
    ene_n  = (ferm43e + ferm43p) *tem *0.5d0
    nemis_n= pair_nrate_base_e * tem**8 *fermi3e*fermi3p!1.06327d26 *teme**8 
    nemis_o= pair_nrate_base_x * tem**8 *fermi3e*fermi3p!2.28337d25 *teme**8 *fermi3e*fermi3p
    !nemis_n= 1.06327d26 *teme**8 *fermi3e*fermi3p
    !nemis_o= 2.28337d25 *teme**8 *fermi3e*fermi3p
    eemis_n= pair_erate_base_e * tem**8 *fermi3e*fermi3p * ene_n
    eemis_o= pair_erate_base_x * tem**8 *fermi3e*fermi3p * ene_n
    
  end subroutine rate_pairann


  subroutine rate_brems(rho,tem,xn,xp,nemis_brems,eemis_brems,ene_brems)
    real(8),intent(in)  :: rho,tem,xn,xp
    real(8),intent(out) :: nemis_brems,eemis_brems,ene_brems
    
    real(8) :: fac_frac
    real(8),parameter :: f283 = 28d0/3d0
    !ene_brems = 4.364d0/2d0*tem
    ene_brems = 2.182d0*tem

    fac_frac = (xn*xn + xp*xp + f283*xn*xp)
    
    eemis_brems = brems_erate_base * tem**5.5d0 * rho**2 * fac_frac
    nemis_brems = eemis_brems/ene_brems * erg_to_mev
    
  end subroutine rate_brems

  subroutine rate_plasm(eta_e,tem,nemis_plasm_e,nemis_plasm_x,eemis_plasm_e,eemis_plasm_x,ene_plasm)
    real(8),intent(in)  :: eta_e,tem
    real(8),intent(out) :: nemis_plasm_e,nemis_plasm_x,eemis_plasm_e,eemis_plasm_x,ene_plasm
    
    real(8) :: gamma_p,gamma_p2,gamma_p6

    gamma_p = 2d0*sqrt(alpha_fine/(9d0*pi)*(pi*pi+3d0*eta_e))
    gamma_p2= gamma_p*gamma_p
    gamma_p6= gamma_p2*gamma_p2*gamma_p2
    
    ene_plasm = 0.5d0*(2d0+gamma_p2/(1d0+gamma_p))*tem
    !ene_plasm = 0.5d0*tem
    
    nemis_plasm_e = plasm_nrate_base_e *tem**8 *gamma_p6 *(1d0+gamma_p) *exp(-gamma_p)
    nemis_plasm_x = plasm_nrate_base_x *tem**8 *gamma_p6 *(1d0+gamma_p) *exp(-gamma_p)

    eemis_plasm_e = nemis_plasm_e * ene_plasm*mev_to_erg
    eemis_plasm_x = nemis_plasm_x * ene_plasm*mev_to_erg
    
  end subroutine rate_plasm

  subroutine init_weak_table(fn,fnn)
    character(*),intent(in) :: fn,fnn
    character(1) :: str1
    integer :: i_temp,i_eta
    real(8) :: hoge

    open(11,file=fn,status="old")
    read(11,*) str1,n_temp,n_eta

    allocate(&
         caprate_logtemp(n_temp),&
         caprate_logeta (n_eta),&
         caprate_num_e(n_temp,n_eta),&
         caprate_num_p(n_temp,n_eta),&
         caprate_ene_e(n_temp,n_eta),&
         caprate_ene_p(n_temp,n_eta))
    
    do i_eta=1,n_eta
       read(11,*)
       do i_temp=1,n_temp
          read(11,*) caprate_logtemp(i_temp),caprate_logeta(i_eta),&
               caprate_num_e(i_temp,i_eta),&
               caprate_num_p(i_temp,i_eta),&
               caprate_ene_e(i_temp,i_eta),&
               caprate_ene_p(i_temp,i_eta)
       enddo
    enddo

    close(11)

    caprate_num_e(:,:)=log10(caprate_num_e(:,:))
    caprate_num_p(:,:)=log10(caprate_num_p(:,:))
    caprate_ene_e(:,:)=log10(caprate_ene_e(:,:))
    caprate_ene_p(:,:)=log10(caprate_ene_p(:,:))

    caprate_logeta(:) = log10(caprate_logeta(:))
    caprate_logtemp(:) = log10(caprate_logtemp(:))

    caprate_logeta_min = caprate_logeta(1)
    caprate_logtemp_min = caprate_logtemp(1)
    
    caprate_dlogeta = caprate_logeta(2) - caprate_logeta(1)
    caprate_dlogtemp = caprate_logtemp(2) - caprate_logtemp(1)

    ! negative eta
    open(11,file=fnn,status="old")
    read(11,*)
    
    allocate(&
         caprate_neg_num_e(n_temp,n_eta),&
         caprate_neg_num_p(n_temp,n_eta),&
         caprate_neg_ene_e(n_temp,n_eta),&
         caprate_neg_ene_p(n_temp,n_eta))

    do i_eta=1,n_eta
       read(11,*)
       do i_temp=1,n_temp
          read(11,*) hoge,hoge, &
               caprate_neg_num_e(i_temp,i_eta),&
               caprate_neg_num_p(i_temp,i_eta),&
               caprate_neg_ene_e(i_temp,i_eta),&
               caprate_neg_ene_p(i_temp,i_eta)
       enddo
    enddo

    close(11)

    caprate_neg_num_e(:,:)=log10(caprate_neg_num_e(:,:))
    caprate_neg_num_p(:,:)=log10(caprate_neg_num_p(:,:))
    caprate_neg_ene_e(:,:)=log10(caprate_neg_ene_e(:,:))
    caprate_neg_ene_p(:,:)=log10(caprate_neg_ene_p(:,:))
    
  end subroutine init_weak_table
  
  subroutine rate_cap_interp_table(temp,eta,caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p)
    real(8),intent(in) :: temp, eta
    real(8),intent(out) :: caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p
    
    integer :: i_eta,i_temp,i_eta1,i_temp1
    real(8) :: ee1,tt1,ee0,tt0
    
    i_eta = max(1,min(n_eta-1,int( (log10(eta)-caprate_logeta_min)/caprate_dlogeta ) + 1))
    i_temp = max(1,min(n_temp-1,int( (log10(temp)-caprate_logtemp_min)/caprate_dlogtemp ) + 1))
    i_eta1=i_eta+1
    i_temp1=i_temp+1
    ee1 = max(0d0,min(1d0,(log10(eta)-caprate_logeta(i_eta))/caprate_dlogeta))
    ee0 = 1d0-ee1
    tt1 = max(0d0,min(1d0,(log10(temp)-caprate_logtemp(i_temp))/caprate_dlogtemp))
    tt0 = 1d0-tt1
    
    caprate1_num_e = tt0*ee0*caprate_num_e(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_num_e(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_num_e(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_num_e(i_temp1,i_eta1)
    caprate1_num_p = tt0*ee0*caprate_num_p(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_num_p(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_num_p(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_num_p(i_temp1,i_eta1)
    caprate1_ene_e = tt0*ee0*caprate_ene_e(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_ene_e(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_ene_e(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_ene_e(i_temp1,i_eta1)
    caprate1_ene_p = tt0*ee0*caprate_ene_p(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_ene_p(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_ene_p(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_ene_p(i_temp1,i_eta1)

    caprate1_num_e = 1d1**caprate1_num_e
    caprate1_num_p = 1d1**caprate1_num_p
    caprate1_ene_e = 1d1**caprate1_ene_e
    caprate1_ene_p = 1d1**caprate1_ene_p

  end subroutine rate_cap_interp_table

  subroutine rate_cap_block_interp_table(temp, eta_e, eta_n, eta_a,caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p)
    real(8),intent(in) :: temp, eta_e, eta_n, eta_a
    real(8),intent(out) :: caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p
    
    integer :: i_eta,i_temp,i_eta1,i_temp1
    real(8) :: ee1,tt1,ee0,tt0

    real(8) :: zeta_n, zeta_np, eta_n_mod, eta_a_mod
    real(8) :: caprateb_num_e,caprateb_num_p,caprateb_ene_e,caprateb_ene_p
    real(8) :: fac_ec, fac_pc

    i_eta = max(1,min(n_eta-1,int( (log10(eta_e)-caprate_logeta_min)/caprate_dlogeta ) + 1))
    i_temp = max(1,min(n_temp-1,int( (log10(temp)-caprate_logtemp_min)/caprate_dlogtemp ) + 1))
    i_eta1=i_eta+1
    i_temp1=i_temp+1
    ee1 = max(0d0,min(1d0,(log10(eta_e)-caprate_logeta(i_eta))/caprate_dlogeta))
    ee0 = 1d0-ee1
    tt1 = max(0d0,min(1d0,(log10(temp)-caprate_logtemp(i_temp))/caprate_dlogtemp))
    tt0 = 1d0-tt1
    
    caprate1_num_e = tt0*ee0*caprate_num_e(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_num_e(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_num_e(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_num_e(i_temp1,i_eta1)
    caprate1_num_p = tt0*ee0*caprate_num_p(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_num_p(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_num_p(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_num_p(i_temp1,i_eta1)
    caprate1_ene_e = tt0*ee0*caprate_ene_e(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_ene_e(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_ene_e(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_ene_e(i_temp1,i_eta1)
    caprate1_ene_p = tt0*ee0*caprate_ene_p(i_temp ,i_eta ) &
                    +tt1*ee0*caprate_ene_p(i_temp1,i_eta ) &
                    +tt0*ee1*caprate_ene_p(i_temp ,i_eta1) &
                    +tt1*ee1*caprate_ene_p(i_temp1,i_eta1)

    caprate1_num_e = 1d1**caprate1_num_e
    caprate1_num_p = 1d1**caprate1_num_p
    caprate1_ene_e = 1d1**caprate1_ene_e
    caprate1_ene_p = 1d1**caprate1_ene_p


    ! Blocking by electron neutrinos
    zeta_n = - dmnp/temp
    eta_n_mod = eta_n - zeta_n
    
    i_eta = max(1,min(n_eta-1,int( (log10(abs(eta_n_mod))-caprate_logeta_min)/caprate_dlogeta ) + 1))
    i_temp = max(1,min(n_temp-1,int( (log10(temp)-caprate_logtemp_min)/caprate_dlogtemp ) + 1))
    i_eta1=i_eta+1
    i_temp1=i_temp+1
    ee1 = max(0d0,min(1d0,(log10(abs(eta_n_mod))-caprate_logeta(i_eta))/caprate_dlogeta))
    ee0 = 1d0-ee1
    tt1 = max(0d0,min(1d0,(log10(temp)-caprate_logtemp(i_temp))/caprate_dlogtemp))
    tt0 = 1d0-tt1
    
    if(eta_n_mod>0d0)then
       
       caprateb_num_e = tt0*ee0*caprate_num_e(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_num_e(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_num_e(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_num_e(i_temp1,i_eta1)
       caprateb_ene_e = tt0*ee0*caprate_ene_e(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_ene_e(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_ene_e(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_ene_e(i_temp1,i_eta1)
    else

       caprateb_num_e = tt0*ee0*caprate_neg_num_e(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_neg_num_e(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_neg_num_e(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_neg_num_e(i_temp1,i_eta1)
       caprateb_ene_e = tt0*ee0*caprate_neg_ene_e(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_neg_ene_e(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_neg_ene_e(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_neg_ene_e(i_temp1,i_eta1)
    endif

    caprateb_num_e = 1d1**caprateb_num_e
    caprateb_ene_e = 1d1**caprateb_ene_e

    ! Blocking by electron antineutrinos
    zeta_np=   dmnp/temp
    eta_a_mod = eta_a - zeta_np

    i_eta = max(1,min(n_eta-1,int( (log10(abs(eta_a_mod))-caprate_logeta_min)/caprate_dlogeta ) + 1))
    i_temp = max(1,min(n_temp-1,int( (log10(temp)-caprate_logtemp_min)/caprate_dlogtemp ) + 1))
    i_eta1=i_eta+1
    i_temp1=i_temp+1
    ee1 = max(0d0,min(1d0,(log10(abs(eta_a_mod))-caprate_logeta(i_eta))/caprate_dlogeta))
    ee0 = 1d0-ee1
    tt1 = max(0d0,min(1d0,(log10(temp)-caprate_logtemp(i_temp))/caprate_dlogtemp))
    tt0 = 1d0-tt1

    if(eta_a_mod>0d0)then

       caprateb_num_p = tt0*ee0*caprate_num_p(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_num_p(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_num_p(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_num_p(i_temp1,i_eta1)
       caprateb_ene_p = tt0*ee0*caprate_ene_p(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_ene_p(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_ene_p(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_ene_p(i_temp1,i_eta1)
    else
       
       !write(6,'(99es12.4)') log10(abs(eta_a_mod)),caprate_logeta(i_eta),caprate_logeta(i_eta1)
       !stop
       caprateb_num_p = tt0*ee0*caprate_neg_num_p(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_neg_num_p(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_neg_num_p(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_neg_num_p(i_temp1,i_eta1)
       caprateb_ene_p = tt0*ee0*caprate_neg_ene_p(i_temp ,i_eta ) &
                       +tt1*ee0*caprate_neg_ene_p(i_temp1,i_eta ) &
                       +tt0*ee1*caprate_neg_ene_p(i_temp ,i_eta1) &
                       +tt1*ee1*caprate_neg_ene_p(i_temp1,i_eta1)

    endif

    caprateb_num_p = 1d1**caprateb_num_p
    caprateb_ene_p = 1d1**caprateb_ene_p
    
    fac_ec = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,eta_n_mod -eta_e)) ))
    fac_pc = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,eta_a_mod +eta_e)) ))

    write(6,'(99es12.4)') fac_ec, fac_pc
    write(6,'(99es12.4)') caprate1_num_p,caprateb_num_p

    caprate1_num_e = fac_ec*max(0d0, caprate1_num_e-caprateb_num_e)
    caprate1_num_p = fac_pc*max(0d0, caprate1_num_p-caprateb_num_p)

  end subroutine rate_cap_block_interp_table

  subroutine nrate_cap_direct(tem,eta,caprate1_num_e,caprate1_num_p)
    real(8),intent(in) :: tem, eta
    real(8),intent(out) :: caprate1_num_e,caprate1_num_p

    real(8) :: zeta_e, zeta_p,tem5,finte,fintp
    
    zeta_e = -dmnp/tem
    zeta_p =  dmnp/tem
    
    tem5=tem*tem*tem*tem*tem
    ! tem6=tem5*tem
    
    call fermicap(zeta_e, eta,tem,0 ,finte)
    call fermicap(zeta_p,-eta,tem,0 ,fintp)
    !call fermicap(zeta_e, eta,tem,1 ,finte_e)
    !call fermicap(zeta_p,-eta,tem,1 ,fintp_e)
    caprate1_num_e=cap_nrate_base*finte*tem5
    caprate1_num_p=cap_nrate_base*fintp*tem5

  end subroutine nrate_cap_direct


  subroutine nrate_cap_block_direct(tem,eta,eta_n,eta_a,caprate1_num_e,caprate1_num_p)
    real(8),intent(in) :: tem, eta, eta_n, eta_a
    real(8),intent(out) :: caprate1_num_e,caprate1_num_p

    real(8) :: zeta_e, zeta_p,tem5,finte,fintp,fintn,finta, fac_ec, fac_pc
    
    zeta_e = -dmnp/tem
    zeta_p =  dmnp/tem
    
    tem5=tem*tem*tem*tem*tem
    ! tem6=tem5*tem
    
    call fermicap(zeta_e, eta,tem,0 ,finte)
    call fermicap(zeta_p,-eta,tem,0 ,fintp)
    !call fermicap(zeta_e, eta,tem,1 ,finte_e)
    !call fermicap(zeta_p,-eta,tem,1 ,fintp_e)

    call fermicap(zeta_e ,eta_n-zeta_e,tem,0,fintn)
    call fermicap(zeta_p ,eta_a-zeta_p,tem,0,finta)

    fac_ec = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,eta_n-zeta_e-eta)) ))
    fac_pc = 1.d0/(1.d0 - exp(max(-1.d2,min(1.d2,eta_a-zeta_p+eta)) ))

    caprate1_num_e=cap_nrate_base*(finte-fintn)*tem5!*fac_ec
    caprate1_num_p=cap_nrate_base*(fintp-finta)*tem5!*fac_pc
    
  end subroutine nrate_cap_block_direct

  ! subroutine fermi_integral(k,eta,fermik)
  !   real(8),intent(in) :: k,eta
  !   real(8),intent(out) :: fermik
    
    
  !   real(8) :: dt,et,xl, xfin,dx
  !   real(8) :: x, x_1,x_2,x_3,x_4, ex_1,ex_2,ex_3,ex_4, f_1,f_2,f_3,f_4, residual

  !   xfin = max(1000.d0,10.d0*(eta-xl))

  !   fint = 0.d0
  !   do i=0,nint(xfin/dx)
  !      x = dx*dble(i)
       
  !      ! x_1 = x;          ex_1 = max(-3d2,min(3d2,x_1+xl-eta))
  !      ! f_1 = (x_1 + xl + zeta)**(2+ih) * ( (x_1 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_1 + xl)**2) *1.d0/(exp(ex_1)+1d0)
  !      ! x_2 = x+0.5d0*dx; ex_2 = max(-3d2,min(3d2,x_2+xl-eta))
  !      ! f_2 = (x_2 + xl + zeta)**(2+ih) * ( (x_2 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_2 + xl)**2) *1.d0/(exp(ex_2)+1d0)
  !      ! x_3 = x+0.5d0*dx; ex_3 = max(-3d2,min(3d2,x_3+xl-eta))
  !      ! f_3 = (x_3 + xl + zeta)**(2+ih) * ( (x_3 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_3 + xl)**2) *1.d0/(exp(ex_3)+1d0)
  !      ! x_4 = x+dx;       ex_4 = max(-3d2,min(3d2,x_4+xl-eta))
  !      ! f_4 = (x_4 + xl + zeta)**(2+ih) * ( (x_4 + xl)**2 - et**2 ) * sqrt(1.d0 - et**2/(x_4 + xl)**2) *1.d0/(exp(ex_4)+1d0)
       
  !      x_1 = x
  !      x_2 = x+0.5d0*dx
  !      x_3 = x+0.5d0*dx
  !      x_4 = x+dx
  !      f_1 = capture_rate(zeta,tem,eta,x_1,ih)
  !      f_2 = capture_rate(zeta,tem,eta,x_2,ih)
  !      f_3 = f_2!capture_rate(zeta,tem,eta,x_3,ih)
  !      f_4 = capture_rate(zeta,tem,eta,x_4,ih)
       
  !      fint = fint + (f_1 + 2d0*f_2 + 2d0*f_3 + f_4)/6d0*dx
  !      !fint = fint + f_1*dx
       
  !      ex_4 = max(-3d2,min(3d2,x_4+xl-eta))
       
  !      residual = exp(-ex_4)* ( (x_4**4 + 4.d0*x_4**3 + 12.d0*x_4**2 + 24.d0*x_4 + 24.d0)  * (5d0*x_4)**ih + 120d0*dble(ih) )
       
  !      !if(x>10.d0*abs(eta-xl).and.x>10.d0*abs(zeta).and.x>10.d0*et.and.x>10.d0*xl.and.fint > 1.d15*residual) goto 10
  !      if(x > 10.d0*max(abs(eta-xl),xl) .and. fint > 1.d15*residual) goto 10
  !   enddo

  ! end subroutine fermi_integral

  subroutine make_epcap_table(eta_min,eta_max,tem_min,tem_max,n_eta,n_tem)
    use omp_lib

    ! A + e- -> B + nu_e
    ! zeta = (m_A - m_B)/T
    ! eta = eta_ele or eta_pos
    ! tem = temperature in MeV
    ! ih = 1 -> energy, 0 -> number
    ! fint = Integral

    real(8),intent(in) :: eta_min,eta_max,tem_min,tem_max
    integer,intent(in) :: n_eta,n_tem
    
    integer :: i,j
    real(8) :: eta,tem,zeta_e,zeta_p,finte,fintp,finte_e,fintp_e
    
    real(8) :: tem5,tem6

    real(8),allocatable :: nrate_ecap(:,:),nrate_pcap(:,:),erate_ecap(:,:),erate_pcap(:,:)

    allocate(nrate_ecap(n_tem,n_eta),nrate_pcap(n_tem,n_eta),erate_ecap(n_tem,n_eta),erate_pcap(n_tem,n_eta))

    !$omp parallel private(i,j,eta ,tem ,zeta_e ,zeta_p ,tem5,tem6,finte,fintp,finte_e,fintp_e)
    !$omp do
    do i=1,n_tem
       do j=1,n_eta
          
          eta = -1d1**(log10(eta_min) + (log10(eta_max)-log10(eta_min))*dble(j-1)/dble(n_eta-1))
          tem = 1d1**(log10(tem_min) + (log10(tem_max)-log10(tem_min))*dble(i-1)/dble(n_tem-1))
          
          ! if(omp_get_thread_num()==0)write(6,*)eta,tem
          
          zeta_e = -dmnp/tem
          zeta_p =  dmnp/tem

          tem5=tem*tem*tem*tem*tem
          tem6=tem5*tem
          
          call fermicap(zeta_e, eta,tem,0 ,finte)
          call fermicap(zeta_p,-eta,tem,0 ,fintp)
          call fermicap(zeta_e, eta,tem,1 ,finte_e)
          call fermicap(zeta_p,-eta,tem,1 ,fintp_e)
          ! write(6,'(99es12.3e3)') tem,eta, &
          !      cap_nrate_base*finte*tem5,&
          !      cap_nrate_base*fintp*tem5,&
          !      cap_erate_base*finte_e*tem6,&
          !      cap_erate_base*fintp_e*tem6
          nrate_ecap(i,j) = cap_nrate_base*finte*tem5
          nrate_pcap(i,j) = cap_nrate_base*fintp*tem5
          erate_ecap(i,j) = cap_erate_base*finte_e*tem6
          erate_pcap(i,j) = cap_erate_base*fintp_e*tem6
          
       enddo
       write(6,*) i
    enddo
    !$omp end do
    !$omp end parallel

    open(11,file="caprate_negeta.dat",status="replace")
    write(11,*) "#",n_tem,n_eta
    do j=1,n_eta
       write(11,*)
       do i=1,n_tem

          eta = -1d1**(log10(eta_min) + (log10(eta_max)-log10(eta_min))*dble(j-1)/dble(n_eta-1))
          tem = 1d1**(log10(tem_min) + (log10(tem_max)-log10(tem_min))*dble(i-1)/dble(n_tem-1))

          write(11,'(99es16.7e3)') tem,eta, &
               nrate_ecap(i,j) , &
               nrate_pcap(i,j) , &
               erate_ecap(i,j) , &
               erate_pcap(i,j)

       enddo
    enddo
    close(11)

  end subroutine make_epcap_table

  subroutine planck_mean_nu_p_scattering(tem,sigma_nup)
    real(8),intent(in) :: tem
    real(8),intent(out) :: sigma_nup
    
    sigma_nup = sigma0*0.25d0*(4d0*wsin4-2d0*wsin2+0.25d0*(1d0+3d0*ga2))&
         * 310d0/147d0*pi*pi *(tem/memev)**2
    
  end subroutine planck_mean_nu_p_scattering
  
  subroutine planck_mean_nu_n_scattering(tem,sigma_nun)
    real(8),intent(in) :: tem
    real(8),intent(out) :: sigma_nun
    
    sigma_nun = sigma0*0.25d0*(1d0+3d0*ga2)*0.25d0 &
         * 310d0/147d0*pi*pi *(tem/memev)**2
    
  end subroutine planck_mean_nu_n_scattering

  subroutine planck_mean_nu_p_scattering_tr(tem,sigma_nup_tr)
    real(8),intent(in) :: tem
    real(8),intent(out) :: sigma_nup_tr
    
    sigma_nup_tr = sigma0/6d0*( (c_v-1d0)**2 + 5d0*ga2*(c_a-1d0)**2) &
         * 310d0/147d0*pi*pi *(tem/memev)**2
    
  end subroutine planck_mean_nu_p_scattering_tr

  subroutine planck_mean_nu_n_scattering_tr(tem,sigma_nun_tr)
    real(8),intent(in) :: tem
    real(8),intent(out) :: sigma_nun_tr
    
    sigma_nun_tr = sigma0/24d0* (1d0 + 5d0*ga2) &
         * 310d0/147d0*pi*pi *(tem/memev)**2
    
  end subroutine planck_mean_nu_n_scattering_tr
  
end module module_weak_interaction
