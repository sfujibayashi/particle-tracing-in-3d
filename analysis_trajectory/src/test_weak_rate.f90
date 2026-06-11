subroutine test_weak_rate
  use module_weak_interaction

  implicit none

  real(8) :: finte,fintp,finte_e,fintp_e

  real(8) :: zeta_e,zeta_p, eta, tem, rho, xp, xn
  integer :: ih

  integer :: i,j
  real(8) :: caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p
  
  real(8) :: nemis_pair_e,nemis_pair_x,eemis_pair_e,eemis_pair_x,ene_pair
  real(8) :: nemis_brems,eemis_brems,ene_brems
  real(8) :: nemis_plasm_e,nemis_plasm_x,eemis_plasm_e,eemis_plasm_x,ene_plasm
  real(8) :: ecap_nrate,ecap_erate,pcap_nrate,pcap_erate

  real(8) :: fermi3_int, fermi4_int
  real(8) :: fmaske, deta, fermi4e, fermi3e

  real(8),parameter :: ppi = pi**2
    
  ! do i=1,100
  !    eta = -1d1**(2d0-4d0*dble(i-1)/dble(100-1))
  !    call fermik( eta,3,fermi3_int)
  !    call fermik( eta,4,fermi4_int)
     
  !    fmaske= (1.d0+sign(1.d0,eta))*0.5d0
  !    deta = exp(-abs(eta))
     
  !    fermi4e = ( eta**5 *0.2d0          &
  !         +eta**3 *2d0/3d0 *ppi       &
  !         +eta    *48.d0          &
  !         +24.d0*deta      ) * fmaske &
  !         +24.d0*deta        *(1.d0-fmaske)
     
  !    fermi3e = ( eta**4 *0.25d0         &
  !         +eta**2 *0.5d0 *ppi     &
  !         +12.d0 -6.d0*deta )     &
  !         * fmaske &
  !           + 6.d0 *deta *(1.d0-fmaske)
     
  !    write(99,'(99es14.6)') eta,fermi3_int,fermi4_int, fermi3e, fermi4e
  ! enddo
  ! do i=1,100
  !    eta = 1d1**(-2d0+4d0*dble(i-1)/dble(100-1))
  !    call fermik( eta,3,fermi3_int)
  !    call fermik( eta,4,fermi4_int)


  !    fmaske= (1.d0+sign(1.d0,eta))*0.5d0
  !    deta = exp(-abs(eta))
     
  !    fermi4e = ( eta**5 *0.2d0          &
  !         +eta**3 *2d0/3d0 *ppi       &
  !         +eta    *48.d0          &
  !         +24.d0*deta      ) * fmaske &
  !         +24.d0*deta        *(1.d0-fmaske)
     
  !    fermi3e = ( eta**4 *0.25d0         &
  !         +eta**2 *0.5d0 *ppi     &
  !         +12.d0 -6.d0*deta )     &
  !         * fmaske &
  !           + 6.d0 *deta *(1.d0-fmaske)
     
  !    write(99,'(99es14.6)') eta,fermi3_int,fermi4_int, fermi3e, fermi4e
  ! enddo
  ! stop
  
  tem=1d11/11.604525006d9
  eta=0d0
  xn=0.5d0
  xp=0.5d0
  rho=1d13

  write(6,*)
  write(6,*) "electron/positron capture"
  call rate_cap_interp_table(tem,eta,caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p)
  ecap_nrate = caprate1_num_e*rho*na*xp
  ecap_erate = caprate1_ene_e*rho*na*xp
  pcap_nrate = caprate1_num_p*rho*na*xn
  pcap_erate = caprate1_ene_p*rho*na*xn
  
  write(6,'("eta, tem, number(ele), number(pos), energy(ele), energy(pos), ene(ele), ene(pos) = ",99es15.7)') eta,tem, &
       ecap_nrate,pcap_nrate, ecap_erate,pcap_erate, ecap_erate/ecap_nrate*erg_to_mev, pcap_erate/pcap_nrate*erg_to_mev

  write(6,'("total energy loss rate = ",99es12.4)') ecap_erate+pcap_erate

  ! zeta_e = -dmnp/tem
  ! zeta_p =  dmnp/tem
  ! call fermicap(zeta_e, eta,tem,0 ,finte)
  ! call fermicap(zeta_p,-eta,tem,0 ,fintp)
  ! call fermicap(zeta_e, eta,tem,1 ,finte_e)
  ! call fermicap(zeta_p,-eta,tem,1 ,fintp_e)
  
  ! write(6,'("eta, tem, number(ele), number(pos), energy(ele), energy(pos), ene(ele), ene(pos) = ",99es15.7)') eta,tem, &
  !            cap_nrate_base*finte*tem**5,&
  !            cap_nrate_base*fintp*tem**5,&
  !            cap_erate_base*finte_e*tem**6,&
  !            cap_erate_base*fintp_e*tem**6,&
  !            cap_erate_base*finte_e*tem**6/(cap_nrate_base*finte*tem**5)*erg_to_mev, &
  !            cap_erate_base*fintp_e*tem**6/(cap_nrate_base*fintp*tem**5)*erg_to_mev
  ! stop

  write(6,*)
  write(6,*) "Pair-annihilation"
  write(6,'("sekig coeff (e): ",99es15.7)') pair_nrate_base_e*memev**8
  write(6,'("sekig coeff (x): ",99es15.7)') pair_nrate_base_x*memev**8

  call rate_pairann(eta,tem,nemis_pair_e,nemis_pair_x,eemis_pair_e,eemis_pair_x,ene_pair)
  write(6,'("eta, tem, number(e), number(x), energy(e), energy(x), ene = ",99es12.4)') eta,tem, nemis_pair_e, nemis_pair_x, eemis_pair_e, eemis_pair_x, ene_pair
  write(6,'("total energy loss rate = ",99es12.4)') eemis_pair_e*2d0+eemis_pair_x*4d0
  write(6,'("total energy loss rate = ",99es12.4)') 

  !
  write(6,*)
  write(6,*) "Plasmon decay"
  write(6,'("sekig coeff (e): ",99es15.7)') c_v**2/(192d0*pi**3*alpha_fine)/memev**2 *sigma0*clight/hbc**6 *memev**8 *mev_to_erg**6
  write(6,'("sekig coeff (x): ",99es15.7)') (c_v-1d0)**2/(192d0*pi**3*alpha_fine)/memev**2 *sigma0*clight/hbc**6 *memev**8 *mev_to_erg**6

  call rate_plasm(eta,tem,nemis_plasm_e,nemis_plasm_x,eemis_plasm_e,eemis_plasm_x,ene_plasm)
  write(6,'("eta, tem, number(e), number(x), energy(e), energy(x), ene = ",99es12.4)') eta,tem, nemis_plasm_e, nemis_plasm_x, eemis_plasm_e, eemis_plasm_x, ene_plasm
  write(6,'("total energy loss rate = ",99es12.4)') eemis_plasm_e*2d0+eemis_plasm_x*4d0

  !
  write(6,*)
  write(6,*) "Bremsstrahlung"
  write(6,'("sekig coeff    : ",99es15.7)')  brems_erate_base * memev**4.5d0/2.182d0/mev_to_erg
  call rate_brems(rho,tem,xn,xp,nemis_brems,eemis_brems,ene_brems)

  write(6,'("rho, tem, number, energy, ene = ",99es12.4)') rho,tem, nemis_brems, eemis_brems, ene_brems
  write(6,'("total energy loss rate = ",99es12.4)') eemis_brems*6d0
  


end subroutine test_weak_rate
