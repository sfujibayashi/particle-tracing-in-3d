module module_equilibrium_ye
end module module_equilibrium_ye

subroutine test_ye_equil
  use module_weak_interaction
  use module_eos

  implicit none
  
  real(8) :: rho,tem,ye,eta,caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p,   ecap_nrate, ecap_erate, pcap_nrate, pcap_erate, xp, xn

  integer :: item, irho, iye, ntem, nrho
  real(8) :: rho_max, rho_min, tem_max, tem_min, ye_max, ye_min

  real(8) :: ss,tt,uu,ss1,tt1,uu1

  real(8) :: item1,irho1,iye1,iye_loop

  integer :: itr,itrlim=30

  real(8) :: ye1,ye2,ye0, f1,f2,f0
  
  real(8) :: rne, rae, abs_n_nrate,abs_a_nrate

  real(8) :: munu, ye_munu0, ye_cap

  ! iye = 11
  ! irho = 326
  ! do item = 51,76
  !    write(6,'(99es15.7)') tem_e(item), rho_e(irho), ye_e(iye), xn_e(item,iye,irho), xp_e(item,iye,irho)
     
  ! enddo
  ! stop


  ! tem = 1.d0
  ! rho = 1d10

  ! call ye_equilibrium_munu0(rho,tem,ye_munu0)
  ! call ye_equilibrium_capture(rho,tem,ye_cap,eta,xn,xp,ecap_nrate,pcap_nrate)

  ! write(6,'(99es12.4)') rho,tem, ye_cap, ye_munu0
  ! stop
  

  
  ! rne = 1d30
  ! rae = 1d30
  ! rho=1d9
  ! tem=1d0
  ! ! ye=0.3d0
  
  ! call ye_equilibrium_abs(rho,tem,rne,rae,ye,xn,xp,abs_n_nrate,abs_a_nrate)
  
  ! write(6,'(99es12.4)') ye,xn,xp,abs_n_nrate,abs_a_nrate
  ! !call nrate_abs(rho,tem,ye,rne,rae,xn,xp,abs_n_nrate,abs_a_nrate)
  ! stop

  rho_min = 1d6
  rho_max = 1d13
  tem_min = 1d-1
  tem_max = 3d1

  ye_min = 0.01d0
  ye_max = 0.60d0
  
  ntem=100
  nrho=100

  open(11,file="ye_equil_table.dat",status="replace",action="write")
  write(11,'("#",99a15)') "rho", "T", "ye(equil)", "eta", "Xn", "Xp", "Rec", "Rpc", "Ye(munu=0)"
  do irho = 1,nrho
     write(11,*)
     do item=1,ntem


        tem = 1.d1**(log10(tem_min)+(log10(tem_max)-log10(tem_min))*dble(item-1)/dble(ntem-1))
        rho = 1.d1**(log10(rho_min)+(log10(rho_max)-log10(rho_min))*dble(irho-1)/dble(nrho-1))
        
        
        ! ye1=ye_min
        ! ye2=ye_max

        ! call nrate(rho,tem,ye1,eta,xn,xp,ecap_nrate,pcap_nrate)
        ! f1 = (ecap_nrate - pcap_nrate)/(ecap_nrate + pcap_nrate)
        ! call nrate(rho,tem,ye2,eta,xn,xp,ecap_nrate,pcap_nrate)
        ! f2 = (ecap_nrate - pcap_nrate)/(ecap_nrate + pcap_nrate)
        
        ! !write(6,*) f1,f2

        ! if(f1*f2>0d0)then
           
        !    if(abs(f1) > abs(f2))then
        !       ye0=ye2
        !    else
        !       ye0=ye1
        !    endif
           
        ! else
        !    do itr=1,itrlim
        !       ye0 = 0.5d0*(ye1+ye2)
        !       call nrate(rho,tem,ye0,eta,xn,xp,ecap_nrate,pcap_nrate)

        !       f0 = (ecap_nrate - pcap_nrate)/(ecap_nrate + pcap_nrate)

        !       if( f0*f1>0d0)then
        !          f1=f0
        !          ye1=ye0
        !       else
        !          f2=f0
        !          ye2=ye0
        !       endif
              
        !       !write(6,'(99es12.4)') ye1,ye0,ye2,f1,f0,f2
              
        !    enddo
        ! endif

        ! call nrate(rho,tem,ye0,eta,xn,xp,ecap_nrate,pcap_nrate)

        call ye_equilibrium_capture(rho,tem,ye0,eta,xn,xp,ecap_nrate,pcap_nrate)
        call ye_equilibrium_munu0(rho,tem,ye_munu0)
        
        write(11,'(" ",99es15.6e3)') rho,tem,ye0,eta,xn,xp, ecap_nrate,pcap_nrate, ye_munu0
        ! write(6 ,'(99es13.4e3)') rho,tem,ye0,eta,xn,xp, ecap_nrate,pcap_nrate, ye_munu0
        !stop
     enddo
  enddo

  close(11)

end subroutine test_ye_equil

subroutine nrate_cap(rho,tem,ye,eta,xn,xp,ecap_nrate, pcap_nrate)
  use module_eos
  use module_weak_interaction

  implicit none

  real(8),intent(in) :: rho,tem,ye
  real(8),intent(out):: ecap_nrate, pcap_nrate,eta,xn,xp
  real(8) :: uu,ss,tt,uu1,ss1,tt1
  integer :: item,irho,iye,item1,irho1,iye1

  real(8) :: caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p

  call interp_coef(rho,tem,ye,irho,item,iye,uu,ss,tt)
  
  item1=item+1
  irho1=irho+1
  iye1=iye+1
  ss1=1d0-ss
  tt1=1d0-tt
  uu1=1d0-uu

  eta =  ss1 *tt1 *uu1 *log10(che_e(item ,iye ,irho ))   &
       + ss  *tt1 *uu1 *log10(che_e(item1,iye ,irho ))   &
       + ss1 *tt  *uu1 *log10(che_e(item ,iye1,irho ))   &
       + ss1 *tt1 *uu  *log10(che_e(item ,iye ,irho1))   &
       + ss  *tt  *uu1 *log10(che_e(item1,iye1,irho ))   &
       + ss  *tt1 *uu  *log10(che_e(item1,iye ,irho1))   &
       + ss1 *tt  *uu  *log10(che_e(item ,iye1,irho1))   &
       + ss  *tt  *uu  *log10(che_e(item1,iye1,irho1))
  eta = 1d1**eta

  xn  =  ss1 *tt1 *uu1 *log10(xn_e(item ,iye ,irho ))   &
       + ss  *tt1 *uu1 *log10(xn_e(item1,iye ,irho ))   &
       + ss1 *tt  *uu1 *log10(xn_e(item ,iye1,irho ))   &
       + ss1 *tt1 *uu  *log10(xn_e(item ,iye ,irho1))   &
       + ss  *tt  *uu1 *log10(xn_e(item1,iye1,irho ))   &
       + ss  *tt1 *uu  *log10(xn_e(item1,iye ,irho1))   &
       + ss1 *tt  *uu  *log10(xn_e(item ,iye1,irho1))   &
       + ss  *tt  *uu  *log10(xn_e(item1,iye1,irho1))
  xn = 1d1**xn
  
  xp  =  ss1 *tt1 *uu1 *log10(xp_e(item ,iye ,irho ))   &
       + ss  *tt1 *uu1 *log10(xp_e(item1,iye ,irho ))   &
       + ss1 *tt  *uu1 *log10(xp_e(item ,iye1,irho ))   &
       + ss1 *tt1 *uu  *log10(xp_e(item ,iye ,irho1))   &
       + ss  *tt  *uu1 *log10(xp_e(item1,iye1,irho ))   &
       + ss  *tt1 *uu  *log10(xp_e(item1,iye ,irho1))   &
       + ss1 *tt  *uu  *log10(xp_e(item ,iye1,irho1))   &
       + ss  *tt  *uu  *log10(xp_e(item1,iye1,irho1))
  xp = 1d1**xp

  ! eta =  ss1 *tt1 *uu1 *che_e(item ,iye ,irho )   &
  !      + ss  *tt1 *uu1 *che_e(item1,iye ,irho )   &
  !      + ss1 *tt  *uu1 *che_e(item ,iye1,irho )   &
  !      + ss1 *tt1 *uu  *che_e(item ,iye ,irho1)   &
  !      + ss  *tt  *uu1 *che_e(item1,iye1,irho )   &
  !      + ss  *tt1 *uu  *che_e(item1,iye ,irho1)   &
  !      + ss1 *tt  *uu  *che_e(item ,iye1,irho1)   &
  !      + ss  *tt  *uu  *che_e(item1,iye1,irho1)

  ! xn  =  ss1 *tt1 *uu1 *xn_e(item ,iye ,irho )   &
  !      + ss  *tt1 *uu1 *xn_e(item1,iye ,irho )   &
  !      + ss1 *tt  *uu1 *xn_e(item ,iye1,irho )   &
  !      + ss1 *tt1 *uu  *xn_e(item ,iye ,irho1)   &
  !      + ss  *tt  *uu1 *xn_e(item1,iye1,irho )   &
  !      + ss  *tt1 *uu  *xn_e(item1,iye ,irho1)   &
  !      + ss1 *tt  *uu  *xn_e(item ,iye1,irho1)   &
  !      + ss  *tt  *uu  *xn_e(item1,iye1,irho1)

  ! xp  =  ss1 *tt1 *uu1 *xp_e(item ,iye ,irho )   &
  !      + ss  *tt1 *uu1 *xp_e(item1,iye ,irho )   &
  !      + ss1 *tt  *uu1 *xp_e(item ,iye1,irho )   &
  !      + ss1 *tt1 *uu  *xp_e(item ,iye ,irho1)   &
  !      + ss  *tt  *uu1 *xp_e(item1,iye1,irho )   &
  !      + ss  *tt1 *uu  *xp_e(item1,iye ,irho1)   &
  !      + ss1 *tt  *uu  *xp_e(item ,iye1,irho1)   &
  !      + ss  *tt  *uu  *xp_e(item1,iye1,irho1)



  call rate_cap_interp_table(tem,eta,caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p)
  ! call nrate_cap_direct(tem,eta,caprate1_num_e,caprate1_num_p)
  
  ! unit of 1/s
  ecap_nrate = caprate1_num_e*xp
  pcap_nrate = caprate1_num_p*xn

end subroutine nrate_cap

subroutine nrate_cap_block(rho,tem,ye,eta_n,eta_a,eta, xn,xp,ecap_nrate, pcap_nrate)
  use module_eos
  use module_weak_interaction

  implicit none

  real(8),intent(in) :: rho,tem,ye,eta_n,eta_a
  real(8),intent(out):: ecap_nrate, pcap_nrate,eta,xn,xp
  real(8) :: uu,ss,tt,uu1,ss1,tt1
  integer :: item,irho,iye,item1,irho1,iye1

  real(8) :: caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p

  call interp_coef(rho,tem,ye,irho,item,iye,uu,ss,tt)
  
  item1=item+1
  irho1=irho+1
  iye1=iye+1
  ss1=1d0-ss
  tt1=1d0-tt
  uu1=1d0-uu

  eta =  ss1 *tt1 *uu1 *log10(che_e(item ,iye ,irho ))   &
       + ss  *tt1 *uu1 *log10(che_e(item1,iye ,irho ))   &
       + ss1 *tt  *uu1 *log10(che_e(item ,iye1,irho ))   &
       + ss1 *tt1 *uu  *log10(che_e(item ,iye ,irho1))   &
       + ss  *tt  *uu1 *log10(che_e(item1,iye1,irho ))   &
       + ss  *tt1 *uu  *log10(che_e(item1,iye ,irho1))   &
       + ss1 *tt  *uu  *log10(che_e(item ,iye1,irho1))   &
       + ss  *tt  *uu  *log10(che_e(item1,iye1,irho1))
  eta = 1d1**eta

  xn  =  ss1 *tt1 *uu1 *log10(xn_e(item ,iye ,irho ))   &
       + ss  *tt1 *uu1 *log10(xn_e(item1,iye ,irho ))   &
       + ss1 *tt  *uu1 *log10(xn_e(item ,iye1,irho ))   &
       + ss1 *tt1 *uu  *log10(xn_e(item ,iye ,irho1))   &
       + ss  *tt  *uu1 *log10(xn_e(item1,iye1,irho ))   &
       + ss  *tt1 *uu  *log10(xn_e(item1,iye ,irho1))   &
       + ss1 *tt  *uu  *log10(xn_e(item ,iye1,irho1))   &
       + ss  *tt  *uu  *log10(xn_e(item1,iye1,irho1))
  xn = 1d1**xn
  
  xp  =  ss1 *tt1 *uu1 *log10(xp_e(item ,iye ,irho ))   &
       + ss  *tt1 *uu1 *log10(xp_e(item1,iye ,irho ))   &
       + ss1 *tt  *uu1 *log10(xp_e(item ,iye1,irho ))   &
       + ss1 *tt1 *uu  *log10(xp_e(item ,iye ,irho1))   &
       + ss  *tt  *uu1 *log10(xp_e(item1,iye1,irho ))   &
       + ss  *tt1 *uu  *log10(xp_e(item1,iye ,irho1))   &
       + ss1 *tt  *uu  *log10(xp_e(item ,iye1,irho1))   &
       + ss  *tt  *uu  *log10(xp_e(item1,iye1,irho1))
  xp = 1d1**xp

  ! eta =  ss1 *tt1 *uu1 *che_e(item ,iye ,irho )   &
  !      + ss  *tt1 *uu1 *che_e(item1,iye ,irho )   &
  !      + ss1 *tt  *uu1 *che_e(item ,iye1,irho )   &
  !      + ss1 *tt1 *uu  *che_e(item ,iye ,irho1)   &
  !      + ss  *tt  *uu1 *che_e(item1,iye1,irho )   &
  !      + ss  *tt1 *uu  *che_e(item1,iye ,irho1)   &
  !      + ss1 *tt  *uu  *che_e(item ,iye1,irho1)   &
  !      + ss  *tt  *uu  *che_e(item1,iye1,irho1)

  ! xn  =  ss1 *tt1 *uu1 *xn_e(item ,iye ,irho )   &
  !      + ss  *tt1 *uu1 *xn_e(item1,iye ,irho )   &
  !      + ss1 *tt  *uu1 *xn_e(item ,iye1,irho )   &
  !      + ss1 *tt1 *uu  *xn_e(item ,iye ,irho1)   &
  !      + ss  *tt  *uu1 *xn_e(item1,iye1,irho )   &
  !      + ss  *tt1 *uu  *xn_e(item1,iye ,irho1)   &
  !      + ss1 *tt  *uu  *xn_e(item ,iye1,irho1)   &
  !      + ss  *tt  *uu  *xn_e(item1,iye1,irho1)

  ! xp  =  ss1 *tt1 *uu1 *xp_e(item ,iye ,irho )   &
  !      + ss  *tt1 *uu1 *xp_e(item1,iye ,irho )   &
  !      + ss1 *tt  *uu1 *xp_e(item ,iye1,irho )   &
  !      + ss1 *tt1 *uu  *xp_e(item ,iye ,irho1)   &
  !      + ss  *tt  *uu1 *xp_e(item1,iye1,irho )   &
  !      + ss  *tt1 *uu  *xp_e(item1,iye ,irho1)   &
  !      + ss1 *tt  *uu  *xp_e(item ,iye1,irho1)   &
  !      + ss  *tt  *uu  *xp_e(item1,iye1,irho1)


  ! call rate_cap_block_interp_table(tem,eta,eta_n,eta_a,caprate1_num_e,caprate1_num_p,caprate1_ene_e,caprate1_ene_p)

  call nrate_cap_block_direct(tem,eta,eta_n,eta_a,caprate1_num_e,caprate1_num_p)
  
  ! call nrate_cap_direct(tem,eta,caprate1_num_e,caprate1_num_p)
  
  ! unit of 1/s
  ecap_nrate = caprate1_num_e*xp
  pcap_nrate = caprate1_num_p*xn

end subroutine nrate_cap_block

subroutine ye_equilibrium_capture(rho,tem,ye,eta,xn,xp,ecap_nrate,pcap_nrate)
  use module_eos
  use module_weak_interaction

  implicit none

  real(8),intent(in)  :: rho,tem
  real(8),intent(out) :: ye,eta,xn,xp,ecap_nrate,pcap_nrate

  real(8) :: ye_min,ye_max
  real(8) :: ye1,ye2,ye0,f1,f2,f0

  integer :: itr,itrlim=30
  
  ye_min = 0.01d0
  ye_max = 0.60d0

  ye1=ye_min
  ye2=ye_max

  call nrate_cap(rho,tem,ye1,eta,xn,xp,ecap_nrate,pcap_nrate)
  f1 = (ecap_nrate - pcap_nrate)/(ecap_nrate + pcap_nrate)
  call nrate_cap(rho,tem,ye2,eta,xn,xp,ecap_nrate,pcap_nrate)
  f2 = (ecap_nrate - pcap_nrate)/(ecap_nrate + pcap_nrate)

  !write(6,*) f1,f2
  !if(tem>4.99d0) write(6,'(99e12.4)')rho,tem,eta
  
  if(f1*f2>0d0)then

     if(abs(f1) > abs(f2))then
        ye0=ye2
     else
        ye0=ye1
     endif

  else
     do itr=1,itrlim
        ye0 = 0.5d0*(ye1+ye2)
        call nrate_cap(rho,tem,ye0,eta,xn,xp,ecap_nrate,pcap_nrate)

        f0 = (ecap_nrate - pcap_nrate)/(ecap_nrate + pcap_nrate)

        if( f0*f1>0d0)then
           f1=f0
           ye1=ye0
        else
           f2=f0
           ye2=ye0
        endif

        !write(6,'(99es12.4)') ye1,ye0,ye2,f1,f0,f2

     enddo
  endif

  call nrate_cap(rho,tem,ye0,eta,xn,xp,ecap_nrate,pcap_nrate)

  ye = ye0

end subroutine ye_equilibrium_capture

subroutine nrate_abs(rho,tem,ye,rne,rae,xn,xp,abs_n_nrate,abs_a_nrate)
  use unit, only : rho_uni, v_uni
  use module_eos
  use module_weak_interaction

  implicit none
  
  real(8),intent(in) :: rho,tem,ye,rne,rae
  real(8),intent(out) :: xn,xp,abs_n_nrate,abs_a_nrate

  real(8) :: rne_comp, rae_comp, tem4, rne_tmp, rae_tmp
  integer :: ien,iea
  real(8) :: renu, rena,vvn,vvnp, vva,vvap, ch_nuf, ch_naf

  real(8) :: uu,ss,tt,uu1,ss1,tt1
  integer :: item,irho,iye,item1,irho1,iye1

  real(8) :: eta,eta_e,eta_nu,eta_na


  real(8),parameter :: f16=1d0/6d0, f56=5d0/6d0, f76=7d0/6d0, f23=2d0/3d0, f13=1d0/3d0, ppi=pi*pi, ppi2=ppi*ppi
  real(8) :: fmask,etax,fd5,fd4,fd3,fd2,ave53,ave43,ave32,ave23,ave54,eblk,blk
  real(8) :: fda3, ave33, dmnpx

  call interp_coef(rho,tem,ye,irho,item,iye,uu,ss,tt)
  
  item1=item+1
  irho1=irho+1
  iye1=iye+1
  ss1=1d0-ss
  tt1=1d0-tt
  uu1=1d0-uu

  eta_e =  ss1 *tt1 *uu1 *che_e(item ,iye ,irho )   &
       + ss  *tt1 *uu1 *che_e(item1,iye ,irho )   &
       + ss1 *tt  *uu1 *che_e(item ,iye1,irho )   &
       + ss1 *tt1 *uu  *che_e(item ,iye ,irho1)   &
       + ss  *tt  *uu1 *che_e(item1,iye1,irho )   &
       + ss  *tt1 *uu  *che_e(item1,iye ,irho1)   &
       + ss1 *tt  *uu  *che_e(item ,iye1,irho1)   &
       + ss  *tt  *uu  *che_e(item1,iye1,irho1)

  xn  =  ss1 *tt1 *uu1 *xn_e(item ,iye ,irho )   &
       + ss  *tt1 *uu1 *xn_e(item1,iye ,irho )   &
       + ss1 *tt  *uu1 *xn_e(item ,iye1,irho )   &
       + ss1 *tt1 *uu  *xn_e(item ,iye ,irho1)   &
       + ss  *tt  *uu1 *xn_e(item1,iye1,irho )   &
       + ss  *tt1 *uu  *xn_e(item1,iye ,irho1)   &
       + ss1 *tt  *uu  *xn_e(item ,iye1,irho1)   &
       + ss  *tt  *uu  *xn_e(item1,iye1,irho1)

  xp  =  ss1 *tt1 *uu1 *xp_e(item ,iye ,irho )   &
       + ss  *tt1 *uu1 *xp_e(item1,iye ,irho )   &
       + ss1 *tt  *uu1 *xp_e(item ,iye1,irho )   &
       + ss1 *tt1 *uu  *xp_e(item ,iye ,irho1)   &
       + ss  *tt  *uu1 *xp_e(item1,iye1,irho )   &
       + ss  *tt1 *uu  *xp_e(item1,iye ,irho1)   &
       + ss1 *tt  *uu  *xp_e(item ,iye1,irho1)   &
       + ss  *tt  *uu  *xp_e(item1,iye1,irho1)

  
!!!
  tem4 = tem**4

  rne_tmp = max(rne,10d0**esce0*tem4*rho_uni*v_uni**2)
  rae_tmp = max(rae,10d0**esce0*tem4*rho_uni*v_uni**2)
  !rne_comp = max(10d0**esce0,rne/(rho_uni*v_uni**2))
  !rae_comp = max(10d0**esce0,rae/(rho_uni*v_uni**2))

  rne_comp = rne_tmp/(rho_uni*v_uni**2)
  rae_comp = rae_tmp/(rho_uni*v_uni**2)
  

  renu = log10(rne_comp/tem4)
  ien  = max(1   , min( int((renu-esce0    )*dscei)+1, ient-1))
  vvn  = max(0.d0, min(     (renu-esce(ien))*dscei   , 1.d0  ))
  vvnp = 1.d0-vvn
  eta_nu = vvnp*ech_ne(ien) + vvn *ech_ne(ien+1)
  
  rena = log10(rae_comp/tem4)
  iea  = max(1   , min( int((rena-esce0    )*dscei)+1, ient-1))
  vva  = max(0.d0, min(     (rena-esce(iea))*dscei   , 1.d0  ))
  vvap = 1.d0-vva
  eta_na = vvap*ech_ne(iea) + vva *ech_ne(iea+1)

  eta = eta_nu
  fmask = (1.d0+dsign(1.d0,(eta-0.d0)))*0.5d0
  etax  = dexp(-eta)
  fd5 = eta**6 *f16        &
       +eta**4 *f56 *ppi   &
       +eta**2 *f76 *ppi2  &
       +240.d0             &
       -120.d0*etax
  fd4 = eta**5 *0.2d0      &
       +eta**3 *f23 *ppi   &
       +eta    *48.d0      &
       + 24.d0*etax
  fd3 = eta**4 *0.25d0     &
       +eta**2 *0.5d0 *ppi &
       +12.d0              &
       - 6.d0*etax
  fd2 = eta**3 *f13        &
       +eta    *4.d0       &
       + 2.d0*etax
  
  ave53 = fd5/fd3 *tem**2 *fmask + 20.d0   *tem**2 *(1.d0-fmask)
  ave43 = fd4/fd3 *tem    *fmask + 4.d0    *tem    *(1.d0-fmask)
  ave32 = fd3/fd2 *tem    *fmask + 3.d0    *tem    *(1.d0-fmask)
  ave23 = fd2/fd3 /tem    *fmask +  f13    /tem    *(1.d0-fmask)
  ave54 = fd5/fd4 *tem    *fmask + 5.d0    *tem    *(1.d0-fmask)

  eblk = max(-200.d0,min(200.d0,-(ave54+dmnp)/tem+eta_e))
  blk = 1.d0/( exp(eblk) + 1.d0 )

  abs_n_nrate = sigma0*clight*(1d0+3d0*ga2)*0.25d0* (ave43 + 2.d0*dmnp       + dmnp**2*ave23)/memev**2 * blk * xn * rne_tmp*erg_to_mev
  ! abs_n_nrate = sigma0*clight*(1d0+3d0*ga2)*0.25d0* (ave43)/memev**2 * rne*erg_to_mev

  !write(6,'(99es12.4)') sigma0*(1d0+3d0*ga2)*0.25d0, xn, blk
  !write(6,'(99es12.4)') eta_nu,ave43,rne,abs_n_nrate
!!!
  eta = eta_na
  if(eta.gt.0.d0) then
     etax = dexp(-eta)
     fda3= eta**4 *0.25d0     &
          +eta**2 *0.5d0 *ppi &
          +12.d0              &
          - 6.d0*etax
  else
     fda3 = 6.d0 *dexp(eta)
  endif
  
  eta = eta_na - dmnp/tem
  if(eta.gt.0.d0) then
     etax = dexp(-eta)
     fd5 = eta**6 *f16        &
          +eta**4 *f56 *ppi   &
          +eta**2 *f76 *ppi2  & 
          +240.d0             &
          -120.d0*etax
     fd4 = eta**5 *0.2d0      &
          +eta**3 *f23 *ppi   &
          +eta    *48.d0      &
          + 24.d0*etax
     fd3 = eta**4 *0.25d0     &
          +eta**2 *0.5d0 *ppi &
          +12.d0              &
          - 6.d0*etax
     fd2 = eta**3 *f13        &
          +eta    *4.d0       &
          + 2.d0*etax
     ave53 = fd5 /fda3 *tem**2 
     ave43 = fd4 /fda3 *tem
     ave33 = fd3 /fda3
     ave23 = fd2 /fda3 /tem
     ave54 = fd5 /fd4  *tem
  else 
     if(eta_na.gt.0.d0) then
        etax = dexp(eta)
        ave53 = 120.d0*etax /fda3 *tem**2
        ave43 =  24.d0*etax /fda3 *tem
        ave33 =   6.d0*etax /fda3
        ave23 =   2.d0*etax /fda3 /tem
     else
        dmnpx = dexp(-dmnp/tem)
        ave53 = 20.d0*dmnpx *tem**2
        ave43 =  4.d0*dmnpx *tem
        ave33 =       dmnpx 
        ave23 =   f13*dmnpx /tem
     end if
     ave54 = 5.d0    *tem
  endif
  
  eblk=max(-200.d0,min(200.d0,-(max(ave54-dmnp,0.d0))/tem-eta_e))
  blk = 1.d0/( exp(eblk) + 1.d0 ) !*etapn

  abs_a_nrate = sigma0*clight*(1d0+3d0*ga2)*0.25d0* (ave43 +2.d0*dmnp +dmnp**2*ave23)/memev**2 * blk * xp * rae_tmp*erg_to_mev
  
  !write(6,'(99es13.5)') rho,tem,ye,eta_nu, eta_na, ave43,  !., abs_n_nrate/xn, abs_a_nrate/xp,xn,xp
  !stop
end subroutine nrate_abs


subroutine ye_equilibrium_abs(rho,tem,rne,rae,ye,xn,xp,abs_n_nrate,abs_a_nrate)
  implicit none
  real(8),intent(in) :: rho,tem,rne,rae
  real(8),intent(out) :: ye,xn,xp,abs_n_nrate,abs_a_nrate

  real(8) :: ye_min,ye_max
  real(8) :: ye1,ye2,ye0,f1,f2,f0

  integer :: itr,itrlim=50

  ye_min = 0.01d0
  ye_max = 0.60d0

  ye1=ye_min
  ye2=ye_max

  call nrate_abs(rho,tem,ye1,rne,rae,xn,xp,abs_n_nrate,abs_a_nrate)
  f1 = (abs_a_nrate - abs_n_nrate)/(abs_a_nrate + abs_n_nrate + 1.d-99)
  call nrate_abs(rho,tem,ye2,rne,rae,xn,xp,abs_n_nrate,abs_a_nrate)
  f2 = (abs_a_nrate - abs_n_nrate)/(abs_a_nrate + abs_n_nrate + 1.d-99)

  ! write(6,*) f1,f2
  
  if(f1*f2>0d0)then

     if(abs(f1) > abs(f2))then
        ye0=ye2
     else
        ye0=ye1
     endif

  else
     do itr=1,itrlim
        ye0 = 0.5d0*(ye1+ye2)
        call nrate_abs(rho,tem,ye0,rne,rae,xn,xp,abs_n_nrate,abs_a_nrate)

        f0 = (abs_a_nrate - abs_n_nrate)/(abs_a_nrate + abs_n_nrate + 1.d-99)

        if( f0*f1>0d0)then
           f1=f0
           ye1=ye0
        else
           f2=f0
           ye2=ye0
        endif

        ! write(6,'(99es12.4)') ye1,ye0,ye2,f1,f0,f2

     enddo
  endif

  call nrate_abs(rho,tem,ye0,rne,rae,xn,xp,abs_n_nrate,abs_a_nrate)

  ye = ye0
  
end subroutine ye_equilibrium_abs


subroutine munu_from_table(rho,tem,ye,munu)
  use module_eos

  implicit none

  real(8),intent(in) :: rho,tem,ye
  real(8),intent(out) :: munu
  
  real(8) :: uu,ss,tt,uu1,ss1,tt1
  integer :: item,irho,iye,item1,irho1,iye1

  real(8) :: eta_e,mu_n,mu_p

  call interp_coef(rho,tem,ye,irho,item,iye,uu,ss,tt)
  
  item1=item+1
  irho1=irho+1
  iye1=iye+1
  ss1=1d0-ss
  tt1=1d0-tt
  uu1=1d0-uu

  eta_e= ss1 *tt1 *uu1 *che_e(item ,iye ,irho )   &
       + ss  *tt1 *uu1 *che_e(item1,iye ,irho )   &
       + ss1 *tt  *uu1 *che_e(item ,iye1,irho )   &
       + ss1 *tt1 *uu  *che_e(item ,iye ,irho1)   &
       + ss  *tt  *uu1 *che_e(item1,iye1,irho )   &
       + ss  *tt1 *uu  *che_e(item1,iye ,irho1)   &
       + ss1 *tt  *uu  *che_e(item ,iye1,irho1)   &
       + ss  *tt  *uu  *che_e(item1,iye1,irho1)

  mu_n = ss1 *tt1 *uu1 *chn_e(item ,iye ,irho )   &
       + ss  *tt1 *uu1 *chn_e(item1,iye ,irho )   &
       + ss1 *tt  *uu1 *chn_e(item ,iye1,irho )   &
       + ss1 *tt1 *uu  *chn_e(item ,iye ,irho1)   &
       + ss  *tt  *uu1 *chn_e(item1,iye1,irho )   &
       + ss  *tt1 *uu  *chn_e(item1,iye ,irho1)   &
       + ss1 *tt  *uu  *chn_e(item ,iye1,irho1)   &
       + ss  *tt  *uu  *chn_e(item1,iye1,irho1)

  mu_p = ss1 *tt1 *uu1 *chp_e(item ,iye ,irho )   &
       + ss  *tt1 *uu1 *chp_e(item1,iye ,irho )   &
       + ss1 *tt  *uu1 *chp_e(item ,iye1,irho )   &
       + ss1 *tt1 *uu  *chp_e(item ,iye ,irho1)   &
       + ss  *tt  *uu1 *chp_e(item1,iye1,irho )   &
       + ss  *tt1 *uu  *chp_e(item1,iye ,irho1)   &
       + ss1 *tt  *uu  *chp_e(item ,iye1,irho1)   &
       + ss  *tt  *uu  *chp_e(item1,iye1,irho1)

  munu = eta_e*tem + mu_p - mu_n
  
end subroutine munu_from_table

subroutine ye_munu0_linear(rho,tem,ye1,ye2,ye_munu0)
  use module_eos

  implicit none

  real(8),intent(in) :: rho,tem,ye1,ye2
  real(8),intent(out) :: ye_munu0
  
  real(8) :: uu,ss,tt,uu1,ss1,tt1
  integer :: item,irho,iye,item1,irho1,iye1,iye_1,iye_2

  real(8) :: eta_e,mu_n,mu_p, munu_1, munu_2

  call interp_coef(rho,tem,ye1,irho,item,iye_1,uu,ss,tt)
  iye_2  = max(jed , min(jeu-1, int((ye2-ye_e_min  )*dyei)+1))

  if(iye_1/=iye_2)then
     stop "not."
  endif

  iye_2=iye_1+1

  item1=item+1
  irho1=irho+1
  ss1=1d0-ss
  uu1=1d0-uu

  eta_e= ss1 *uu1 *che_e(item ,iye_1 ,irho )   &
       + ss  *uu1 *che_e(item1,iye_1 ,irho )   &
       + ss1 *uu  *che_e(item ,iye_1 ,irho1)   &
       + ss  *uu  *che_e(item1,iye_1 ,irho1) 
  mu_n = ss1 *uu1 *chn_e(item ,iye_1 ,irho )   &
       + ss  *uu1 *chn_e(item1,iye_1 ,irho )   &
       + ss1 *uu  *chn_e(item ,iye_1 ,irho1)   &
       + ss  *uu  *chn_e(item1,iye_1 ,irho1) 
  mu_p = ss1 *uu1 *chp_e(item ,iye_1 ,irho )   &
       + ss  *uu1 *chp_e(item1,iye_1 ,irho )   &
       + ss1 *uu  *chp_e(item ,iye_1 ,irho1)   &
       + ss  *uu  *chp_e(item1,iye_1 ,irho1)
  munu_1 = eta_e*tem + mu_p - mu_n

  eta_e= ss1 *uu1 *che_e(item ,iye_2 ,irho )   &
       + ss  *uu1 *che_e(item1,iye_2 ,irho )   &
       + ss1 *uu  *che_e(item ,iye_2 ,irho1)   &
       + ss  *uu  *che_e(item1,iye_2 ,irho1) 
  mu_n = ss1 *uu1 *chn_e(item ,iye_2 ,irho )   &
       + ss  *uu1 *chn_e(item1,iye_2 ,irho )   &
       + ss1 *uu  *chn_e(item ,iye_2 ,irho1)   &
       + ss  *uu  *chn_e(item1,iye_2 ,irho1) 
  mu_p = ss1 *uu1 *chp_e(item ,iye_2 ,irho )   &
       + ss  *uu1 *chp_e(item1,iye_2 ,irho )   &
       + ss1 *uu  *chp_e(item ,iye_2 ,irho1)   &
       + ss  *uu  *chp_e(item1,iye_2 ,irho1)
  munu_2 = eta_e*tem + mu_p - mu_n
  
  !write(6,*) munu_1,munu_2
   
  ye_munu0 = ye_e(iye_1) + (ye_e(iye_2)-ye_e(iye_1))*(0d0-munu_1)/(munu_2-munu_1)
end subroutine ye_munu0_linear

subroutine ye_equilibrium_munu0(rho,tem,ye_munu0)
  use module_eos
  implicit none

  real(8),intent(in) :: rho,tem
  real(8),intent(out) :: ye_munu0
  
  real(8) :: ye,ye_min,ye_max,ye0,ye1,ye2,munu1,munu2,munu0,f1,f2,f0

  integer :: itr, itrlim=30

  ye_min=0.01d0
  ye_max=0.60d0

  ye1=ye_min
  ye2=ye_max
  
  call munu_from_table(rho,tem,ye1,munu1)
  f1=munu1/tem
  call munu_from_table(rho,tem,ye2,munu2)
  f2=munu2/tem

  if(f1*f2>0d0)then

     if(abs(f1) > abs(f2))then
        ye0=ye2
     else
        ye0=ye1
     endif
  else
     do itr=1,itrlim
        ye0 = 0.5d0*(ye1+ye2)
        call munu_from_table(rho,tem,ye0,munu0)
        
        f0 = munu0/tem
        
        if( f0*f1>0d0)then
           f1=f0
           ye1=ye0
        else
           f2=f0
           ye2=ye0
        endif
        
        !write(6,'(i5,99es12.4)') itr,ye1,ye0,ye2,f1,f0,f2
        
        if(int(ye1*100d0)==int(ye2*100d0))then
           exit
        endif
     
     enddo

     call ye_munu0_linear(rho,tem,ye1,ye2,ye0)
     
  endif

  
  ye_munu0 = ye0
end subroutine ye_equilibrium_munu0
