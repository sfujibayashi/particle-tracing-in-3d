module module_eos
  implicit none
  
  integer :: ntemp, nrho, nye
  
  real(8) :: rho_e_min,tem_e_min,ye_e_min
  real(8) :: drho_e, dye_e, dtem_e
  real(8) :: drhoi  , dyei, dtemi

  real(8),allocatable :: rho_e(:),  ye_e(:), tem_e (:) &
        ,pres_e (:,:,:), eps_e (:,:,:)  &
        ,sen_e  (:,:,:), cs_e  (:,:,:)  &
        ,chn_e  (:,:,:), chp_e (:,:,:), che_e (:,:,:)  &
        ,aa_e   (:,:,:), zz_e  (:,:,:)  &
        ,xn_e   (:,:,:), xp_e  (:,:,:)  &
        ,xA_e   (:,:,:), xy_e  (:,:,:)  &
        ,xd_e   (:,:,:), xt_e  (:,:,:), xh_e  (:,:,:)

  real(8) :: hhh_min

  real(8),allocatable :: h_min_tab(:)

contains

  subroutine readeos(fn,nrho_in,ntemp_in,nye_in)
    use unit
    character(*),intent(in) :: fn
    integer,intent(in) :: nrho_in,ntemp_in,nye_in
    integer :: irho,itemp,iye
    real(8) :: ePmin , eEmin, eCmin

    real(8) :: h

    nrho=nrho_in
    ntemp=ntemp_in
    nye=nye_in

    allocate(  rho_e(nrho),  ye_e(0:nye), tem_e (ntemp)                  &
        ,pres_e (ntemp,0:nye,nrho), eps_e (ntemp,0:nye,nrho)  &
        ,sen_e  (ntemp,0:nye,nrho), cs_e  (ntemp,0:nye,nrho)  &
        ,chn_e  (ntemp,0:nye,nrho), chp_e (ntemp,0:nye,nrho), che_e (ntemp,0:nye,nrho)  &
        ,aa_e   (ntemp,0:nye,nrho), zz_e  (ntemp,0:nye,nrho)  &
        ,xn_e   (ntemp,0:nye,nrho), xp_e  (ntemp,0:nye,nrho)  &
        ,xA_e   (ntemp,0:nye,nrho), xy_e  (ntemp,0:nye,nrho)  &
        ,xd_e   (ntemp,0:nye,nrho), xt_e  (ntemp,0:nye,nrho), xh_e  (ntemp,0:nye,nrho) )
    
    
    open(20,file=trim(adjustl(fn)) ,status="old", action="read", form="binary")

    read(20) tem_e, ye_e, rho_e         &
         ,pres_e, eps_e, sen_e, cs_e &
         ,chn_e , chp_e, che_e       &
         ,aa_e  , zz_e               &
         ,xn_e  , xp_e , xA_e        &
         ,xd_e  , xt_e , xh_e , xy_e &
         ,ePmin , eEmin, eCmin
    close(20)
    
    
    do itemp = 1 ,ntemp
       tem_e(itemp) = log10(tem_e(itemp))
    enddo
    do irho = 1,nrho
       rho_e(irho) = log10(rho_e(irho))
    enddo
    rho_e_min = rho_e(1)
    tem_e_min = tem_e(1)
    
    do itemp = 1 ,ntemp
       do iye = 1 ,nye
          do irho = 1 ,nrho
             pres_e (itemp,iye,irho) = log10(pres_e (itemp,iye,irho))
             eps_e  (itemp,iye,irho) = log10( eps_e (itemp,iye,irho)/v_uni**2 + 1)
             cs_e   (itemp,iye,irho) = log10(  cs_e (itemp,iye,irho))
          enddo
       enddo
    enddo
    
    ye_e_min = 0d0
    rho_e_min = rho_e(1)
    tem_e_min = tem_e(1)
    
    drho_e = rho_e(2) - rho_e(1)
    dtem_e = tem_e(2) - tem_e(1)
    dye_e  = ye_e (2) - ye_e (1)
    
    drhoi = 1d0/drho_e
    dtemi = 1d0/dtem_e
    dyei  = 1d0/dye_e

    ! write(6,*) drho_e,dtem_e, dye_e
    hhh_min = 1d99
    do itemp = 1 ,ntemp
       do iye = 1 ,nye
          do irho = 1 ,nrho
             if(.not.isnan(eps_e(itemp,iye,irho)))then
                h = 10d0**eps_e(itemp,iye,irho) + 10d0**pres_e(itemp,iye,irho)/10d0**rho_e(irho)/v_uni**2
                hhh_min = min(hhh_min,h)
             endif
          enddo
       enddo
    enddo

    ! iye=50
    ! itemp=1
    ! irho=1
    ! write(6,*) rho_e(irho), tem_e(itemp), ye_e(iye)
    ! h = 10d0**eps_e(itemp,iye,irho) + 10d0**pres_e(itemp,iye,irho)/10d0**rho_e(irho)/v_uni**2
    ! hhh_min = h
    
    write(6,*) "Minimum enthalpy : ",  hhh_min

    ! Ye-dependent h-min
    allocate(h_min_tab(nye))
    block
      real(8) :: hm
      do iye=1,nye
         hm=1d99
         
         do itemp = 1 ,ntemp
            do irho = 1 ,nrho
               if(.not.isnan(eps_e(itemp,iye,irho)))then
                  h = 10d0**eps_e(itemp,iye,irho) + 10d0**pres_e(itemp,iye,irho)/10d0**rho_e(irho)/v_uni**2
                  hm = min(hm,h)
               endif
            enddo
         enddo
         h_min_tab(iye) = hm
         ! write(6,*) ye_e(iye), hm
      enddo
    end block
    
  end subroutine readeos


  subroutine get_h_min_ye(ye,h)

    real(8),intent(in) :: ye
    real(8),intent(out):: h
    
    integer :: iye, iye1
    real(8) :: tt, ttp

    iye  = max(1 , min(nye-1, int((ye-ye_e_min )*dyei)+1))
    iye1 = iye +1
    tt   = max(0.d0, min(1.d0 ,     (ye-ye_e(iye))*dyei)   )
    ttp  = 1.d0-tt

    h = ttp*h_min_tab(iye) + tt *h_min_tab(iye1)

    ! write(6,*) tt,ttp
    ! write(6,*) hmintab_ye(iye),ye,hmintab_ye(iye1)
    ! write(6,*) hmintab_hhh(iye),hhh_min,hmintab_hhh(iye1)
    
  end subroutine get_h_min_ye

  
end module module_eos
