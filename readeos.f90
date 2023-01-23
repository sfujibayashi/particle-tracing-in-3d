subroutine readeos(ceos,ceosb,cynu,cenu)
!  use mpipara
  use eostab
  use eos00
  use unit
  implicit real(8) (a-h,o-z) 
  character(*) :: ceos,ceosb,cynu,cenu
  character(200) :: fn

  open(20,file=trim(adjustl(ceos)) ,status='old',form='binary')
  open(21,file=trim(adjustl(ceosb)),status='old',form='binary')
  open(22,file=trim(adjustl(cynu)) ,status='old',form='binary')
  open(23,file=trim(adjustl(cenu)) ,status='old',form='binary')
  
!      do inu = 1,inut
!         read(22,*) escn(inu), ech_nn(inu), eprnn(inu)
!      enddo
  read(22) escn, ech_nn, eprnn
  close(22)

!      do inu = 1,ient
!         read(23,*) esce(inu), ech_ne(inu), eprne(inu)  &
!                   ,ynue(inu), eave  (inu)
!      enddo
  read(23)esce, ech_ne, eprne, ynue, eave
  close(23)

  read(20) tem_e, ye_e, rho_e         &
          ,pres_e, eps_e, sen_e, cs_e &
          ,chn_e , chp_e, che_e       &
          ,aa_e  , zz_e               &
          ,xn_e  , xp_e , xA_e        &
          ,xd_e  , xt_e , xh_e , xy_e &
          ,ePmin , eEmin, eCmin
  close(20)
      
!       ePmin  = 1.d50
!       eEmin  = 1.d50
!       eCmin  = 1.d50
!       do ie = ied,ieu
!       do je = jed,jeu
!       do ke = ked,keu
!          read(20,'(19e28.18e3)') tem_e(ie), ye_e(je), rho_e(ke)      &
!               ,pres_e(ie,je,ke), eps_e(ie,je,ke), sen_e(ie,je,ke), cs_e(ie,je,ke)  &
!               ,chn_e (ie,je,ke), chp_e(ie,je,ke), che_e(ie,je,ke)  &
!               ,aa_e  (ie,je,ke), zz_e (ie,je,ke)                   &
!               ,xn_e  (ie,je,ke), xp_e (ie,je,ke), xA_e (ie,je,ke)  &
!               ,xd_e  (ie,je,ke), xt_e (ie,je,ke), xh_e (ie,je,ke), xy_e(ie,je,ke)

!          eps_e(ie,je,ke) = eps_e(ie,je,ke) !+ 1.045d0/1.660539274d-24/624150.6479963235d0
!          che_e(ie,je,ke) = che_e(ie,je,ke)/tem_e(ie)

!          if(pres_e(ie,je,ke).gt.0.d0) then
!             ePmin = min(ePmin, pres_e(ie,je,ke))
!          endif
!          if(eps_e(ie,je,ke).gt.0.d0) then
!             if(eps_e(ie,je,ke) .lt. eEmin) then
!                eEmin = eps_e(ie,je,ke)
!                ie_eEmin = ie
!                je_eEmin = je
!                ke_eEmin = ke
!             end if
!          endif
!          if(cs_e(ie,je,ke).gt.0.d0) then
!             eCmin = min(eCmin, cs_e(ie,je,ke))
!          endif
!       enddo
!       enddo
!       enddo
!       close(20)
!----------
      ! --- set region where relativistic beta-equil. EOS is valid
  iflag_beos = 0
  do ie = ierd,ieru
  do je = jerd,jeru
  do ke = kerd,keru
     iflag_beos(ie,je,ke)=1
  enddo
  enddo
  enddo
  ! ---

  read(21)yl_eb                          &
         ,pres_eb, eps_eb, sen_eb, cs_eb &
         ,chn_eb , chp_eb, che_eb        &
         ,aa_eb  , zz_eb                 &
         ,xn_eb  , xp_eb , xA_eb         &
         ,xd_eb  , xt_eb , xh_eb , xy_eb &
         ,ye_eb  , yn_eb , ya_eb         &
         ,ePnmin , eEnmin, eCnmin
  
!   if(myrank==0)then
!    do ie = ied,ieu
!      write(6,'(i4,99es13.4)') ie,tem_e(ie)
!    enddo
!    do je = jed,jeu
!      write(6,'(i4,99es13.4)') je, ye_e(je)
!    enddo
!    do ke = ked,keu
!      write(6,'(i4,99es13.4)') ke,rho_e(ke)
!    enddo
!    do je = jed,jeu
!      write(6,'(i4,99es13.4)') je,yl_eb(je)
!    enddo
!   endif
  ! call mpi_finalize(ierr)
  ! stop      
  
!   ePnmin = 1.d50
!   eEnmin = 1.d50
!   eCnmin = 1.d50
!   do ie = ied,ieu
!   do je = jed,jeu
!   do ke = ked,keu
!      read(21,'(22e28.18e3)')  hoge,  yl_eb(je),  hoge      &
!           ,pres_eb(ie,je,ke), eps_eb(ie,je,ke), sen_eb(ie,je,ke), cs_eb(ie,je,ke)  &
!           ,chn_eb (ie,je,ke), chp_eb(ie,je,ke), che_eb(ie,je,ke)  &
!           ,aa_eb  (ie,je,ke), zz_eb (ie,je,ke)                    &
!           ,xn_eb  (ie,je,ke), xp_eb (ie,je,ke), xA_eb (ie,je,ke)  &
!           ,xd_eb  (ie,je,ke), xt_eb (ie,je,ke), xh_eb (ie,je,ke), xy_eb(ie,je,ke)  &
!           ,ye_eb  (ie,je,ke), yn_eb (ie,je,ke), ya_eb (ie,je,ke)

!      eps_eb(ie,je,ke) = eps_eb(ie,je,ke) !+ 1.045d0/1.660539274d-24/624150.6479963235d0
!      che_eb(ie,je,ke) = che_eb(ie,je,ke)/tem_e(ie)

!      if(pres_eb(ie,je,ke).gt.0.d0) then
!         ePnmin = min(ePnmin, pres_eb(ie,je,ke))
!      endif
!      if(eps_eb(ie,je,ke).gt.0.d0) then
!         if(eps_eb(ie,je,ke) .lt. eEnmin) then
!            eEnmin = eps_eb(ie,je,ke)
!            ie_eEnmin = ie
!            je_eEnmin = je
!            ke_eEnmin = ke
!         end if
!      endif
!      if(cs_eb(ie,je,ke).gt.0.d0) then
!         eCnmin = min(eCnmin, cs_eb(ie,je,ke))
!      endif
!   enddo
!   enddo
!   enddo

  do ie = ied ,ieu
  do je = jed ,jeu
  do ke = ked ,keu
     ye_eb(ie,je,ke) = yl_eb(je) - yn_eb(ie,je,ke) + ya_eb(ie,je,ke) 
  enddo
  enddo
  enddo
  close(21)

!----------

  do ie = ied ,ieu
     tem_e(ie) = log10(tem_e(ie))
  enddo
  do ke = ked,keu
     rho_e(ke) = log10(rho_e(ke)/rho_uni)
  enddo
  rho_e_min = rho_e(1)
  tem_e_min = tem_e(1)

  do ie = ied ,ieu
  do je = jed ,jeu
  do ke = ked ,keu
!     pres_e (ie,je,ke) = max( ePmin, pres_e (ie,je,ke))
!     eps_e  (ie,je,ke) = max( eEmin,  eps_e (ie,je,ke))
!     cs_e   (ie,je,ke) = min(0.999d0, max( eCmin,  cs_e (ie,je,ke)))
!     pres_eb(ie,je,ke) = max(ePnmin, pres_eb(ie,je,ke))
!     eps_eb (ie,je,ke) = max(eEnmin,  eps_eb(ie,je,ke))
!     cs_eb  (ie,je,ke) = min(0.999d0, max(eCnmin,  cs_eb(ie,je,ke)))
     pres_e (ie,je,ke) = log10(pres_e (ie,je,ke)/rho_uni/v_uni**2)
     eps_e  (ie,je,ke) = log10( eps_e (ie,je,ke)/v_uni**2 + 1.d0)
     cs_e   (ie,je,ke) = log10(  cs_e (ie,je,ke))
     pres_eb(ie,je,ke) = log10(pres_eb(ie,je,ke)/rho_uni/v_uni**2)
     eps_eb (ie,je,ke) = log10( eps_eb(ie,je,ke)/v_uni**2 + 1.d0)
     cs_eb  (ie,je,ke) = log10(  cs_eb(ie,je,ke))
  enddo
  enddo
  enddo

!---------------------------------!
!   EOS at atmosphere
!---------------------------------!
  
  item_at = 1
  krho_at = 2
  rho_at = 10.d0**rho_e(krho_at) *1.00000001d0
  tem_at = 10.d0**tem_e(item_at) *1.00000001d0
   ye_at = 0.5d0
  jye_at = int((ye_at - ye_e_min    )*dyei) + 1
  fye_at =     (ye_at - ye_e(jye_at))*dyei
  peos_at = (1.d0-fye_at)*pres_e(item_at, jye_at  , krho_at)  &
           +      fye_at *pres_e(item_at, jye_at+1, krho_at)
  eeos_at = (1.d0-fye_at)* eps_e(item_at, jye_at  , krho_at)  &
           +      fye_at * eps_e(item_at, jye_at+1, krho_at)
  ceos_at = (1.d0-fye_at)*  cs_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  cs_e(item_at, jye_at+1, krho_at)
   eps_at = 10.d0**eeos_at - 1.d0
  pres_at = 10.d0**peos_at
    cs_at = 10.d0**ceos_at
  ehat_at = 1.d0 +eps_at
  eps2_at = pres_at/rho_at
   hhh_at = 1.d0 +eps_at +pres_at/rho_at
  ch_n_at = (1.d0-fye_at)* chn_e(item_at, jye_at  , krho_at)  &
           +      fye_at * chn_e(item_at, jye_at+1, krho_at)
  ch_p_at = (1.d0-fye_at)* chp_e(item_at, jye_at  , krho_at)  &
           +      fye_at * chp_e(item_at, jye_at+1, krho_at)
  ch_e_at = (1.d0-fye_at)* che_e(item_at, jye_at  , krho_at)  &
           +      fye_at * che_e(item_at, jye_at+1, krho_at)
   x_n_at = (1.d0-fye_at)*  xn_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xn_e(item_at, jye_at+1, krho_at)
   x_p_at = (1.d0-fye_at)*  xp_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xp_e(item_at, jye_at+1, krho_at)
   x_A_at = (1.d0-fye_at)*  xA_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xA_e(item_at, jye_at+1, krho_at)
   x_d_at = (1.d0-fye_at)*  xd_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xd_e(item_at, jye_at+1, krho_at)
   x_t_at = (1.d0-fye_at)*  xd_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xd_e(item_at, jye_at+1, krho_at)
   x_h_at = (1.d0-fye_at)*  xd_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xd_e(item_at, jye_at+1, krho_at)
   x_y_at = (1.d0-fye_at)*  xd_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  xd_e(item_at, jye_at+1, krho_at)
    aa_at = (1.d0-fye_at)*  aa_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  aa_e(item_at, jye_at+1, krho_at) 
    zz_at = (1.d0-fye_at)*  zz_e(item_at, jye_at  , krho_at)  &
           +      fye_at *  zz_e(item_at, jye_at+1, krho_at)
  fbeta = (ch_p_at-ch_n_at)/tem_at
  ch_nu_at = ch_e_at + fbeta
  ch_na_at =-ch_nu_at
  sen_at = 0.d0
  yn_at  = 1.d-99
  ya_at  = 1.d-99
  yo_at  = 1.d-99

!---------------------------------!
!   EOS miscellaneous
!---------------------------------!

  krho0 = 1
  rho00 = 10.d0**rho_e(krho0) *1.00000001d0
  jye0  = 1
  yemin =         ye_e(jye0 ) *1.00000001d0
  item0 = 1
  tem_at= 10.d0**tem_e(item0) *1.00000001d0
  tem01 = 0.1d0
  
  yemax = 0.50d0
  
  rho_at1 = 1.d2 
  
  rho6    = 1.d6 
  rho7    = 1.d7 
  rho8    = 1.d8 
  rho10   = 1.d10
  rho13   = 1.d13
  rho14   = 1.d14

  rho_core = rho8
  rho_nb   = rho10
  rhocut   = rho_at*1.1d0 

  qcrit0 = rho_at *2.d0
  rcrit0 = 1.d-40
  fqcrt  = 1.d-12
  frcrt  = 1.d-25
  ffcr   = 2.d0

  rhomin = qcrit0

  rnemin = rcrit0
  raemin = rcrit0
  roemin = rcrit0

  rnemax = rcrit0
  raemax = rcrit0
  roemax = rcrit0

!------------------------------------------------------!
!            Read Weak Rate Table
!------------------------------------------------------!


!      open(30,file='Data/Rate/nse_rate_v0.6.10.dat')
!      open(31,file='Data/Rate/PAIRrate_Itoh.dat')
!      open(32,file='Data/Rate/PLASMArate_Itoh2.dat')
!      open(33,file='Data/Rate/BREMSrate_Itoh.dat')

!      do jw = jwd,jwu
!      do iw = iwd,iwu
!      do kw = kwd,kwu
!         read(30,*)   &
!               yeNSE(jw), t9NSE(iw), rhoNSE(kw)    &
!             , zzNSE(iw,jw,kw),  aaNSE(iw,jw,kw)   &
!             ,qijNSE(iw,jw,kw)                     & 
!             ,dqnNSE(iw,jw,kw), dqaNSE(iw,jw,kw)   &
!             ,dynNSE(iw,jw,kw), dyaNSE(iw,jw,kw)   &
!             ,xnNSE (iw,jw,kw), xpNSE (iw,jw,kw)   &
!             ,xHeNSE(iw,jw,kw), xANSE (iw,jw,kw)
!      enddo
!      enddo
!      enddo
!      do iw = iwd,iwu
!      do jw = jwd,jwu
!      do kw = kwd,kwu
!         dqnNSE(iw,jw,kw) = dlog10(dmin1(dmax1(1.d-99, dqnNSE(iw,jw,kw)), 1.d60))
!         dqaNSE(iw,jw,kw) = dlog10(dmin1(dmax1(1.d-99, dqaNSE(iw,jw,kw)), 1.d60))
!         dynNSE(iw,jw,kw) = dlog10(dmin1(dmax1(1.d-99,-dynNSE(iw,jw,kw)), 1.d60))
!         dyaNSE(iw,jw,kw) = dlog10(dmin1(dmax1(1.d-99, dyaNSE(iw,jw,kw)), 1.d60))
!      end do
!      end do
!      end do

!      do i=iprd,ipru
!      do j=jprd,jpru
!         read(31,*)   &
!              t9ipr(i), roeipr(j)  &
!             ,dqPAIR1(i,j), dqPAIR2(i,j) 
!      end do
!      t9ipr(i) = t9ipr(i) - 9.d0
!      end do 
!      do i=ipld,iplu
!      do j=jpld,jplu
!         read(32,*)   &
!              t9ipl(i), roeipl(j)  &
!             ,dqPLAS1(i,j), dqPLAS2(i,j) 
!      end do
!      end do 
!      do i=ibrd,ibru
!      do j=jbrd,jbru
!         read(33,*)   &
!              t9ibr(i), rhoibr(j)  &
!             ,dqBRMS1(i,j), dqBRMS2(i,j) 
!      end do
!      end do 

  return
end subroutine readeos
