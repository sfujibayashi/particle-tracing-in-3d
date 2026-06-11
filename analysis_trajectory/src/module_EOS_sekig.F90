module module_EOS_sekig
#include "macro.h"
  implicit none

  integer :: ied,ieu,ked,keu,jed,jeu
#ifdef Timmes
  ! --- Nuc. EOS + Timmes EOS
  ! -- Full Range
  parameter(ied =1, ieu =131 )
  parameter(jed =1, jeu =61 )
#ifdef SFHo
  parameter(ked =1, keu =408) ! --- SFHo+Timmes
#else
  parameter(ked =1, keu =426) ! --- DD2+Timmes
#endif

#else
  ! -- Full Range
  parameter(ied =1, ieu =81 )
  parameter(jed =1, jeu =61 )
#ifdef SFHo
  parameter(ked =1, keu =308) ! --- SFHo+Timmes
#else
  parameter(ked =1, keu =326) ! --- DD2+Timmes
#endif
  
#endif

  real(8) :: rho_e_min,tem_e_min,ye_e_min,drho_e,dye_e,dtem_e,drhoi,dyei,dtemi
  parameter(drho_e =  0.04d0, drhoi=1.d0/drho_e)
  parameter( dye_e =  0.01d0, dyei =1.d0/ dye_e)
  parameter(dtem_e =  0.04d0, dtemi=1.d0/dtem_e)

  ! --- EOS for Stiner : T>=0.1MeV  NO-BETA-eq
  real(8) :: &
         rho_e(ked:keu),  ye_e(jed:jeu), tem_e (ied:ieu)                  &
        ,pres_e (ied:ieu,jed:jeu,ked:keu), eps_e (ied:ieu,jed:jeu,ked:keu)  &
        ,sen_e  (ied:ieu,jed:jeu,ked:keu), cs_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,chn_e  (ied:ieu,jed:jeu,ked:keu), chp_e (ied:ieu,jed:jeu,ked:keu), che_e (ied:ieu,jed:jeu,ked:keu)  &
        ,aa_e   (ied:ieu,jed:jeu,ked:keu), zz_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,xn_e   (ied:ieu,jed:jeu,ked:keu), xp_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,xA_e   (ied:ieu,jed:jeu,ked:keu), xy_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,xd_e   (ied:ieu,jed:jeu,ked:keu), xt_e  (ied:ieu,jed:jeu,ked:keu), xh_e  (ied:ieu,jed:jeu,ked:keu) 
  ! --- EOS for Steiner : T>=0.1MeV  BETA-eq
  real(8) :: &
         yl_eb(jed:jeu)                , ye_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,yn_eb  (ied:ieu,jed:jeu,ked:keu), ya_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,pres_eb(ied:ieu,jed:jeu,ked:keu), eps_eb(ied:ieu,jed:jeu,ked:keu)  &
        ,sen_eb (ied:ieu,jed:jeu,ked:keu), cs_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,chn_eb (ied:ieu,jed:jeu,ked:keu), chp_eb(ied:ieu,jed:jeu,ked:keu), che_eb(ied:ieu,jed:jeu,ked:keu)  &
        ,aa_eb  (ied:ieu,jed:jeu,ked:keu), zz_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,xn_eb  (ied:ieu,jed:jeu,ked:keu), xp_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,xA_eb  (ied:ieu,jed:jeu,ked:keu), xy_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,xd_eb  (ied:ieu,jed:jeu,ked:keu), xt_eb (ied:ieu,jed:jeu,ked:keu), xh_eb (ied:ieu,jed:jeu,ked:keu) 

#ifdef Timmes
  integer ::ierd,ieru,jerd,jeru,kerd,keru,knuc
  ! --- High-density EOS region
  parameter(ierd =51    , ieru =ieu )
  parameter(jerd =jed, jeru =jeu )
  parameter(kerd =101   , keru =keu )

  ! --- above this density only Nuc. EOS is used in detehat,
  ! --- i.e., the temperature search is performed above 0.1 MeV.
  ! parameter(knuc = 276)
  parameter(knuc = 326)
  ! --- if Timmes EOS is not used, set knuc=0 and ierd=jerd=kerd=1.

  ! --- 1 for the region where beta-EOS is available, and 0 otherwise.
  integer :: iflag_beos(ied:ieu,jed:jeu,ked:keu)
#endif

  real(8) :: escn0,dscn,dscni,escn_min,esce0, dsce,dscei,esce_min
  integer :: inut,ient
  ! --- Range neutrino data
  parameter(escn0=-50.d0, dscn=0.003d0)
  parameter(inut=15001)
  parameter(dscni =1.d0/dscn, escn_min=1.d-50)
  real(8) escn(1:inut), ech_nn(1:inut), eprnn(1:inut)

  parameter(esce0=-47.d0, dsce=0.002d0)
  parameter(ient=20001)
  parameter(dscei = 1.d0/dsce, esce_min=1.d-47)
  real(8) esce(1:ient), ech_ne(1:ient), eprne(1:ient)  &
         ,ynue(1:ient), eave  (1:ient)

contains
  
  subroutine readeos(fn_eos,fn_eosb,fn_ynu,fn_enu)
    use unit

    character(*) :: fn_eos,fn_eosb,fn_ynu,fn_enu
    real(8) :: ePmin, eEmin, eCmin, ePnmin, eEnmin, eCnmin
    integer :: ie,je,ke

    integer :: inu

    open(22,file=trim(adjustl(fn_ynu)) ,status='old',form='unformatted',access='stream')
    read(22) escn, ech_nn, eprnn
    close(22)

    do inu=1,inut
       escn(inu)=escn(inu)
    enddo
    
    open(23,file=trim(adjustl(fn_enu)) ,status='old',form='unformatted',access='stream')
    read(23) esce, ech_ne, eprne, ynue, eave
    close(23)
    
    open(20,file=trim(adjustl(fn_eos)) ,status='old',form='unformatted',access='stream')
    read(20) tem_e, ye_e, rho_e         &
         ,pres_e, eps_e, sen_e, cs_e &
         ,chn_e , chp_e, che_e       &
         ,aa_e  , zz_e               &
         ,xn_e  , xp_e , xA_e        &
         ,xd_e  , xt_e , xh_e , xy_e &
         ,ePmin , eEmin, eCmin
    close(20)

    
    open(21,file=trim(adjustl(fn_eosb)),status='old',form='unformatted',access='stream')
    read(21)yl_eb                        &
         ,pres_eb, eps_eb, sen_eb, cs_eb &
         ,chn_eb , chp_eb, che_eb        &
         ,aa_eb  , zz_eb                 &
         ,xn_eb  , xp_eb , xA_eb         &
         ,xd_eb  , xt_eb , xh_eb , xy_eb &
         ,ye_eb  , yn_eb , ya_eb         &
         ,ePnmin , eEnmin, eCnmin
    close(21)

    do ke=ked,keu
       rho_e(ke) = log10(rho_e(ke))
    enddo
    do ie=ied,ieu
       tem_e(ie) = log10(tem_e(ie))
    enddo

    pres_e(:,:,:) = log10(pres_e(:,:,:))
    eps_e(:,:,:)  = log10(eps_e(:,:,:) + v_uni**2)
    pres_eb(:,:,:) = log10(pres_eb(:,:,:))
    eps_eb(:,:,:)  = log10(eps_eb(:,:,:) + v_uni**2)

    rho_e_min=rho_e(ked)
    tem_e_min=tem_e(ied)
    ye_e_min = 0.d0

  end subroutine readeos

  subroutine interp_coef(rho,tem,ye,irho,item,iye,uu,ss,tt)
    real(8),intent(in)  :: rho,tem,ye
    real(8),intent(out) :: uu,ss,tt
    integer,intent(out) :: irho,item,iye
    
    item = max(ied , min(ieu-1, int((log10(tem)-tem_e_min  )*dtemi)+1))
    ss   = max(0.d0, min(1.d0 ,     (log10(tem)-tem_e(item))*dtemi))
    
    irho = max(ked , min(keu-1, int((log10(rho)-rho_e_min  )*drhoi)+1))
    uu   = max(0.d0, min(1.d0,      (log10(rho)-rho_e(irho))*drhoi))
    
    iye  = max(jed , min(jeu-1, int((ye-ye_e_min  )*dyei)+1))
    tt   = max(0.d0, min(1.d0 ,     (ye-ye_e(iye) )*dyei))
    
  end subroutine interp_coef

  subroutine interp_val(rho,tem,ye,table_val,val)
    
    real(8),intent(in)  :: table_val(ied:ieu,jed:jeu,ked:keu)
    real(8),intent(in)  :: rho,tem,ye
    real(8),intent(out) :: val
    integer :: irho,item,iye
    integer :: irho1,item1,iye1
    real(8) :: uu,ss,tt,uu1,ss1,tt1

    item = max(ied , min(ieu-1, int((log10(tem)-tem_e_min  )*dtemi)+1))
    ss   = max(0.d0, min(1.d0 ,     (log10(tem)-tem_e(item))*dtemi))
    
    irho = max(ked , min(keu-1, int((log10(rho)-rho_e_min  )*drhoi)+1))
    uu   = max(0.d0, min(1.d0,      (log10(rho)-rho_e(irho))*drhoi))
    
    iye  = max(jed , min(jeu-1, int((ye-ye_e_min  )*dyei)+1))
    tt   = max(0.d0, min(1.d0 ,     (ye-ye_e(iye) )*dyei))

    item1=item+1
    irho1=irho+1
    iye1=iye+1
    ss1=1d0-ss
    tt1=1d0-tt
    uu1=1d0-uu

    val =  ss1 *tt1 *uu1 *table_val(item ,iye ,irho )   &
         + ss  *tt1 *uu1 *table_val(item1,iye ,irho )   &
         + ss1 *tt  *uu1 *table_val(item ,iye1,irho )   &
         + ss1 *tt1 *uu  *table_val(item ,iye ,irho1)   &
         + ss  *tt  *uu1 *table_val(item1,iye1,irho )   &
         + ss  *tt1 *uu  *table_val(item1,iye ,irho1)   &
         + ss1 *tt  *uu  *table_val(item ,iye1,irho1)   &
         + ss  *tt  *uu  *table_val(item1,iye1,irho1)
    
  end subroutine interp_val
  
end module module_EOS_sekig
