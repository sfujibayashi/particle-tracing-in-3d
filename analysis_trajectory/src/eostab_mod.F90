!------------------------------------------------------!
!                  EOS  Table
!------------------------------------------------------!
module eostab
#include "macro.h"
  implicit real(8) (a-h,o-z)
  ! --- Nuc. EOS + Timmes EOS
  parameter(ied =1, ieu =131 )
  parameter(jed =1, jeu =61 )
#ifdef DD2
  parameter(ked =1, keu =426)   ! --- DD2+Timmes
#endif
#ifdef SFHo
  parameter(ked =1, keu =408)  ! --- SFHo+Timmes
#endif
  
!   !---- DD2 only
!   parameter(ied =1, ieu =81 )
!   parameter(jed =1, jeu =61 )
!   parameter(ked =1, keu =326)
  ! --- for SFHo
  ! parameter(ked =1, keu =308)

  real(8) rho_e_min,tem_e_min
  parameter( ye_e_min = 0.00d0)
  parameter(drho_e =  0.04d0, drhoi=1.d0/drho_e)
  parameter( dye_e =  0.01d0, dyei =1.d0/ dye_e)
  parameter(dtem_e =  0.04d0, dtemi=1.d0/dtem_e)
  ! --- EOS for Stiner : T>=0.1MeV  NO-BETA-eq
  real(8)  rho_e(ked:keu),  ye_e(jed:jeu), tem_e (ied:ieu)                  &
        ,pres_e (ied:ieu,jed:jeu,ked:keu), eps_e (ied:ieu,jed:jeu,ked:keu)  &
        ,sen_e  (ied:ieu,jed:jeu,ked:keu), cs_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,chn_e  (ied:ieu,jed:jeu,ked:keu), chp_e (ied:ieu,jed:jeu,ked:keu), che_e (ied:ieu,jed:jeu,ked:keu)  &
        ,aa_e   (ied:ieu,jed:jeu,ked:keu), zz_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,xn_e   (ied:ieu,jed:jeu,ked:keu), xp_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,xA_e   (ied:ieu,jed:jeu,ked:keu), xy_e  (ied:ieu,jed:jeu,ked:keu)  &
        ,xd_e   (ied:ieu,jed:jeu,ked:keu), xt_e  (ied:ieu,jed:jeu,ked:keu), xh_e  (ied:ieu,jed:jeu,ked:keu) 
  ! --- EOS for Steiner : T>=0.1MeV  BETA-eq
  !   real(8) &
  !   parameter(yl_eb_min = )
  !   parameter(yl_eb_max = )
  !   parameter(dyl_eb=, dyli=1.d0/dyl_eb)
  real(8)  yl_eb(jed:jeu)                , ye_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,yn_eb  (ied:ieu,jed:jeu,ked:keu), ya_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,pres_eb(ied:ieu,jed:jeu,ked:keu), eps_eb(ied:ieu,jed:jeu,ked:keu)  &
        ,sen_eb (ied:ieu,jed:jeu,ked:keu), cs_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,chn_eb (ied:ieu,jed:jeu,ked:keu), chp_eb(ied:ieu,jed:jeu,ked:keu), che_eb(ied:ieu,jed:jeu,ked:keu)  &
        ,aa_eb  (ied:ieu,jed:jeu,ked:keu), zz_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,xn_eb  (ied:ieu,jed:jeu,ked:keu), xp_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,xA_eb  (ied:ieu,jed:jeu,ked:keu), xy_eb (ied:ieu,jed:jeu,ked:keu)  &
        ,xd_eb  (ied:ieu,jed:jeu,ked:keu), xt_eb (ied:ieu,jed:jeu,ked:keu), xh_eb (ied:ieu,jed:jeu,ked:keu) 
  real(8) ye_0b(ked:keu)
  
  ! --- Nuc. EOS region
  parameter(ierd =51 , ieru =ieu )
  parameter(jerd =jed, jeru =jeu )
  parameter(kerd =101, keru =keu )
!   parameter(ierd =1  , ieru =ieu )
!   parameter(jerd =jed, jeru =jeu )
!   parameter(kerd =1  , keru =keu )

  ! --- above this density only Nuc. EOS is used in detehat,
  ! --- i.e., the temperature search is performed above 0.1 MeV.
  parameter(knuc = 326)
  ! --- if Timmes EOS is not used, set knuc=0 and ierd=jerd=kerd=1.

  ! --- 1 for the region where beta-EOS is available, and 0 otherwise.
  integer :: iflag_beos(ied:ieu,jed:jeu,ked:keu)

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
end module eostab
