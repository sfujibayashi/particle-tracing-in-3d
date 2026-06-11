module eos00
   implicit none
   real(8)  rho00,  tem01, tem00, yemin, yemax
   real(8)  rho_at,   ye_at,   yn_at,  ya_at,   yo_at,  tem_at  &
          , eps_at, pres_at,   cs_at, hhh_at, ehat_at, eps2_at  &
          ,ch_n_at, ch_p_at, ch_e_at, x_n_at,  x_p_at,  x_A_at  &
          ,  aa_at,   zz_at,  sen_at, x_L_at
   real(8)  rho_nb, rho_tt, rrr_at, rho_at1, rhocut
   integer  krho_at, jye_at, item_at, krho0, jye0, item0
end module eos00
