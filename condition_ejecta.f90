logical function condition_ejecta(j,k,l,lv,rfl,rin,hhh_crit)
  use simdata3D
  use particle_data
  
  implicit none
  integer,intent(in) :: j,k,l,lv
  real(8),intent(in) :: rfl,rin,hhh_crit

  condition_ejecta = .false.

  if(-ut(j,k,l,lv)*hhh(j,k,l,lv)-hhh_crit> 0.d0  &
       .and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
       .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
       .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>rin )then
     condition_ejecta=.true.
  endif

  return
end function condition_ejecta

logical function condition_ejecta_geo(j,k,l,lv,rfl,rin)
  use simdata3D
  use particle_data
  
  implicit none
  integer,intent(in) :: j,k,l,lv
  real(8),intent(in) :: rfl,rin

  condition_ejecta_geo = .false.

  if(ut(j,k,l,lv)+1d0 < 0.d0  &
       .and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
       .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
       .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>rin )then
     condition_ejecta_geo =.true.
  endif
  
  return
end function condition_ejecta_geo
