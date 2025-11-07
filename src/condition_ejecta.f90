module ejecta_mod
  implicit none
contains

  logical function condition_ejecta(j,k,l,lv,rfl,rin,hhh_crit)
    use simdata3D
    use particle_data

    implicit none
    integer,intent(in) :: j,k,l,lv
    real(8),intent(in) :: rfl,rin,hhh_crit

    ! logical :: condition_ejecta_geo, condition_ejecta_bernoulli
    ! condition_ejecta = .false.

    ! if(-ut(j,k,l,lv)*hhh(j,k,l,lv)-hhh_crit> 0.d0  &
    !      .and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
    !      .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
    !      !.and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>rin &
    !      )then
    !    condition_ejecta=.true.
    ! endif

    !condition_ejecta = condition_ejecta_geo(j,k,l,lv,rfl,rin)
    condition_ejecta = condition_ejecta_bernoulli(j,k,l,lv,rfl,rin,hhh_crit)
    !condition_ejecta = condition_bound_bernoulli(j,k,l,lv,rfl,rin,hhh_crit)
    
    return
  end function condition_ejecta


  logical function condition_ejecta_bernoulli(j,k,l,lv,rfl,rin,hhh_crit)
    use simdata3D
    use particle_data

    implicit none
    integer,intent(in) :: j,k,l,lv
    real(8),intent(in) :: rfl,rin,hhh_crit

    real(8) :: r2

    condition_ejecta_bernoulli = .false.
    r2 = x(j,lv)**2+y(k,lv)**2+z(l,lv)**2
    if(-ut(j,k,l,lv)*hhh(j,k,l,lv)-hhh_crit> 0.d0  &
         ! .and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
         .and. r2<rfl**2 &
         .and. r2>=rin**2 &
         )then
       condition_ejecta_bernoulli = .true.
    endif

    return
  end function condition_ejecta_bernoulli


  logical function condition_ejecta_hut1(j,k,l,lv,rfl,rin,hhh_crit)
    use simdata3D
    use particle_data

    implicit none
    integer,intent(in) :: j,k,l,lv
    real(8),intent(in) :: rfl,rin,hhh_crit

    condition_ejecta_hut1 = .false.

    if(-ut(j,k,l,lv)*hhh(j,k,l,lv)-1d0> 0.d0  &
                                ! .and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
         .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
         .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>=rin &
         )then
       condition_ejecta_hut1 = .true.
    endif

    return
  end function condition_ejecta_hut1


  logical function condition_ejecta_bernoulli_positive_velocity(j,k,l,lv,rfl,rin,hhh_crit)
    use simdata3D
    use particle_data

    implicit none
    integer,intent(in) :: j,k,l,lv
    real(8),intent(in) :: rfl,rin,hhh_crit

    condition_ejecta_bernoulli_positive_velocity = .false.

    if(-ut(j,k,l,lv)*hhh(j,k,l,lv)-hhh_crit> 0.d0  &
         .and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
         .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
         .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>=rin &
         )then
       condition_ejecta_bernoulli_positive_velocity = .true.
    endif

    return
  end function condition_ejecta_bernoulli_positive_velocity

  logical function condition_ejecta_geo(j,k,l,lv,rfl,rin)
    use simdata3D
    use particle_data

    implicit none
    integer,intent(in) :: j,k,l,lv
    real(8),intent(in) :: rfl,rin

    condition_ejecta_geo = .false.

    if(ut(j,k,l,lv)+1d0 < 0.d0  &
                                !.and. vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv) > 0.d0 &
         .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)<rfl &
         .and. sqrt(x(j,lv)**2+y(k,lv)**2+z(l,lv)**2)>=rin &
         )then
       condition_ejecta_geo =.true.
    endif

    return
  end function condition_ejecta_geo

  logical function condition_bound_bernoulli(j,k,l,lv,rfl,rin,hhh_crit)
    use simdata3D
    use particle_data

    implicit none
    integer,intent(in) :: j,k,l,lv
    real(8),intent(in) :: rfl,rin,hhh_crit

    real(8) :: r2

    condition_bound_bernoulli = .false.
    r2 = x(j,lv)**2+y(k,lv)**2+z(l,lv)**2
    if(-ut(j,k,l,lv)*hhh(j,k,l,lv)-hhh_crit < 0.d0  &
         .and. r2<rfl**2 &
         .and. r2>=rin**2 &
         .and. qrho(j,k,l,lv) < 1d14 &
         )then
       condition_bound_bernoulli = .true.
    endif

    return
  end function condition_bound_bernoulli


end module ejecta_mod
