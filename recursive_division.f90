recursive subroutine recursive_division(ib,lv,jjd,jju,kkd,kku,lld,llu,ip,rfl,rin,hhh_crit,mass_crit,mass_min,mass_traj)
  use simdata3D
  use particle_data
  use divide
  implicit none
  integer,intent(in) :: lv,ib
  integer,intent(inout) :: jjd,jju,kkd,kku,lld,llu,ip
  real(8),intent(in) :: rfl,rin,hhh_crit,mass_crit,mass_min
  real(8),intent(inout) :: mass_traj

  integer :: j,k,l
  logical :: condition_ejecta

  integer :: jd_divided8(8),ju_divided8(8),kd_divided8(8),ku_divided8(8),ld_divided8(8),lu_divided8(8)
  integer :: i8

  real(8) :: mass,comx,comy,comz,vol,vr2_av,vr_av,vr,sigmav,rho_av

  mass = 0d0
  comx = 0.d0
  comy = 0.d0
  comz = 0.d0
  vol  = 0.d0
  
  vr_av = 0.d0
  vr2_av= 0.d0
  
  do l=lld,llu
     do k=kkd,kku
        do j=jjd,jju
           if(condition_ejecta(j,k,l,lv,rfl,rin,hhh_crit))then
              mass = mass + vol3D(j,k,l,lv)*qb(j,k,l,lv)
              comx = comx + vol3D(j,k,l,lv)*qb(j,k,l,lv)*x(j,lv)
              comy = comy + vol3D(j,k,l,lv)*qb(j,k,l,lv)*y(k,lv)
              comz = comz + vol3D(j,k,l,lv)*qb(j,k,l,lv)*z(l,lv)

              vr = (vlx(j,k,l,lv)*x(j,lv) + vly(j,k,l,lv)*y(k,lv) + vlz(j,k,l,lv)*z(l,lv))/sqrt(x(j,lv)**2 + y(k,lv)**2 + z(l,lv)**2 + 1.d0)
              
              vr_av  = vr_av  + vr   *vol3D(j,k,l,lv)*qb(j,k,l,lv)
              vr2_av = vr2_av + vr**2*vol3D(j,k,l,lv)*qb(j,k,l,lv)
              vol = vol + vol3D(j,k,l,lv)
           endif
        enddo
     enddo
  enddo

  if(mass > 0.d0)then
     vr_av  = vr_av /mass
     vr2_av = vr2_av/mass
     sigmav = sqrt( vr2_av - vr_av**2 )
     rho_av = mass/vol
  else
     vr_av  = 0.d0
     vr2_av = 0.d0
     sigmav = 0.d0
     rho_av = 0.d0
  endif
  ! write(6,'(6i5,99es12.4)') jjd,jju,kkd,kku,lld,llu,mass

  if( ( mass > mass_crit ) &
       .and. jju-jjd+1>=2 .and. kku-kkd+1>=2 .and. llu-lld+1>=2)then
     call divide8(jjd,jju,kkd,kku,lld,llu,jd_divided8,ju_divided8,kd_divided8,ku_divided8,ld_divided8,lu_divided8)
     do i8=1,8
        jjd=jd_divided8(i8)
        jju=ju_divided8(i8)
        kkd=kd_divided8(i8)
        kku=ku_divided8(i8)
        lld=ld_divided8(i8)
        llu=lu_divided8(i8)
        call recursive_division(ib,lv,jjd,jju,kkd,kku,lld,llu,ip,rfl,rin,hhh_crit,mass_crit,mass_min,mass_traj)
     enddo
  elseif(mass > mass_min.and.rho_av > 1d0)then
  !elseif(mass > mass_min)then
     ip = ip + 1
     mass_traj = mass_traj + mass
     
     if(ib==1)then
        x_p(ip) = comx/mass
        y_p(ip) = comy/mass
        z_p(ip) = comz/mass
        dm_p(ip)= mass
     endif
     
  endif

end subroutine recursive_division
