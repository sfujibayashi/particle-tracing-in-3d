subroutine set_ejecta_inside_3D_divide(ib,rfl,rin,mass_crit,mass_min,npv)
#include "macro.h"
  use simdata3D
  use particle_data
  use module_eos
  use divide
  use io
  implicit none
!!! functions
  logical :: condition_ejecta, condition_ejecta_geo


  integer,intent(in) :: ib
  real(8),intent(in) :: rfl
  real(8),intent(in) :: mass_crit, mass_min, rin
  integer,intent(inout) :: npv

  ! real(8) :: rfl
  ! real(8) :: mass_crit, mass_min, rin
  ! integer :: npv

  
  ! real(8),parameter :: v_uni = 2.99792458d10
  integer :: ip,j,k,l,lv,lv0,j1,j0,k1,k0,l1,l0,jj,kk,ll,jjd,jju,kkd,kku,lld,llu
  real(8) :: dx,dy,dz,vr,xx,yy,zz,x0,y0,z0,x1,y1,z1
  integer :: i_wing,i_skip
  ! integer,parameter :: i_skip = 17
  ! integer,parameter :: i_wing = i_skip/2
  
  real(8) :: mass,mye,msen,comx,comy,comz,mass_total,mass_total_therm,ye_av,ye_av_therm,sen_av,sen_av_therm, hhh_i,ut_i
  integer :: flag,lv_min_pset
  
  !real(8) :: rin
  real(8) :: hhh_crit, floss, hhh_r, gam_inf_r
  real(8) :: rp_max
  ! real(8),parameter :: v_uni = 2.99792458d10
   ! integer :: ip,j,k,l,j1,j0,k1,k0,l1,l0,jj,kk,ll,jjd,jju,kkd,kku,lld,llu
  ! real(8) :: dx,dy,dz,rfl,vr,xx,yy,zz,x0,y0,z0,x1,y1,z1
  ! integer :: i_wing,i_skip
  
  ! real(8) :: mass,mye,msen,comx,comy,comz,mass_total,mass_total_therm,ye_av,ye_av_therm,sen_av,sen_av_therm,facqb
  ! integer :: flag,np_test, npv

  ! real(8),parameter :: mass_crit=5.d-4 * 1.989d33, mass_min=1.d-12*1.989d33
  ! 
  !integer :: jd_divided4(4),ju_divided4(4),kd_divided4(4),ku_divided4(4),ld_divided4(4),lu_divided4(4)
  integer,allocatable :: jd_divided(:),ju_divided(:),kd_divided(:),ku_divided(:),ld_divided(:),lu_divided(:)
  integer :: nzone,izone
  real(8) :: mass_traj

  integer :: ld_read

  real(8) :: vx_i,vy_i,vz_i,vr_i,ye_i,sen_i

  real(8) :: mass_p_min,mass_p_max,mass_p_tot

  
  real(8) :: dm, vel
  
  integer :: unit_num


  ld_read=ld
#ifdef STAGGERED
#ifndef FULL
  ld_read=ld+1
#endif
#endif
  
  ! hhh_crit = 1.d0 + 7.4785d-3
  hhh_crit = hhh_min
  !rin = 5d7

  hhh_r = 1d0

  ! write(6,*) "test",ld_read,rin,rfl,hhh_crit

  lv=lv_max
  setlevel: do
     ! write(6,'(i5,99es12.4)') lv,x(ju,lv),rfl
     if(x(ju,lv)>rfl.or.lv==lv_min) exit setlevel
     lv=lv-1
  enddo setlevel
  lv_min_pset = lv
  
  if(ib==0)then

     mass_total = 0.d0
     do lv=lv_min_pset,lv_max
        !$omp parallel default(none) &
        !$omp shared(lv,ld,lu,kd,ku,jd,ju,rfl,rin,hhh_crit,vol3D,qb) &
        !$omp private(dm) &
        !$omp reduction(+:mass_total)
        !$omp do
        do l=ld,lu
           do k=kd,ku
              do j=jd,ju
                 if(condition_ejecta(j,k,l,lv,rfl,rin,hhh_crit))then
                    dm = vol3D(j,k,l,lv)*qb(j,k,l,lv)
                    mass_total = mass_total + dm
                 endif
              enddo
           enddo
        enddo
        !$omp end do
        !$omp end parallel
     enddo
     
     write(*,*) "ejecta mass (total integration) =", mass_total
  endif
!  return
#ifdef FULL
  nzone = 8
#else
  nzone = 4
#endif
  allocate(jd_divided(nzone),ju_divided(nzone),kd_divided(nzone),ku_divided(nzone),ld_divided(nzone),lu_divided(nzone))

  ip_ejecta_vol(:,:,:,:) = 0
  mass_traj = 0.d0
  rp_max = 0d0
  ip = 0
  do lv=lv_min_pset,lv_max

#ifdef FULL
     call divide8(jd,ju,kd,ku,ld_read,lu,jd_divided,ju_divided,kd_divided,ku_divided,ld_divided,lu_divided)
#else
     call divide4(jd,ju,kd,ku,ld_read,lu,jd_divided,ju_divided,kd_divided,ku_divided,ld_divided,lu_divided)
#endif
     
     do izone=1,nzone
        ! jjd=jd_divided4(izone)
        ! jju=ju_divided4(izone)
        ! kkd=kd_divided4(izone)
        ! kku=ku_divided4(izone)
        ! lld=ld_divided4(izone)
        ! llu=lu_divided4(izone)
        jjd=jd_divided(izone)
        jju=ju_divided(izone)
        kkd=kd_divided(izone)
        kku=ku_divided(izone)
        lld=ld_divided(izone)
        llu=lu_divided(izone)
        ! write(6,'(99i5)') i4,jjd,jju,kkd,kku,lld,llu

        call recursive_division(ib,lv,jjd,jju,kkd,kku,lld,llu,ip,rfl,rin,hhh_crit,mass_crit,mass_min,mass_traj)
        ! mass_total = 0d0
        ! do l=lld,llu
        !    do k=kkd,kku
        !       do j=jjd,jju
        !          if(condition_ejecta(j,k,l,lv,rfl,rin,hhh_crit))then
        !             mass_total = mass_total + vol3D(j,k,l,lv)*qb(j,k,l,lv)
        !          endif
        !       enddo
        !    enddo
        ! enddo
        ! write(6,*) mass_total
        
        ! if(mass_total > mass_crit)then
        !    call divide8(jjd,jju,kkd,kku,lld,llu,jd_divided8,ju_divided8,kd_divided8,ku_divided8,ld_divided8,lu_divided8)
        !    do i8=1,8
        !       jjd=jd_divided8(i8)
        !       jju=ju_divided8(i8)
        !       kkd=kd_divided8(i8)
        !       kku=ku_divided8(i8)
        !       lld=ld_divided8(i8)
        !       llu=lu_divided8(i8)
        !       write(6,'(99i5)') i8,jjd,jju,kkd,kku,lld,llu
              
        !    enddo
        ! elseif(mass_total > mass_min)then
        !    ip = ip + 1
        ! endif

     enddo
     
  enddo

  deallocate(jd_divided,ju_divided,kd_divided,ku_divided,ld_divided,lu_divided)

  if(ib == 0)then
     npv = ip
     write(*,*) "ejecta mass (sum of particle mass) =", mass_traj
     write(*,*) "missing mass =", mass_total - mass_traj
     write(*,*) "# of particles =", npv
     ! write(*,*) "maximum radius of particle:",rp_max,rfl
  endif

  if(ib == 1)then

     open(101,file=trim(dir_out)//'/p_inside.dat',status='replace')
     write(101,*) "#", npv
     write(101,'(a1,99a14)') "#", &
          "ip","m","ut-1","hut","x","y","z","vx","vy","vz","ye","Gamma_inf(r)"
     do ip=1,npv
        
        flag_evol(ip) = 1
        xx=var_p(index_x,ip)
        yy=var_p(index_y,ip)
        zz=var_p(index_z,ip)
  
        call coorindex3D(xx,yy,zz,j1,k1,l1,lv0)
        j0=j1-1
        k0=k1-1
        l0=l1-1
        
        x1 = (xx-x(j0,lv0))/(x(j1,lv0)-x(j0,lv0))
        x0 = 1.d0-x1
        y1 = (yy-y(k0,lv0))/(y(k1,lv0)-y(k0,lv0))
        y0 = 1.d0-y1
        z1 = (zz-z(l0,lv0))/(z(l1,lv0)-z(l0,lv0))
        z0 = 1.d0-z1
  
        hhh_i = x1*y1*z1* hhh   (j1,k1,l1,lv0) &
              + x0*y1*z1* hhh   (j0,k1,l1,lv0) &
              + x1*y0*z1* hhh   (j1,k0,l1,lv0) &
              + x0*y0*z1* hhh   (j0,k0,l1,lv0) &
              + x1*y1*z0* hhh   (j1,k1,l0,lv0) &
              + x0*y1*z0* hhh   (j0,k1,l0,lv0) &
              + x1*y0*z0* hhh   (j1,k0,l0,lv0) &
              + x0*y0*z0* hhh   (j0,k0,l0,lv0)
        ut_i  = x1*y1*z1* ut    (j1,k1,l1,lv0) &
              + x0*y1*z1* ut    (j0,k1,l1,lv0) &
              + x1*y0*z1* ut    (j1,k0,l1,lv0) &
              + x0*y0*z1* ut    (j0,k0,l1,lv0) &
              + x1*y1*z0* ut    (j1,k1,l0,lv0) &
              + x0*y1*z0* ut    (j0,k1,l0,lv0) &
              + x1*y0*z0* ut    (j1,k0,l0,lv0) &
              + x0*y0*z0* ut    (j0,k0,l0,lv0)
        ye_i  = x1*y1*z1* ye    (j1,k1,l1,lv0) &
              + x0*y1*z1* ye    (j0,k1,l1,lv0) &
              + x1*y0*z1* ye    (j1,k0,l1,lv0) &
              + x0*y0*z1* ye    (j0,k0,l1,lv0) &
              + x1*y1*z0* ye    (j1,k1,l0,lv0) &
              + x0*y1*z0* ye    (j0,k1,l0,lv0) &
              + x1*y0*z0* ye    (j1,k0,l0,lv0) &
              + x0*y0*z0* ye    (j0,k0,l0,lv0)

        vx_i  = x1*y1*z1* vlx    (j1,k1,l1,lv0) &
              + x0*y1*z1* vlx    (j0,k1,l1,lv0) &
              + x1*y0*z1* vlx    (j1,k0,l1,lv0) &
              + x0*y0*z1* vlx    (j0,k0,l1,lv0) &
              + x1*y1*z0* vlx    (j1,k1,l0,lv0) &
              + x0*y1*z0* vlx    (j0,k1,l0,lv0) &
              + x1*y0*z0* vlx    (j1,k0,l0,lv0) &
              + x0*y0*z0* vlx    (j0,k0,l0,lv0)
        vy_i  = x1*y1*z1* vly    (j1,k1,l1,lv0) &
              + x0*y1*z1* vly    (j0,k1,l1,lv0) &
              + x1*y0*z1* vly    (j1,k0,l1,lv0) &
              + x0*y0*z1* vly    (j0,k0,l1,lv0) &
              + x1*y1*z0* vly    (j1,k1,l0,lv0) &
              + x0*y1*z0* vly    (j0,k1,l0,lv0) &
              + x1*y0*z0* vly    (j1,k0,l0,lv0) &
              + x0*y0*z0* vly    (j0,k0,l0,lv0)
        vz_i  = x1*y1*z1* vlz    (j1,k1,l1,lv0) &
              + x0*y1*z1* vlz    (j0,k1,l1,lv0) &
              + x1*y0*z1* vlz    (j1,k0,l1,lv0) &
              + x0*y0*z1* vlz    (j0,k0,l1,lv0) &
              + x1*y1*z0* vlz    (j1,k1,l0,lv0) &
              + x0*y1*z0* vlz    (j0,k1,l0,lv0) &
              + x1*y0*z0* vlz    (j1,k0,l0,lv0) &
              + x0*y0*z0* vlz    (j0,k0,l0,lv0)

        floss = max(0d0,0.0032d0 - 0.0085d0*ye_i)
        
        gam_inf_r = - hhh_i*ut_i/hhh_r * (1d0-floss)
        
        var_p(index_ut1,ip) = ut_i + 1.d0
        var_p(index_hut,ip) = ut_i*hhh_i + hhh_min
        
        write(101,'(a1,i14,99es14.6)') " ",ip, var_p(index_dm,ip), var_p(index_ut1,ip), var_p(index_hut,ip), var_p(index_x,ip), var_p(index_y,ip), var_p(index_z,ip), vx_i, vy_i, vz_i, ye_i, gam_inf_r
     enddo
     
     mass_p_min = 1d99
     mass_p_max = 0d0
     mass_p_tot = 0d0
     do ip=1,npv
        mass_p_tot = mass_p_tot + var_p(index_dm,ip)
        mass_p_min = min(mass_p_min, var_p(index_dm,ip))
        mass_p_max = max(mass_p_max, var_p(index_dm,ip))
     enddo
     write(6,'("Max, min, average mass of particles = ",99es15.7)') mass_p_max,mass_p_min,mass_p_tot/dble(npv)

     write(101,'("# Max, min, average mass of particles = ",99es15.7)') mass_p_max,mass_p_min,mass_p_tot/dble(npv)

     close(101)


!!!  histogram
     block
       !histogram
       real(8),parameter :: dye=0.01d0
       real(8),parameter :: ye_max = 0.60d0+0.5d0*dye, ye_min = 0.01d0 - 0.5*dye
       integer,parameter :: n_ye = nint((ye_max-ye_min)/dye)
       
       real(8),parameter :: dvel=0.01d0
       real(8),parameter :: vel_max = 1.d0, vel_min = 0.0d0
       integer,parameter :: n_vel = nint((vel_max-vel_min)/dvel)
       
       real(8),parameter :: logentr_max = log10(1d3), logentr_min = log10(0.1d0)
       integer,parameter :: n_entr = 100
       real(8),parameter :: dlogentr=(logentr_max-logentr_min)/dble(n_entr)
       
       real(8) :: histogram_ye_entr_total(n_ye,n_entr), histogram_ye_entr_vel_ejecta(n_ye,n_entr,n_vel), histogram_ye_entr_vel_ejectah(n_ye,n_entr,n_vel), histogram_ye_entr_vel_particles(n_ye,n_entr,n_vel), histogram_ye_entr_vel_particlesh(n_ye,n_entr,n_vel)
       real(8) :: hist_v_ye(n_ye), hist_v_entr(n_entr), hist_v_vel(n_vel)
       
       integer :: i_ye, i_vel, i_entr

       histogram_ye_entr_vel_particlesh(:,:,:) = 0d0
       histogram_ye_entr_vel_particles(:,:,:) = 0d0
       histogram_ye_entr_vel_ejectah(:,:,:)=0d0
       histogram_ye_entr_vel_ejecta(:,:,:)=0d0
       histogram_ye_entr_total(:,:)=0d0
       
       do ip=1,npv

          xx=var_p(index_x,ip)
          yy=var_p(index_y,ip)
          zz=var_p(index_z,ip)
    
          call coorindex3D(xx,yy,zz,j1,k1,l1,lv0)
          j0=j1-1
          k0=k1-1
          l0=l1-1
          
          x1 = (xx-x(j0,lv0))/(x(j1,lv0)-x(j0,lv0))
          x0 = 1.d0-x1
          y1 = (yy-y(k0,lv0))/(y(k1,lv0)-y(k0,lv0))
          y0 = 1.d0-y1
          z1 = (zz-z(l0,lv0))/(z(l1,lv0)-z(l0,lv0))
          z0 = 1.d0-z1
    
          hhh_i = x1*y1*z1* hhh   (j1,k1,l1,lv0) &
                + x0*y1*z1* hhh   (j0,k1,l1,lv0) &
                + x1*y0*z1* hhh   (j1,k0,l1,lv0) &
                + x0*y0*z1* hhh   (j0,k0,l1,lv0) &
                + x1*y1*z0* hhh   (j1,k1,l0,lv0) &
                + x0*y1*z0* hhh   (j0,k1,l0,lv0) &
                + x1*y0*z0* hhh   (j1,k0,l0,lv0) &
                + x0*y0*z0* hhh   (j0,k0,l0,lv0)
          ut_i  = x1*y1*z1* ut    (j1,k1,l1,lv0) &
                + x0*y1*z1* ut    (j0,k1,l1,lv0) &
                + x1*y0*z1* ut    (j1,k0,l1,lv0) &
                + x0*y0*z1* ut    (j0,k0,l1,lv0) &
                + x1*y1*z0* ut    (j1,k1,l0,lv0) &
                + x0*y1*z0* ut    (j0,k1,l0,lv0) &
                + x1*y0*z0* ut    (j1,k0,l0,lv0) &
                + x0*y0*z0* ut    (j0,k0,l0,lv0)
          ye_i  = x1*y1*z1* ye    (j1,k1,l1,lv0) &
                + x0*y1*z1* ye    (j0,k1,l1,lv0) &
                + x1*y0*z1* ye    (j1,k0,l1,lv0) &
                + x0*y0*z1* ye    (j0,k0,l1,lv0) &
                + x1*y1*z0* ye    (j1,k1,l0,lv0) &
                + x0*y1*z0* ye    (j0,k1,l0,lv0) &
                + x1*y0*z0* ye    (j1,k0,l0,lv0) &
                + x0*y0*z0* ye    (j0,k0,l0,lv0)
          
          vx_i  = x1*y1*z1* vlx    (j1,k1,l1,lv0) &
                + x0*y1*z1* vlx    (j0,k1,l1,lv0) &
                + x1*y0*z1* vlx    (j1,k0,l1,lv0) &
                + x0*y0*z1* vlx    (j0,k0,l1,lv0) &
                + x1*y1*z0* vlx    (j1,k1,l0,lv0) &
                + x0*y1*z0* vlx    (j0,k1,l0,lv0) &
                + x1*y0*z0* vlx    (j1,k0,l0,lv0) &
                + x0*y0*z0* vlx    (j0,k0,l0,lv0)
          vy_i  = x1*y1*z1* vly    (j1,k1,l1,lv0) &
                + x0*y1*z1* vly    (j0,k1,l1,lv0) &
                + x1*y0*z1* vly    (j1,k0,l1,lv0) &
                + x0*y0*z1* vly    (j0,k0,l1,lv0) &
                + x1*y1*z0* vly    (j1,k1,l0,lv0) &
                + x0*y1*z0* vly    (j0,k1,l0,lv0) &
                + x1*y0*z0* vly    (j1,k0,l0,lv0) &
                + x0*y0*z0* vly    (j0,k0,l0,lv0)
          vz_i  = x1*y1*z1* vlz    (j1,k1,l1,lv0) &
                + x0*y1*z1* vlz    (j0,k1,l1,lv0) &
                + x1*y0*z1* vlz    (j1,k0,l1,lv0) &
                + x0*y0*z1* vlz    (j0,k0,l1,lv0) &
                + x1*y1*z0* vlz    (j1,k1,l0,lv0) &
                + x0*y1*z0* vlz    (j0,k1,l0,lv0) &
                + x1*y0*z0* vlz    (j1,k0,l0,lv0) &
                + x0*y0*z0* vlz    (j0,k0,l0,lv0)
          
          vr_i = (vx_i*xx + vy_i*yy + vz_i*zz)/sqrt(xx**2+yy**2+zz**2)
          
          ye_i  = x1*y1*z1* ye  (j1,k1,l1,lv0) &
                + x0*y1*z1* ye  (j0,k1,l1,lv0) &
                + x1*y0*z1* ye  (j1,k0,l1,lv0) &
                + x0*y0*z1* ye  (j0,k0,l1,lv0) &
                + x1*y1*z0* ye  (j1,k1,l0,lv0) &
                + x0*y1*z0* ye  (j0,k1,l0,lv0) &
                + x1*y0*z0* ye  (j1,k0,l0,lv0) &
                + x0*y0*z0* ye  (j0,k0,l0,lv0)
          sen_i = x1*y1*z1* sen  (j1,k1,l1,lv0) &
                + x0*y1*z1* sen  (j0,k1,l1,lv0) &
                + x1*y0*z1* sen  (j1,k0,l1,lv0) &
                + x0*y0*z1* sen  (j0,k0,l1,lv0) &
                + x1*y1*z0* sen  (j1,k1,l0,lv0) &
                + x0*y1*z0* sen  (j0,k1,l0,lv0) &
                + x1*y0*z0* sen  (j1,k0,l0,lv0) &
                + x0*y0*z0* sen  (j0,k0,l0,lv0)
          
          i_ye = max(1,min(n_ye  ,int((ye_i-ye_min)/dye)+1))
          i_entr=max(1,min(n_entr,int((log10(sen_i)-logentr_min)/dlogentr)+1))
          
          if( ut_i*hhh_i + hhh_min < 0d0)then
             vel = sqrt(1d0-1d0/(-ut_i*hhh_i/hhh_min)**2)        
             i_vel= max(1,min(n_vel,int((vel-vel_min)/dvel)+1))
             histogram_ye_entr_vel_particlesh(i_ye,i_entr,i_vel) = histogram_ye_entr_vel_particlesh(i_ye,i_entr,i_vel) + var_p(index_dm,ip)
          endif
          
          if( ut_i+1d0<0d0)then
             vel = sqrt(1d0-1d0/(-ut_i)**2)
             i_vel= max(1,min(n_vel,int((vel-vel_min)/dvel)+1))
             histogram_ye_entr_vel_particles(i_ye,i_entr,i_vel) = histogram_ye_entr_vel_particles(i_ye,i_entr,i_vel) + var_p(index_dm,ip)
          endif
          
       enddo

       do lv=lv_min_pset,lv_max
!          !$omp parallel default(none) &
!          !$omp num_threads(1) &
!          !$omp shared(lv,ld,lu,kd,ku,jd,ju,rfl,rin,hhh_crit,hhh_min,vol3D,qb,ye,ut,hhh,sen) &
!          !$omp private(dm,vel,i_ye,i_entr,i_vel) &
!          !$omp reduction(+: histogram_ye_entr_total, histogram_ye_entr_vel_ejectah, histogram_ye_entr_vel_ejecta)
!          !$omp do
          do l=ld,lu
             do k=kd,ku
                do j=jd,ju
                   dm = vol3D(j,k,l,lv)*qb(j,k,l,lv)
                   i_ye = max(1, min(n_ye  , int((ye(j,k,l,lv)-ye_min)/dye)+1))
                   i_entr=max(1, min(n_entr, int((log10(sen(j,k,l,lv))-logentr_min)/dlogentr)+1))
                   histogram_ye_entr_total(i_ye,i_entr) = histogram_ye_entr_total(i_ye,i_entr) + dm

                   if(condition_ejecta(j,k,l,lv,rfl,rin,hhh_crit))then
                      vel = sqrt(1d0-1d0/(-ut(j,k,l,lv)*hhh(j,k,l,lv)/hhh_min)**2)
                      i_vel= max(1,min(n_vel,int((vel-vel_min)/dvel)+1))
                      histogram_ye_entr_vel_ejectah(i_ye,i_entr,i_vel) = histogram_ye_entr_vel_ejectah(i_ye,i_entr,i_vel) + dm
                   endif

                   if(condition_ejecta_geo(j,k,l,lv,rfl,rin))then
                      vel = sqrt(1d0-1d0/(-ut(j,k,l,lv))**2)
                      i_vel= max(1,min(n_vel,int((vel-vel_min)/dvel)+1))
                      histogram_ye_entr_vel_ejecta(i_ye,i_entr,i_vel) = histogram_ye_entr_vel_ejecta(i_ye,i_entr,i_vel) + dm
                   endif

                enddo
             enddo
          enddo
!          !$omp end do
!          !$omp end parallel
       enddo


       do i_ye=1,n_ye
          hist_v_ye(i_ye) = ye_min + dye*dble(i_ye-1)
          !write(6,*) i_ye, hist_v_ye(i_ye)
       end do

       do i_vel=1,n_vel
          hist_v_vel(i_vel) = vel_min + dvel*dble(i_vel-1)
          !write(6,*) i_vel, hist_v_vel(i_vel)
       end do

       do i_entr=1,n_entr
          hist_v_entr(i_entr) = 10d0**(logentr_min + dlogentr*dble(i_entr-1))
          !write(6,*) i_entr, hist_v_entr(i_entr)
       end do

       open(newunit=unit_num,file=trim(dir_out)//"/ye_hist_inside.dat",status="replace")
       write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       write(unit_num,'(a1,99a15)') "#", "Ye", "particle", "particles(geo)", "ejecta", "ejecta(geo)", "total"
       do i_ye = 1,n_ye
          write(unit_num,'(a1,99es15.7)') " ",hist_v_ye(i_ye), sum(histogram_ye_entr_vel_particlesh(i_ye,:,:)), sum(histogram_ye_entr_vel_particles(i_ye,:,:)), sum(histogram_ye_entr_vel_ejectah(i_ye,:,:)), sum(histogram_ye_entr_vel_ejecta(i_ye,:,:)), sum(histogram_ye_entr_total(i_ye,:))
       enddo
       close(unit_num)

       open(newunit=unit_num,file=trim(dir_out)//"/sen_hist_inside.dat",status="replace")
       write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       write(unit_num,'(a1,99a15)') "#", "entropy", "particle", "particles(geo)", "ejecta", "ejecta(geo)", "total"
       do i_entr = 1,n_entr
          write(unit_num,'(a1,99es15.7)') " ",hist_v_entr(i_entr), sum(histogram_ye_entr_vel_particlesh(:,i_entr,:)), sum(histogram_ye_entr_vel_particles(:,i_entr,:)), sum(histogram_ye_entr_vel_ejectah(:,i_entr,:)), sum(histogram_ye_entr_vel_ejecta(:,i_entr,:)), sum(histogram_ye_entr_total(:,i_entr))
       enddo
       close(unit_num)

       open(newunit=unit_num,file=trim(dir_out)//"/vel_hist_inside.dat",status="replace")
       write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       write(unit_num,'(a1,99a15)') "#", "vel_inf", "particle", "particles(geo)", "ejecta", "ejecta(geo)", "total"
       do i_vel = 1,n_vel
          write(unit_num,'(a1,99es15.7)') " ",hist_v_vel(i_vel), sum(histogram_ye_entr_vel_particlesh(:,:,i_vel)), sum(histogram_ye_entr_vel_particles(:,:,i_vel)), sum(histogram_ye_entr_vel_ejectah(:,:,i_vel)), sum(histogram_ye_entr_vel_ejecta(:,:,i_vel))
       enddo
       close(unit_num)

       open(newunit=unit_num,file=trim(dir_out)//"/ye_sen_hist_inside.dat",status="replace")
       write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       write(unit_num,'(a1,99a15)') "#", "Ye", "entropy", "particle", "particles(geo)", "ejecta", "ejecta(geo)", "total"
       do i_entr = 1,n_entr
          write(unit_num,*)
          do i_ye = 1,n_ye
             write(unit_num,'(a1,99es15.7)') " ",hist_v_ye(i_ye), hist_v_entr(i_entr), sum(histogram_ye_entr_vel_particlesh(i_ye,i_entr,:)), sum(histogram_ye_entr_vel_particles(i_ye,i_entr,:)), sum(histogram_ye_entr_vel_ejectah(i_ye,i_entr,:)), sum(histogram_ye_entr_vel_ejecta(i_ye,i_entr,:)), histogram_ye_entr_total(i_ye,i_entr)
          enddo
       enddo
       close(unit_num)

       open(newunit=unit_num,file=trim(dir_out)//"/ye_vel_hist_inside.dat",status="replace")
       write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       write(unit_num,'(a1,99a15)') "#", "Ye", "vel_inf", "particle", "particles(geo)", "ejecta", "ejecta(geo)"
       do i_vel = 1,n_vel
          write(unit_num,*)
          do i_ye = 1,n_ye
             write(unit_num,'(a1,99es15.7)') " ",hist_v_ye(i_ye), hist_v_vel(i_vel), sum(histogram_ye_entr_vel_particlesh(i_ye,:,i_vel)), sum(histogram_ye_entr_vel_particles(i_ye,:,i_vel)), sum(histogram_ye_entr_vel_ejectah(i_ye,:,i_vel)), sum(histogram_ye_entr_vel_ejecta(i_ye,:,i_vel))
          enddo
       enddo
       close(unit_num)

       open(newunit=unit_num,file=trim(dir_out)//"/sen_vel_hist_inside.dat",status="replace")
       write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       write(unit_num,'(a1,99a15)') "#", "Ye", "entropy", "particle", "particles(geo)", "ejecta", "ejecta(geo)"
       do i_vel = 1,n_vel
          write(unit_num,*)
          do i_entr = 1,n_entr
             write(unit_num,'(a1,99es15.7)') " ",hist_v_entr(i_entr), hist_v_vel(i_vel), sum(histogram_ye_entr_vel_particlesh(:,i_entr,i_vel)), sum(histogram_ye_entr_vel_particles(:,i_entr,i_vel)), sum(histogram_ye_entr_vel_ejectah(:,i_entr,i_vel)), sum(histogram_ye_entr_vel_ejecta(:,i_entr,i_vel))
          enddo
       enddo
       close(unit_num)

       ! open(newunit=unit_num,file=trim(dir_out)//"/ye_sen_vel_hist_inside.dat",status="replace")
       ! write(unit_num,'(a1,15x,99es15.7)') "#",sum(histogram_ye_entr_vel_particlesh(:,:,:)), sum(histogram_ye_entr_vel_particles(:,:,:)), sum(histogram_ye_entr_vel_ejectah(:,:,:)), sum(histogram_ye_entr_vel_ejecta(:,:,:)),sum(histogram_ye_entr_total(:,:))
       ! write(unit_num,'(a1,99a15)') "#", "Ye", "entropy", "vel_inf", "particle", "particles(geo)", "ejecta", "ejecta(geo)"
       ! do i_vel = 1,n_vel
       !    write(unit_num,*)
       !    do i_entr = 1,n_entr
       !       write(unit_num,*)
       !       do i_ye = 1,n_ye
       !          write(unit_num,'(a1,99es15.7)') " ",hist_v_ye(i_ye), hist_v_entr(i_entr), hist_v_vel(i_vel), sum(histogram_ye_entr_vel_particlesh(i_ye,:,i_vel)), sum(histogram_ye_entr_vel_particles(i_ye,:,i_vel)), sum(histogram_ye_entr_vel_ejectah(i_ye,:,i_vel)), sum(histogram_ye_entr_vel_ejecta(i_ye,:,i_vel))
       !       enddo
       !    enddo
       ! enddo
       ! close(unit_num)

     end block
  endif

  return
end subroutine set_ejecta_inside_3D_divide
