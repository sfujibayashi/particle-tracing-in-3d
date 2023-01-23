subroutine tr_analysis(model, dir_read)
  !$ use omp_lib
  implicit none

  character(*),intent(in) :: model,dir_read
  integer :: np,itt_max,itt_min  
  
  real(8),parameter :: tem_uni = 1.160445d10
  real(8),parameter :: v_uni = 2.99792458d10

  integer :: access
  integer :: nunit, punit

  real(8),allocatable :: time(:)
  real(8),allocatable :: x_p(:),y_p(:),z_p(:),&
       qrho_p(:),&
       ye_p  (:),&
       tem_p (:),&
       ut_p  (:),&
       qb_p  (:),&
       sen_p (:),&
       vlx_p (:),&
       vly_p (:),&
       vlz_p (:),&
       hhh_p (:),&
       dt_p  (:),&
       rne_p (:),&
       rae_p (:),&
       deptn_p(:),&
       depta_p(:)
  character(200) :: fn, str1,dir_out

  integer,parameter :: nflag=5
  integer :: iflag_cond(nflag)

  integer :: i,ip,n,ibuf,it,it_max,it_5gk,it_10gk,it_1gk,it_3gk,it_tem_max, &
       iflag_tem05,iflag_tem07,iflag_tem10,iflag_sc70,iflag_sc74, &
       it_5gk_after_10gk

  real(8) :: smax_traj_7_0,smax_traj_7_4,smin_traj_7_0,smin_traj_7_4

  integer,allocatable :: n_cond(:)
  real(8),allocatable :: &
       mass_traj(:), &
       t_fin_traj(:), x_fin_traj(:), y_fin_traj(:), z_fin_traj(:), vr_fin_traj(:), s_fin_traj(:), ye_fin_traj(:), ut1_fin_traj(:), hut_fin_traj(:), &
       tem_max_traj(:), tem_max_af3gk_traj(:),t_tem_max_traj(:), t_tem_max_af3gk_traj(:), &
       t_5gk_traj(:), s_5gk_traj(:), ye_5gk_traj(:), texp_5gk_traj(:), &
       t_3gk_traj(:), t_1gk_traj(:), ye_10gk_traj(:), t_10gk_traj(:), &
       time_50gk_25gk_traj(:), rho_fin_traj(:), vx_fin_traj(:), vy_fin_traj(:), vz_fin_traj(:), r_ini_traj(:)

  real(8) :: &
       t_fin, x_fin, y_fin, z_fin, vr_fin, s_fin, ye_fin, &
       tem_max, tem_max_af3gk, t_tem_max_af3gk, t_tem_max, &
       t_5gk, s_5gk, ye_5gk, texp_5gk, &
       t_3gk, t_1gk, ye_10gk, t_10gk, &
       time_50gk_25gk, rho_fin, vx_fin, vy_fin, vz_fin, r_ini
  real(8) :: temp_gk, x_ini, y_ini, z_ini,r_fin, x_5gk, y_5gk, z_5gk, r_5gk, vx_5gk, vy_5gk, vz_5gk, vr_5gk, s_10gk, s_3gk, ye_3gk, s_1gk,ye_1gk

  real(8) :: s1,s0

  integer :: my_thr, max_thr

  integer :: job, job_min, job_max, nsteps_job

  ! additional values
  real(8),allocatable :: &
       v_max_traj(:)
  !
  
  dir_out = dir_read

  job_min = 0
  job_max = 0
  job = 1
  jobs:do
     write(str1,'(i10)') job
     
     fn = trim(dir_read) // "/steps_"//trim(adjustl(str1))//".dat"
     if(access(fn," ")==0)then
        if(job_min==0) job_min = job
        job_max=job
     else
        if(job_min/=0)exit jobs
     endif
     job=job+1
  enddo jobs
  write(6,'("job:",i5,"--",i5)') job_min,job_max
  
  
  itt_min = 1
  itt_max = 1
  do job=job_min,job_max
     write(str1,'(i10)') job
     fn = trim(dir_read) // "/steps_"//trim(adjustl(str1))//".dat"
     open(10,file=fn,status="old",action="read")
     read(10,*) nsteps_job
     close(10)
     itt_max = itt_max + nsteps_job
  enddo
  
  write(*,*) "model = ", trim(model)
  write(*,'(a)') trim(dir_read)//"/report_ptr.dat"
  open(10,file=trim(dir_read)//"/report_ptr.dat",status="old")
  read(10,*)
  read(10,*) np
  
  itt_max = itt_max + 1
  write(*,'("np, itt_min, itt_max=",3i7)') np,itt_min,itt_max

  ! allocate
  allocate( &
       time     (itt_min:itt_max),&
       x_p      (itt_min:itt_max),&
       y_p      (itt_min:itt_max),&
       z_p      (itt_min:itt_max),&
       qrho_p   (itt_min:itt_max),&
       ye_p     (itt_min:itt_max),&
       tem_p    (itt_min:itt_max),&
       ut_p     (itt_min:itt_max),&
       qb_p     (itt_min:itt_max),&
       sen_p    (itt_min:itt_max),&
       vlx_p    (itt_min:itt_max),&
       vly_p    (itt_min:itt_max),&
       vlz_p    (itt_min:itt_max),&
       hhh_p    (itt_min:itt_max),&
       dt_p     (itt_min:itt_max),&
       rne_p    (itt_min:itt_max),&
       rae_p    (itt_min:itt_max),&
       deptn_p  (itt_min:itt_max),&
       depta_p  (itt_min:itt_max) )
  allocate (n_cond(np) )
  allocate ( &
       mass_traj(np), &
       t_fin_traj(np), x_fin_traj(np), y_fin_traj(np), z_fin_traj(np), vr_fin_traj(np), s_fin_traj(np), ye_fin_traj(np), ut1_fin_traj(np), hut_fin_traj(np), &
       tem_max_traj(np), tem_max_af3gk_traj(np), t_tem_max_traj(np), t_tem_max_af3gk_traj(np), &
       t_5gk_traj(np), s_5gk_traj(np), ye_5gk_traj(np), texp_5gk_traj(np), &
       t_3gk_traj(np), t_1gk_traj(np), ye_10gk_traj(np), t_10gk_traj(np), &
       time_50gk_25gk_traj(np), rho_fin_traj(np), vx_fin_traj(np), vy_fin_traj(np), vz_fin_traj(np), r_ini_traj(np) )

  allocate ( &
       v_max_traj(np) )

  open(newunit=nunit,file=trim(dir_out)//"/condition_number.dat",status="replace")
  write(nunit,'("# number - condition correspondence")')
  write(nunit,'("# number, Tmax > 5 GK?, Tmax > 7 GK?, Tmax >10 GK?, dS/<S> < 0.2 in T=7-0 GK?, dS/<S> < 0.2 in T=7-4 GK?")')
  do n=0,2**(nflag) - 1
     iflag_cond(:) = 0
     i = 1
     ibuf = n
     do while(ibuf>0)
        iflag_cond(i) = mod(ibuf,2)
        ibuf = ibuf/2
        i = i + 1
     enddo
     write(nunit,'(i8,3i14,2i27)') n, iflag_cond(1:5)
  enddo
  close(nunit)

  !$omp parallel default(none) &
  !$omp shared(np,dir_out,n_cond,mass_traj, &
  !$omp   t_fin_traj, x_fin_traj, y_fin_traj, z_fin_traj, vr_fin_traj, s_fin_traj, ye_fin_traj, ut1_fin_traj, hut_fin_traj, &
  !$omp   tem_max_traj, tem_max_af3gk_traj, t_tem_max_traj, t_tem_max_af3gk_traj, &
  !$omp   t_5gk_traj, s_5gk_traj, ye_5gk_traj, texp_5gk_traj, &
  !$omp   t_3gk_traj, t_1gk_traj, ye_10gk_traj, t_10gk_traj, &
  !$omp   time_50gk_25gk_traj, rho_fin_traj, vx_fin_traj, vy_fin_traj, vz_fin_traj, r_ini_traj, &
  !$omp   v_max_traj) &
  !$omp private( fn,str1,nunit,it,time, x_p,y_p, z_p, vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, sen_p, rne_p, rae_p, &
  !$omp   it_max,it_tem_max,it_5gk,it_10gk,it_3gk,it_1gk, temp_gk, x_ini,y_ini,z_ini, &
  !$omp   smax_traj_7_0,smax_traj_7_4,smin_traj_7_0,smin_traj_7_4, r_fin, &
  !$omp   s1,s0,x_5gk, y_5gk, z_5gk, r_5gk, vx_5gk, vy_5gk, vz_5gk, vr_5gk, s_10gk, s_3gk, ye_3gk, s_1gk,ye_1gk, &
  !$omp   iflag_tem05,iflag_tem07,iflag_tem10,iflag_sc70,iflag_sc74, it_5gk_after_10gk, iflag_cond, &
  !$omp   t_fin, x_fin, y_fin, z_fin, vr_fin, s_fin, ye_fin, &
  !$omp   tem_max, tem_max_af3gk, t_tem_max_af3gk, t_tem_max, &
  !$omp   t_5gk, s_5gk, ye_5gk, texp_5gk, &
  !$omp   t_3gk, t_1gk, ye_10gk, t_10gk, &
  !$omp   time_50gk_25gk, rho_fin, vx_fin, vy_fin, vz_fin, r_ini, &
  !$omp   my_thr,max_thr)
  my_thr  = omp_get_thread_num()
  max_thr = omp_get_max_threads()
  if(my_thr==0)write(6,*) "max thread = ",max_thr
  !$omp do
  do ip =1,np

     write(str1,'(i8.8)') ip
     fn=trim(dir_out)//"/traj_"//trim(str1)//".dat"

     if(access(fn," ") == 0)then
        open(newunit=nunit,file=fn,status="old")
        read(nunit,*)
        ! read(nunit,*)
        read(nunit,*)
        read(nunit,'(16x,es13.5,20x,2es13.5)') mass_traj(ip), ut1_fin_traj(ip), hut_fin_traj(ip)
        read(nunit,*)

        it=1
        do
           read(nunit,*,end=99) &
                time(it), &
                x_p(it), &
                y_p(it), &
                z_p(it), &
                vlx_p(it), &
                vly_p(it), &
                vlz_p(it), &
                qrho_p(it), &
                tem_p(it), &
                ye_p(it), &
                sen_p(it), &
                rne_p(it), &
                rae_p(it)
           
           it = it + 1
        enddo
99      continue
        close(nunit)
        it_max = it - 1

        ! analysis of the trajectories
        smax_traj_7_0=0.d0
        smin_traj_7_0=1.d99
        smax_traj_7_4=0.d0
        smin_traj_7_4=1.d99
        tem_max=0.d0
        tem_max_af3gk = 0.d0
        it_tem_max=0
        it_5gk=0
        it_10gk=0
        it_3gk = 0
        it_1gk = 0
        t_tem_max = 0.d0
        t_tem_max_af3gk = 0.d0

        time_50gk_25gk = 0.d0

        v_max = 0d0
        do it = 1,it_max
           temp_gk = tem_p(it)
           if(tem_max < temp_gk)then
              tem_max = temp_gk
              it_tem_max = it
           endif
           
           v_max = max(v_max, sqrt(sqrt(vx_p(it)**2 + vy_p(it)**2 + vz_p(it)**2)))
           
           ! max, minimum entropy values in 0 < T < 7
           if(0.d0 < temp_gk .and. temp_gk < 7.d9)then
              smax_traj_7_0  = max(smax_traj_7_0 ,sen_p(it))
              smin_traj_7_0  = min(smin_traj_7_0 ,sen_p(it))
           endif
           ! max, minimum entropy values in 4 < T < 7
           if(4.d9 < temp_gk .and. temp_gk < 7.d9)then
              smax_traj_7_4  = max(smax_traj_7_4 ,sen_p(it))
              smin_traj_7_4  = min(smin_traj_7_4 ,sen_p(it))
           endif
           ! time at which the particle crosses T = 5 GK
           if( it < it_max )then
              if( 5.d9 <= tem_p(it) .and. 5.d9 > tem_p(it+1) ) then
                 if(it_5gk==0)it_5gk=it
              endif
              if( 3.d9 <= tem_p(it) .and. 3.d9 > tem_p(it+1) ) then
                 if(it_3gk==0)it_3gk=it
              endif
              if( 1.d9 <= tem_p(it) .and. 1.d9 > tem_p(it+1) .and. it_3gk > 0 ) then
                 it_1gk=it
              endif

              ! 10GK
              if( 10.d9 <= tem_p(it) .and. 10.d9 > tem_p(it+1) ) then
                 it_10gk = it
                 it_5gk=0
                 it_3gk = 0
                 it_1gk = 0

                 tem_max_af3gk = 0.d0
                 time_50gk_25gk = 0.d0
              endif

              if( it_3gk > 0 .and. tem_max_af3gk < temp_gk)then
                 tem_max_af3gk = temp_gk
                 t_tem_max_af3gk = time(it)
              endif

              ! integrated time of T<10 GK
              if( 2.5d9 < tem_p(it) .and. tem_p(it) <5.0d9 ) then
                 time_50gk_25gk = time_50gk_25gk + (time(it+1)-time(it))
              endif

           endif

        enddo

        ! values at initial time
        x_ini  = x_p(1)
        y_ini  = y_p(1)
        z_ini  = z_p(1)
        r_ini  = sqrt(x_ini**2 + y_ini**2 + z_ini**2)

        ! values at final time
        t_fin  = time(it_max)
        x_fin  = x_p(it_max)
        y_fin  = y_p(it_max)
        z_fin  = z_p(it_max)
        r_fin  = sqrt(x_fin**2 + y_fin**2 + z_fin**2)

        vx_fin = vlx_p(it_max)
        vy_fin = vly_p(it_max)
        vz_fin = vlz_p(it_max)
        vr_fin = vx_fin *x_fin/r_fin + vy_fin*y_fin/r_fin + vz_fin *z_fin/r_fin
        s_fin  = sen_p(it_max)
        ye_fin = ye_p(it_max)
        rho_fin= qrho_p(it_max)
        
        ! important value for nuc. reaction
        it = it_5gk
        if(it>0)then
           s1 = (5.d9-tem_p(it))/(tem_p(it+1)-tem_p(it))
           s0 = 1.d0-s1
           t_5gk = s1*time (it+1) + s0*time (it)
           s_5gk = s1*sen_p(it+1) + s0*sen_p(it)
           ye_5gk= s1* ye_p(it+1) + s0* ye_p(it)
           x_5gk = s1*  x_p(it+1) + s0*  x_p(it)
           y_5gk = s1*  y_p(it+1) + s0*  y_p(it)
           z_5gk = s1*  z_p(it+1) + s0*  z_p(it)
           vx_5gk= s1*vlx_p(it+1) + s0*vlx_p(it)
           vy_5gk= s1*vly_p(it+1) + s0*vly_p(it)
           vz_5gk= s1*vlz_p(it+1) + s0*vlz_p(it)
           r_5gk = sqrt(x_5gk**2 +y_5gk**2 +z_5gk**2) + 1.d-20
           vr_5gk= vx_5gk *x_5gk/r_5gk + vy_5gk *y_5gk/r_5gk + vz_5gk *z_5gk/r_5gk
           texp_5gk = r_5gk/vr_5gk
        else
           s_5gk    = 0.d0
           ye_5gk   = 0.d0
           texp_5gk = 0.d0
           t_5gk    = 0.d0
        endif

        it = it_10gk
        if(it>0)then
           s1 = (10.d9-tem_p(it))/(tem_p(it+1)-tem_p(it))
           s0 = 1.d0-s1
           t_10gk = s1*time (it+1) + s0*time (it)
           s_10gk = s1*sen_p(it+1) + s0*sen_p(it)
           ye_10gk= s1* ye_p(it+1) + s0* ye_p(it)
           ! x_10gk = s1*  x_p(it+1) + s0*  x_p(it)
           ! y_10gk = s1*  y_p(it+1) + s0*  y_p(it)
           ! z_10gk = s1*  z_p(it+1) + s0*  z_p(it)
           ! vx_10gk= s1*vlx_p(it+1) + s0*vlx_p(it)
           ! vy_10gk= s1*vly_p(it+1) + s0*vly_p(it)
           ! vz_10gk= s1*vlz_p(it+1) + s0*vlz_p(it)
           ! r_10gk = sqrt(x_10gk**2 +y_10gk**2 +z_10gk**2) + 1.d-20
           ! vr_10gk= vx_10gk *x_10gk/r_10gk + vy_10gk *y_10gk/r_10gk + vz_10gk *z_10gk/r_10gk
           ! texp_10gk = r_10gk/vr_10gk
        else
           s_10gk    = sen_p(it_tem_max)
           ye_10gk   = ye_p(it_tem_max)
           ! texp_10gk = 0.d0
           t_10gk    = 0.d0
        endif

        t_tem_max = time(it_tem_max)

        it = it_3gk
        if(it>0)then
           s1 = (3.d9-tem_p(it))/(tem_p(it+1)-tem_p(it))
           s0 = 1.d0-s1
           t_3gk = s1*time (it+1) + s0*time (it)
           s_3gk = s1*sen_p(it+1) + s0*sen_p(it)
           ye_3gk= s1* ye_p(it+1) + s0* ye_p(it)
        else
           t_3gk    = 0.d0
           s_3gk    = 0.d0
           ye_3gk   = 0.d0
        endif
        
        it = it_1gk
        if(it>0)then
           s1 = (1.d9-tem_p(it))/(tem_p(it+1)-tem_p(it))
           s0 = 1.d0-s1
           t_1gk = s1*time (it+1) + s0*time (it)
           s_1gk = s1*sen_p(it+1) + s0*sen_p(it)
           ye_1gk= s1* ye_p(it+1) + s0* ye_p(it)
        else
           t_1gk    = 0.d0
           s_1gk    = 0.d0
           ye_1gk   = 0.d0
        endif

        ! making flags
        if( tem_max > 5.d9 )then
           iflag_tem05 = 1
        else
           iflag_tem05 = 0
        endif

        if( tem_max > 7.d9 )then
           iflag_tem07 = 1
        else
           iflag_tem07 = 0
        endif

        if( tem_max > 10.d9 )then
           iflag_tem10 = 1
        else
           iflag_tem10 = 0
        endif

        if(abs((smax_traj_7_0-smin_traj_7_0)/(smax_traj_7_0+smin_traj_7_0)) < 0.2d0)then
           iflag_sc70 = 1
        else
           iflag_sc70 = 0
        endif

        if(abs((smax_traj_7_4-smin_traj_7_4)/(smax_traj_7_4+smin_traj_7_4)) < 0.2d0)then
           iflag_sc74 = 1
        else
           iflag_sc74 = 0
        endif

        iflag_cond(:) = 0
        iflag_cond(1) = iflag_tem05
        iflag_cond(2) = iflag_tem07
        iflag_cond(3) = iflag_tem10
        iflag_cond(4) = iflag_sc70
        iflag_cond(5) = iflag_sc74
        n_cond(ip) = 0
        do i=1,nflag
           if(iflag_cond(i)==1) n_cond(ip) = n_cond(ip) + 2**(i-1)
        enddo
        
        
        t_fin_traj(ip)           = t_fin
        x_fin_traj(ip)           = x_fin
        y_fin_traj(ip)           = y_fin
        z_fin_traj(ip)           = z_fin
        vr_fin_traj(ip)          = vr_fin
        s_fin_traj(ip)           = s_fin
        ye_fin_traj(ip)          = ye_fin
        tem_max_traj(ip)         = tem_max
        tem_max_af3gk_traj(ip)   = tem_max_af3gk
        t_tem_max_traj(ip)       = t_tem_max
        t_tem_max_af3gk_traj(ip) = t_tem_max_af3gk
        t_5gk_traj(ip)           = t_5gk
        s_5gk_traj(ip)           = s_5gk
        ye_5gk_traj(ip)          = ye_5gk
        texp_5gk_traj(ip)        = texp_5gk
        t_3gk_traj(ip)           = t_3gk
        t_1gk_traj(ip)           = t_1gk
        ye_10gk_traj(ip)         = ye_10gk
        t_10gk_traj(ip)          = t_10gk
        time_50gk_25gk_traj(ip)  = time_50gk_25gk
        rho_fin_traj(ip)         = rho_fin
        vx_fin_traj(ip)          = vx_fin
        vy_fin_traj(ip)          = vy_fin
        vz_fin_traj(ip)          = vz_fin
        r_ini_traj(ip)           = r_ini
        v_max_traj(ip) = v_max

     endif

     !if(mod(ip,100)==0) write(6,'("ip = ", i8)')  ip
     if(my_thr==0) write(6,*) "ip=",ip,np/max_thr
  enddo
  !$omp end do
  !$omp end parallel
  
  open(newunit=nunit,file=trim(dir_out)//"/ana_traj.dat",status="replace")  
  write(nunit,'("# model: ",a)') trim(model)
  write(nunit,'("# 1:id  2:condition num.  3:mass [g]  4:time(final) [s]   5:x(final) [cm]   6:y(final) [cm]  7:z(final) [cm]   8:Vrad(final) [cm/s]   9:S(final) [k_b/nuc]   10:Ye(final)  11:Time(5GK) [s]  12:S(5GK) [k_b/nuc]  13:Ye(5GK)  14:t_exp(5GK) [s]  15:Tmax [K]   16:ut+1(final)  17:hut+h_atm(final)   18:Ye(10GK)   19:Time(3GK)   20:Time(1GK)  21:Tmax(after3GK)  22:Time(after3GK)  23:Time(5-2.5GK)  24:rho(final)  25:vx(final)  26:vy(final)  27:vz(final)  28:time(10GK)  29:time(Tmax) 30:r(init) [cm] 31:time(first 5GK after 10GK)")')
  
  !write(nunit,'("#     id condition num.       mass [g]       Time [s]         x [cm]         y [cm]         z [cm]    Vrad [cm/s]    S [k_b/nuc]             Ye       Time [s]    S [k_b/nuc]             Ye      t_exp [s]       Tmax [K]           ut+1      hut+h_atm       Ye(10GK)      Time(3GK)      Time(1GK) Tmax(after3GK) Time(after3GK)  Time(5-2.5GK) rho (boundary)             vx             vy             vz     Time(10GK)     Time(Tmax)")')
  do ip=1,np
     write(nunit,'(i8,i15,99es15.7)') &
          ip, n_cond(ip), &
          mass_traj(ip), &
          t_fin_traj(ip),x_fin_traj(ip),y_fin_traj(ip),z_fin_traj(ip),vr_fin_traj(ip),s_fin_traj(ip),ye_fin_traj(ip),t_5gk_traj(ip),s_5gk_traj(ip),ye_5gk_traj(ip),texp_5gk_traj(ip), &
          tem_max_traj(ip),ut1_fin_traj(ip),hut_fin_traj(ip),ye_10gk_traj(ip),t_3gk_traj(ip),t_1gk_traj(ip),tem_max_af3gk_traj(ip),t_tem_max_af3gk_traj(ip),time_50gk_25gk_traj(ip),&
          rho_fin_traj(ip),vx_fin_traj(ip),vy_fin_traj(ip),vz_fin_traj(ip),t_10gk_traj(ip),t_tem_max_traj(ip),r_ini_traj(ip)
  enddo
  
  close(nunit)

  open(newunit=nunit,file=trim(dir_out)//"/other_traj.dat",status="replace")  
  write(nunit,'("# model: ",a)') trim(model)
  write(nunit,'("# 1:id, 2:mass, 3:v_max")')
  
  !write(nunit,'("#     id condition num.       mass [g]       Time [s]         x [cm]         y [cm]         z [cm]    Vrad [cm/s]    S [k_b/nuc]             Ye       Time [s]    S [k_b/nuc]             Ye      t_exp [s]       Tmax [K]           ut+1      hut+h_atm       Ye(10GK)      Time(3GK)      Time(1GK) Tmax(after3GK) Time(after3GK)  Time(5-2.5GK) rho (boundary)             vx             vy             vz     Time(10GK)     Time(Tmax)")')
  do ip=1,np
     write(nunit,'(i8,99es15.7)') &
          ip,mass_traj(ip),v_max_traj(ip)
  enddo
  
  close(nunit)


  write(6,*) "ana_traj.dat done."
  
end subroutine tr_analysis
