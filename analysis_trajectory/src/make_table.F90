program make_table

  use unit
  use module_weak_interaction
  use module_etanu

  use, intrinsic :: ieee_arithmetic

  implicit none

  logical,parameter :: output_once = .false.

  real(8),parameter :: emev=0.51099996d0
  real(8),parameter :: mev2t9= 1.160445d1

  character(256) :: fn_eos,fn_eosb,fn_enu,fn_ynu,dir_read,fn,fn_out,label
  character(10) :: str1
  integer :: ntemp_in, nrho_in, nye_in

  character(256) :: fn_caprate, fn_caprate_neg, fn_etanu

  integer :: itt_min,itt_max,np,ip,it,nt,np_skip,np_start
  real(8),allocatable :: time(:),mass_p(:)
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


  real(8),allocatable :: yn_p(:),ya_p(:)

  integer :: unit_traj
  integer :: itemp,ntemp_list
  real(8) :: temp_min,temp_max
  integer,allocatable :: np_temp_list(:)
  real(8),allocatable :: temp_list(:),ye_av_list(:),mass_temp_list(:),ye2_av_list(:),eta_av_list(:),eta2_av_list(:),sen_av_list(:),sen2_av_list(:),texp_av_list(:),texp2_av_list(:)
  real(8),allocatable :: caprate_av_list(:),caprate2_av_list(:),absrate_av_list(:),absrate2_av_list(:),yecap_av_list(:),yecap2_av_list(:),yeabs_av_list(:),yeabs2_av_list(:)
  integer :: idx
  
  real(8) :: hoge
  integer :: ip_current
  integer,allocatable :: np_active(:)
  real(8),allocatable :: time_temp_traj_list(:,:),ye_temp_traj_list(:,:),mass_temp_traj_list(:,:),eta_temp_traj_list(:,:),sen_temp_traj_list(:,:),texp_temp_traj_list(:,:), vel_temp_traj_list(:,:), r_temp_traj_list(:,:)
  real(8),allocatable :: caprate_temp_traj_list(:,:),absrate_temp_traj_list(:,:),yecap_temp_traj_list(:,:),yeabs_temp_traj_list(:,:),yemu0_temp_traj_list(:,:)
  integer,allocatable :: ip_list(:)
  real(8),allocatable :: mass_list(:)
  real(8),allocatable :: dummy_list(:)

  real(8),parameter :: temp_weak = 1d0*tem_uni
  integer :: it_last_hot, it_search_start
  real(8) :: temp_target 
  ! integer :: ip_lw, ip_med, ip_rw
  ! real(8) :: fac_lw=0.15d0, fac_med=0.5d0, fac_rw=0.85d0
  ! real(8) :: ss_lw, ss_med, ss_rw
  ! real(8),allocatable :: ye_lw_temp_list(:),ye_med_temp_list(:),ye_rw_temp_list(:)
  ! real(8),allocatable :: eta_lw_temp_list(:),eta_med_temp_list(:),eta_rw_temp_list(:)
  ! real(8),allocatable :: sen_lw_temp_list(:),sen_med_temp_list(:),sen_rw_temp_list(:)
  ! real(8),allocatable :: texp_lw_temp_list(:),texp_med_temp_list(:),texp_rw_temp_list(:)
  ! real(8),allocatable :: caprate_lw_temp_list(:),caprate_med_temp_list(:),caprate_rw_temp_list(:)
  ! real(8),allocatable :: absrate_lw_temp_list(:),absrate_med_temp_list(:),absrate_rw_temp_list(:)
  ! real(8),allocatable :: yecap_lw_temp_list(:),yecap_med_temp_list(:),yecap_rw_temp_list(:)
  ! real(8),allocatable :: yeabs_lw_temp_list(:),yeabs_med_temp_list(:),yeabs_rw_temp_list(:)
  ! real(8),allocatable :: yemu0_lw_temp_list(:),yemu0_med_temp_list(:),yemu0_rw_temp_list(:)

  real(8) :: tt,tt1,ye_tmp,den_tmp,tem_tmp, yn_tmp, ya_tmp
  
  integer :: ke,ke1,ie,ie1,je,je1
  real(8) :: uu,uup,ss,ssp,ttp,che_tmp,eta_tmp,sen_tmp,vr_tmp,r_tmp,texp_tmp,t_tmp, dlntemp_dt_tmp, dlntemp_dt_long

  real(8) :: ye_equil_cap,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate
  real(8) :: ye_equil_abs,abs_n_nrate,abs_a_nrate,rne_tmp,rae_tmp,taun_tmp,taua_tmp
  real(8) :: ye_mu0

  logical :: mass_weighted
  real(8) :: my_time

  integer :: access
  real(8) :: eta_n, eta_a, rnetrap_tmp, raetrap_tmp, abs_nt_nrate,abs_at_nrate

  real(8) :: ecap_nrate_block,pcap_nrate_block

  !!! ana_traj
  real(8) :: buf(100)
  real(8),allocatable :: hut(:), ut1(:)

  integer :: nunit, nunit_time, nunit_ye, nunit_eta, nunit_entr, nunit_texp, nunit_rcap, nunit_yecap, nunit_vel, nunit_r
  integer :: i
    
  !call make_epcap_table(1d-2,1d3,1d-2,1d2,401,401)
  !stop

  block
    use inputparser
    call get_string_parameter("parameters","dir_read",dir_read)
    call get_string_parameter("parameters","label",label)
    call get_logical_parameter("parameters","mass_weighted",mass_weighted)
    call get_string_parameter("parameters","fn_eos",fn_eos)
    call get_string_parameter("parameters","fn_eosb",fn_eosb)
    call get_integer_parameter("parameters","nrho",nrho_in)
    call get_integer_parameter("parameters","ntemp",ntemp_in)
    call get_integer_parameter("parameters","nye",nye_in)
    write(6,*) fn_eos
    write(6,*) nrho_in, ntemp_in, nye_in
    call get_integer_parameter("parameters","ntemp_list",ntemp_list)
    call get_double_parameter("parameters","temp_min",temp_min)
    call get_double_parameter("parameters","temp_max",temp_max)
    write(6,*) ntemp_list, temp_min, temp_max
    itt_min = 1
    call get_integer_parameter("parameters","itt_max",itt_max)
    call get_integer_parameter("parameters","np",np)
    call get_integer_parameter("parameters","np_start",np_start)
    call get_integer_parameter("parameters","np_skip",np_skip)
    write(6,*) itt_min, itt_max
    write(6,*) np_start, np, np_skip
    call get_string_parameter("parameters","fn_caprate",fn_caprate)
    write(6,*)trim( fn_caprate)
    call get_string_parameter("parameters","fn_caprate_neg",fn_caprate_neg)
    call get_string_parameter("parameters","fn_etanu",fn_etanu)

    
  end block
  fn_enu="/data/scratch/sfujibayashi/EOS/enu_to_chnu2"
  fn_ynu="/data/scratch/sfujibayashi/EOS/ynu_to_chnu"
  
  call init_weak_table(fn_caprate, fn_caprate_neg)
  
  call etanu_init(fn_etanu)
  ! rhoynu = 1d11*2d-2

  ! call etanu_eta_f3(rhoynu*6d23*1d-39/5d0**3*197d0**3*2d0*pi**2, eta_nu, f3_nu)
  ! write(*,*)eta_nu, f3_nu
  ! stop
  
!!! 3D models
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6/data_post"; fn_out="ye_ave_q6.dat"
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q4/data_post"; label="q4"; mass_weighted=.true.
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_3D_SFHo_135_135-150mstg/data"; label="sfho135135dyn_ut1"; mass_weighted=.true.

!!! 2D models  
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_SFHo_120_150-150mstg/data";label="sfho1215_trap_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_SFHo_130_140-150mstg/data";label="sfho1314_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_SFHo_135_135-150mstg/data";label="sfho135135_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_DD2_125M/data";label="DD2-125M_number"; mass_weighted=.false.

!!! collapsar in 2D
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_T20/data_bind_1ms128r2.0e+09"; label="T20_reduced"; mass_weighted=.true.
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_AD20x2/data_bind_1ms128r2.0e+09"; label="AD20x2_reduced"; mass_weighted=.true.
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_AD20x1/data_bind_1ms128r2.0e+09"; label="AD20x1_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_AD09x1/data_bind_1ms128r2.0e+09"; label="AD09x1_reduced"; mass_weighted=.true.
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_BHdisk/data_bind_1ms64r2.0e+09"; label="BHdisk"; mass_weighted=.true.

  
  ! ntemp_list=20
  ! temp_min = 0.4d0
  ! temp_max = 3.d0
  allocate(temp_list(ntemp_list),np_temp_list(ntemp_list),ye_av_list(ntemp_list),ye2_av_list(ntemp_list),mass_temp_list(ntemp_list),eta_av_list(ntemp_list),eta2_av_list(ntemp_list),sen_av_list(ntemp_list),sen2_av_list(ntemp_list), &
       texp_av_list(ntemp_list),texp2_av_list(ntemp_list))

  allocate(caprate_av_list(ntemp_list),caprate2_av_list(ntemp_list),absrate_av_list(ntemp_list),absrate2_av_list(ntemp_list),yecap_av_list(ntemp_list),yecap2_av_list(ntemp_list),yeabs_av_list(ntemp_list),yeabs2_av_list(ntemp_list))

  do itemp=1,ntemp_list
     temp_list(itemp) = 10d0**(log10(temp_min) + (log10(temp_max)-log10(temp_min))*dble(itemp-1)/dble(ntemp_list-1))
     ! write(6,*) itemp,temp_list(itemp)
  enddo

  block
    use module_eos
    ! call readeos(fn_eos,fn_eosb,fn_ynu,fn_enu)
    call readeos(fn_eos,nrho_in, ntemp_in, nye_in)
  end block

  ! call test_ye_equil

  allocate(ut1(np),hut(np))
  open(11,file=trim(dir_read)//"/ana_traj.dat",status="old",action="read")
  read(11,*);read(11,*)
  do ip=1,np
     read(11,*) buf(1:30)
     ut1(ip) = buf(16)
     hut(ip) = buf(17)
  enddo
  close(11)

  allocate( ip_list(np), mass_list(np))
  allocate( time_temp_traj_list(ntemp_list,np),ye_temp_traj_list(ntemp_list,np),mass_temp_traj_list(ntemp_list,np),eta_temp_traj_list(ntemp_list,np),sen_temp_traj_list(ntemp_list,np),texp_temp_traj_list(ntemp_list,np),vel_temp_traj_list(ntemp_list,np),r_temp_traj_list(ntemp_list,np) )
  allocate( caprate_temp_traj_list(ntemp_list,np),absrate_temp_traj_list(ntemp_list,np),yecap_temp_traj_list(ntemp_list,np),yeabs_temp_traj_list(ntemp_list,np),yemu0_temp_traj_list(ntemp_list,np) )

  itt_max = itt_max + 1
  write(*,'("np, itt_min, itt_max=",3i7)') np,itt_min,itt_max

  allocate( &
       mass_p(np), &
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

  allocate(yn_p    (itt_min:itt_max),&
           ya_p    (itt_min:itt_max) )

  
  
  if(.not.output_once)then
     open(newunit=nunit_time, file="traj_time", status="replace", action="write")
     write(nunit_time,'("#", 99a15)') "id", "mass", "T"
     write(nunit_time,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_ye, file="traj_Ye", status="replace", action="write")
     write(nunit_ye,'("#", 99a15)') "id", "mass", "T"
     write(nunit_ye,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_eta, file="traj_eta", status="replace", action="write")
     write(nunit_eta,'("#", 99a15)') "id", "mass", "T"
     write(nunit_eta,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_entr, file="traj_sen", status="replace", action="write")
     write(nunit_entr,'("#", 99a15)') "id", "mass", "T"
     write(nunit_entr,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_texp, file="traj_texp", status="replace", action="write")
     write(nunit_texp,'("#", 99a15)') "id", "mass", "T"
     write(nunit_texp,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_vel, file="traj_vel", status="replace", action="write")
     write(nunit_vel,'("#", 99a15)') "id", "mass", "T"
     write(nunit_vel,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_r, file="traj_r", status="replace", action="write")
     write(nunit_r,'("#", 99a15)') "id", "mass", "T"
     write(nunit_r,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_rcap, file="traj_caprate", status="replace", action="write")
     write(nunit_rcap,'("#", 99a15)') "id", "mass", "T"
     write(nunit_rcap,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     
     open(newunit=nunit_yecap, file="traj_yecap", status="replace", action="write")
     write(nunit_yecap,'("#", 99a15)') "id", "mass", "T"
     write(nunit_yecap,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
  endif

  idx = 0
  do ip=np_start,np,np_skip
     
     write(str1,'(i8.8)') ip
     fn = trim(dir_read)//"/traj_"//trim(str1)//".dat"
     write(6,'(a)') trim(fn)
     open(newunit=unit_traj,file=fn,status="old",action="read")
     read(unit_traj,*)
     read(unit_traj,*)
     read(unit_traj,'(16x,es13.5)') mass_p(ip)
     read(unit_traj,*)
     it = 0
     do
        nt = it
        it = it + 1
        read(unit_traj,*,end=99) &
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
             rae_p(it), &
             deptn_p(it), &
             depta_p(it)

     enddo
99   continue
     close(unit_traj)

!!! misc. file
     fn = trim(dir_read)//"/misc_"//trim(str1)//".dat"
     if(access(fn," ")==0)then
        
        open(newunit=unit_traj,file=fn,status="old")
        read(unit_traj,*)
        read(unit_traj,*)
        read(unit_traj,*)
        read(unit_traj,*)
        it = 0
        do
           it = it + 1
           read(unit_traj,*,end=991) &
                hoge, &
                yn_p(it), &
                ya_p(it)
           
        enddo
991     continue
        close(unit_traj)
     else
        yn_p(:) = 0d0
        ya_p(:) = 0d0
        
     endif

     if(.not.mass_weighted) mass_p(ip) = 1d0
     
     if( &
          ! ut1(ip)<0d0.and. &
          hut(ip)<0d0.and. &
          !sen_p(nt)<30.d0 .and. &
          !maxval(tem_p(1:nt))>=temp_list(ntemp_list)*tem_uni.and. &
          !minval(tem_p(1:nt))<=temp_list(1)*tem_uni.and. &
          maxval(tem_p(1:nt))>1d0*tem_uni.and. &
          ! tem_p(nt)<=temp_list(1)*tem_uni.and. &
          tem_p(nt)<=temp_weak.and. &
          !ye_p(nt)<0.15d0.and.  &
          !ye_p(nt)>0.15d0.and.  &
          .true. &
          )then
        
        idx = idx + 1
        
        ip_list(idx) = ip
        mass_list(idx) = mass_p(ip)

        !it=nt
        ! do while(tem_p(it) > minval(tem_p(1:nt)))
           !    it = it - 1
        ! enddo
        
        it_last_hot = 0
        do it = 1, nt
           if (tem_p(it) >= temp_weak) it_last_hot = it
        enddo
        
        if (it_last_hot == 0) then
           ! この tracer は一度も 1 MeV 以上になっていない
           ! Ye を決める高温 phase なし
        endif
        
        if (it_last_hot >= nt) then
           ! 最後まで 1 MeV 以上。まだ freeze-out していない可能性あり
           it_search_start = nt - 1
        else
           do it=it_last_hot,nt
              if(tem_p(it) < tem_p(it-1))then
                 it_search_start = it
              endif
           enddo
           ! it_search_start = it_last_hot
        endif

        it = it_search_start
        if(it==nt) it = nt-1

        write(6,'("it,t,T=",i10,99es12.4)') it_search_start, time(it_search_start), tem_p(it_search_start)/tem_uni
        loop_itemp:do itemp=1,ntemp_list
           
           temp_target = temp_list(itemp)*tem_uni
        
           !write(6,'(99es12.4)') temp_list(itemp)*tem_uni, minval(tem_p(1:nt))
           !if( temp_list(itemp)*tem_uni < minval(tem_p(1:nt)) )then
           ! if( temp_list(itemp)*tem_uni < tem_p(nt) )then
           !    !write(6,*) "skipped"
           !    it=0
           ! endif

           if (temp_target < tem_p(it_search_start)) then
              time_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              ye_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              eta_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              sen_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              texp_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              vel_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              r_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)

              caprate_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              absrate_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              yecap_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              yeabs_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              yemu0_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              cycle loop_itemp
           endif
   
           !write(99,*)
           search_it:do
              if (it <= 0) exit search_it
              !write(99,'(i5,99es12.4)') it, tem_p(it),  temp_list(itemp)*tem_uni, tem_p(it+1)
              !write(6,'(i5,99es12.4)') it, sqrt(x_p(it)**2 + y_p(it)**2 + z_p(it)**2), tem_p(it)-temp_list(itemp)*tem_uni, tem_p(it+1)-temp_list(itemp)*tem_uni
              if(tem_p(it) > temp_target .and. temp_target >= tem_p(it+1)) exit search_it
              !if( (tem_p(it)-temp_list(itemp)*tem_uni)*(tem_p(it+1)-temp_list(itemp)*tem_uni) < 0d0)exit search_it
              it = it - 1
              !if( abs(sen_p(it+1)-sen_p(it))/sen_p(it+1) >0.5d0 ) exit loop_itemp
           enddo search_it
           
           if(it>0)then
              tt1= (temp_list(itemp)*tem_uni-tem_p(it))/(tem_p(it+1)-tem_p(it))
              tt = 1.d0-tt1
              ! write(6,'(i6,99es14.6)') it,tt1, tem_p(it), temp_list(itemp)*tem_uni, tem_p(it+1)
              t_tmp = tt*time(it) + tt1*time(it+1)
              den_tmp = tt*qrho_p(it) + tt1*qrho_p(it+1)
              ye_tmp  = tt*ye_p(it) + tt1*ye_p(it+1)
              !tem_tmp  = tt*tem_p(it) + tt1*tem_p(it+1)
              !tem_tmp = tem_tmp / (mev2t9*1.d9)
              tem_tmp = temp_list(itemp)
              sen_tmp = tt*sen_p(it) + tt1*sen_p(it+1)
              rne_tmp = tt*rne_p(it) + tt1*rne_p(it+1)
              rae_tmp = tt*rae_p(it) + tt1*rae_p(it+1)
              taun_tmp = tt*deptn_p(it) + tt1*deptn_p(it+1)
              taua_tmp = tt*depta_p(it) + tt1*depta_p(it+1)
              yn_tmp = tt*yn_p(it) + tt1*yn_p(it+1)
              ya_tmp = tt*ya_p(it) + tt1*ya_p(it+1)

              vr_tmp = tt *(x_p(it  )*vlx_p(it  ) + y_p(it  )*vly_p(it  ) + z_p(it  )*vlz_p(it  ))/sqrt(x_p(it  )**2+y_p(it  )**2+z_p(it  )**2) &
                      +tt1*(x_p(it+1)*vlx_p(it+1) + y_p(it+1)*vly_p(it+1) + z_p(it+1)*vlz_p(it+1))/sqrt(x_p(it+1)**2+y_p(it+1)**2+z_p(it+1)**2)
              r_tmp =  tt *sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2) &
                      +tt1*sqrt(x_p(it+1)**2+y_p(it+1)**2+z_p(it+1)**2)
              texp_tmp=r_tmp/abs(vr_tmp)

              ! dlntemp_dt_tmp = (tem_p(it+1) - tem_p(it))/(time(it+1) - time(it)) / ((tem_p(it+1) + tem_p(it))/2d0)
              dlntemp_dt_tmp = (log(tem_p(it+1)) - log(tem_p(it)))/(time(it+1) - time(it))
              dlntemp_dt_long = (log(tem_p(it+1000)) - log(tem_p(it-1000)))/(time(it+1000) - time(it-1000))

              ! write(6,'(99es15.7)') temp_list(itemp), vr_tmp, r_tmp, texp_tmp, 1d0/dlntemp_dt_tmp, 1d0/dlntemp_dt_long
              ! write(99,'(99es14.6)') t_tmp, tem_tmp*tem_uni
              block
                use module_eos
                ke  = max(1 , min(nrho-1, int((log10(den_tmp)-rho_e_min  )*drhoi)+1))
                ke1 = ke+1
                uu    = max(0.d0, min(1.d0 ,     (log10(den_tmp)-rho_e(ke))*drhoi)   )
                uup   = 1.d0-uu

                ie = max(1 , min(ntemp-1, int((log10(tem_tmp)-tem_e_min  )*dtemi)+1))
                ie1=ie+1
                ss    = max(0.d0, min(1.d0 ,     (log10(tem_tmp)-tem_e(ie))*dtemi))
                ssp   = 1.d0-ss

                je  = max(1 , min(nye-1, int((ye_tmp-ye_e_min )*dyei )+1))
                je1 = je +1
                tt   = max(0.d0, min(1.d0 ,     (ye_tmp-ye_e(je))*dyei))
                ttp  = 1.d0-tt

                eta_tmp = ssp *ttp *uup * che_e(ie ,je ,ke )   &
                     + ss  *ttp *uup * che_e(ie1,je ,ke )   &
                     + ssp *tt  *uup * che_e(ie ,je1,ke )   &
                     + ssp *ttp *uu  * che_e(ie ,je ,ke1)   &
                     + ss  *tt  *uup * che_e(ie1,je1,ke )   &
                     + ss  *ttp *uu  * che_e(ie1,je ,ke1)   &
                     + ssp *tt  *uu  * che_e(ie ,je1,ke1)   &
                     + ss  *tt  *uu  * che_e(ie1,je1,ke1)
              end block
              
              call ye_equilibrium_capture(den_tmp,tem_tmp,ye_equil_cap,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)
              !call ye_equilibrium_abs(den_tmp,tem_tmp,rne_tmp,rae_tmp,ye_equil_abs,xn_dummy,xp_dummy,abs_n_nrate,abs_a_nrate)
              ! call ye_equilibrium_munu0(den_tmp,tem_tmp,ye_mu0)
              ye_mu0 = 0d0

              call nrate_cap(den_tmp,tem_tmp,ye_tmp,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)
              !call nrate_abs(den_tmp,tem_tmp,ye_tmp,rne_tmp,rae_tmp,xn_dummy,xp_dummy,abs_n_nrate,abs_a_nrate)

              !
              !write(6,'(99es12.4)') 
              !call etanu_eta_ene(den_tmp*yn_tmp,tem_tmp, eta_n, rnetrap_tmp)
              !call etanu_eta_ene(den_tmp*ya_tmp,tem_tmp, eta_a, raetrap_tmp)
              !call nrate_abs(den_tmp,tem_tmp,ye_tmp,rnetrap_tmp,raetrap_tmp,xn_dummy,xp_dummy,abs_nt_nrate,abs_at_nrate)

              ! call nrate_cap_block(den_tmp,tem_tmp,ye_tmp,eta_n, eta_a, eta_dummy,xn_dummy,xp_dummy,ecap_nrate_block,pcap_nrate_block)
              !write(6,'(99es12.4)')  tem_tmp,den_tmp, taun_tmp, taua_tmp, yn_tmp, ya_tmp, eta_n, eta_a, abs_nt_nrate, abs_n_nrate, ecap_nrate_block, ecap_nrate, pcap_nrate_block, pcap_nrate
              !abs_n_nrate = abs_n_nrate *exp(-2d0*taun_tmp)
              !abs_a_nrate = abs_a_nrate *exp(-2d0*taua_tmp)
              abs_n_nrate = abs_n_nrate + abs_nt_nrate
              abs_a_nrate = abs_a_nrate + abs_at_nrate

              ecap_nrate = ecap_nrate/xp_dummy
              pcap_nrate = pcap_nrate/xn_dummy
              abs_n_nrate = abs_n_nrate/xn_dummy
              abs_a_nrate = abs_a_nrate/xp_dummy

              time_temp_traj_list(itemp,idx) = t_tmp
              ye_temp_traj_list(itemp,idx) = ye_tmp
              eta_temp_traj_list(itemp,idx) = eta_tmp
              sen_temp_traj_list(itemp,idx) = sen_tmp
              texp_temp_traj_list(itemp,idx) = texp_tmp
              vel_temp_traj_list(itemp,idx) = vr_tmp
              r_temp_traj_list(itemp,idx) = r_tmp

              caprate_temp_traj_list(itemp,idx) = max(ecap_nrate,pcap_nrate)
              absrate_temp_traj_list(itemp,idx) = max(abs_n_nrate,abs_a_nrate)
              yecap_temp_traj_list(itemp,idx) = ye_equil_cap
              yeabs_temp_traj_list(itemp,idx) = ye_equil_abs
              yemu0_temp_traj_list(itemp,idx) = ye_mu0

           else
              
              time_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              ye_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              eta_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              sen_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              texp_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              vel_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              r_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)

              caprate_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              absrate_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              yecap_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              yeabs_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              yemu0_temp_traj_list(itemp,idx) = ieee_value(hoge, ieee_quiet_nan)
              
           endif
999        continue
           
        enddo loop_itemp
        
        if(.not.output_once)then
           i = idx
           write(nunit_time,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (time_temp_traj_list(itemp,i),itemp=1,ntemp_list)
           
           write(nunit_ye,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (ye_temp_traj_list(itemp,i),itemp=1,ntemp_list)
           
           write(nunit_eta,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (eta_temp_traj_list(itemp,i),itemp=1,ntemp_list)
           
           write(nunit_entr,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (sen_temp_traj_list(itemp,i),itemp=1,ntemp_list)
           
           write(nunit_texp,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (texp_temp_traj_list(itemp,i),itemp=1,ntemp_list)

           write(nunit_vel,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (vel_temp_traj_list(itemp,i),itemp=1,ntemp_list)
           
           write(nunit_r,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (r_temp_traj_list(itemp,i),itemp=1,ntemp_list)

           write(nunit_rcap,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (caprate_temp_traj_list(itemp,i),itemp=1,ntemp_list)
           
           write(nunit_yecap,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (yecap_temp_traj_list(itemp,i),itemp=1,ntemp_list)


           ! do itemp=1,ntemp_list
           !    write(6,'(99es15.7)') temp_list(itemp), r_temp_traj_list(itemp,i), vel_temp_traj_list(itemp,i), texp_temp_traj_list(itemp,i)
           ! end do

        endif
           
     endif
     
  enddo
  
  write(6,*) "# of tracers considered: "
  write(6,*) idx

  if(output_once)then
     open(newunit=nunit_time, file="traj_time", status="replace", action="write")
     write(nunit_time,'("#", 99a15)') "id", "mass", "T"
     write(nunit_time,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_time,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (time_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_time)

     open(newunit=nunit_ye, file="traj_Ye", status="replace", action="write")
     write(nunit_ye,'("#", 99a15)') "id", "mass", "T"
     write(nunit_ye,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_ye,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (ye_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_ye)

     open(newunit=nunit_eta, file="traj_eta", status="replace", action="write")
     write(nunit_eta,'("#", 99a15)') "id", "mass", "T"
     write(nunit_eta,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_eta,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (eta_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_eta)

     open(newunit=nunit_entr, file="traj_sen", status="replace", action="write")
     write(nunit_entr,'("#", 99a15)') "id", "mass", "T"
     write(nunit_entr,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_entr,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (sen_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_entr)

     open(newunit=nunit_texp, file="traj_texp", status="replace", action="write")
     write(nunit_texp,'("#", 99a15)') "id", "mass", "T"
     write(nunit_texp,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_texp,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (texp_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_texp)

     open(newunit=nunit_vel, file="traj_vel", status="replace", action="write")
     write(nunit_vel,'("#", 99a15)') "id", "mass", "T"
     write(nunit_vel,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_vel,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (vel_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_vel)

     open(newunit=nunit_r, file="traj_r", status="replace", action="write")
     write(nunit_r,'("#", 99a15)') "id", "mass", "T"
     write(nunit_r,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_r,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (r_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_r)

     open(newunit=nunit_rcap, file="traj_caprate", status="replace", action="write")
     write(nunit_rcap,'("#", 99a15)') "id", "mass", "T"
     write(nunit_rcap,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_rcap,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (caprate_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_rcap)

     open(newunit=nunit_yecap, file="traj_yecap", status="replace", action="write")
     write(nunit_yecap,'("#", 99a15)') "id", "mass", "T"
     write(nunit_yecap,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
     do i=1,idx
        write(nunit_yecap,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (yecap_temp_traj_list(itemp,i),itemp=1,ntemp_list)
     enddo
     close(nunit_yecap)
  else

     close(nunit_time)
     close(nunit_ye)
     close(nunit_eta)
     close(nunit_entr)
     close(nunit_texp)
     close(nunit_vel)
     close(nunit_r)
     close(nunit_rcap)
     close(nunit_yecap)

  endif
    ! open(newunit=nunit, file="traj_yemu0", status="replace", action="write")
    ! write(nunit,'("#", 99a15)') "id", "mass", "T"
    ! write(nunit,'(" ", 15x,15x,99es15.7)') (temp_list(itemp),itemp=1,ntemp_list)
    ! do i=1,idx
    !    write(nunit,'(" ", i15,99es15.7)') ip_list(i), mass_list(i), (yemu0_temp_traj_list(itemp,i),itemp=1,ntemp_list)
    ! enddo
    ! close(nunit)

end program make_table
