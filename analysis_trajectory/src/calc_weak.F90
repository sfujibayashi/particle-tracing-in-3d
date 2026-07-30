program weak

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
  real(8),allocatable :: caprate_temp_traj_list(:,:),absrate_temp_traj_list(:,:),yecap_temp_traj_list(:,:),yeabs_temp_traj_list(:,:),yemu0_temp_traj_list(:,:), avtexp_temp_traj_list(:,:)
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
  real(8) :: uu,uup,ss,ssp,ttp,che_tmp,eta_tmp,sen_tmp,vr_tmp,r_tmp,texp_tmp,t_tmp, dlntemp_dt_tmp, dlntemp_dt_long, t_exp

  real(8) :: ye_equil_cap,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate
  real(8) :: ye_equil_abs,abs_n_nrate,abs_a_nrate,rne_tmp,rae_tmp,taun_tmp,taua_tmp
  real(8) :: ye_mu0

  logical :: mass_weighted
  real(8) :: my_time

  logical :: file_exists
  real(8) :: eta_n, eta_a, rnetrap_tmp, raetrap_tmp, abs_nt_nrate,abs_at_nrate

  real(8) :: ecap_nrate_block,pcap_nrate_block

  !!! ana_traj
  real(8) :: buf(100)
  real(8),allocatable :: hut(:), ut1(:)

  integer :: nunit, nunit_time, nunit_ye, nunit_eta, nunit_entr, nunit_texp, nunit_rcap, nunit_yecap, nunit_vel, nunit_r, nunit_avtexp
  integer :: i

  integer :: nunit_out
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
    write(6,'(a)') trim(fn_eos)
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
    write(6,'(a)') trim(fn_caprate)
    call get_string_parameter("parameters","fn_caprate_neg",fn_caprate_neg)
    call get_string_parameter("parameters","fn_etanu",fn_etanu)

    
  end block
  fn_enu="/data/scratch/sfujibayashi/EOS/enu_to_chnu2"
  fn_ynu="/data/scratch/sfujibayashi/EOS/ynu_to_chnu"
  
  call init_weak_table(fn_caprate, fn_caprate_neg)
  
  call etanu_init(fn_etanu)

  block
    use module_eos
    ! call readeos(fn_eos,fn_eosb,fn_ynu,fn_enu)
    call readeos(fn_eos,nrho_in, ntemp_in, nye_in)
  end block


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
  allocate( caprate_temp_traj_list(ntemp_list,np),absrate_temp_traj_list(ntemp_list,np),yecap_temp_traj_list(ntemp_list,np),yeabs_temp_traj_list(ntemp_list,np),yemu0_temp_traj_list(ntemp_list,np),avtexp_temp_traj_list(ntemp_list,np))


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
     yn_p(:) = 0d0
     ya_p(:) = 0d0
     fn = trim(dir_read)//"/misc_"//trim(str1)//".dat"
     inquire(file=fn, exist=file_exists)
     if(file_exists)then
        
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
     endif

     fn = "./weak/weak_"//trim(str1)//".dat"
     open(newunit=nunit_out, file=fn, status="replace", action="write")
     write(nunit_out, '("#",99a15)') "t", "rho", "T", "Ye", "eta", "Xn", "Xp", "lambda_ec*Xp", "lambda_pc*Xn", "r", "v^r", "t_exp", "Ye(eq,cap)"
     do it=1,nt
        den_tmp = qrho_p(it)
        tem_tmp = tem_p(it)/tem_uni
        ye_tmp = ye_p(it)
        
        call ye_equilibrium_capture(den_tmp,tem_tmp,ye_equil_cap,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)
        call nrate_cap(den_tmp,tem_tmp,ye_tmp,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)

        r_tmp =  sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2)
        vr_tmp= (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/r_tmp
        t_exp = r_tmp/abs(vr_tmp)

        write(nunit_out, '(" ",99es15.6e3)') time(it), den_tmp, tem_tmp, ye_tmp, eta_dummy, xn_dummy, xp_dummy, ecap_nrate, pcap_nrate, r_tmp, vr_tmp, t_exp, ye_equil_cap
     enddo
     close(nunit_out)
  enddo
  
end program weak
