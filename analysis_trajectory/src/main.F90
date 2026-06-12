program main
#include "macro.h"
  use unit
  use module_EOS_sekig
  use module_weak_interaction
  use module_etanu

  implicit none

  real(8),parameter :: emev=0.51099996d0
  real(8),parameter :: mev2t9= 1.160445d1

  character(256) :: fn_eos,fn_eosb,fn_enu,fn_ynu,dir_read,fn,fn_out,label
  character(10) :: str1

  character(256) :: fn_caprate, fn_caprate_neg, fn_etanu

  integer :: itt_min,itt_max,np,ip,it,nt,np_skip
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

  integer :: itemp,ntemp
  real(8) :: temp_min,temp_max
  integer,allocatable :: np_temp_list(:)
  real(8),allocatable :: temp_list(:),ye_av_list(:),mass_temp_list(:),ye2_av_list(:),eta_av_list(:),eta2_av_list(:),sen_av_list(:),sen2_av_list(:),texp_av_list(:),texp2_av_list(:)
  real(8),allocatable :: caprate_av_list(:),caprate2_av_list(:),absrate_av_list(:),absrate2_av_list(:),yecap_av_list(:),yecap2_av_list(:),yeabs_av_list(:),yeabs2_av_list(:)
  
  real(8) :: hoge
  integer :: ip_current
  integer,allocatable :: np_active(:)
  real(8),allocatable :: ye_temp_traj_list(:,:),mass_temp_traj_list(:,:),eta_temp_traj_list(:,:),sen_temp_traj_list(:,:),texp_temp_traj_list(:,:)
  real(8),allocatable :: caprate_temp_traj_list(:,:),absrate_temp_traj_list(:,:),yecap_temp_traj_list(:,:),yeabs_temp_traj_list(:,:),yemu0_temp_traj_list(:,:)
  integer,allocatable :: ip_list(:)
  real(8),allocatable :: dummy_list(:)

  integer :: ip_lw, ip_med, ip_rw
  real(8) :: fac_lw=0.15d0, fac_med=0.5d0, fac_rw=0.85d0
  real(8) :: ss_lw, ss_med, ss_rw
  real(8),allocatable :: ye_lw_temp_list(:),ye_med_temp_list(:),ye_rw_temp_list(:)
  real(8),allocatable :: eta_lw_temp_list(:),eta_med_temp_list(:),eta_rw_temp_list(:)
  real(8),allocatable :: sen_lw_temp_list(:),sen_med_temp_list(:),sen_rw_temp_list(:)
  real(8),allocatable :: texp_lw_temp_list(:),texp_med_temp_list(:),texp_rw_temp_list(:)
  real(8),allocatable :: caprate_lw_temp_list(:),caprate_med_temp_list(:),caprate_rw_temp_list(:)
  real(8),allocatable :: absrate_lw_temp_list(:),absrate_med_temp_list(:),absrate_rw_temp_list(:)
  real(8),allocatable :: yecap_lw_temp_list(:),yecap_med_temp_list(:),yecap_rw_temp_list(:)
  real(8),allocatable :: yeabs_lw_temp_list(:),yeabs_med_temp_list(:),yeabs_rw_temp_list(:)
  real(8),allocatable :: yemu0_lw_temp_list(:),yemu0_med_temp_list(:),yemu0_rw_temp_list(:)

  real(8) :: tt,tt1,ye_tmp,den_tmp,tem_tmp, yn_tmp, ya_tmp
  
  integer :: ke,ke1,ie,ie1,je,je1
  real(8) :: uu,uup,ss,ssp,ttp,che_tmp,eta_tmp,sen_tmp,vr_tmp,r_tmp,texp_tmp,t_tmp

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

  integer :: np_analized

  !call make_epcap_table(1d-2,1d3,1d-2,1d2,401,401)
  !stop

  block
    use inputparser
    call get_string_parameter("parameters","dir_read",dir_read)
    call get_string_parameter("parameters","label",label)
    call get_logical_parameter("parameters","mass_weighted",mass_weighted)
    call get_string_parameter("parameters","fn_eos",fn_eos)
    call get_string_parameter("parameters","fn_eosb",fn_eosb)
    write(6,*) fn_eosb
    call get_integer_parameter("parameters","ntemp",ntemp)
    call get_double_parameter("parameters","temp_min",temp_min)
    call get_double_parameter("parameters","temp_max",temp_max)
    write(6,*) ntemp, temp_min, temp_max
    itt_min = 1
    call get_integer_parameter("parameters","itt_max",itt_max)
    call get_integer_parameter("parameters","np",np)
    call get_integer_parameter("parameters","np_skip",np_skip)
    write(6,*) itt_min, itt_max, np, np_skip
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

  
  ! ntemp=20
  ! temp_min = 0.4d0
  ! temp_max = 3.d0
  allocate(temp_list(ntemp),np_temp_list(ntemp),ye_av_list(ntemp),ye2_av_list(ntemp),mass_temp_list(ntemp),eta_av_list(ntemp),eta2_av_list(ntemp),sen_av_list(ntemp),sen2_av_list(ntemp), &
       texp_av_list(ntemp),texp2_av_list(ntemp))

  allocate(caprate_av_list(ntemp),caprate2_av_list(ntemp),absrate_av_list(ntemp),absrate2_av_list(ntemp),yecap_av_list(ntemp),yecap2_av_list(ntemp),yeabs_av_list(ntemp),yeabs2_av_list(ntemp))

  do itemp=1,ntemp
     temp_list(itemp) = 10d0**(log10(temp_min) + (log10(temp_max)-log10(temp_min))*dble(itemp-1)/dble(ntemp-1))
     ! write(6,*) itemp,temp_list(itemp)
  enddo

  call readeos(fn_eos,fn_eosb,fn_ynu,fn_enu)

  ! ! write(*,*) "model = ", trim(cmodel)
  ! open(10,file=trim(dir_read)//"/report_ptr.dat",status="old",action="read")
  ! read(10,*); read(10,*) itt_min
  ! read(10,*); read(10,*) itt_max
  ! read(10,*); read(10,*)
  ! read(10,*); read(10,*)
  ! !read(10,*); read(10,*);read(10,*); read(10,*)
  ! read(10,*);read(10,*) np
  ! close(10)

  allocate(ut1(np),hut(np))
  open(11,file=trim(dir_read)//"/ana_traj.dat",status="old",action="read")
  read(11,*);read(11,*)
  do ip=1,np
     read(11,*) buf(1:30)
     ut1(ip) = buf(16)
     hut(ip) = buf(17)
  enddo
  close(11)

  allocate(np_active(ntemp))
  allocate( ye_temp_traj_list(ntemp,np),mass_temp_traj_list(ntemp,np),eta_temp_traj_list(ntemp,np),sen_temp_traj_list(ntemp,np),texp_temp_traj_list(ntemp,np) )
  allocate( caprate_temp_traj_list(ntemp,np),absrate_temp_traj_list(ntemp,np),yecap_temp_traj_list(ntemp,np),yeabs_temp_traj_list(ntemp,np),yemu0_temp_traj_list(ntemp,np) )

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

  np_temp_list(:) = 0
  mass_temp_list(:) = 0.d0
  ye_av_list(:) = 0.d0
  ye2_av_list(:) = 0.d0
  eta_av_list(:) = 0.d0
  eta2_av_list(:) = 0.d0
  sen_av_list(:) = 0.d0
  sen2_av_list(:) = 0.d0

  ! ip_current = 1
  np_active(:) = 0
  
  do ip=1,np,np_skip
     !write(6,*) ip
     !if(ut1(ip)>1d0) continue
     
     write(str1,'(i8.8)') ip
     fn = trim(dir_read)//"/traj_"//trim(str1)//".dat"
     write(6,'(a)') trim(fn)
     open(12,file=fn,status="old",action="read")
     read(12,*)
     read(12,*)
     read(12,'(16x,es13.5)') mass_p(ip)
     read(12,*)
     it = 0
     do
        nt = it
        it = it + 1
        read(12,*,end=99) &
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

        !write(6,'(i8,99es12.4)') it, time(it), qrho_p(it), tem_p(it), ye_p(it), vlx_p(it)
     enddo
99   continue
     close(12)

     !write(6,*) "file read done"
     
!!! misc. file
     fn = trim(dir_read)//"/misc_"//trim(str1)//".dat"
     if(access(fn," ")==0)then
        
        open(12,file=fn,status="old")
        read(12,*)
        read(12,*)
        read(12,*)
        read(12,*)
        it = 0
        do
           it = it + 1
           read(12,*,end=991) &
                hoge, &
                yn_p(it), &
                ya_p(it)
           
        enddo
991     continue
        close(12)
     else
        yn_p(:) = 0d0
        ya_p(:) = 0d0
        
     endif

     if(.not.mass_weighted) mass_p(ip) = 1d0

     !write(6,*) ip,maxval(tem_p(1:nt))/tem_uni, minval(tem_p(1:nt))/tem_uni
     !stop
     ! if(time(nt)>0.5d0.and. &
     if( &
          !ut1(ip)<0d0.and. &
          !sen_p(nt)<30.d0 .and. &
          !maxval(tem_p(1:nt))>=temp_list(ntemp)*tem_uni.and. &
          !minval(tem_p(1:nt))<=temp_list(1)*tem_uni.and. &
          maxval(tem_p(1:nt))>1d0*tem_uni.and. &
          ! tem_p(nt)<=temp_list(1)*tem_uni.and. &
          !ye_p(nt)<0.15d0.and.  &
          !ye_p(nt)>0.15d0.and.  &
          .true. &
        )then

        it=nt
        ! do while(tem_p(it) > minval(tem_p(1:nt)))
        !    it = it - 1
        ! enddo
        if(it==nt) it = nt-1

        loop_itemp:do itemp=1,ntemp

           !write(6,'(99es12.4)') temp_list(itemp)*tem_uni, minval(tem_p(1:nt))
           !if( temp_list(itemp)*tem_uni < minval(tem_p(1:nt)) )then
           if( temp_list(itemp)*tem_uni < tem_p(nt) )then
              !write(6,*) "skipped"
              goto 999
           endif
           
           np_active(itemp) = np_active(itemp) + 1
           
           search_it:do
              !write(6,'(i5,99es12.4)') it, sqrt(x_p(it)**2 + y_p(it)**2 + z_p(it)**2), tem_p(it)-temp_list(itemp)*tem_uni, tem_p(it+1)-temp_list(itemp)*tem_uni
              if(tem_p(it) > temp_list(itemp)*tem_uni .and. temp_list(itemp)*tem_uni >= tem_p(it+1)) exit search_it
              !if( (tem_p(it)-temp_list(itemp)*tem_uni)*(tem_p(it+1)-temp_list(itemp)*tem_uni) < 0d0)exit search_it
              it = it - 1
              if(it==0) exit loop_itemp
              !if( abs(sen_p(it+1)-sen_p(it))/sen_p(it+1) >0.5d0 ) exit loop_itemp
           enddo search_it
           tt1= (temp_list(itemp)*tem_uni-tem_p(it))/(tem_p(it+1)-tem_p(it))
           tt = 1.d0-tt1
           
           !write(6,*) itemp, it, nt, tt1, tem_p(it+1)/tem_uni, temp_list(itemp), tem_p(it)/tem_uni
           !stop
! !!! specify time
!            my_time = 1d0
!            it=1
!            do while(time(it+1)<my_time.and.it+1<nt)
!               it = it + 1
!            enddo
!            tt1= (my_time - time(it))/(time(it+1) - time(it))
!            tt = 1.d0-tt1

           np_temp_list(itemp) = np_temp_list(itemp) + 1
           mass_temp_list(itemp) = mass_temp_list(itemp) + mass_p(ip)
           
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

           vr_tmp = tt *(x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2) &
                   +tt1*(x_p(it+1)*vlx_p(it+1) + y_p(it+1)*vly_p(it+1) + z_p(it+1)*vlz_p(it+1))/sqrt(x_p(it+1)**2+y_p(it+1)**2+z_p(it+1)**2)
           r_tmp =  tt *sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2) &
                   +tt1*sqrt(x_p(it+1)**2+y_p(it+1)**2+z_p(it+1)**2)
           texp_tmp=r_tmp/abs(vr_tmp)
           

           ke  = max(ked , min(keu-1, int((log10(den_tmp)-rho_e_min  )*drhoi)+1))
           ke1 = ke+1
           uu    = max(0.d0, min(1.d0 ,     (log10(den_tmp)-rho_e(ke))*drhoi)   )
           uup   = 1.d0-uu
           
           ie = max(ied , min(ieu-1, int((log10(tem_tmp)-tem_e_min  )*dtemi)+1))
           ie1=ie+1
           ss    = max(0.d0, min(1.d0 ,     (log10(tem_tmp)-tem_e(ie))*dtemi))
           ssp   = 1.d0-ss
           
           je  = max(jed+1 , min(jeu-1, int((ye_tmp-ye_e_min )*dyei )+1))
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

           
           ye_av_list(itemp) = ye_av_list(itemp) + mass_p(ip)*ye_tmp
           ye2_av_list(itemp) = ye2_av_list(itemp) + mass_p(ip)*ye_tmp**2
           eta_av_list(itemp) = eta_av_list(itemp) + mass_p(ip)*eta_tmp
           eta2_av_list(itemp) = eta2_av_list(itemp) + mass_p(ip)*eta_tmp**2
           sen_av_list(itemp) = sen_av_list(itemp) + mass_p(ip)*sen_tmp
           sen2_av_list(itemp) = sen2_av_list(itemp) + mass_p(ip)*sen_tmp**2
           !ye_av_list(itemp) = ye_av_list(itemp) + ye_tmp
           !ye2_av_list(itemp) = ye2_av_list(itemp) + ye_tmp**2
           
           call ye_equilibrium_capture(den_tmp,tem_tmp,ye_equil_cap,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)
           call ye_equilibrium_abs(den_tmp,tem_tmp,rne_tmp,rae_tmp,ye_equil_abs,xn_dummy,xp_dummy,abs_n_nrate,abs_a_nrate)
           call ye_equilibrium_munu0(den_tmp,tem_tmp,ye_mu0)
        
           call nrate_cap(den_tmp,tem_tmp,ye_tmp,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)
           call nrate_abs(den_tmp,tem_tmp,ye_tmp,rne_tmp,rae_tmp,xn_dummy,xp_dummy,abs_n_nrate,abs_a_nrate)

           !
           !write(6,'(99es12.4)') 
           call etanu_eta_ene(den_tmp*yn_tmp,tem_tmp, eta_n, rnetrap_tmp)
           call etanu_eta_ene(den_tmp*ya_tmp,tem_tmp, eta_a, raetrap_tmp)
           call nrate_abs(den_tmp,tem_tmp,ye_tmp,rnetrap_tmp,raetrap_tmp,xn_dummy,xp_dummy,abs_nt_nrate,abs_at_nrate)
           
           call nrate_cap_block(den_tmp,tem_tmp,ye_tmp,eta_n, eta_a, eta_dummy,xn_dummy,xp_dummy,ecap_nrate_block,pcap_nrate_block)
           !write(6,'(99es12.4)')  tem_tmp,den_tmp, taun_tmp, taua_tmp, yn_tmp, ya_tmp, eta_n, eta_a, abs_nt_nrate, abs_n_nrate, ecap_nrate_block, ecap_nrate, pcap_nrate_block, pcap_nrate
           !abs_n_nrate = abs_n_nrate *exp(-2d0*taun_tmp)
           !abs_a_nrate = abs_a_nrate *exp(-2d0*taua_tmp)
           abs_n_nrate = abs_n_nrate + abs_nt_nrate
           abs_a_nrate = abs_a_nrate + abs_at_nrate

           ecap_nrate = ecap_nrate/xp_dummy
           pcap_nrate = pcap_nrate/xn_dummy
           abs_n_nrate = abs_n_nrate/xn_dummy
           abs_a_nrate = abs_a_nrate/xp_dummy

           ! if(itemp==1.and. (ecap_nrate >1d0 .or.pcap_nrate >1d0 ))then
           ! !if(itemp==1.and. (ecap_nrate < 1d-5 .or.pcap_nrate < 1d-5 ))then
           !    write(6,'(99es12.4)') ecap_nrate, pcap_nrate, den_tmp, ye_tmp, tem_tmp,eta_tmp
           !    stop
           ! endif
           
           caprate_av_list(itemp)  = caprate_av_list(itemp)  + mass_p(ip)*log10(max(ecap_nrate,pcap_nrate))
           caprate2_av_list(itemp) = caprate2_av_list(itemp) + mass_p(ip)*log10(max(ecap_nrate,pcap_nrate))**2
           absrate_av_list(itemp)  = absrate_av_list(itemp)  + mass_p(ip)*log10(max(abs_n_nrate,abs_a_nrate))
           absrate2_av_list(itemp) = absrate2_av_list(itemp) + mass_p(ip)*log10(max(abs_n_nrate,abs_a_nrate))**2

           yecap_av_list(itemp) = yecap_av_list(itemp) + mass_p(ip)*ye_equil_cap
           yecap2_av_list(itemp) = yecap2_av_list(itemp) + mass_p(ip)*ye_equil_cap**2
           yeabs_av_list(itemp) = yeabs_av_list(itemp) + mass_p(ip)*ye_equil_abs
           yeabs2_av_list(itemp) = yeabs2_av_list(itemp) + mass_p(ip)*ye_equil_abs**2

           texp_av_list(itemp) = texp_av_list(itemp) + mass_p(ip)*texp_tmp
           texp2_av_list(itemp) = texp2_av_list(itemp) + mass_p(ip)*texp_tmp**2
           !
           !if(itemp==ntemp) write(6,'(99es12.4)') den_tmp, tem_tmp, ye_tmp, eta_tmp, taun_tmp, taua_tmp
           !write(6,'(99es12.4)') r_tmp, rne_tmp, rae_tmp, abs_n_nrate, abs_a_nrate, tem_tmp, eta_tmp, ye_tmp, ecap_nrate, pcap_nrate

           ip_current = np_active(itemp)
           
           ye_temp_traj_list(itemp,ip_current) = ye_tmp
           mass_temp_traj_list(itemp,ip_current) = mass_p(ip)
           eta_temp_traj_list(itemp,ip_current) = eta_tmp
           sen_temp_traj_list(itemp,ip_current) = sen_tmp
           texp_temp_traj_list(itemp,ip_current) = texp_tmp
           
           caprate_temp_traj_list(itemp,ip_current) = max(ecap_nrate,pcap_nrate)
           absrate_temp_traj_list(itemp,ip_current) = max(abs_n_nrate,abs_a_nrate)
           yecap_temp_traj_list(itemp,ip_current) = ye_equil_cap
           yeabs_temp_traj_list(itemp,ip_current) = ye_equil_abs
           yemu0_temp_traj_list(itemp,ip_current) = ye_mu0

999        continue
           
        enddo loop_itemp
        ! exit
        
        write(6,'(99i6)') ip,np,np_active(:)
        !ip_current = ip_current + 1
        !if(np_active(ntemp)>0)stop
        !stop
     endif
     
  enddo
  
  ! sort and obtain median
  ! np_active = ip_current - 1
  write(6,*) "# of tracers considered: "
  do itemp=1,ntemp
     write(6,*) temp_list(itemp), np_active(itemp)
  enddo
  !stop

  ! output

  allocate( ye_lw_temp_list(ntemp),ye_med_temp_list(ntemp),ye_rw_temp_list(ntemp))
  allocate( eta_lw_temp_list(ntemp),eta_med_temp_list(ntemp),eta_rw_temp_list(ntemp))
  allocate( sen_lw_temp_list(ntemp),sen_med_temp_list(ntemp),sen_rw_temp_list(ntemp))
  allocate( texp_lw_temp_list(ntemp),texp_med_temp_list(ntemp),texp_rw_temp_list(ntemp))
  allocate( caprate_lw_temp_list(ntemp),caprate_med_temp_list(ntemp),caprate_rw_temp_list(ntemp))
  allocate( absrate_lw_temp_list(ntemp),absrate_med_temp_list(ntemp),absrate_rw_temp_list(ntemp))
  allocate( yecap_lw_temp_list(ntemp),yecap_med_temp_list(ntemp),yecap_rw_temp_list(ntemp))
  allocate( yeabs_lw_temp_list(ntemp),yeabs_med_temp_list(ntemp),yeabs_rw_temp_list(ntemp))
  allocate( yemu0_lw_temp_list(ntemp),yemu0_med_temp_list(ntemp),yemu0_rw_temp_list(ntemp))
  
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,ye_temp_traj_list,mass_temp_traj_list,ye_lw_temp_list,ye_med_temp_list,ye_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,eta_temp_traj_list,mass_temp_traj_list,eta_lw_temp_list,eta_med_temp_list,eta_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,sen_temp_traj_list,mass_temp_traj_list,sen_lw_temp_list,sen_med_temp_list,sen_rw_temp_list)

  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,texp_temp_traj_list,mass_temp_traj_list,texp_lw_temp_list,texp_med_temp_list,texp_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,caprate_temp_traj_list,mass_temp_traj_list,caprate_lw_temp_list,caprate_med_temp_list,caprate_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,absrate_temp_traj_list,mass_temp_traj_list,absrate_lw_temp_list,absrate_med_temp_list,absrate_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,yecap_temp_traj_list,mass_temp_traj_list,yecap_lw_temp_list,yecap_med_temp_list,yecap_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,yeabs_temp_traj_list,mass_temp_traj_list,yeabs_lw_temp_list,yeabs_med_temp_list,yeabs_rw_temp_list)
  call sort_median(ntemp,np,np_active,fac_lw,fac_med,fac_rw,yemu0_temp_traj_list,mass_temp_traj_list,yemu0_lw_temp_list,yemu0_med_temp_list,yemu0_rw_temp_list)

  open(12,file="median_"//trim(label)//".dat",status="replace")
  write(12,'("#",99es12.4)') fac_lw,fac_med,fac_rw
  write(12,'("#",99a15)')  "T",&
       "Ye(L)", "Ye(med)", "Ye(R)",&
       "texp(L)", "texp(med)", "texp(R)",&
       "Rcap(L)", "Rcap(med)", "Rcap(R)",&
       "Rabs(L)", "Rabs(med)", "Rabs(R)",&
       "Ye(eq,cap)(L)", "Ye(eq,cap)(med)", "Ye(eq,cap)(R)",&
       "Ye(eq,abs)(L)", "Ye(eq,abs)(med)", "Ye(eq,abs)(R)",&
       "Ye(mu=0)(L)", "Ye(mu=0)(med)", "Ye(mu=0)(R)",&
       "eta(L)", "eta(med)", "eta(R)",&
       "entr(L)", "entr(med)", "entr(R)"
  do itemp=1,ntemp
     write(12,'(" ",99es15.7)') temp_list(itemp), &
          ye_lw_temp_list(itemp),ye_med_temp_list(itemp),ye_rw_temp_list(itemp), &
          texp_lw_temp_list(itemp),texp_med_temp_list(itemp),texp_rw_temp_list(itemp), &
          caprate_lw_temp_list(itemp),caprate_med_temp_list(itemp),caprate_rw_temp_list(itemp), &
          absrate_lw_temp_list(itemp),absrate_med_temp_list(itemp),absrate_rw_temp_list(itemp), &
          yecap_lw_temp_list(itemp),yecap_med_temp_list(itemp),yecap_rw_temp_list(itemp), &
          yeabs_lw_temp_list(itemp),yeabs_med_temp_list(itemp),yeabs_rw_temp_list(itemp), &
          yemu0_lw_temp_list(itemp),yemu0_med_temp_list(itemp),yemu0_rw_temp_list(itemp), &
          eta_lw_temp_list(itemp),eta_med_temp_list(itemp),eta_rw_temp_list(itemp), &
          sen_lw_temp_list(itemp),sen_med_temp_list(itemp),sen_rw_temp_list(itemp)
  enddo
  close(12)
  stop
  
     
  do itemp=1,ntemp
     ye_av_list(itemp) =ye_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     ye2_av_list(itemp)=ye2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     eta_av_list(itemp) =eta_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     eta2_av_list(itemp)=eta2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     sen_av_list(itemp) =sen_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     sen2_av_list(itemp)=sen2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     caprate_av_list(itemp) =caprate_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     caprate2_av_list(itemp)=caprate2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     absrate_av_list(itemp) =absrate_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     absrate2_av_list(itemp)=absrate2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     texp_av_list(itemp) =texp_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     texp2_av_list(itemp)=texp2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     yecap_av_list(itemp) =yecap_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     yecap2_av_list(itemp)=yecap2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     yeabs_av_list(itemp) =yeabs_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
     yeabs2_av_list(itemp)=yeabs2_av_list(itemp)/(mass_temp_list(itemp)+1.d-99)
  enddo

  ! caprate_av_list(:) =1d1**caprate_av_list(:)
  ! caprate2_av_list(:)=1d1**caprate2_av_list(:)
  ! absrate_av_list(:) =1d1**absrate_av_list(:)
  ! absrate2_av_list(:)=1d1**absrate2_av_list(:)
     
  open(13,file="ye_average_"//trim(label)//".dat",status="replace")
  do itemp=1,ntemp
     ! write(6,'(es14.6,i5,99es14.6)') temp_list(itemp), np_temp_list(itemp), &
     !      ye_av_list(itemp), ye2_av_list(itemp) - ye_av_list(itemp)**2, &
     !      eta_av_list(itemp), eta2_av_list(itemp) - eta_av_list(itemp)**2,&
     !      sen_av_list(itemp), sen2_av_list(itemp) - sen_av_list(itemp)**2,&
     !      caprate_av_list(itemp), caprate2_av_list(itemp) - caprate_av_list(itemp)**2,&
     !      absrate_av_list(itemp), absrate2_av_list(itemp) - absrate_av_list(itemp)**2,&
     !      texp_av_list(itemp), texp2_av_list(itemp) - texp_av_list(itemp)**2,&
     !      mass_temp_list(itemp)

     write(13,'(es14.6,i5,99es14.6)') temp_list(itemp), np_temp_list(itemp), &
          ye_av_list(itemp), sqrt(ye2_av_list(itemp) - ye_av_list(itemp)**2), &
          eta_av_list(itemp), sqrt(eta2_av_list(itemp) - eta_av_list(itemp)**2),&
          sen_av_list(itemp), sqrt(sen2_av_list(itemp) - sen_av_list(itemp)**2),&
          1.d1**caprate_av_list(itemp), 1.d1**sqrt(caprate2_av_list(itemp) - caprate_av_list(itemp)**2),&
          1.d1**absrate_av_list(itemp), 1.d1**sqrt(absrate2_av_list(itemp) - absrate_av_list(itemp)**2),&
          texp_av_list(itemp), sqrt(texp2_av_list(itemp) - texp_av_list(itemp)**2),&
          yecap_av_list(itemp), sqrt(yecap2_av_list(itemp) - yecap_av_list(itemp)**2),&
          yeabs_av_list(itemp), sqrt(yeabs2_av_list(itemp) - yeabs_av_list(itemp)**2),&
          mass_temp_list(itemp)
  enddo
  close(13)

end program main
