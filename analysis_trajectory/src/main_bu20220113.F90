program main
#include "macro.h"
  use unit
  use module_EOS_sekig
  use module_weak_interaction
  implicit none

  real(8),parameter :: emev=0.51099996d0
  real(8),parameter :: mev2t9= 1.160445d1

  character(256) :: fn_eos,fn_eosb,fn_enu,fn_ynu,dir_read,fn,fn_out,label
  character(10) :: str1

  integer :: itt_min,itt_max,np,ip,it,nt
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


  integer :: itemp,ntemp
  real(8) :: temp_min,temp_max
  integer,allocatable :: np_temp_list(:)
  real(8),allocatable :: temp_list(:),ye_av_list(:),mass_temp_list(:),ye2_av_list(:),eta_av_list(:),eta2_av_list(:),sen_av_list(:),sen2_av_list(:),texp_av_list(:),texp2_av_list(:)
  real(8),allocatable :: caprate_av_list(:),caprate2_av_list(:),absrate_av_list(:),absrate2_av_list(:),yecap_av_list(:),yecap2_av_list(:),yeabs_av_list(:),yeabs2_av_list(:)
  
  real(8) :: hoge
  integer :: ip_current, np_active, ipp, ip_dummy
  !integer,allocatable :: np_active(:)
  real(8),allocatable :: ye_temp_traj_list(:,:),mass_temp_traj_list(:,:),eta_temp_traj_list(:,:),sen_temp_traj_list(:,:),texp_temp_traj_list(:,:)
  real(8),allocatable :: caprate_temp_traj_list(:,:),absrate_temp_traj_list(:,:),yecap_temp_traj_list(:,:),yeabs_temp_traj_list(:,:)
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

  real(8) :: tt,tt1,ye_tmp,den_tmp,tem_tmp
  
  integer :: ke,ke1,ie,ie1,je,je1
  real(8) :: uu,uup,ss,ssp,ttp,che_tmp,eta_tmp,sen_tmp,vr_tmp,r_tmp,texp_tmp,t_tmp

  real(8) :: ye_equil_cap,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate
  real(8) :: ye_equil_abs,abs_n_nrate,abs_a_nrate,rne_tmp,rae_tmp,taun_tmp,taua_tmp

  logical :: mass_weighted


  real(8) :: my_time
  
!!! 3D models
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6/data_post"; fn_out="ye_ave_q6.dat"
  ! dir_read="/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q4/data_post"; label="q4"; mass_weighted=.true.
  dir_read="/scratch/sfujibayashi/Particle_trace_data/data_3D_SFHo_135_135-150mstg/data"; label="sfho135135dyn"; mass_weighted=.true.
  
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_SFHo_12_15_150stg/data";label="sfho1215_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_SFHo_130_140-150mstg/data";label="sfho1314_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_SFHo_135_135-150mstg/data";label="sfho135135_reduced"; mass_weighted=.true.
  !dir_read="/scratch/sfujibayashi/Particle_trace_data/data_2D_DD2_125M/data";label="DD2-125M_number"; mass_weighted=.false.

  call init_weak_table("caprate.dat")


  ntemp=50
  temp_min = 0.1d0
  temp_max = 5.d0
  allocate(temp_list(ntemp),np_temp_list(ntemp),ye_av_list(ntemp),ye2_av_list(ntemp),mass_temp_list(ntemp),eta_av_list(ntemp),eta2_av_list(ntemp),sen_av_list(ntemp),sen2_av_list(ntemp), &
       texp_av_list(ntemp),texp2_av_list(ntemp))

  allocate(caprate_av_list(ntemp),caprate2_av_list(ntemp),absrate_av_list(ntemp),absrate2_av_list(ntemp),yecap_av_list(ntemp),yecap2_av_list(ntemp),yeabs_av_list(ntemp),yeabs2_av_list(ntemp))

  do itemp=1,ntemp
     temp_list(itemp) = 10d0**(log10(temp_min) + (log10(temp_max)-log10(temp_min))*dble(itemp-1)/dble(ntemp-1))
     ! write(6,*) itemp,temp_list(itemp)
  enddo

#ifdef DD2
  write(6,*) "DD2 EOS used"
  fn_eos ="/scratch/sfujibayashi/EOS/EOS_Hempel_DD2Tim_TF_326"
  fn_eosb="/scratch/sfujibayashi/EOS/EOS_Hempel_DD2Tim_TFB"
#endif
#ifdef SFHo
  write(6,*) "SFHo EOS used"
  fn_eos ="/scratch/sfujibayashi/EOS/EOS_Hempel_SFHoTim_TF_326"
  fn_eosb="/scratch/sfujibayashi/EOS/EOS_Hempel_SFHoTim_TFB"
#endif
  fn_enu="/scratch/sfujibayashi/EOS/enu_to_chnu2"
  fn_ynu="/scratch/sfujibayashi/EOS/ynu_to_chnu"

  call readeos(fn_eos,fn_eosb,fn_ynu,fn_enu)

  ! write(*,*) "model = ", trim(cmodel)
  open(10,file=trim(dir_read)//"/report_ptr.dat",status="old")
  read(10,*)
  read(10,*) itt_min
  read(10,*)
  read(10,*) itt_max
  read(10,*)
  !read(10,*); read(10,*);read(10,*); read(10,*)
  read(10,*) np


  ! allocate(np_active(ntemp))
  allocate( ye_temp_traj_list(ntemp,np),mass_temp_traj_list(ntemp,np),eta_temp_traj_list(ntemp,np),sen_temp_traj_list(ntemp,np),texp_temp_traj_list(ntemp,np) )
  allocate( caprate_temp_traj_list(ntemp,np),absrate_temp_traj_list(ntemp,np),yecap_temp_traj_list(ntemp,np),yeabs_temp_traj_list(ntemp,np) )


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

  np_temp_list(:) = 0
  mass_temp_list(:) = 0.d0
  ye_av_list(:) = 0.d0
  ye2_av_list(:) = 0.d0
  eta_av_list(:) = 0.d0
  eta2_av_list(:) = 0.d0
  sen_av_list(:) = 0.d0
  sen2_av_list(:) = 0.d0

  ip_current = 1
  !np_active(:) = 0

  do ip=1,np,10
     
     write(str1,'(i8.8)') ip
     fn = trim(dir_read)//"/traj_"//trim(str1)//".dat"
     write(6,'(a)') trim(fn)
     open(12,file=fn,status="old")
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
        !stop
     enddo
99   continue
     close(12)

     if(.not.mass_weighted) mass_p(ip) = 1d0
     ! write(6,*)nt
     ! do it=1,nt
     !    write(6,*) it,time(it),tem_p(it)/tem_uni
     ! enddo
     ! stop

     
     !stop
     ! if(time(nt)>0.5d0.and. &
     if( &
          !sen_p(nt)<30.d0 .and. &
          maxval(tem_p(1:nt))/tem_uni>temp_list(ntemp).and. &
          minval(tem_p(1:nt))/tem_uni<temp_list(1).and. &
          !maxval(tem_p(1:nt))>1d10.and. &
          !ye_p(nt)>0.1d0 &
          1==1 &
        )then

        it=nt
        do while(tem_p(it) > minval(tem_p(1:nt)))
           it = it - 1
        enddo
        if(it==nt) it = nt-1

        loop_itemp:do itemp=1,ntemp

           ! if( temp_list(itemp) < minval(tem_p(1:nt)) )then
           !    goto 999
           ! endif
           
           do while(tem_p(it)/tem_uni < temp_list(itemp))
              it = it - 1
              if(it==0) exit loop_itemp
           enddo
           tt1= (temp_list(itemp)*tem_uni-tem_p(it))/(tem_p(it+1)-tem_p(it))
           tt = 1.d0-tt1
           
           ! write(6,*) it, nt, temp_list(itemp), tem_p(it)/tem_uni,  tem_p(it+1)/tem_uni
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

           call nrate_cap(den_tmp,tem_tmp,ye_tmp,eta_dummy,xn_dummy,xp_dummy,ecap_nrate,pcap_nrate)
           call nrate_abs(den_tmp,tem_tmp,ye_tmp,rne_tmp,rae_tmp,xn_dummy,xp_dummy,abs_n_nrate,abs_a_nrate)

           !abs_n_nrate = abs_n_nrate *exp(-2d0*taun_tmp)
           !abs_a_nrate = abs_a_nrate *exp(-2d0*taua_tmp)

           ecap_nrate = ecap_nrate/xp_dummy
           pcap_nrate = pcap_nrate/xn_dummy
           abs_n_nrate = abs_n_nrate/xn_dummy
           abs_a_nrate = abs_a_nrate/xp_dummy
           
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

           ye_temp_traj_list(itemp,ip_current) = ye_tmp
           mass_temp_traj_list(itemp,ip_current) = mass_p(ip)
           eta_temp_traj_list(itemp,ip_current) = eta_tmp
           sen_temp_traj_list(itemp,ip_current) = sen_tmp
           texp_temp_traj_list(itemp,ip_current) = texp_tmp
           
           caprate_temp_traj_list(itemp,ip_current) = max(ecap_nrate,pcap_nrate)
           absrate_temp_traj_list(itemp,ip_current) = max(abs_n_nrate,abs_a_nrate)
           yecap_temp_traj_list(itemp,ip_current) = ye_equil_cap
           yeabs_temp_traj_list(itemp,ip_current) = ye_equil_abs

        enddo loop_itemp
        ! exit
        
        write(6,*) ip,np,ip_current
        ip_current = ip_current + 1
        !stop
     endif
  

  enddo
  
  ! sort and obtain median
  np_active = ip_current - 1
  write(6,*) "# of tracers considered: ", np_active
  allocate( ye_lw_temp_list(ntemp),ye_med_temp_list(ntemp),ye_rw_temp_list(ntemp))
  allocate( eta_lw_temp_list(ntemp),eta_med_temp_list(ntemp),eta_rw_temp_list(ntemp))
  allocate( sen_lw_temp_list(ntemp),sen_med_temp_list(ntemp),sen_rw_temp_list(ntemp))
  allocate( texp_lw_temp_list(ntemp),texp_med_temp_list(ntemp),texp_rw_temp_list(ntemp))
  allocate( caprate_lw_temp_list(ntemp),caprate_med_temp_list(ntemp),caprate_rw_temp_list(ntemp))
  allocate( absrate_lw_temp_list(ntemp),absrate_med_temp_list(ntemp),absrate_rw_temp_list(ntemp))
  allocate( yecap_lw_temp_list(ntemp),yecap_med_temp_list(ntemp),yecap_rw_temp_list(ntemp))
  allocate( yeabs_lw_temp_list(ntemp),yeabs_med_temp_list(ntemp),yeabs_rw_temp_list(ntemp))
  
  call sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,ye_temp_traj_list,mass_temp_traj_list,ye_lw_temp_list,ye_med_temp_list,ye_rw_temp_list)
  call sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,texp_temp_traj_list,mass_temp_traj_list,texp_lw_temp_list,texp_med_temp_list,texp_rw_temp_list)
  call sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,caprate_temp_traj_list,mass_temp_traj_list,caprate_lw_temp_list,caprate_med_temp_list,caprate_rw_temp_list)
  call sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,absrate_temp_traj_list,mass_temp_traj_list,absrate_lw_temp_list,absrate_med_temp_list,absrate_rw_temp_list)
  call sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,yecap_temp_traj_list,mass_temp_traj_list,yecap_lw_temp_list,yecap_med_temp_list,yecap_rw_temp_list)
  call sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,yeabs_temp_traj_list,mass_temp_traj_list,yeabs_lw_temp_list,yeabs_med_temp_list,yeabs_rw_temp_list)

  open(12,file="median_"//trim(label)//".dat",status="replace")
  write(12,'("# ",99es12.4)') fac_lw,fac_med,fac_rw
  do itemp=1,ntemp
     write(12,'(99es12.4)') temp_list(itemp), &
          ye_lw_temp_list(itemp),ye_med_temp_list(itemp),ye_rw_temp_list(itemp), &
          texp_lw_temp_list(itemp),texp_med_temp_list(itemp),texp_rw_temp_list(itemp), &
          caprate_lw_temp_list(itemp),caprate_med_temp_list(itemp),caprate_rw_temp_list(itemp), &
          absrate_lw_temp_list(itemp),absrate_med_temp_list(itemp),absrate_rw_temp_list(itemp), &
          yecap_lw_temp_list(itemp),yecap_med_temp_list(itemp),yecap_rw_temp_list(itemp), &
          yeabs_lw_temp_list(itemp),yeabs_med_temp_list(itemp),yeabs_rw_temp_list(itemp)
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
