module analysis

  private
  public :: tr_analysis

contains

  subroutine tr_analysis(model, dir_read, old_format, format_2d, time_map)
    !$ use omp_lib
    implicit none

    character(*),intent(in) :: model,dir_read
    logical,intent(in) :: old_format, format_2d
    real(8),intent(in),optional :: time_map

    integer :: np,np_inside,itt_max,itt_min  

    real(8),parameter :: tem_uni = 1.160445d10
    real(8),parameter :: v_uni = 2.99792458d10
    real(8),parameter :: rho_drip = 4d11

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

    integer :: i,ip,n,ibuf,it,it_max,it_init,it_5gk,it_10gk,it_1gk,it_3gk,it_tem_max, &
         iflag_tem05,iflag_tem07,iflag_tem10,iflag_sc70,iflag_sc74, &
         it_5gk_after_10gk,it_rho_max,it_drip

    real(8) :: smax_traj_7_0,smax_traj_7_4,smin_traj_7_0,smin_traj_7_4

    integer,allocatable :: n_cond(:)
    real(8),allocatable :: &
         mass_traj(:), &
         t_fin_traj(:), x_fin_traj(:), y_fin_traj(:), z_fin_traj(:), vr_fin_traj(:), s_fin_traj(:), ye_fin_traj(:), ut1_fin_traj(:), hut_fin_traj(:), &
         tem_max_traj(:), tem_max_af3gk_traj(:),t_tem_max_traj(:), t_tem_max_af3gk_traj(:), &
         t_5gk_traj(:), s_5gk_traj(:), ye_5gk_traj(:), texp_5gk_traj(:), &
         t_3gk_traj(:), t_1gk_traj(:), ye_10gk_traj(:), t_10gk_traj(:), &
         time_50gk_25gk_traj(:), rho_fin_traj(:), vx_fin_traj(:), vy_fin_traj(:), vz_fin_traj(:), r_ini_traj(:)
    real(8),allocatable :: x_ini_traj(:), y_ini_traj(:), z_ini_traj(:), rho_ini_traj(:), t_ini_traj(:), rho_max_traj(:), &
         s_rho_max_traj(:), texp_rho_max_traj(:), ye_rho_max_traj(:), &
         s_drip_traj(:), texp_drip_traj(:), ye_drip_traj(:), &
         s_tem_max_traj(:), texp_tem_max_traj(:), ye_tem_max_traj(:)

    real(8) :: &
         t_fin, x_fin, y_fin, z_fin, vr_fin, s_fin, ye_fin, &
         tem_max, tem_max_af3gk, t_tem_max_af3gk, t_tem_max, &
         t_5gk, s_5gk, ye_5gk, texp_5gk, &
         t_3gk, t_1gk, ye_10gk, t_10gk, &
         time_50gk_25gk, rho_fin, vx_fin, vy_fin, vz_fin, r_ini, v_max
    real(8) :: temp_gk, x_ini, y_ini, z_ini,r_fin, x_5gk, y_5gk, z_5gk, r_5gk, vx_5gk, vy_5gk, vz_5gk, vr_5gk, s_10gk, s_3gk, ye_3gk, s_1gk,ye_1gk, rho_ini, t_ini, rho_max
    real(8) :: x_drip,y_drip,z_drip,r_drip,vx_drip,vy_drip,vz_drip,vr_drip,texp_drip,s_drip,ye_drip
    real(8) :: x_dmax,y_dmax,z_dmax,r_dmax,vx_dmax,vy_dmax,vz_dmax,vr_dmax,texp_rho_max,s_rho_max,ye_rho_max
    real(8) :: x_tmax,y_tmax,z_tmax,r_tmax,vx_tmax,vy_tmax,vz_tmax,vr_tmax,texp_tem_max,s_tem_max,ye_tem_max

    real(8) :: s1,s0

    integer :: my_thr, max_thr

    integer :: job, job_min, job_max, nsteps_job

    ! logical :: old_format = .true.
    ! logical :: format_2d = .true.

    real(8) :: time_init

    ! additional values
    real(8),allocatable :: &
         v_max_traj(:)
    !

    dir_out = dir_read
    !dir_out = "."
    !old_format = .false.

    if(old_format)then
       write(6,*) "itt_max from report file."
    else
       write(6,*) "counting itt_max..."

       job_min = 0
       job_max = 0
       job = 1
       jobs:do
          write(str1,'(i10)') job

          fn = trim(dir_read) // "/steps_"//trim(adjustl(str1))//".dat"
          ! write(6,'(a)') fn
          if(access(fn," ")==0)then
             if(job_min==0) job_min = job
             job_max=job
          else
             if(job_min/=0)exit jobs
          endif
          job=job+1
          if(job>1000)then
             write(6,*) "something is wrong"
             stop
          endif
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

       itt_max = itt_max + 1

    endif

    write(*,*) "model = ", trim(model)
    write(*,'(a)') trim(dir_read)//"/report_ptr.dat"
    open(10,file=trim(dir_read)//"/report_ptr.dat",status="old")

    if(old_format)then

       read(10,*); read(10,*) itt_min
       read(10,*); read(10,*) itt_max
       if(format_2d)then
          read(10,*); read(10,*)
          read(10,*); read(10,*)
       endif
       read(10,*); read(10,*) np
       read(10,*); read(10,*) np_inside
    else
       read(10,*); read(10,*) np
       read(10,*); read(10,*) np_inside
    endif

    write(*,'("np, np_inside, itt_min, itt_max=",4i7)') np,np_inside,itt_min,itt_max

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
         v_max_traj(np),x_ini_traj(np),y_ini_traj(np),z_ini_traj(np),rho_ini_traj(np),t_ini_traj(np),rho_max_traj(np), &
         s_rho_max_traj(np), texp_rho_max_traj(np), ye_rho_max_traj(np), &
         s_drip_traj(np), texp_drip_traj(np), ye_drip_traj(np), &
         s_tem_max_traj(np), texp_tem_max_traj(np), ye_tem_max_traj(np) )

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


    if(present(time_map))then
       time_init = time_map
    else
       time_init = 0d0
    endif

    !$omp parallel default(none) &
    !$omp shared(np,dir_out,dir_read,n_cond,mass_traj, &
    !$omp   t_fin_traj, x_fin_traj, y_fin_traj, z_fin_traj, vr_fin_traj, s_fin_traj, ye_fin_traj, ut1_fin_traj, hut_fin_traj, &
    !$omp   tem_max_traj, tem_max_af3gk_traj, t_tem_max_traj, t_tem_max_af3gk_traj, &
    !$omp   t_5gk_traj, s_5gk_traj, ye_5gk_traj, texp_5gk_traj, &
    !$omp   t_3gk_traj, t_1gk_traj, ye_10gk_traj, t_10gk_traj, &
    !$omp   time_50gk_25gk_traj, rho_fin_traj, vx_fin_traj, vy_fin_traj, vz_fin_traj, r_ini_traj, &
    !$omp   v_max_traj,x_ini_traj, y_ini_traj, z_ini_traj, rho_ini_traj, t_ini_traj, time_init, rho_max_traj, &
    !$omp   s_rho_max_traj, texp_rho_max_traj, ye_rho_max_traj, &
    !$omp   s_drip_traj, texp_drip_traj, ye_drip_traj, &
    !$omp   s_tem_max_traj, texp_tem_max_traj, ye_tem_max_traj) &
    !$omp private( fn,str1,nunit,it,time, x_p,y_p, z_p, vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, sen_p, rne_p, rae_p, &
    !$omp   it_max,it_tem_max,it_5gk,it_10gk,it_3gk,it_1gk, temp_gk, x_ini,y_ini,z_ini, &
    !$omp   smax_traj_7_0,smax_traj_7_4,smin_traj_7_0,smin_traj_7_4, r_fin, &
    !$omp   s1,s0,x_5gk, y_5gk, z_5gk, r_5gk, vx_5gk, vy_5gk, vz_5gk, vr_5gk, s_10gk, s_3gk, ye_3gk, s_1gk,ye_1gk, &
    !$omp   iflag_tem05,iflag_tem07,iflag_tem10,iflag_sc70,iflag_sc74, it_5gk_after_10gk, iflag_cond, &
    !$omp   t_fin, x_fin, y_fin, z_fin, vr_fin, s_fin, ye_fin, &
    !$omp   tem_max, tem_max_af3gk, t_tem_max_af3gk, t_tem_max, &
    !$omp   t_5gk, s_5gk, ye_5gk, texp_5gk, &
    !$omp   t_3gk, t_1gk, ye_10gk, t_10gk, &
    !$omp   time_50gk_25gk, rho_fin, vx_fin, vy_fin, vz_fin, r_ini,v_max, rho_ini, t_ini, it_init,  rho_max, &
    !$omp   x_drip,y_drip,z_drip,r_drip,vx_drip,vy_drip,vz_drip,vr_drip,texp_drip, s_drip, ye_drip, it_drip, &
    !$omp   x_dmax,y_dmax,z_dmax,r_dmax,vx_dmax,vy_dmax,vz_dmax,vr_dmax,texp_rho_max,s_rho_max,ye_rho_max, it_rho_max, &
    !$omp   x_tmax,y_tmax,z_tmax,r_tmax,vx_tmax,vy_tmax,vz_tmax,vr_tmax,texp_tem_max,s_tem_max,ye_tem_max, &
    !$omp   my_thr,max_thr)
    my_thr  = omp_get_thread_num()
    max_thr = omp_get_max_threads()
    if(my_thr==0)write(6,*) "max thread = ",max_thr
    !$omp do
    do ip =1,np

       write(str1,'(i8.8)') ip
       fn=trim(dir_read)//"/traj_"//trim(str1)//".dat"

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
99        continue
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

          rho_max = 0d0
          it_rho_max = 0
          it_drip    = 0

          v_max = 0d0
          do it = 1,it_max
             temp_gk = tem_p(it)
             if(tem_max < temp_gk)then
                tem_max = temp_gk
                it_tem_max = it
             endif

             v_max = max(v_max, sqrt(sqrt(vlx_p(it)**2 + vly_p(it)**2 + vlz_p(it)**2)))

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
                   if(it_5gk/=it_10gk)it_5gk = 0
                   if(it_3gk/=it_10gk)it_3gk = 0
                   if(it_1gk/=it_10gk)it_1gk = 0

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

                if( rho_max < qrho_p(it) )then
                   rho_max = qrho_p(it)
                   it_rho_max = it
                endif
                
                ! drip density
                if(rho_drip <= qrho_p(it) .and. rho_drip > qrho_p(it+1))then
                   it_drip = it
                endif
                
             endif

             ! write(6,'(i7,3es12.4,i7)') it, qrho_p(it),tem_p(it),time(it), it_5gk
          enddo

          it_init = 1
          do it=1,it_max
             if( abs(time_init - time(it)) < (time(2)-time(1))*0.5d0 )then
                it_init = it
                exit
             endif
          enddo
          
          ! values at initial time
          t_ini = time(it_init)
          x_ini  = x_p(it_init)
          y_ini  = y_p(it_init)
          z_ini  = z_p(it_init)
          r_ini  = sqrt(x_ini**2 + y_ini**2 + z_ini**2)
          
          rho_ini = qrho_p(it_init)

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

          ! at maximum temperature
          it = it_tem_max
          
          t_tem_max = time(it)
          
          s_tem_max = sen_p(it)
          ye_tem_max = ye_p(it)
          
          x_tmax =   x_p(it)
          y_tmax =   y_p(it)
          z_tmax =   z_p(it)
          vx_tmax= vlx_p(it)
          vy_tmax= vly_p(it)
          vz_tmax= vlz_p(it)
          
          r_tmax = sqrt(x_tmax**2 +y_tmax**2 +z_tmax**2) + 1.d-20
          vr_tmax= vx_tmax *x_tmax/r_tmax + vy_tmax *y_tmax/r_tmax + vz_tmax *z_tmax/r_tmax
          texp_tem_max = r_tmax/vr_tmax


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

          ! at drip density
          it = it_drip
          if(it>0)then
             s1 = (rho_drip-qrho_p(it))/(qrho_p(it+1)-qrho_p(it))
             s0 = 1.d0-s1

             s_drip = s1*sen_p(it+1) + s0*sen_p(it)
             ye_drip = s1*ye_p(it+1) + s0*ye_p(it)

             x_drip = s1*  x_p(it+1) + s0*  x_p(it)
             y_drip = s1*  y_p(it+1) + s0*  y_p(it)
             z_drip = s1*  z_p(it+1) + s0*  z_p(it)
             vx_drip= s1*vlx_p(it+1) + s0*vlx_p(it)
             vy_drip= s1*vly_p(it+1) + s0*vly_p(it)
             vz_drip= s1*vlz_p(it+1) + s0*vlz_p(it)
             r_drip = sqrt(x_drip**2 +y_drip**2 +z_drip**2) + 1.d-20
             vr_drip= vx_drip *x_drip/r_drip + vy_drip *y_drip/r_drip + vz_drip *z_drip/r_drip
             texp_drip = r_drip/vr_drip
          else
             s_drip    = 0.d0
             texp_drip = 0.d0
          endif
          
          ! at maximum density
          it = it_rho_max
          if(it>0)then
             s_rho_max = sen_p(it)
             ye_rho_max = ye_p(it)
             
             x_dmax =   x_p(it)
             y_dmax =   y_p(it)
             z_dmax =   z_p(it)
             vx_dmax= vlx_p(it)
             vy_dmax= vly_p(it)
             vz_dmax= vlz_p(it)

             r_dmax = sqrt(x_dmax**2 +y_dmax**2 +z_dmax**2) + 1.d-20
             vr_dmax= vx_dmax *x_dmax/r_dmax + vy_dmax *y_dmax/r_dmax + vz_dmax *z_dmax/r_dmax
             texp_rho_max = r_dmax/vr_dmax
          else
             s_rho_max    = 0.d0
             texp_rho_max = 0.d0
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
          x_ini_traj(ip) = x_ini
          y_ini_traj(ip) = y_ini
          z_ini_traj(ip) = z_ini
          rho_ini_traj(ip) = rho_ini
          t_ini_traj(ip)   = t_ini

          rho_max_traj(ip) = rho_max
          
          s_drip_traj(ip) = s_drip
          texp_drip_traj(ip) = texp_drip
          ye_drip_traj(ip) = ye_drip

          s_rho_max_traj(ip) = s_rho_max
          texp_rho_max_traj(ip) = texp_rho_max
          ye_rho_max_traj(ip) = ye_rho_max

          s_tem_max_traj(ip) = s_tem_max
          texp_tem_max_traj(ip) = texp_tem_max
          ye_tem_max_traj(ip) = ye_tem_max

          
       endif

       !if(mod(ip,100)==0) write(6,'("ip = ", i8)')  ip
       if(my_thr==0) write(6,*) "ip=",ip,np/max_thr
    enddo
    !$omp end do
    !$omp end parallel

    ! more reader-friendly file
    open(newunit=nunit,file=trim(dir_out)//"/stat_traj.dat",status="replace")  
    write(nunit,'("# model: ",a)') trim(model)
    write(nunit,'("# ntraj,ntraj_inside: ",99i15)') np, np_inside
    write(nunit,'("#",99i20)') (i,i=1,33)
    write(nunit,'("#",99a20)') "id",  "mass [g]",  "time_fin [s]",  "x_fin [cm]", "y_fin [cm]", "z_fin [cm]", "v^r_fin [cm/s]", "v^x_fin", "v^y_fin", "v^z_fin", &
                          "s_fin [k_b/nuc]", "Ye_fin", "time_5GK [s]", "s_5GK [k_b/nuc]", "Ye_5GK", "t_exp_5GK [s]", "Tmax [K]", "ut+1", "hut+h_at", "Ye(10GK)", &
                          "t(T=10GK)", "t(T=Tmax)", "(r/v)_5GK [s]", "rho_max [g/cm^3]", "s_drip [k_b/nuc]", "(r/v)_drip [s]", "Ye_drip", "s_dmax [k_b/nuc]", "(r/v)_dmax [s]", "Ye_dmax", &
                          "s_Tmax [kb/nuc]", "(r/v)_Tmax [s]", "Ye_Tmax"
    do ip=1,np
       write(nunit,'(" ",i20,99es20.7)') &
            ip, &
            mass_traj(ip), t_fin_traj(ip),x_fin_traj(ip),y_fin_traj(ip),z_fin_traj(ip),vr_fin_traj(ip),vx_fin_traj(ip),vy_fin_traj(ip),vz_fin_traj(ip),s_fin_traj(ip),ye_fin_traj(ip),t_5gk_traj(ip),s_5gk_traj(ip),ye_5gk_traj(ip),texp_5gk_traj(ip), &
            tem_max_traj(ip),ut1_fin_traj(ip),hut_fin_traj(ip),ye_10gk_traj(ip), t_10gk_traj(ip), t_tem_max_traj(ip), texp_5gk_traj(ip), rho_max_traj(ip), &
            s_drip_traj(ip),texp_drip_traj(ip),ye_drip_traj(ip), &
            s_rho_max_traj(ip), texp_rho_max_traj(ip), ye_rho_max_traj(ip), &
            s_tem_max_traj(ip), texp_tem_max_traj(ip), ye_tem_max_traj(ip)
    enddo

    close(nunit)
    stop


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
    if(present(time_map))then
       write(nunit,'("# model: ",a," time(map):", es12.4)') trim(model), time_map
    else
       write(nunit,'("# model: ",a)') trim(model)
    endif
    write(nunit,'("#",99i15)') np, np_inside
    write(nunit,'("#",99i15)') (i,i=1,16)
    write(nunit,'("#",99a15)') "id","mass","v_max", "x_initial", "y_initial", "z_initial", "x_final", "y_final", "z_final", "vx_final", "vy_final", "vz_final", "rho_initial", "rho_final", "time_initial", "time_final"

    do ip=1,np
       write(nunit,'(" ",i15,99es15.7)') &
            ip,mass_traj(ip),v_max_traj(ip),x_ini_traj(ip),y_ini_traj(ip),z_ini_traj(ip),x_fin_traj(ip),y_fin_traj(ip),z_fin_traj(ip),vx_fin_traj(ip),vy_fin_traj(ip),vz_fin_traj(ip),rho_ini_traj(ip),rho_fin_traj(ip),t_ini_traj(ip), t_fin_traj(ip)
    enddo

    close(nunit)

    write(6,*) "ana_traj.dat done."

    

  end subroutine tr_analysis

end module analysis
