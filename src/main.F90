program main
#include "macro.h"

  use const
  use simdata3D
  use particle_data
  use particle_set
  use io
  use hdf5
  use h5lt
  
  use unit

  use module_rho_ye
  use module_restart_hdf
  use module_eos
  
  use module_flux_sample
  implicit none

  real(8),parameter :: msun = 1.988d33

  !logical,parameter :: mode_backward = .true.
  !logical,parameter :: mode_volbased = .true.

  ! skip the timestep
  integer :: it_skip
  ! # of skip for data
  integer:: it_skip_out
  ! angular resolution
  integer :: n_theta
  ! for flux calculation
  real(8) :: rfl, rin
  ! for mass
  real(8) :: mass_crit, mass_min
  ! backward -> True, forward -> False
  logical :: mode_backward
  ! volume-based or flux-based. if flux-based (mode_volbased=False), particles are distributed in volme-based way at the last snapshot.
  logical :: mode_volbased
  ! incrementation in the next job
  integer :: incr_next
  ! whether reading file for pset timing
  logical :: read_pset_file = .false.

  integer :: step

  ! quantities for particles
  integer :: np,npv,np_flux,np_evolve,np_set

  character(200) :: model,fn,str1,str2,str3,dir_base
  
  character(200),allocatable :: filename(:)
  
  ! for hdf5 I/O
  INTEGER        :: error         ! Error flag
  INTEGER(HID_T) :: file_id,file2_id       ! File identifier
  integer(HSIZE_T) :: dims1(1),dims2(2),dims3(3)
  real(4),allocatable :: buf1(:), buf2(:,:), buf3(:,:,:)

  ! computational time
  real(8) :: time0,time1

  ! job-number and timestep related
  integer :: job_max,job_min,job_start,job_prv,job,job1,job2,it1,it2,it_start,nsteps,it0,it_save
  integer :: it,substep_max
  integer :: count_skip

  ! time
  real(8) :: time, time_prv, dt

  ! particle
  integer :: ip,ips,ipu

  ! setting particle
  integer :: it_skip_pset
  integer :: count_pset
  real(8) :: v_average, v_max,v_min,m_average,m_max,m_min, dt_pset
  
  ! EOS
  character(200) :: fn_eos
  integer :: nrho_in, ntemp_in, nye_in

  ! output
  integer :: count_out

  integer,allocatable :: nstep_job(:)
  
  integer :: job_restart, it_restart, it_prv
  character(1) :: restart
  logical :: first
  character(200) :: fn_read

  !!! 3D
  real(8) :: tms_start,tms_end,t_min,t_max!,tms_glo_min,tms_glo_max
  integer :: n_pset
  logical :: link_exists
  !real(8),allocatable :: tsta_job(:),tend_job(:),dt_job(:)
  !real(8) :: dtms_pset, dtms_snap
    
  integer :: unum, unum2
  integer :: access

  integer,parameter :: incl_next=20
  
  call h5open_f (error)


  ! fn = "/sakura/ptmp/shofu/SFHoTim276_13_14_0025_150mstg_B0_HLLC/Analysis_ptr/data_ana/map.h5"
  ! call all_proc(fn)
  ! stop
  dir_out = "."

  ! open(10,file="para.dat",status="old",action="read")
  ! read(10,*);read(10,'(a)') dir_out
  ! close(10)

  open(10,file=trim(dir_out)//"/parameters.dat",status="old",action="read")
  read(10,*);read(10,'(a)') model
  read(10,*);read(10,*) job_min
  read(10,*);read(10,*) job_max
  read(10,*);read(10,'(a)') dir_read
  read(10,*);read(10,'(a)') dir_out
  read(10,*);read(10,*) mode_backward
  read(10,*);read(10,*) mode_volbased
  
  read(10,*);read(10,*) it_start
  read(10,*);read(10,*) it_skip
  read(10,*);read(10,*) it_skip_out
  read(10,*);read(10,*) n_theta
  read(10,*);read(10,*) rfl
  read(10,*);read(10,*) rin
  read(10,*);read(10,*) mass_crit; mass_crit = mass_crit*msun
  read(10,*);read(10,*) mass_min; mass_min = mass_min*msun
  read(10,*);read(10,'(a)') fn_eos
  read(10,*);read(10,*) nrho_in, nye_in, ntemp_in
  read(10,*);read(10,*) incr_next

  read(10,*);read(10,*) restart
  if(restart=="Y")then
     read(10,*);read(10,*) job_restart
     read(10,*);read(10,*) it_restart
     if(it_restart/=0)then
        write(6,*) "finite it_restart is not supported yet. sorry!"
        stop
     endif
  endif
  close(10)

  fn = trim(dir_out)//"/restart_info.dat"
  if(access(fn," ")==0)then
     open(10,file=fn,status="old",action="read")
     read(10,*) job_min
     read(10,*) job_max
     read(10,*) restart
     close(10)
     if(restart=="Y")then
        if(mode_backward)then
           job_restart = job_max+1
        else
           job_restart = job_min-1
        endif
     endif
  endif
  
  
  write(*,'("model name      : ",a)') trim(model)
  write(*,'("job             : ",2i5)') job_min,job_max
  write(*,'("restart flag    : ",a)') restart
  write(*,'("data read from  : ",a)') trim(dir_read)
  write(*,'("result saved in : ",a)') trim(dir_out)
  
  ! call ascii(model,dir_out,it_skip_out)
  ! call tr_analysis(model,dir_out)
  ! stop

  
  write(6,'("it_skip, it_skip_out      : ",2i5)') it_skip,it_skip_out
  write(6,'("r_out, r_in (cm)          : ",2es12.4)') rfl,rin
  write(6,'("mass_crit, mass_min (Msun): ",2es12.4)') mass_crit/msun,mass_min/msun

  if(mode_backward)then
     write(6,'("### Evolution backward in time")')
     step=-1
     if(it_start==0)then
        write(6,'("trace from the last snapshot")')
     else
        write(6,'("it_start : ",2i5)') it_start
     endif

  else
     write(6,'("### Evolution forward in time")')
     step=1
     if(it_start==0)then
        write(6,'("trace from the first snapshot")')
     else
        write(6,'("it_start : ",2i5)') it_start
     endif

  endif
  
  if(mode_volbased)then
     write(6,'("### Flux-based particle setting : OFF")')
  else
     write(6,'("### Flux-based particle setting : ON")')
     write(6,'("    n_theta                     : ",i5)') n_theta
  endif

  ! fn_eos = "/sakura/ptmp/shofu/EOS/EOS_Hempel_SFHoTim326_TF"; nrho=408; nye=60; ntemp=131

  write(6,*)
  write(6,'("EOS table        : ",a)') trim(fn_eos)
  write(6,'("nrho, nye, ntemp : ",3i5)') nrho_in,nye_in,ntemp_in

  tms_end = 5.d0
  
! !!!
  call readeos(fn_eos,nrho_in,ntemp_in,nye_in)

!!! give division number of polar angle, then it calculates total number on a sphere
  call init_angle(n_theta,n_pset)
  call init_angle_sample(3*n_theta)

  write(6,*)
  write(6,*) "### particle set based on the flux"
  write(6,*) "n_theta, # of particles set in single step = ", n_theta,n_pset
  
  allocate(nstep_job(job_min:job_max))
  !allocate(it1_job(job_min:job_max),it2_job(job_min:job_max) )
  !allocate(tsta_job(job_min:job_max),tend_job(job_min:job_max),dt_job(job_min:job_max))
  
  write(6,*)

  ! set file name
  allocate(filename(job_min:job_max))
  do job = job_min,job_max
     write(str1,'(i10)') job
     fn = trim(dir_read)//trim(adjustl(str1))//"/raw3d.h5"
     write(*,'(i5,": ",a)') job,trim(fn)
     filename(job) = fn
  enddo
  
  do job = job_min,job_max
     
     write(str1,'(i10)') job
     fn = filename(job)
     
     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)

     it = 0
     setnum: do
        it = it + 1
        write(str1,'(i10)') it
        call h5lexists_f(file_id,"/level1/data"//trim(adjustl(str1)),link_exists,error)
        if(.not. link_exists)then
           nstep_job(job) = it - 1
           exit setnum
        endif
     enddo setnum

     tms(:)=0.d0
     write(str2,'(i10)') 1
     call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str2))//"/time",tms,dims1,error)
     t_min = tms(1)*time_unit_h5
     write(str2,'(i10)') nstep_job(job)
     call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str2))//"/time",tms,dims1,error)
     t_max = tms(1)*time_unit_h5

     write(6,'("job, steps in job, time(min,max) = ",2i4,2es12.4)') job,nstep_job(job),t_min,t_max
     !tsta_job(job) = t_min
     !tend_job(job) = t_max

     call h5fclose_f(file_id, error)

     !if(job==job_min) tms_glo_min = t_min
     !if(job==job_max) tms_glo_max = t_max
     
     write(str1,'(i10)') job
     open(11,file=trim(dir_out)//"/steps_"//trim(adjustl(str1))//".dat",status="replace",action="write")
     write(11,*) nstep_job(job)
     write(11,*) t_min
     write(11,*) t_max
     close(11)
  enddo

  ! dtms_snap = (tms_glo_max-tms_glo_min)/dble(sum(nstep_job(:)))
  ! do job=job_min,job_max
     
  !    ! if(tend_job(job)<rfl/(0.1d0*v_uni)*1d3)then
  !    !    it_skip_pset(job) = int(rfl*pi/2d0/dble(n_th)/(0.1d0*v_uni)*1d3/dtms_snap)
  !    ! else
  !    !    it_skip_pset(job) = (itt_max-itt_min+1)*n_pset/8000
  !    ! endif
  !    ! it_skip_pset(job) = 10

  !    it_skip_pset(job) = int(rfl*pi/2d0/dble(n_th)/(0.05d0*v_uni)*1d3/dtms_snap)
     
  !    ! if(job<=10)then
  !    !    it_skip_pset(job) = -1
  !    ! endif
  !    !it_skip_pset(job) = -1
     
  !    write(6,'("job, skip for particle set, dt_pset, dt_snap (ms): ",i4,i7,2es12.4,i5)') job, it_skip_pset(job), dble(it_skip_pset(job))*dtms_snap,dtms_snap!,nstep_job(job)/it_skip_pset(job)*n_pset
  ! enddo
  
  np_flux=0
  npv = 0
  !call calc_particle_number(nstep_job,it_skip_pset,job_min,job_max,n_pset,np_flux)
  !write(6,*) "# of particles set during evolution", np_flux

  if    (mode_backward)then
     job1 = job_max
     job2 = job_min

     job_start = job_max
     
     if(it_start==0.or.restart=="Y")it_start = nstep_job(job_start)
     
  elseif(.not.mode_backward)then
     job1 = job_min
     job2 = job_max
     
     job_start = job_min
     
     if(it_start==0)it_start = 1

  endif
  
!!! if you want to specify the start time...
  !t_start = 2.926d-2
  ! write(6,'("Searching for initial time....")')
  !call it_from_time(job_min,job_max,filename,nstep_job,t_start,step,job_start,it_start)

  ! if(job_start == job_min)then
  !    itt_start = it_start
  ! else
  !    itt_start = sum(nstep_job(job_min:job_start-1)) + it_start     
  ! endif
  
  ! if(mode_backward)then
  !    if(mod(itt_start,it_skip)==0)then
  !       itt_min = it_skip
  !    else
  !       itt_min = mod(itt_start,it_skip)
  !    endif
  !    itt_max = itt_start
  ! elseif(.not.mode_backward)then
  !    itt_min = itt_start
  !    itt_max = sum(nstep_job(job_min:job_max))
  ! endif
  ! nsteps = itt_max-itt_min+1
  
  write(6,'("# job_min,job_max,job_start,it_start = ",4i5)') job_min,job_max,job_start,it_start
  !write(6,'("# itt_min,itt_max,itt_start = ",3i10)') itt_min,itt_max,itt_start

!!! set particle data to be stored

!!! set start and end step in each job
  ! it1_job(:)=0
  ! it2_job(:)=0
  
  ! if(mode_backward)then

  !    job = job_start! job_max
  !    ittot=itt_start! sum(nstep_job(job_min:job-1)) + it_start
  !    do while(job>=job_min)
        
  !       do while(ittot > sum(nstep_job(job_min:job+step)))
  !          it = ittot-sum(nstep_job(job_min:job+step))
           
  !          if(it2_job(job)==0) it1_job(job) = it
  !          it2_job(job) = min(it1_job(job),it)
  !          ittot = ittot - it_skip
  !       enddo
  !       if(it2_job(job) == 0) it2_job(job) = +step
  !       write(6,'("job, Nstep, step, it1, it2",99i5)') job,nstep_job(job),it_skip*step,it1_job(job),it2_job(job)
  !       job = job - 1
  !    enddo

  ! else

  !    job  = job_start
  !    ittot= itt_start
  !    do while(job<=job_max)
        
  !       do while(ittot <= sum(nstep_job(job_min:job)))
  !          if(job==job_min)then
  !             it = ittot
  !          else
  !             it = ittot-sum(nstep_job(job_min:job-1))
  !          endif

  !          if(it2_job(job)==0) it1_job(job) = it
  !          it2_job(job) = max(it1_job(job),it)
  !          ittot = ittot + it_skip*step
  !       enddo
  !       if(it2_job(job) == 0) it2_job(job) = -1
  !       write(6,'("job, Nstep, step, it1, it2",99i5)') job,nstep_job(job),it_skip*step,it1_job(job),it2_job(job)
  !       job = job + step
  !    enddo
     
  ! endif

!!! show initial time step and time  

  fn = filename(job1)
  call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
  
  write(str2,'(i10)') it_start
  call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str2))//"/time",tms,dims1,error)
  call h5fclose_f(file_id, error)
  
  write(6,'("Initial time step: job, it = ",i3,i10, ", t = ",es13.5," s")') job1,it_start,tms(1)*time_unit_h5
  
  fn = filename(job2)
  call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
  
  write(str2,'(i10)') 1
  call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str2))//"/time",tms,dims1,error)
  call h5fclose_f(file_id, error)
  write(6,'("Last    time step: job, it = ",i3,i10, ", t = ",es13.5," s")') job2,1,tms(1)*time_unit_h5


  fn = filename(job_start)
  call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)

!!! obtain information on the grid
  call get_ngrid_info(file_id)
!!! allocate simulation data variables
  call allocate_simdata
!!! set coordinate (x,z)
  call get_coor(file_id)
  
  call h5fclose_f(file_id, error)

  if (restart=='N')then
!!! ips: accumulated number of particle. ipu: 
     first = .true.
     time = 0d0
     ips = 0
     ipu = 0
     !itt_pset_next = itt_start-it_skip
     count_pset = 1
     count_out  = 0
     count_skip = 0
     
     if(mode_backward)then
        job1 = job_start
        job2 = job_min
     else
        job1 = job_start
        job2 = job_max
     endif

! !!! save total steps job, it
!      open(11,file=trim(dir_out)//"/steps.dat",status="replace")
!      do job = job1,job2, step
        
!         it1=it1_job(job)
!         it2=it2_job(job)
        
!         ! write(11,'("# job, it1,it2,step = ",i3,2i5,i3)') job,it1,it2,step*it_skip
    
!         if(job==job_min)then
!            it0 = 0
!         else
!            it0 = sum(nstep_job(job_min:job-1))
!         endif
    
!         do it = it1, it2, step*it_skip
!            ittot = it0 + it
      
!            write(11,'(99i8)') job, it, ittot
!         enddo

!      enddo
!      close(11)
! !!!


  elseif(restart=='Y')then

     first = .false.

     write(str1,'(i3.3)') job_restart
     !write(str2,'(i6.6)') it_restart
     !fn = trim(dir_out) // "/data_"//trim(str1)//"_"//trim(str2)//".h5"
     fn = trim(dir_out) // "/res_"//trim(str1)//".h5"
     
!!! set ips and ipu
!!! allocate particle var.
!!! read all particle data
     ipu = 0
     ips = 0
     call read_checkpoint_hdf(fn,job_prv,it_prv,np,ipu,time,count_pset,count_out,count_skip,npv,fn_read)
     
     write(6,*)
     write(6,*) " -- Checkpoint file info -- "
     write(6,'("job,it  = ",2i5)') job_prv,it_prv
     write(6,'("np, npv = ",2i7)') np,npv
     write(6,'("time    = ",es12.4)') time
     write(6,'("counts  = ",3i5)') count_pset,count_out,count_skip
     write(6,'(a,a)') "file read: ",trim(fn_read)
     write(6,*) " -------------------------- "
     write(6,*)

     ips = ipu
     
     fn = fn_read
     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
!!! read sim data for the previous-step velocity
     call read_simdata(file_id,it_prv,time)
     call h5fclose_f(file_id, error)
     
     vlx_b(:,:,:,:) = vlx(:,:,:,:)
     vly_b(:,:,:,:) = vly(:,:,:,:)
     vlz_b(:,:,:,:) = vlz(:,:,:,:)

     ! if(mode_backward)then
     !    if(it_restart==it1_job(job_restart))then
     !       job1 = job_restart-1
     !    else
     !       job1 = job_restart
     !       it1_job(job_restart)=it_restart + step*it_skip
     !    endif

     !    job2 = max(job_min,job_end)
     ! elseif(.not.mode_backward)then
     !    job1 = job_restart
     !    job2 = min(job_max,job_end)
     ! endif
      
  else
     write(6,*) "restart flag is not set correctly"
     stop
  endif

  write(*,'("Tracing start...")')
  call cpu_time(time0)

  
!!! LOOP !!!!
  do job = job1,job2, step

     ! write(str1,'(i3.3)') job
     ! fn = trim(dir_out)//"/pset_"//trim(str1)//".h5"

     ! if(read_pset_file)then
     !    call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_pset_id, error)
        
     !    count_pset_max = 0
     !    count_pset = 0
     !    set_count_pset: do
     !       count_pset = count_pset + 1
     !       write(str1,'(i10)') count_pset
     !       call h5lexists_f(file_pset_id,"/"//trim(adjustl(str1)),link_exists,error)
     !       if(.not. link_exists)then
     !          count_pset_max = count_pset - 1
     !          exit set_count_pset
     !       endif
     !    enddo set_count_pset
        
     !    count_pset = 1
        
     !    if(count_pset <= count_pset_max)then
     !       write(str1,'(i10)') count_pset
     !       dims1(1) = 1
     !       allocate(ibuf1(1))
     !       call H5LTread_dataset_int_f(file_pset_id,"/"//trim(adjustl(str1))//"/it",ibuf1,dims1,error)
     !       it_pset_next=ibuf1(1)
     !       deallocate(ibuf1)
     !    else
     !       it_pset_next=-1
     !    endif

     !    write(6,'("particles are set ",i5," times in this job.")') count_pset_max
     !    write(6,*) "next : ",it_pset_next

     ! else
     !    count_pset = 0
     !    call h5fcreate_f(fn, H5F_ACC_TRUNC_F, file_pset_id, error)
     ! endif


     fn = filename(job)
     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
     
     if(job==job_min)then
        it0 = 0
     else
        it0 = sum(nstep_job(job_min:job-1))
     endif
     
     ! it1=it1_job(job)
     ! it2=it2_job(job)

     if(mode_backward)then

        if    (restart=="N".and.job==job1)then
           it1 = it_start
        else
           it1 = nstep_job(job)-count_skip
        endif
        it2 = 1

     else

        if    (restart=="N".and.job==job1)then
           it1 = it_start
        else
           it1 = count_skip
        endif
        it2 = nstep_job(job)
        
     endif

     write(6,'("job, it1,it2,step = ",i3,2i5,i3)') job,it1,it2,step*it_skip

     
     write(str1,'(i10)') job
     fn = trim(dir_out)//"/flux_"//trim(adjustl(str1))//".dat"
     open(newunit=unum,file=fn,status="replace",action="write")
     
     do it = it1, it2, step*it_skip
        
        time_prv = time
        ! itt_prv = ittot
        ! ittot = it0 + it
        
!!! read profile
        ! write(6,*) "read 3D data"
        call read_simdata(file_id,it,time)
        ! write(6,*) "set secondary"
        call set_secondary
        if(first)then
           dt = 0.d0
        else
           dt = time - time_prv
        endif
        
!!! evolve particles
        ! write(6,*) "evolve particles"
        call evolution_particle_3D(ipu,time,time_prv,substep_max)
        
!!! set particle
!!! set max number of particle at the first step
        if(first)then
           write(6,*) "first-time task"
           call analysis(dir_out,time)
           !call print_data(time,job,it)
           call partial_output_hdf(dir_out, job, it, time)
           !stop

           ! call set_ejecta_uniform(rfl,rin,mass_crit,mass_min,npv)
           ! np = npv
           call set_ejecta_inside_3D_divide(0,rfl,rin,mass_crit,mass_min,npv)
           np = npv

           write(*,*) "# of particles set in volume-based way : ",npv
           
           call allocate_particle_data(np)

           call set_ejecta_inside_3D_divide(1,rfl,rin,mass_crit,mass_min,npv)
           
           ips=npv
           ipu=ips

           !call partial_output_hdf(dir_out, job, it, time)
           !stop

        endif ! procedure only in the first step end
        
        if(.not.first)then
           
           if(.not.mode_volbased .and. count_pset == 0 )then
              
              if(ips+n_points>np)then
                 call reallocate_particle_data(np,ips+n_points)
                 np = ips+n_points
                 write(6,*) "np reset to", np
              endif
                 
              call set_new_particle(ips,rfl,dt,v_average,it_skip,it_skip_out,it_skip_pset,m_max,m_min,m_average,v_max,v_min,np_set)
              
              !itt_pset_next = ittot - it_skip_pset
              count_pset = it_skip_pset
              write(6,'("# of particles set = ",i5,", v/c(max,min,ave) = ",3es12.4,", m(max,min,ave) = ",3es12.4,". Next: dt (s), skip = ",es12.4,i5)') np_set, v_max,v_min,v_average, m_max,m_min,m_average, dt*dble(it_skip_pset), it_skip_pset
              
              write(unum,'(2i10,es15.7,i10,99es15.7)') job, it, time, np_set, m_average*dble(np_set)/(abs(dt)*dble(it_skip_pset)), sum(dm_p(1:ips)), m_average*dble(np_set), m_average, m_max, m_min, v_average, v_max,v_min
              
              ipu = ips
              
           endif
           
        endif ! particle set end
        
        ! output
        if( count_out==0)then
           ! procedure before output
           call set_particle_data(ipu)
           
           if(count_pset == it_skip_pset)then
              call sample_compare(rfl, ipu, np_set, job, it,time)
           endif
           
        endif
        
        vlx_b(:,:,:,:) = vlx(:,:,:,:)
        vly_b(:,:,:,:) = vly(:,:,:,:)
        vlz_b(:,:,:,:) = vlz(:,:,:,:)
        
        np_evolve = sum(flag_evol(:))
        write(6,'("job,it =",2i6,", Time, dt (s) = ",2es12.4, ", # of particle evolving = ",i8 "/",i8,i6,es12.4,2i6 )') job,it,time,dt, np_evolve,ipu,substep_max,sum(dm_p(1:ipu))/1.989d33,count_pset,count_out
        
        if((     mode_backward .and. (count_out==0 .or. first))   .or. &
           (.not.mode_backward .and. (count_out==0 .or. first))   )then
           
           ! write(str1,'(i6.6)') ittot
           ! open (11,file=trim(dir_out) // "/position_"//trim(str1)//".dat",status="replace")
           ! write(11,'("# ",99es12.4)') time, sum(dm_p(1:ipu))
           ! do ip = 1,ipu
           !    if(flag_evol(ip)==1)then
           !       !write(11,'(i10,99es11.3)') ip,dm_p(ip) ,x_p(ip),y_p(ip),z_p(ip),qrho_p(ip),ye_p(ip),tem_p(ip),sen_p(ip)
           !       write(11,'(i10,99es25.17)') ip, x_p(ip),y_p(ip),z_p(ip),&
           !            qrho_p(ip),&
           !            ye_p  (ip),&
           !            tem_p (ip),&
           !            ut_p  (ip),&
           !            qb_p  (ip),&
           !            sen_p (ip),&
           !            vlx_p (ip),&
           !            vly_p (ip),&
           !            vlz_p (ip),&
           !            hhh_p (ip),&
           !            rne_p (ip),&
           !            rae_p (ip),&
           !            deptn_p(ip),&
           !            depta_p(ip),&
           !            dm_p  (ip), &
           !            ut1_p (ip), &
           !            hut_p(ip)
           !    endif
           ! enddo
           ! close(11)
           
           if    (it-it_skip==0)then
              count_skip = 0
           elseif(it-it_skip< 0)then
              count_skip = it_skip - it
           else
              count_skip = nstep_job(job)-it+it_skip
           endif
           !if(nstep_job(job)-count_skip < 0) count_skip = count_skip - it_skip
           write(str1,'(i3.3)') job
           write(str2,'(i6.6)') it
           fn = trim(dir_out) // "/data_"//trim(str1)//"_"//trim(str2)//".h5"
           call save_checkpoint_hdf(fn,job,it,np,ipu,time,count_pset,count_out,count_skip,npv,filename(job))
           
        endif
        
        if( mode_backward.and. count_out==0 )then

           do ip=1,ipu
              if(flag_evol(ip)==1.and. &
                   ( &
                   ! tem_p(ip,ittot)>5.d0.or. &
                   ! ittot==itt_min.or. &
                   time*1.d3<tms_end) )then
                 
                 flag_evol(ip) = 0
                 
              endif
           enddo
           
        endif

        
        if(.not.mode_backward .and. count_out==0)then

           do ip = 1,ipu
              
              if(flag_evol(ip)==1.and. &
                   ( &
                   ! tem_p(ip,ittot)>5.d0.or. &
                   ! ittot==itt_max.or. &
                   sqrt(x_p(ip)**2+y_p(ip)**2+z_p(ip)**2)>rfl&
                   ! time*1.d3>tms_end &
                   ))then
                 
                 flag_evol(ip) = 0
                 
              endif
           enddo
           
        endif
        
        if(sum(flag_evol(:))==0)goto 100

        if(first) first = .false.
        if(count_out == 0) count_out = it_skip_out
        
        count_out = count_out - 1
        
        if(.not.mode_volbased) count_pset = count_pset - 1
     
        it_save = it
     enddo !end of the iteration of this job

     close(unum)

     call h5fclose_f(file_id, error)
     
     write(6,'("Output restart data")')
     write(str1,'(i3.3)') job
     fn = trim(dir_out) // "/res_"//trim(str1)//".h5"
     call save_checkpoint_hdf(fn,job,it_save,np,ipu,time,count_pset,count_out,count_skip,npv,filename(job))
     
     
  enddo !end of this job
100 continue
  write(6,'("Tracing finished.")')

  ! write(6,*) "Output particle data"
  ! do ip = 1,ipu
     
  !    write(str1,'(i8.8)') ip
  !    open (11,file=trim(dir_out)//"/traj_"//trim(str1)//".dat",status="replace")
  !    write(11,'("# particle id:",i8)') ip
  !    write(11,'("# model:",a)') trim(model)
  !    write(11,'("# particle mass:",es13.5," g, ut+1, hut+h_atm:",2es13.5)') dm_p(ip), ut1_p(ip), hut_p(ip)
  !    write(11,'("#     Time [s]        x [cm]        y [cm]        z [cm]     Vx [cm/s]     Vy [cm/s]     Vz [cm/s]  rho [g/cm^3]         T [K]            Ye   S [k_b/nuc] Ee [erg/cm^3] Ea [erg/cm^3]         tau_e         tau_a")')
  !    do itp_tmp=itp_sta(ip),itp_end(ip)
  !       write(11,'(99es14.6)') &
  !            time_p(itp_tmp), &
  !            x_p(ip,itp_tmp), &
  !            y_p(ip,itp_tmp), &
  !            z_p(ip,itp_tmp), &
  !            vlx_p(ip,itp_tmp)*v_uni, &
  !            vly_p(ip,itp_tmp)*v_uni, &
  !            vlz_p(ip,itp_tmp)*v_uni, &
  !            qrho_p(ip,itp_tmp), &
  !            tem_p(ip,itp_tmp)*tem_uni, &
  !            ye_p(ip,itp_tmp), &
  !            sen_p(ip,itp_tmp), &
  !            rne_p(ip,itp_tmp), &
  !            rae_p(ip,itp_tmp), &
  !            deptn_p(ip,itp_tmp), &
  !            deptn_p(ip,itp_tmp)
  !    enddo
  !    close(11)
     
  ! enddo
  
  call cpu_time(time1)
  write(6,*) time1-time0
  
  ! if((mode_backward.and.job2/=job_min).or.(.not.mode_backward.and.job2/=job_max))then
  ! ip=ip_test
  ! do it=1,it_max
  !    write(100,'(i8,99es13.5)') it,time,x_p(ip,ittot),y_p(ip,ittot),z_p(ip,ittot),qrho_p(ip,ittot),ye_p(ip,ittot),sen_p(ip,ittot)
  ! enddo
  
  write(6,'("Output summary")')
  !write(str1,'(i10)') job2
  open(11,file=trim(dir_out)//"/report_ptr.dat",status="replace")
  ! write(11,'("possible min:")')
  ! write(11,'(i6,es15.7)') !itt_min!,time_p(itp_min)
  ! write(11,'("Analysed to:")')
  ! write(11,'(i6,es15.7)') !itt_max!,time_p(itp_max)
  write(11,'("Number of trajectories evolved:")')
  write(11,'(i6)') ipu
  write(11,'("Number of trajectories inside the extraction:")')
  write(11,'(i6)') npv
  
  ! write(11,'("Analysed from:")')
  ! write(11,'(i6,es15.7)') itt_min!,time(itt_min)
  ! write(11,'("Analysed to:")')
  ! write(11,'(i6,es15.7)') itt_max!,time(itt_max)
  if(mode_volbased)then
     write(11,'("Volume-based particle set")')
  else
     write(11,'("flux-based particle set")')
  endif
  ! write(11,'("Mass crit, min of vol-based part:")')
  ! write(11,'(99es15.7)') mass_crit/1.989d33, mass_min/1.989d33
  close(11)

  fn = trim(dir_out)//"/restart_info.dat"
  open(10,file=fn,status="replace",action="write")
  if(mode_backward)then
     job_max = job2-1
     job_min = max(1,job_max-incl_next)
  else
     job_min = job2+1
     job_max = job_min+incl_next
  endif
  write(10,*) job_min
  write(10,*) job_max
  write(10,*) "Y"
  close(10)
  
  call ascii(model,dir_out,it_skip_out)
  call tr_analysis(model,dir_out)
  
end program main
