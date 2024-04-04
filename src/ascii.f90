subroutine ascii(model, dir_read, it_skip_out)
  ! use module_binary_to_ascii
  !use module_hdf_to_ascii
  use particle_data
  use module_store
  use module_restart_hdf
  use hdf5
  use h5lt
  implicit none

  character(*),intent(in) :: model, dir_read
  integer,intent(in) :: it_skip_out

  real(8) :: time
  integer :: ipu,np,npv,np_buf,count_pset,count_out,count_skip,job,it,itt
  character(256) :: fn, dir_out, str1,str2,fn_hdf,fn_vel

  integer :: access

  logical :: first
  real(8) :: mem_max = 1d9/2d0 ! byte
  integer :: it_out

  integer :: hdferr
  integer :: i, nunit

  integer :: job_min, job_max
  integer :: nsteps_job

  integer :: job_buf,it_buf
  logical :: flag_amr

  call h5open_f (hdferr)
  
  call deallocate_particle_data

  dir_out = trim(dir_read)

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

  first=.true.
  it_out=0
  itt=0
  do job=job_min,job_max
     write(str1,'(i10)') job
     fn = trim(dir_read) // "/steps_"//trim(adjustl(str1))//".dat"
     open(newunit=nunit,file=fn,status="old",action="read")
     read(nunit,*) nsteps_job
     close(nunit)

     do it=1,nsteps_job
        write(str1,'(i3.3)') job
        write(str2,'(i6.6)') it
        fn = trim(dir_read) // "/data_"//trim(str1)//"_"//trim(str2)//".h5"
        if(access(fn," ")==0)then
           itt = itt + 1
           ! call read_checkpoint_hdf(fn,job_buf,it_buf,np_buf,ipu,time,count_pset,count_out,count_skip,npv,fn_hdf,fn_vel,flag_amr)
           call read_checkpoint_hdf(fn,job_buf,it_buf,np_buf,ipu,time,count_pset,count_out,count_skip,npv,fn_hdf)
           write(6,'(4i8,es12.4)') job,it,itt,it_out,time

           if(first)then
              np=ipu
              write(6,*) "# of particles:",np
              call allocate_store(np,mem_max)
              call output_first(np,model,dir_out)

              first=.false.

              write(str1,'(i6.6)') itt
              fn = trim(dir_out) // "/position_"//trim(str1)//".dat"
              call output_position(fn,ipu,time)
              
           endif
           
           it_out=it_out+1
           call store_data(np,it_out,time)
           if(it_out==nt)then
              call output_stored_data(np,dir_out,it_out)
              it_out=0
           endif

        endif

     enddo

  enddo

  if(it_out>0)then
     call output_stored_data(np,dir_out,it_out)
  endif

  write(6,*) "files of the trajectories done."
  
end subroutine ascii
