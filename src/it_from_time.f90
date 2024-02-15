subroutine it_from_time(job_min,job_max,filename,it_num,t_start,step,job_start,it_start)
  use hdf5
  use h5lt

  implicit none
  real(8),intent(in) :: t_start
  integer,intent(in) :: job_min, job_max, step
  integer,intent(in) :: it_num(job_min:job_max)
  character(*),intent(in) :: filename(job_min:job_max)
  integer,intent(out) :: job_start, it_start
  
  character(256) :: fn,str1
  integer :: job, it, job_
  real(8) :: t

  ! hdf5 I/O
  INTEGER        :: error         ! Error flag
  INTEGER(HID_T) :: file_id       ! File identifier
  integer(HSIZE_T) :: dims1(1)
  real(4),allocatable :: buf1(:)

  allocate(buf1(1))
  dims1(1) = 1

  search_job:do job=job_min,job_max
     fn = filename(job)
     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
     write(str1,'(i10)') it_num(job)
     call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str1))//"/time",buf1,dims1,error)
     call h5fclose_f(file_id, error)

     t = buf1(1)/1d3
     write(6,*) job,t,t_start
     if(t>t_start)then
        job_start=job
        exit search_job
     endif
  enddo search_job
  
  job=job_start
  fn = filename(job)
  call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
  search_it:do it=1,it_num(job)
     write(str1,'(i10)') it
     call H5LTread_dataset_float_f(file_id,"/level1/data"//trim(adjustl(str1))//"/time",buf1,dims1,error)

     t = buf1(1)/1d3
     write(6,*) it,t,t_start
     if(t>t_start)then
        it_start=it
        exit search_it
     endif
  enddo search_it
  call h5fclose_f(file_id, error)
  
end subroutine it_from_time
