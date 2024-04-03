subroutine set_particle_from_hdf(fn_hdf,np,ipu)

  use particle_data
  use hdf5
  use h5lt

  implicit none

  character(*),intent(in) :: fn_hdf

  integer,intent(out) :: np, ipu

  INTEGER        :: error         ! Error flag
  INTEGER(HID_T) :: file_id       ! File identifier
  integer(HSIZE_T) :: dims1(1)
  real(8),allocatable :: buf1(:)
  integer,allocatable :: ibuf1(:)

  write(6,*) "particles are set from the file"
  write(6,'(a)') trim(fn_hdf)

  call h5fopen_f(fn_hdf, H5F_ACC_RDONLY_F, file_id, error)
  
  dims1(1)=1
  allocate(ibuf1(1))

  call H5LTread_dataset_int_f(file_id,"/np",ibuf1,dims1,error)
  np=ibuf1(1)
  call H5LTread_dataset_int_f(file_id,"/ipu",ibuf1,dims1,error)
  ipu=ibuf1(1)
  
  deallocate(ibuf1)
  
  call allocate_particle_data(np)
  
  dims1(1) = ipu
  call H5LTread_dataset_double_f(file_id,"/dm_p",dm_p(1:ipu),dims1,error)
  call H5LTread_dataset_double_f(file_id,"/x_p",x_p(1:ipu),dims1,error)
  call H5LTread_dataset_double_f(file_id,"/y_p",y_p(1:ipu),dims1,error)
  call H5LTread_dataset_double_f(file_id,"/z_p",z_p(1:ipu),dims1,error)
  call H5LTread_dataset_double_f(file_id,"/ut1_p",ut1_p(1:ipu),dims1,error)
  call H5LTread_dataset_double_f(file_id,"/hut_p",hut_p(1:ipu),dims1,error)
  call H5LTread_dataset_double_f(file_id,"/t_p",t_p(1:ipu),dims1,error)


  flag_evol(1:ipu) = 1
  
  call h5fclose_f(file_id,error)
  
end subroutine set_particle_from_hdf
