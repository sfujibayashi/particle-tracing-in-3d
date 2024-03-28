module module_restart_hdf
  
  implicit none

contains
  subroutine save_checkpoint_hdf(fn,job,it,np,ipu,time,count_pset,count_out,count_skip,npv,fn_3d_read,fn_vel_read,flag_amr)
    use particle_data
    use hdf5
    use h5lt

    integer,intent(in) :: ipu,np,count_pset,count_out,count_skip,npv,job,it
    real(8),intent(in) :: time
    character(*),intent(in) :: fn
    character(*),intent(in),optional :: fn_3d_read,fn_vel_read
    logical,intent(in) :: flag_amr
    
    integer :: ip

    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HSIZE_T) :: dims1(1)
    real(8),allocatable :: buf1(:)
    integer,allocatable :: ibuf1(:)

    integer :: iflag
    INTEGER(HID_T) :: dspace_id, dset_id
    ! write(6,*) fn,np,ipu,itt_tmp,itt_max
    
    !INTEGER(HID_T)  :: filetype, memtype
    !integer :: sdim

    call h5fcreate_f(fn, H5F_ACC_TRUNC_F, file_id, error)
    
    ! dim-1 integers
    dims1(1) = 1
    allocate(ibuf1(1))
    ibuf1(1) = job
    call h5ltmake_dataset_int_f(file_id, "/job", 1, dims1, ibuf1, error)
    ibuf1(1) = it
    call h5ltmake_dataset_int_f(file_id, "/it" , 1, dims1, ibuf1, error)
    ibuf1(1) = np
    call h5ltmake_dataset_int_f(file_id, "/np", 1, dims1, ibuf1, error)
    ibuf1(1) = ipu
    call h5ltmake_dataset_int_f(file_id, "/ipu", 1, dims1, ibuf1, error)
    ibuf1(1) = count_pset
    call h5ltmake_dataset_int_f(file_id, "/count_pset", 1, dims1, ibuf1, error)
    ibuf1(1) = count_out
    call h5ltmake_dataset_int_f(file_id, "/count_out", 1, dims1, ibuf1, error)
    ibuf1(1) = count_skip
    call h5ltmake_dataset_int_f(file_id, "/count_skip", 1, dims1, ibuf1, error)
    ibuf1(1) = npv
    call h5ltmake_dataset_int_f(file_id, "/npv", 1, dims1, ibuf1, error)
    deallocate(ibuf1)
    
    if(flag_amr)then
       iflag = 1
    else
       iflag = 0
    endif
    CALL h5screate_simple_f(1, dims1, dspace_id, error)
    CALL h5dcreate_f(file_id, "flag_AMR", H5T_NATIVE_INTEGER, dspace_id, &
         dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, iflag, dims1, error)
    CALL h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
    
    
    ! dim-1 double
    allocate(buf1(1))
    buf1(1) = time
    call h5ltmake_dataset_double_f(file_id, "/time", 1, dims1, buf1, error)
    deallocate(buf1)


    ! particle data
    dims1(1) = ipu

    call h5ltmake_dataset_double_f(file_id, "/dm_p", 1, dims1, dm_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ut1_p", 1, dims1, ut1_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/hut_p", 1, dims1, hut_p(1:ipu), error)
    !call h5ltmake_dataset_double_f(file_id, "/ebind_p", 1, dims1, ebind_p(1:ipu), error)
    call h5ltmake_dataset_int_f   (file_id, "/flag_evol", 1, dims1, flag_evol(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/t_p", 1, dims1, t_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/x_p", 1, dims1, x_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/y_p", 1, dims1, y_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/z_p", 1, dims1, z_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/qrho_p", 1, dims1, qrho_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/tem_p", 1, dims1, tem_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ye_p", 1, dims1, ye_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/sen_p", 1, dims1, sen_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/hhh_p", 1, dims1, hhh_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ut_p", 1, dims1, ut_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/vlx_p", 1, dims1, vlx_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/vly_p", 1, dims1, vly_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/vlz_p", 1, dims1, vlz_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/rne_p", 1, dims1, rne_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/rae_p", 1, dims1, rae_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/deptn_p", 1, dims1, deptn_p(1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/depta_p", 1, dims1, depta_p(1:ipu), error)

    if(present(fn_3d_read))then

       ! dims1(1) = 1
       ! sdim = len(trim(fn_3d_read))
       
       ! CALL H5Tcopy_f(H5T_C_S1, filetype, error)
       ! CALL H5Tset_size_f(filetype, sdim+1, error)

       ! CALL H5Tcopy_f( H5T_FORTRAN_S1, memtype, error)
       ! CALL H5Tset_size_f(memtype, sdim, error)
       ! CALL h5screate_simple_f(1, dims1, dspace_id, error)
       ! CALL h5dcreate_f(file_id, "/hdf_read", filetype, dspace_id, dset_id, error)
       ! f_ptr = C_LOC(wdata(1)(1:1))
       ! CALL H5Dwrite_f(dset, memtype, f_ptr, hdferr)

       ! allocate(ibuf1(1))
       ! ibuf1(1) = len(trim(fn_3d_read))
       ! call h5ltmake_dataset_int_f(file_id, "/nstring" , 1, dims1, ibuf1, error)
       ! deallocate(ibuf1)
       call h5ltmake_dataset_string_f(file_id, "/hdf_read", trim(fn_3d_read), error)
    endif

    if(present(fn_vel_read))then
       call h5ltmake_dataset_string_f(file_id, "/vel_read", trim(fn_vel_read), error)
    endif


    call h5fclose_f(file_id,error)

  end subroutine save_checkpoint_hdf
  

  subroutine read_checkpoint_hdf(fn,job,it,np,ipu,time,count_pset,count_out,count_skip,npv,fn_3d_read,fn_vel_read,flag_amr)
    use particle_data
    use hdf5
    use h5lt

    character(*),intent(in) :: fn
    integer,intent(out) :: ipu,np,count_pset,count_out,count_skip,job,it,npv
    real(8),intent(out) :: time
    character(*),intent(out),optional :: fn_3d_read, fn_vel_read
    logical,intent(out) :: flag_amr

    integer :: ip

    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HSIZE_T) :: dims1(1)
    real(8),allocatable :: buf1(:)
    integer,allocatable :: ibuf1(:)
    
    !integer :: nstring
    !character(200) :: str1
    INTEGER(HID_T) :: dset_id,filetype,dspace_id,memtype
    INTEGER(SIZE_T) :: size
    INTEGER(SIZE_T),parameter :: sdim=256
    INTEGER(HSIZE_T), DIMENSION(1:1) :: maxdims
    character(256),target :: str_read(1)
    TYPE(c_ptr) :: f_ptr
    
    logical :: link_exists
    integer :: iflag

    write(6,'(a)') fn

    call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)
    
    dims1(1) = 1
    allocate(ibuf1(1),buf1(1))

    call H5LTread_dataset_int_f(file_id,"/job",ibuf1,dims1,error)
    job=ibuf1(1)
    call H5LTread_dataset_int_f(file_id,"/it",ibuf1,dims1,error)
    it=ibuf1(1)
    call H5LTread_dataset_int_f(file_id,"/np",ibuf1,dims1,error)
    np=ibuf1(1)
    call H5LTread_dataset_int_f(file_id,"/ipu",ibuf1,dims1,error)
    ipu=ibuf1(1)
    call H5LTread_dataset_double_f(file_id,"/time",buf1,dims1,error)
    time=buf1(1)
    call H5LTread_dataset_int_f(file_id,"/count_pset",ibuf1,dims1,error)
    count_pset=ibuf1(1)
    call H5LTread_dataset_int_f(file_id,"/count_out",ibuf1,dims1,error)
    count_out=ibuf1(1)
    call H5LTread_dataset_int_f(file_id,"/count_skip",ibuf1,dims1,error)
    count_skip=ibuf1(1)
    call H5LTread_dataset_int_f(file_id,"/npv",ibuf1,dims1,error)
    npv=ibuf1(1)

    call h5lexists_f(file_id,"/flag_AMR",link_exists,error)
    if(.not.link_exists)then
       write(6,*) "flag_AMR not exists. STOP"
       stop
    endif

    call H5LTread_dataset_int_f(file_id,"/flag_AMR",ibuf1,dims1,error)
    iflag=ibuf1(1)
    if(iflag==1)then
       flag_amr=.true.
    else
       flag_amr=.false.
    endif

    deallocate(ibuf1,buf1)

    call allocate_particle_data(np)

    dims1(1) = ipu
    call H5LTread_dataset_double_f(file_id,"/dm_p",dm_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ut1_p",ut1_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/hut_p",hut_p(1:ipu),dims1,error)
    !call H5LTread_dataset_double_f(file_id,"/ebind_p",ebind_p(1:ipu),dims1,error)
    call H5LTread_dataset_int_f   (file_id,"/flag_evol",flag_evol(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/t_p",t_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/x_p",x_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/y_p",y_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/z_p",z_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/qrho_p",qrho_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/tem_p" ,tem_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ye_p"  ,ye_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/sen_p" ,sen_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/hhh_p" ,hhh_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ut_p" ,ut_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/vlx_p" ,vlx_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/vly_p" ,vly_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/vlz_p" ,vlz_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/rne_p" ,rne_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/rae_p" ,rae_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/deptn_p",deptn_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/depta_p",depta_p(1:ipu),dims1,error)

    if(present(fn_3d_read))then
       ! dims1(1) = 1
       ! allocate(ibuf1(1))
       ! call H5LTread_dataset_int_f(file_id,"/nstring",ibuf1,dims1,error)
       ! nstring=ibuf1(1)
       ! deallocate(ibuf1)

       ! call H5LTread_dataset_string_f(file_id,"/hdf_read",str1,error)
       ! fn_3d_read = str1(:nstring)

       dims1(1) = 1
       
       CALL h5dopen_f(file_id, "/hdf_read", dset_id, error)
       CALL H5Dget_type_f(dset_id, filetype, error)
       CALL H5Tget_size_f(filetype, size, error)
       
       IF(size.GT.sdim+1)THEN
          PRINT*,'ERROR: Character LEN is too small'
          STOP
       ENDIF

       CALL H5Dget_space_f(dset_id, dspace_id, error)
       CALL H5Sget_simple_extent_dims_f(dspace_id, dims1, maxdims, error)

       CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, error)
       ! CALL H5Tcopy_f(H5T_C_S1, memtype, error)
       CALL H5Tset_size_f(memtype, sdim, error)
       
       f_ptr = C_LOC(str_read(1))
       CALL H5Dread_f(dset_id, memtype, f_ptr, error, dspace_id)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)
       call h5tclose_f(filetype,error)
       call h5tclose_f(memtype,error)

       fn_3d_read = trim(str_read(1))
       ! write(6,'(a)') str_read(1)
       ! write(6,'(a)') fn_3d_read
       ! stop
       
    endif

    if(present(fn_vel_read))then
       ! dims1(1) = 1
       ! allocate(ibuf1(1))
       ! call H5LTread_dataset_int_f(file_id,"/nstring",ibuf1,dims1,error)
       ! nstring=ibuf1(1)
       ! deallocate(ibuf1)

       ! call H5LTread_dataset_string_f(file_id,"/hdf_read",str1,error)
       ! fn_3d_read = str1(:nstring)

       dims1(1) = 1
       
       CALL h5dopen_f(file_id, "/vel_read", dset_id, error)
       CALL H5Dget_type_f(dset_id, filetype, error)
       CALL H5Tget_size_f(filetype, size, error)
       
       IF(size.GT.sdim+1)THEN
          PRINT*,'ERROR: Character LEN is too small'
          STOP
       ENDIF

       CALL H5Dget_space_f(dset_id, dspace_id, error)
       CALL H5Sget_simple_extent_dims_f(dspace_id, dims1, maxdims, error)

       CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, error)
       ! CALL H5Tcopy_f(H5T_C_S1, memtype, error)
       CALL H5Tset_size_f(memtype, sdim, error)
       
       f_ptr = C_LOC(str_read(1))
       CALL H5Dread_f(dset_id, memtype, f_ptr, error, dspace_id)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)
       call h5tclose_f(filetype,error)
       call h5tclose_f(memtype,error)

       fn_vel_read = trim(str_read(1))
       !write(6,'(a)') str_read
       !stop
       
    endif

    call h5fclose_f(file_id,error)
    
  end subroutine read_checkpoint_hdf
  
end module module_restart_hdf
