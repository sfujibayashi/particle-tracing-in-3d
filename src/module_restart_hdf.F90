module module_restart_hdf
  
  implicit none

contains
  subroutine save_checkpoint_hdf(fn,job,it,np,ipu,time,count_pset,count_out,count_skip,npv,fn_read)
    use particle_data
    use hdf5
    use h5lt

    integer,intent(in) :: ipu,np,count_pset,count_out,count_skip,npv,job,it
    real(8),intent(in) :: time
    character(*),intent(in) :: fn
    character(*),intent(in),optional :: fn_read
    
    integer :: ip

    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HSIZE_T) :: dims1(1)
    real(8),allocatable :: buf1(:)
    integer,allocatable :: ibuf1(:)

    ! write(6,*) fn,np,ipu,itt_tmp,itt_max
    
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
    
    ! dim-1 double
    allocate(buf1(1))
    buf1(1) = time
    call h5ltmake_dataset_double_f(file_id, "/time", 1, dims1, buf1, error)
    deallocate(buf1)


    ! particle data
    dims1(1) = ipu

    call h5ltmake_dataset_int_f   (file_id, "/flag_evol", 1, dims1, flag_evol(1:ipu), error)

    call h5ltmake_dataset_double_f(file_id, "/dm_p", 1, dims1, var_p(index_dm,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ut1_p", 1, dims1, var_p(index_ut1,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/hut_p", 1, dims1, var_p(index_hut,1:ipu), error)
    !call h5ltmake_dataset_double_f(file_id, "/ebind_p", 1, dims1, ebind_p(index_,1:ipu), error)

    call h5ltmake_dataset_double_f(file_id, "/x_p", 1, dims1, var_p(index_x,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/y_p", 1, dims1, var_p(index_y,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/z_p", 1, dims1, var_p(index_z,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/vlx_p", 1, dims1, var_p(index_vlx,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/vly_p", 1, dims1, var_p(index_vly,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/vlz_p", 1, dims1, var_p(index_vlz,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/qrho_p", 1, dims1, var_p(index_rho,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/tem_p", 1, dims1, var_p(index_tem,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ye_p", 1, dims1, var_p(index_ye,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/sen_p", 1, dims1, var_p(index_sen,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/rne_p", 1, dims1, var_p(index_rne,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/rae_p", 1, dims1, var_p(index_rae,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/deptn_p", 1, dims1, var_p(index_deptn,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/depta_p", 1, dims1, var_p(index_depta,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ut_p", 1, dims1, var_p(index_ut,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/hhh_p", 1, dims1, var_p(index_hhh,1:ipu), error)

    call h5ltmake_dataset_double_f(file_id, "/ch_nuf_p", 1, dims1, var_p(index_ch_nuf,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/ch_naf_p", 1, dims1, var_p(index_ch_naf,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/b2_p", 1, dims1, var_p(index_b2,1:ipu), error)
    call h5ltmake_dataset_double_f(file_id, "/pres_p", 1, dims1, var_p(index_pres,1:ipu), error)

    if(present(fn_read))then
       dims1(1) = 1
       allocate(ibuf1(1))
       ibuf1(1) = len(trim(fn_read))
       call h5ltmake_dataset_int_f(file_id, "/nstring" , 1, dims1, ibuf1, error)
       deallocate(ibuf1)
       call h5ltmake_dataset_string_f(file_id, "/hdf_read", trim(fn_read), error)
    endif

    call h5fclose_f(file_id,error)

  end subroutine save_checkpoint_hdf
  

  subroutine read_checkpoint_hdf(fn,job,it,np,ipu,time,count_pset,count_out,count_skip,npv,fn_read)
    use particle_data
    use hdf5
    use h5lt

    character(*),intent(in) :: fn
    integer,intent(out) :: ipu,np,count_pset,count_out,count_skip,job,it,npv
    real(8),intent(out) :: time
    character(*),intent(out),optional :: fn_read

    integer :: ip

    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HSIZE_T) :: dims1(1)
    real(8),allocatable :: buf1(:)
    integer,allocatable :: ibuf1(:)
    
    integer :: nstring
    character(200) :: str1
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
    
    deallocate(ibuf1,buf1)

    call allocate_particle_data(np)

    dims1(1) = ipu

    call H5LTread_dataset_int_f   (file_id,"/flag_evol",flag_evol(1:ipu),dims1,error)

    call H5LTread_dataset_double_f(file_id,"/dm_p",var_p(index_dm,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ut1_p",var_p(index_ut1,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/hut_p",var_p(index_hut,1:ipu),dims1,error)
    !call H5LTread_dataset_double_f(file_id,"/ebind_p",ebind_p(index_,1:ipu),dims1,error)

    call H5LTread_dataset_double_f(file_id,"/x_p"    , var_p(index_x,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/y_p"    , var_p(index_y,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/z_p"    , var_p(index_z,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/vlx_p"  , var_p(index_vlx,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/vly_p"  , var_p(index_vly,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/vlz_p"  , var_p(index_vlz,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/qrho_p" , var_p(index_rho,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/tem_p"  , var_p(index_tem,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ye_p"   , var_p(index_ye,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/sen_p"  , var_p(index_sen,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/rne_p"  , var_p(index_rne,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/rae_p"  , var_p(index_rae,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/deptn_p", var_p(index_deptn,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/depta_p", var_p(index_depta,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ut_p"   , var_p(index_ut,1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/hhh_p"  , var_p(index_hhh,1:ipu),dims1,error)

    call h5ltread_dataset_double_f(file_id, "/ch_nuf_p", var_p(index_ch_nuf,1:ipu),dims1, error)
    call h5ltread_dataset_double_f(file_id, "/ch_naf_p", var_p(index_ch_naf,1:ipu),dims1, error)
    call h5ltread_dataset_double_f(file_id, "/b2_p"    , var_p(index_b2,1:ipu),dims1, error)
    call h5ltread_dataset_double_f(file_id, "/pres_p"  , var_p(index_pres,1:ipu),dims1, error)


    if(present(fn_read))then
       dims1(1) = 1
       allocate(ibuf1(1))
       call H5LTread_dataset_int_f(file_id,"/nstring",ibuf1,dims1,error)
       nstring=ibuf1(1)
       deallocate(ibuf1)

       call H5LTread_dataset_string_f(file_id,"/hdf_read",str1,error)
       fn_read = str1(:nstring)
    endif

    call h5fclose_f(file_id,error)
    
  end subroutine read_checkpoint_hdf
  
end module module_restart_hdf
