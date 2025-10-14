program interpolation
  
  implicit none
  
  
  character(256) :: fn_orig, fn_interp, fn_ref
  
  integer :: it_interp

  block
    use inputparser
    character(256) :: fn_para
    fn_para = "interp.para"
    call get_string_parameter(fn_para, "fn_orig", fn_orig)
    call get_string_parameter(fn_para, "fn_interp", fn_interp)
    call get_string_parameter(fn_para, "fn_ref", fn_ref)
    call get_integer_parameter(fn_para, "it_interp", it_interp)
  end block

  block
    use interp_data3D
    use simdata3D

    use hdf5
    
    ! for hdf5 I/O
    INTEGER        :: hdf_err         ! Hdf_Err flag
    INTEGER(HID_T) :: file_id_orig, file_id_ref       ! File identifier
    ! integer(HSIZE_T) :: dims1(1),dims2(2),dims3(3)

    real(8) :: time

    call h5open_f (hdf_err)
    
    ! original data
    call h5fopen_f(fn_orig, H5F_ACC_RDONLY_F, file_id_orig, hdf_err)
    ! obtain information on the grid
    call get_ngrid_info(file_id_orig)
    ! allocate simulation data variables
    call allocate_simdata
    ! set coordinate (x,y,z)
    call get_coor(file_id_orig)
    ! read data
    call read_simdata(file_id_orig,it_interp,time)
    
    call h5fclose_f(file_id_orig, hdf_err)



    ! original data
    call h5fopen_f(fn_ref, H5F_ACC_RDONLY_F, file_id_ref, hdf_err)
    ! obtain information on the grid
    call interp_get_grid_info(file_id_ref)
    ! allocate simulation data variables
    call interp_allocate_simdata
    ! set coordinate (x,y,z)
    call interp_get_coor(file_id_ref)
    
    call h5fclose_f(file_id_ref, hdf_err)

    ! interpolate
    call interp_data()
    
    
    call interp_output_hdf(fn_interp, time, it_interp)
    
    
    call h5close_f(hdf_err)
     
  end block
end program interpolation
