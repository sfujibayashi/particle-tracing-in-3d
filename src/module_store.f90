module module_store

  implicit none

  integer :: nt
  real(8),allocatable :: time_store(:)

  integer,allocatable :: flag_evol_store(:,:)
  real(8),allocatable :: var_p_store(:,:,:)

contains

  subroutine allocate_store(np,mem_max)
    use particle_data, only: nvar
    integer,intent(in) :: np
    real(8),intent(in) :: mem_max
    integer :: ip,it
    
    nt = int(mem_max/ dble( (20*8 + 1*8)*np ))
    write(6,*) "number of snapshots that can be stored: ",nt
    
    allocate( &
         time_store(nt), &
         flag_evol_store(np,nt), &
         var_p_store(nvar,np,nt))
    
    ! initialize
    flag_evol_store(:,:) = 0
    
  end subroutine allocate_store
  

  subroutine store_data(np,it_out,time)
    ! use module_binary_to_ascii
    ! use module_hdf_to_ascii
    use particle_data
    integer,intent(in) :: np,it_out
    real(8),intent(in) :: time
    integer :: ip

    ! write(6,*) np,it_out,time

    !ip=1000
    !write(6,*) time,flag_evol(ip)

    time_store(it_out) = time

    do ip=1,np
       flag_evol_store(ip,it_out) = flag_evol(ip)
       var_p_store(:,ip,it_out) = var_p(:,ip)
    enddo
    
  end subroutine store_data

  subroutine output_first(np,model,dir_out)
    ! use module_binary_to_ascii
    ! use module_hdf_to_ascii
    use particle_data
    integer,intent(in) :: np
    character(*),intent(in) :: model, dir_out
    integer :: ip,it_out
    integer :: unit
    character(10) :: str1
    
    !$omp parallel default(none) &
    !$omp   shared(np,flag_evol,dir_out,model,var_p) &
    !$omp   private(str1,unit)
    !$omp do
    do ip=1,np
       !if(flag_evol(ip))then
       
       write(str1,'(i8.8)') ip
       open (newunit=unit,file=trim(dir_out)//"/traj_"//trim(str1)//".dat",status="replace")
       write(unit,'("# particle id:",i8)') ip
       write(unit,'("# model: ",a)') trim(model)
       write(unit,'("# particle mass:",es13.5," g, ut+1, hut+h_atm:",2es13.5)') var_p(index_dm,ip), var_p(index_ut1,ip), var_p(index_hut,ip)
       ! write(unit,'("#     Time [s]        x [cm]        y [cm]        z [cm]     Vx [cm/s]     Vy [cm/s]     Vz [cm/s]  rho [g/cm^3]         T [K]            Ye   S [k_b/nuc] Ee [erg/cm^3] Ea [erg/cm^3]         tau_e         tau_a           u_t         h/c^2")')
       write(unit,'("#",99a14)') "Time [s]", label_var(1:nvar_out)
       close(unit)
       
       !endif

    enddo
    !$omp end do
    !$omp end parallel

  end subroutine output_first

  subroutine output_stored_data(np,dir_out,nt_output)
    integer,intent(in) :: np,nt_output
    character(*),intent(in) :: dir_out

    integer :: ip,it_out
    integer :: unit
    character(10) :: str1

    real(8),parameter :: mev_to_kelvin  = 1.160445d10, clight = 2.99792458d10
    
    !$omp parallel default(none) &
    !$omp   shared(np,nt_output,flag_evol_store,dir_out,time_store,var_p_store) &
    !$omp   private(str1,unit)
    !$omp do
    do ip=1,np
       write(str1,'(i8.8)') ip
       open(newunit=unit,file=trim(dir_out)//"/traj_"//trim(str1)//".dat",status="old",position="append")
       do it_out = 1, nt_output
          if(flag_evol_store(ip,it_out)==1)then
             write(unit,'(99es14.6)') &
                  time_store(it_out), &
                  var_p_store(1:20,ip,it_out)
          endif
          
       enddo
       close(unit)
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine output_stored_data

  
  ! subroutine output_hdf_first(np,model,dir_out,)
  !   ! use module_binary_to_ascii
  !   use module_hdf_to_ascii
  !   use hdf5
  !   use h5lt

  !   integer,intent(in) :: np
  !   character(*),intent(in) :: model, dir_out
  !   integer :: ip,it_out
  !   integer :: unit
  !   character(10) :: str1
  !   character(256) :: fn

  !   INTEGER        :: hdferr       ! Error flag
  !   INTEGER(HID_T) :: file_id      ! File identifier
  !   integer(HSIZE_T) :: dims1(1)
  !   real(8),allocatable :: buf1(:)
  !   character,allocatable :: sbuf1(:)
  !   integer,allocatable :: ibuf1(:)
  !   INTEGER(HSIZE_T) :: maxdims(2) = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
    
  !   INTEGER(HID_T) :: dataspace     ! Dataspace identifier
  !   INTEGER(HID_T) :: memspace      ! Memory dataspace identifier
  !   INTEGER(HID_T) :: crp_list      ! Dataset creation property identifier

  !   call h5open_f(hdferr)

  !   fn = trim(dir_out)//"/traj.h5"
    
  !   call h5fcreate_f(fn, H5F_ACC_TRUNC_F, file_id, hdferr)

  !   ! model name
  !   call h5ltmake_dataset_string_f(file_id, "/model", trim(model), hdferr)

  !   ! dim-1 integers
  !   dims1(1) = 1
  !   allocate(ibuf1(1))
  !   ibuf1(1) = np
  !   call h5ltmake_dataset_int_f(file_id, "/np", 1, dims1, ibuf1, hdferr)
  !   deallocate(ibuf1)
    
  !   ! ! dim-1 double
  !   ! allocate(buf1(1))
  !   ! buf1(1) = time
  !   ! call h5ltmake_dataset_double_f(file_id, "/time", 1, dims1, buf1, hdferr)
  !   ! deallocate(buf1)

  !   !Create the data space with unlimited dimensions.
  !   dims1(1) = nt
  !   call h5screate_simple_f(1, dims1, dataspace, error, maxdims)
    
  !   !Modify dataset creation properties, i.e. enable chunking
  !   call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)

  !   call h5pset_chunk_f(crp_list, 1, dimsc, error)

  !   !Create a dataset with  dimensions dims1 using cparms creation properties.
  !   CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
  !        dset_id, error, crp_list )

    
  !   call h5fclose_f(file_id,hdferr)
  !   stop

    
  ! end subroutine output_hdf_first

end module module_store
