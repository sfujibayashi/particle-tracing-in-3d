module module_store

  implicit none

  integer :: nt
  real(8),allocatable :: time_store(:)

  integer,allocatable :: flag_evol_store(:,:)
  real(8),allocatable :: x_p_store(:,:),y_p_store(:,:),z_p_store(:,:),&
       qrho_p_store(:,:),&
       ye_p_store(:,:),&
       tem_p_store(:,:),&
       ut_p_store(:,:),&
       qb_p_store(:,:),&
       sen_p_store(:,:),&
       vlx_p_store(:,:),&
       vly_p_store(:,:),&
       vlz_p_store(:,:),&
       hhh_p_store(:,:),&
       rne_p_store(:,:),&
       rae_p_store(:,:),&
       deptn_p_store(:,:),&
       depta_p_store(:,:),&
       dm_p_store(:,:), &
       ut1_p_store(:,:), &
       hut_p_store(:,:)


contains

  subroutine allocate_store(np,mem_max)
    integer,intent(in) :: np
    real(8),intent(in) :: mem_max
    integer :: ip,it
    
    nt = int(mem_max/ dble( (20*8 + 1*8)*np ))
    write(6,*) "number of snapshots that can be stored: ",nt
    
    allocate( &
         time_store(nt), &
         flag_evol_store(np,nt), &
         x_p_store(np,nt),y_p_store(np,nt),z_p_store(np,nt),&
         qrho_p_store(np,nt),&
         ye_p_store(np,nt),&
         tem_p_store(np,nt),&
         ut_p_store(np,nt),&
         qb_p_store(np,nt),&
         sen_p_store(np,nt),&
         vlx_p_store(np,nt),&
         vly_p_store(np,nt),&
         vlz_p_store(np,nt),&
         hhh_p_store(np,nt),&
         rne_p_store(np,nt),&
         rae_p_store(np,nt),&
         deptn_p_store(np,nt),&
         depta_p_store(np,nt),&
         dm_p_store(np,nt), &
         ut1_p_store(np,nt), &
         hut_p_store(np,nt) )
    
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
       x_p_store(ip,it_out) = x_p(ip)
       y_p_store(ip,it_out) = y_p(ip)
       z_p_store(ip,it_out) = z_p(ip)
       qrho_p_store(ip,it_out) = qrho_p(ip)
       ye_p_store(ip,it_out) = ye_p(ip)
       tem_p_store(ip,it_out) = tem_p(ip)
       ut_p_store(ip,it_out) = ut_p(ip)
       qb_p_store(ip,it_out) = qb_p(ip)
       sen_p_store(ip,it_out) = sen_p(ip)
       vlx_p_store(ip,it_out) = vlx_p(ip)
       vly_p_store(ip,it_out) = vly_p(ip)
       vlz_p_store(ip,it_out) = vlz_p(ip)
       hhh_p_store(ip,it_out) = hhh_p(ip)
       rne_p_store(ip,it_out) = rne_p(ip)
       rae_p_store(ip,it_out) = rae_p(ip)
       deptn_p_store(ip,it_out) = deptn_p(ip)
       depta_p_store(ip,it_out) = depta_p(ip)
       dm_p_store(ip,it_out) = dm_p(ip)
       ut1_p_store(ip,it_out) = ut1_p(ip)
       hut_p_store(ip,it_out) = hut_p(ip)

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
    !$omp   shared(np,flag_evol,dir_out,model,dm_p,ut1_p,hut_p) &
    !$omp   private(str1,unit)
    !$omp do
    do ip=1,np
       
       write(str1,'(i8.8)') ip
       open (newunit=unit,file=trim(dir_out)//"/traj_"//trim(str1)//".dat",status="replace")
       write(unit,'("# particle id:",i8)') ip
       write(unit,'("# model: ",a)') trim(model)
       write(unit,'("# particle mass:",es13.5," g, ut+1, hut+h_atm:",2es13.5)') dm_p(ip), ut1_p(ip), hut_p(ip)
       ! write(unit,'("#     Time [s]        x [cm]        y [cm]        z [cm]     Vx [cm/s]     Vy [cm/s]     Vz [cm/s]  rho [g/cm^3]         T [K]            Ye   S [k_b/nuc] Ee [erg/cm^3] Ea [erg/cm^3]         tau_e         tau_a")')
       write(unit,'("#",99a14)') "Time [s]","x [cm]","y [cm]","z [cm]","Vx [cm/s]","Vy [cm/s]","Vz [cm/s]","rho [g/cm^3]","T [K]","Ye","S [k_b/nuc]","Ee [erg/cm^3]"," Ea [erg/cm^3]","tau_e","tau_a","u_t+1","h-1"
       
       close(unit)

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
    !$omp   shared(np,nt_output,flag_evol_store,dir_out,time_store,x_p_store,y_p_store,z_p_store,vlx_p_store,vly_p_store,vlz_p_store, &
    !$omp          qrho_p_store, tem_p_store, ye_p_store,sen_p_store,rne_p_store,rae_p_store,deptn_p_store,depta_p_store,ut_p_store,hhh_p_store) &
    !$omp   private(str1,unit)
    !$omp do
    do ip=1,np
       write(str1,'(i8.8)') ip
       open(newunit=unit,file=trim(dir_out)//"/traj_"//trim(str1)//".dat",status="old",position="append")
       do it_out = 1, nt_output
          if(flag_evol_store(ip,it_out)==1)then
             write(unit,'(" ",99es14.6)') &
                  time_store(it_out), &
                  x_p_store(ip,it_out), &
                  y_p_store(ip,it_out), &
                  z_p_store(ip,it_out), &
                  vlx_p_store(ip,it_out)*clight, &
                  vly_p_store(ip,it_out)*clight, &
                  vlz_p_store(ip,it_out)*clight, &
                  qrho_p_store(ip,it_out), &
                  tem_p_store(ip,it_out)*mev_to_kelvin, &
                  ye_p_store(ip,it_out), &
                  sen_p_store(ip,it_out), &
                  rne_p_store(ip,it_out), &
                  rae_p_store(ip,it_out), &
                  deptn_p_store(ip,it_out), &
                  depta_p_store(ip,it_out), &
                  ut_p_store(ip,it_out)+1d0, &
                  hhh_p_store(ip,it_out)-1d0
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
