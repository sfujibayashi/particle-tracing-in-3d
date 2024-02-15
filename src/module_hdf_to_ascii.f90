module module_hdf_to_ascii
  implicit none

  ! They have only np
  integer,allocatable :: flag_evol(:)
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
       rne_p (:),&
       rae_p (:),&
       deptn_p(:),&
       depta_p(:),&
       dm_p  (:), &
       ut1_p (:), &
       hut_p(:)

contains

  subroutine allocate_particle_data(np)
    integer,intent(in) :: np
    integer :: ip

    if(.not.allocated(flag_evol))then
       
       allocate( &
            flag_evol(np), &
            x_p   (np),&
            y_p   (np),&
            z_p   (np),&
            qrho_p(np),&
            ye_p  (np),&
            tem_p (np),&
            ut_p  (np),&
            qb_p  (np),&
            sen_p (np),&
            vlx_p (np),&
            vly_p (np),&
            vlz_p (np),&
            hhh_p (np),&
            rne_p (np),&
            rae_p (np),&
            deptn_p(np),&
            depta_p(np),&
            dm_p  (np), &
            ut1_p (np), &
            hut_p (np)  )
    endif

    ! initialize
    flag_evol(:) = 0
    
  end subroutine allocate_particle_data
  
  subroutine read_checkpoint(fn,job,it,np,ipu,time,itt_pset_next)
    use hdf5
    use h5lt

    character(*),intent(in) :: fn
    integer,intent(out) :: ipu,np,itt_pset_next,job,it
    real(8),intent(out) :: time
    
    integer :: ip

    INTEGER        :: error         ! Error flag
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HSIZE_T) :: dims1(1)
    real(8),allocatable :: buf1(:)
    integer,allocatable :: ibuf1(:)


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
    call H5LTread_dataset_int_f(file_id,"/itt_pset_next",ibuf1,dims1,error)
    itt_pset_next=ibuf1(1)
    
    deallocate(ibuf1,buf1)
    
    ! write(6,*) "job,it = ",job,it
    ! write(6,*) "np = ",np
    ! write(6,*) "time = ",time
    call allocate_particle_data(np)

    dims1(1) = ipu
    call H5LTread_dataset_double_f(file_id,"/dm_p",dm_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/ut1_p",ut1_p(1:ipu),dims1,error)
    call H5LTread_dataset_double_f(file_id,"/hut_p",hut_p(1:ipu),dims1,error)
    call H5LTread_dataset_int_f   (file_id,"/flag_evol",flag_evol(1:ipu),dims1,error)
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
    
    call h5fclose_f(file_id,error)
    
  end subroutine read_checkpoint

  subroutine output_position(fn,ipu,time)
    integer,intent(in) :: ipu
    character(*),intent(in) :: fn
    real(8),intent(in) :: time
    
    character(10) :: str1
    integer :: unit,ip
    
    open (newunit=unit,file=fn,status="replace")
    write(unit,'("# ",99es12.4)') time
    do ip = 1,ipu
       if(flag_evol(ip)==1)then
          ! write(unit,'(i10,99es10.2)') ip, x_p(ip),z_p(ip),qrho_p(ip),ye_p(ip),tem_p(ip),sen_p(ip),dm_p(ip)
          write(unit,'(i10,99es25.17)') ip, x_p(ip),z_p(ip),qrho_p(ip),ye_p(ip),tem_p(ip),sen_p(ip),dm_p(ip)
       endif
    enddo
    close(unit)
    
  end subroutine output_position

  
end module module_hdf_to_ascii
