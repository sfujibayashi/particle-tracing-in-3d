recursive subroutine recursive_evolution(lv,lv_min,lv_max,lvf1,lvf2,it_p,fvel_id,mode_backward,substep_max,ipu,famr_id,lvf2_limit)
  use hdf5
  use simdata3D, only: read_velocity, time_level, vlx, vly, vlz, vlx_b, vly_b, vlz_b, jd,ju,kd,ku,ld,lu
  use particle_data
  use io

  implicit none
  integer,intent(in) :: lv_min,lv_max,lvf1,lvf2
  integer,intent(in) :: lvf2_limit
  ! level evolved, timeslice of its parent level toward which particles are evolved.
  integer,intent(in) :: lv,it_p,ipu
  logical,intent(in) :: mode_backward
  integer,intent(out) :: substep_max
  
  INTEGER(HID_T),intent(in) :: fvel_id ! File identifier of velocity data
  INTEGER(HID_T),intent(in) :: famr_id ! File identifier of output
  
  character(256) :: str1
  integer :: it_min, it_max, it1, it2, step, it, it_skip

  integer :: ip, np_evolved
  logical,allocatable :: evolution_finished(:)
  real(8) :: time_prv

  allocate(evolution_finished(ipu))

  if(lv<=lvf1.or.lv>lvf2)then
     it_min = it_p
     it_max = it_p
  else
     it_min = it_p*2
     it_max = it_p*2 + 1
  endif

  it_skip = 1
  ! write(6,*) "lv,it_p,it_min,it_max,it_skip=", lv,it_p,it_min,it_max,it_skip
  
  if(mode_backward)then
     it1 = it_max
     it2 = it_min
     step = -1
  else
     it1 = it_min
     it2 = it_max
     step = +1  
  endif


  if(lv>lvf2_limit)then
     it1=it2
  endif


  do it=it1,it2,step*it_skip

     vlx_b(:,:,:,lv) = vlx(:,:,:,lv)
     vly_b(:,:,:,lv) = vly(:,:,:,lv)
     vlz_b(:,:,:,lv) = vlz(:,:,:,lv)
     time_prv = time_level(lv)
     
     call read_velocity(fvel_id,it,lv)
     write(nunit_timestep,'("read data",i4,i8,es15.7,99es12.4)') lv, it, time_level(lv),  time_level(lv_min:lv)-time_level(lv_min)
     ! write(6,'("read data",i4,i8,es15.7,99es12.4)') lv, it, time_level(lv),  time_level(lv_min:lv_max)-time_level(lv_min)
     
     ! block
     !   integer :: j,k,l
     !   if( lv==lv_max )then
     !      write(99,*) "#", it, time_prv, time_level(lv)
     !      l=ld+1
     !      do k=kd,ku,2
     !         do j=jd,ju,2
     !            write(99,'(3i4,99es12.4)') j,k,l,vlx(j,k,l,lv),vly(j,k,l,lv),vlz(j,k,l,lv),vlx_b(j,k,l,lv),vly_b(j,k,l,lv),vlz_b(j,k,l,lv)
     !         enddo
     !      enddo
     !      stop
     !   endif

     ! end block
     
     if(lv<lv_max)then
        call recursive_evolution(lv+1,lv_min,lv_max,lvf1,lvf2,it,fvel_id, mode_backward, substep_max, ipu, famr_id, lvf2_limit)
     endif
     
     ! np_evolved = 0; substep_max=0; remain(:)=0
     call evolution_particle_3D_levels(ipu,substep_max,lv,lv,np_evolved,evolution_finished)
     write(nunit_timestep,'("evolution done. lv=", i7, ", time=",es15.7, ", evolved=", i7, ", substep=", i5, ", remainings=", i7)') lv,time_level(lv), np_evolved,substep_max,sum(remain(:))

     
     output:block
       use hdf5
       character(256) :: str1,str2, str_group
       integer :: hdf_err
       INTEGER(HID_T) :: group_id, dset_id
       INTEGER(HID_T) :: dspace_id
       
       integer(HSIZE_T) :: dims1(1)
       integer :: ndat, idat
       real(8),allocatable :: fbuf1_1(:), fbuf1_2(:), fbuf1_3(:), fbuf1_4(:)
       integer,allocatable :: ibuf1(:)

       ndat=0
       do ip = 1,ipu
          if(evolution_finished(ip)) ndat = ndat+1
          !write(1000+ip,'(i5,es17.9,99es12.4)') lv,t_p(ip),x_p(ip),y_p(ip),z_p(ip)
       enddo
      
       write(str1,'(i10)') lv
       write(str2,'(i10)') it
       str_group = "level"//trim(adjustl(str1))//"/data"//trim(adjustl(str2))
       call h5gcreate_f(famr_id, trim(str_group), group_id, hdf_err)

       dims1(1) = 1
       CALL h5screate_simple_f(1, dims1, dspace_id, hdf_err)
       CALL h5dcreate_f(group_id, "Nevolved", H5T_NATIVE_INTEGER, dspace_id, &
            dset_id, hdf_err)
       CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ndat, dims1, hdf_err)
       CALL h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       CALL h5screate_simple_f(1, dims1, dspace_id, hdf_err)
       CALL h5dcreate_f(group_id, "time", H5T_NATIVE_DOUBLE, dspace_id, &
            dset_id, hdf_err)
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_level(lv), dims1, hdf_err)
       CALL h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)
       
       if(ndat>0)then

          dims1(1) = ndat
          allocate(ibuf1(ndat))
          allocate(fbuf1_1(ndat), fbuf1_2(ndat), fbuf1_3(ndat), fbuf1_4(ndat))
          
          ! list ip evolved
          idat = 1
          do ip = 1,ipu
             if(evolution_finished(ip))then
                ibuf1(idat) = ip
                fbuf1_1(idat) = x_p(ip)
                fbuf1_2(idat) = y_p(ip)
                fbuf1_3(idat) = z_p(ip)
                fbuf1_4(idat) = t_p(ip)

                idat = idat+1
             endif
          enddo
          call h5screate_simple_f(1, dims1, dspace_id, hdf_err)
          call h5dcreate_f(group_id, "list_evolved", H5T_NATIVE_INTEGER, &
               dspace_id, dset_id, hdf_err)
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,ibuf1, dims1, hdf_err)
          call h5dclose_f(dset_id, hdf_err)
          call h5sclose_f(dspace_id, hdf_err)

          ! position: x
          call h5screate_simple_f(1, dims1, dspace_id, hdf_err)
          call h5dcreate_f(group_id, "x_evolved", H5T_NATIVE_DOUBLE, &
               dspace_id, dset_id, hdf_err)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fbuf1_1, dims1, hdf_err)
          call h5dclose_f(dset_id, hdf_err)
          call h5sclose_f(dspace_id, hdf_err)

          ! position: y
          call h5screate_simple_f(1, dims1, dspace_id, hdf_err)
          call h5dcreate_f(group_id, "y_evolved", H5T_NATIVE_DOUBLE, &
               dspace_id, dset_id, hdf_err)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fbuf1_2, dims1, hdf_err)
          call h5dclose_f(dset_id, hdf_err)
          call h5sclose_f(dspace_id, hdf_err)

          ! position: z
          call h5screate_simple_f(1, dims1, dspace_id, hdf_err)
          call h5dcreate_f(group_id, "z_evolved", H5T_NATIVE_DOUBLE, &
               dspace_id, dset_id, hdf_err)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fbuf1_3, dims1, hdf_err)
          call h5dclose_f(dset_id, hdf_err)
          call h5sclose_f(dspace_id, hdf_err)

          ! position: t
          call h5screate_simple_f(1, dims1, dspace_id, hdf_err)
          call h5dcreate_f(group_id, "t_evolved", H5T_NATIVE_DOUBLE, &
               dspace_id, dset_id, hdf_err)
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fbuf1_4, dims1, hdf_err)
          call h5dclose_f(dset_id, hdf_err)
          call h5sclose_f(dspace_id, hdf_err)
          
          deallocate(ibuf1)
          deallocate(fbuf1_1, fbuf1_2, fbuf1_3, fbuf1_4)
       endif

       call h5gclose_f(group_id, hdf_err)
     end block output
     
  enddo
  
  deallocate(evolution_finished)
  
end subroutine recursive_evolution

