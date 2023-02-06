subroutine partial_output_hdf(dir_out,job,it,t)
#include "macro.h"
  use hdf5
  use h5lt
  use simdata3D
  implicit none
  
  character(*),intent(in) :: dir_out
  real(8),intent(in)  :: t
  integer,intent(in)  :: job,it
  
  character(256) :: fn, str1,str2,str3

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: group_id, group2_id      ! Group identifier
  integer :: hdf_err
  integer(HSIZE_T) :: dims1(1),dims2(2),dims3(3)
  integer :: jdat,kdat,ldat, ld_write

  integer :: j,k,l,lv

  fn = trim(dir_out)//"/test.h5"
  call h5fcreate_f(fn, H5F_ACC_TRUNC_F, file_id, hdf_err)
  
#ifdef STAGGERED
#ifndef FULL
  ld_write=ld+1
#endif
#endif
    
  jdat = ju-jd+1
  kdat = ku-kd+1
  ldat = lu-ld_write+1
  
  do lv=lv_min,lv_max
     
     write(str1,'(i10)') lv
     call h5gcreate_f(file_id, "level"//trim(adjustl(str1)) , group_id, hdf_err)
     
     dims1(1) = jdat
     call h5ltmake_dataset_float_f(group_id, "x", 1, dims1, x(:,lv), hdf_err)
     dims1(1) = kdat
     call h5ltmake_dataset_float_f(group_id, "y", 1, dims1, y(:,lv), hdf_err)
     dims1(1) = ldat
     call h5ltmake_dataset_float_f(group_id, "z", 1, dims1, z(ld_write:lu,lv), hdf_err)
     
     ! t-dependent data

     write(str2,'(i10)') it
     call h5gcreate_f(group_id, "data"//trim(adjustl(str2)) , group2_id, hdf_err)
     
     dims1(1) = 1
     tms(1) = t*1d3
     call h5ltmake_dataset_float_f(group2_id,"time", 1, dims1, tms, hdf_err)

     dims3(1) = jdat
     dims3(2) = kdat
     dims3(3) = ldat

     call h5ltmake_dataset_float_f(group2_id, "density", 3, dims3, qrho(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "u_t", 3, dims3, ut(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "ye", 3, dims3, ye(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "entropy", 3, dims3, sen(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "temperature", 3, dims3, tem(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "vx", 3, dims3, vlx(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "vy", 3, dims3, vly(:,:,ld_write:lu,lv), hdf_err)
     call h5ltmake_dataset_float_f(group2_id, "vz", 3, dims3, vlz(:,:,ld_write:lu,lv), hdf_err)
     
     call h5gclose_f(group_id, hdf_err)
     call h5gclose_f(group2_id, hdf_err)

  enddo

  call h5fclose_f(file_id, hdf_err)

  
end subroutine partial_output_hdf
