recursive subroutine recursive_evolution_to_end(lv,lv_min,lv_max,lvf1,lvf2,it_p,fvel_id,mode_backward,substep_max,ipu,famr_id,lvf2_limit)
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

  !write(6,*) "lv,it=",lv,it_p

  if(lv<lvf1)then
     it = it_p
     !write(6,*) "since", lv, "<", lvf1, ", go to the next level with it=",it
     call recursive_evolution_to_end(lv+1,lv_min,lv_max,lvf1,lvf2,it,fvel_id, mode_backward, substep_max, ipu, famr_id, lvf2_limit)
     return
  endif

  if(lv<lvf1)then
     it = it_p
  else
     it = it_p*2 - 1
  endif
  
  if(lv<lv_max.and.lv<lvf2_limit)then
     ! write(6,*) "call recursive evolution for lv, it_p = ",lv+1, it
     call recursive_evolution(lv+1,lv_min,lv_max,lvf1,lvf2,it,it,fvel_id, mode_backward, substep_max, ipu, famr_id, lvf2_limit)

     ! write(6,*) "call itself for a higher level",lv+1, it
     call recursive_evolution_to_end(lv+1,lv_min,lv_max,lvf1,lvf2,it,fvel_id, mode_backward, substep_max, ipu, famr_id, lvf2_limit)
  endif
  
  deallocate(evolution_finished)
  
end subroutine recursive_evolution_to_end

