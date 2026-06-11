subroutine sort_median(ntemp,np_active,fac_lw,fac_med,fac_rw,val_temp_traj_list,mass_temp_traj_list,val_lw_temp_list,val_med_temp_list,val_rw_temp_list)
  implicit none

  integer,intent(in) :: ntemp, np_active(ntemp)
  real(8),intent(in) :: fac_lw,fac_med,fac_rw
  real(8),intent(in) :: val_temp_traj_list(ntemp,np_active),mass_temp_traj_list(ntemp,np_active)
  real(8),intent(out) :: val_lw_temp_list(ntemp),val_med_temp_list(ntemp),val_rw_temp_list(ntemp)
  
  integer,allocatable :: ip_list(:)
  real(8) :: dummy_list(np_active)
  
  integer :: itemp, ip, ipp,ip_dummy
  integer :: ip_lw,ip_med,ip_rw
  integer :: ss_lw,ss_med,ss_rw
  
  ! median
  do itemp=1,ntemp

!     write(6,*) itemp
     np_active_itemp = np_active(itemp)
     allocate(ip_list(np_active_itemp))

     do ip=1,np_active_itemp
        ip_list(ip) = ip
     enddo
     do ip=1,np_active_itemp
        do ipp=ip,np_active_itemp
           if( val_temp_traj_list(itemp,ip_list(ip)) > val_temp_traj_list(itemp,ip_list(ipp)) )then
              ip_dummy = ip_list(ip)
              ip_list(ip)=ip_list(ipp)
              ip_list(ipp) = ip_dummy
           endif
        enddo
     enddo

     do ip=1,np_active_itemp
        dummy_list(ip) = 0d0
        do ipp=1,ip
           dummy_list(ip) = dummy_list(ip) + mass_temp_traj_list(itemp,ip_list(ipp))/sum(mass_temp_traj_list(itemp,:))
        enddo
     enddo
     
     ! do ip=1,np_active_itemp
     !    write(6,*) dummy_list(ip), val_temp_traj_list(itemp,ip_list(ip))
     ! enddo

     do ip=2,np_active_itemp
        if(dummy_list(ip) > fac_lw .and. fac_lw > dummy_list(ip-1) )then
           ip_lw= ip
           ss_lw = (fac_lw-dummy_list(ip-1))/(dummy_list(ip)-dummy_list(ip-1))
        endif
        if(dummy_list(ip) > fac_med .and. fac_med > dummy_list(ip-1) )then
           ip_med= ip
           ss_med = (fac_med-dummy_list(ip-1))/(dummy_list(ip)-dummy_list(ip-1))
        endif
        if(dummy_list(ip) > fac_rw .and. fac_rw > dummy_list(ip-1) )then
           ip_rw= ip
           ss_rw = (fac_rw-dummy_list(ip-1))/(dummy_list(ip)-dummy_list(ip-1))
        endif
     enddo

     val_lw_temp_list(itemp) = ss_lw*val_temp_traj_list(itemp,ip_list(ip_lw)) + (1d0-ss_lw)*val_temp_traj_list(itemp,ip_list(ip_lw-1))
     val_med_temp_list(itemp) = ss_med*val_temp_traj_list(itemp,ip_list(ip_med)) + (1d0-ss_med)*val_temp_traj_list(itemp,ip_list(ip_med-1))
     val_rw_temp_list(itemp) = ss_rw*val_temp_traj_list(itemp,ip_list(ip_rw)) + (1d0-ss_rw)*val_temp_traj_list(itemp,ip_list(ip_rw-1))

     deallocate(ip_list)
     
  enddo

end subroutine sort_median
