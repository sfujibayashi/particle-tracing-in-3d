module particle_data

  implicit none

  ! Their size is only np
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
         x_p   (np), y_p   (np), z_p   (np),&
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


  subroutine reallocate_particle_data(np_old,np_new)
    integer,intent(in) :: np_old,np_new
    integer :: ip
    
    integer,allocatable :: flag_evol_buf(:)
    real(8),allocatable :: x_buf(:),y_buf(:),z_buf(:),&
         qrho_buf(:),&
         ye_buf  (:),&
         tem_buf (:),&
         ut_buf  (:),&
         qb_buf  (:),&
         sen_buf (:),&
         vlx_buf (:),&
         vly_buf (:),&
         vlz_buf (:),&
         hhh_buf (:),&
         rne_buf (:),&
         rae_buf (:),&
         deptn_buf(:),&
         depta_buf(:),&
         dm_buf  (:), &
         ut1_buf (:), &
         hut_buf(:)
    
    
    allocate( &
         flag_evol_buf(np_old), &
         x_buf   (np_old),&
         y_buf   (np_old),&
         z_buf   (np_old),&
         qrho_buf(np_old),&
         ye_buf  (np_old),&
         tem_buf (np_old),&
         ut_buf  (np_old),&
         qb_buf  (np_old),&
         sen_buf (np_old),&
         vlx_buf (np_old),&
         vly_buf (np_old),&
         vlz_buf (np_old),&
         hhh_buf (np_old),&
         rne_buf (np_old),&
         rae_buf (np_old),&
         deptn_buf(np_old),&
         depta_buf(np_old),&
         dm_buf  (np_old), &
         ut1_buf (np_old), &
         hut_buf (np_old)  )
    
    do ip=1,np_old
       flag_evol_buf(ip) = flag_evol(ip)
       x_buf        (ip) = x_p      (ip)
       y_buf        (ip) = y_p      (ip)
       z_buf        (ip) = z_p      (ip)
       qrho_buf     (ip) = qrho_p   (ip)
       ye_buf       (ip) = ye_p     (ip)
       tem_buf      (ip) = tem_p    (ip)
       ut_buf       (ip) = ut_p     (ip)
       qb_buf       (ip) = qb_p     (ip)
       sen_buf      (ip) = sen_p    (ip)
       vlx_buf      (ip) = vlx_p    (ip)
       vly_buf      (ip) = vly_p    (ip)
       vlz_buf      (ip) = vlz_p    (ip)
       hhh_buf      (ip) = hhh_p    (ip)
       rne_buf      (ip) = rne_p    (ip)
       rae_buf      (ip) = rae_p    (ip)
       deptn_buf    (ip) = deptn_p  (ip)
       depta_buf    (ip) = depta_p  (ip)
       dm_buf       (ip) = dm_p     (ip)
       ut1_buf      (ip) = ut1_p    (ip)
       hut_buf      (ip) = hut_p    (ip)
    enddo

    deallocate( &
         flag_evol, &
         x_p   ,&
         y_p   ,&
         z_p   ,&
         qrho_p,&
         ye_p  ,&
         tem_p ,&
         ut_p  ,&
         qb_p  ,&
         sen_p ,&
         vlx_p ,&
         vly_p ,&
         vlz_p ,&
         hhh_p ,&
         rne_p ,&
         rae_p ,&
         deptn_p,&
         depta_p,&
         dm_p  , &
         ut1_p , &
         hut_p   )

    call allocate_particle_data(np_new)

    do ip=1,np_old
       flag_evol(ip) = flag_evol_buf(ip)
       x_p      (ip) = x_buf        (ip)
       y_p      (ip) = y_buf        (ip)
       z_p      (ip) = z_buf        (ip)
       qrho_p   (ip) = qrho_buf     (ip)
       ye_p     (ip) = ye_buf       (ip)
       tem_p    (ip) = tem_buf      (ip)
       ut_p     (ip) = ut_buf       (ip)
       qb_p     (ip) = qb_buf       (ip)
       sen_p    (ip) = sen_buf      (ip)
       vlx_p    (ip) = vlx_buf      (ip)
       vly_p    (ip) = vly_buf      (ip)
       vlz_p    (ip) = vlz_buf      (ip)
       hhh_p    (ip) = hhh_buf      (ip)
       rne_p    (ip) = rne_buf      (ip)
       rae_p    (ip) = rae_buf      (ip)
       deptn_p  (ip) = deptn_buf    (ip)
       depta_p  (ip) = depta_buf    (ip)
       dm_p     (ip) = dm_buf       (ip)
       ut1_p    (ip) = ut1_buf      (ip)
       hut_p    (ip) = hut_buf      (ip)
    enddo
    
    deallocate( &
         flag_evol_buf, &
         x_buf   ,&
         y_buf   ,&
         z_buf   ,&
         qrho_buf,&
         ye_buf  ,&
         tem_buf ,&
         ut_buf  ,&
         qb_buf  ,&
         sen_buf ,&
         vlx_buf ,&
         vly_buf ,&
         vlz_buf ,&
         hhh_buf ,&
         rne_buf ,&
         rae_buf ,&
         deptn_buf,&
         depta_buf,&
         dm_buf  , &
         ut1_buf , &
         hut_buf   )
    
  end subroutine reallocate_particle_data

  subroutine deallocate_particle_data
    
    if(allocated(flag_evol))then
       deallocate( &
         flag_evol, &
         x_p   , y_p   , z_p   ,&
         qrho_p,&
         ye_p  ,&
         tem_p ,&
         ut_p  ,&
         qb_p  ,&
         sen_p ,&
         vlx_p ,&
         vly_p ,&
         vlz_p ,&
         hhh_p ,&
         rne_p ,&
         rae_p ,&
         deptn_p,&
         depta_p,&
         dm_p  , &
         ut1_p , &
         hut_p   )
    else
       write(6,*) "particle data not allocated yet."
    endif
    
  end subroutine deallocate_particle_data
  
  subroutine output_position(fn,ipu,time)
    integer,intent(in) :: ipu
    character(*),intent(in) :: fn
    real(8),intent(in) :: time
    
    character(10) :: str1
    integer :: unit,ip
    
    open (newunit=unit,file=fn,status="replace",action="write")
    write(unit,'("# ",99es12.4)') time
    do ip = 1,ipu
       if(flag_evol(ip)==1)then
          ! write(unit,'(i10,99es10.2)') ip, x_p(ip),z_p(ip),qrho_p(ip),ye_p(ip),tem_p(ip),sen_p(ip),dm_p(ip)
          write(unit,'(i10,99es25.17)') ip, x_p(ip),z_p(ip),qrho_p(ip),ye_p(ip),tem_p(ip),sen_p(ip),dm_p(ip)
       endif
    enddo
    close(unit)
    
  end subroutine output_position

end module particle_data
