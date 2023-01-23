module module_rho_ye
  implicit none

contains
  subroutine make_histogram_rho_ye(t,it)
    use simdata3D
    use io
    use module_eos
    integer,intent(in) :: it
    real(8),intent(in) :: t

    real(8),allocatable :: histogram_rho_ye_total(:,:),histogram_rho_ye_outside_ah(:,:),histogram_rho_ye_ejecta(:,:),histogram_vel_ye_ejecta(:,:)

    real(8),allocatable :: hist_v_ye(:), histogram_ye_ejecta(:), histogram_ye_ejectah(:), histogram_ye_outside_ah(:)
    real(8),parameter :: dye=0.0010d0
    real(8),parameter :: ye_max = 0.55d0+0.5d0*dye, ye_min = 0.01d0 - 0.5*dye
    integer,parameter :: n_ye = nint((ye_max-ye_min)/dye)

    real(8),allocatable :: hist_v_vel(:), histogram_vel_ejecta(:), histogram_vel_ejectah(:)
    real(8),parameter :: dvel=0.005d0
    real(8),parameter :: vel_max = 1.d0, vel_min = 0.0d0
    integer,parameter :: n_vel = nint((vel_max-vel_min)/dvel)
    
    integer,parameter :: n_rho = 130
    real(8),allocatable :: hist_v_lrho(:)
    real(8),parameter :: lrho_max = log10(1d15), lrho_min = log10(1d2)
    real(8),parameter :: dlrho=(lrho_max-lrho_min)/dble(n_rho)
    

    real(8),allocatable :: histogram_tem_ye_total(:,:),histogram_tem_ye_outside_ah(:,:),histogram_tem_ye_bound(:,:),hist_v_ltem(:)
    real(8),allocatable :: histogram_rho_tem_total(:,:),histogram_rho_tem_outside_ah(:,:),histogram_rho_tem_bound(:,:)
    real(8),parameter :: ltem_max = log10(3d1), ltem_min = log10(1d-2)
    integer,parameter :: n_tem = 100

    real(8),parameter :: lsen_max = log10(1d2), lsen_min = log10(1d0)
    integer,parameter :: n_sen = 100
    real(8),allocatable :: histogram_rho_sen_total(:,:),histogram_rho_sen_outside_ah(:,:),histogram_rho_sen_bound(:,:),hist_v_lsen(:)


    real(8),parameter :: dltem=(ltem_max-ltem_min)/dble(n_tem)
    real(8),parameter :: dlsen=(lsen_max-lsen_min)/dble(n_sen)
    integer :: i_ye, i_vel, i_rho, i_tem, i_sen
    real(8) :: dm,vel

    integer :: j,k,l,lv

    integer :: unit_num
    character(256) :: str1
    

    allocate(histogram_rho_ye_total(n_rho,n_ye),histogram_rho_ye_outside_ah(n_rho,n_ye),histogram_rho_ye_ejecta(n_rho,n_ye),hist_v_lrho(n_rho))
    allocate(histogram_rho_sen_total(n_rho,n_sen),histogram_rho_sen_outside_ah(n_rho,n_sen),histogram_rho_sen_bound(n_rho,n_sen),hist_v_lsen(n_sen))
    allocate(histogram_rho_tem_total(n_rho,n_tem),histogram_rho_tem_outside_ah(n_rho,n_tem),histogram_rho_tem_bound(n_rho,n_tem))
    allocate(histogram_tem_ye_total(n_tem,n_ye),histogram_tem_ye_outside_ah(n_tem,n_ye),histogram_tem_ye_bound(n_tem,n_ye),hist_v_ltem(n_tem))

    allocate(histogram_vel_ye_ejecta(n_vel,n_ye))
    allocate(hist_v_ye(n_ye),histogram_ye_ejecta(n_ye),histogram_ye_ejectah(n_ye),histogram_ye_outside_ah(n_ye))
    allocate(hist_v_vel(n_vel),histogram_vel_ejecta(n_vel),histogram_vel_ejectah(n_vel))

    histogram_rho_ye_total(:,:)=0.d0
    histogram_rho_ye_outside_ah(:,:)=0.d0
    histogram_rho_ye_ejecta(:,:)=0.d0

    histogram_rho_sen_total(:,:)=0.d0
    histogram_rho_sen_outside_ah(:,:)=0.d0
    histogram_rho_sen_bound(:,:)=0.d0

    histogram_rho_tem_total(:,:)=0.d0
    histogram_rho_tem_outside_ah(:,:)=0.d0
    histogram_rho_tem_bound(:,:)=0.d0

    histogram_tem_ye_total(:,:)=0.d0
    histogram_tem_ye_outside_ah(:,:)=0.d0
    histogram_tem_ye_bound(:,:)=0.d0

    do i_rho=1,n_rho
       hist_v_lrho(i_rho) = lrho_min + dlrho*dble(i_rho-1)
       !write(6,*) i_rho, hist_v_lrho(i_rho)
    end do
    !write(6,*)
    do i_ye=1,n_ye
       hist_v_ye(i_ye) = ye_min + dye*dble(i_ye-1)
       !write(6,*) i_ye, hist_v_ye(i_ye)
    end do
    !write(6,*)
    do i_vel=1,n_vel
       hist_v_vel(i_vel) = vel_min + dvel*dble(i_vel-1)
       !write(6,*) i_vel, hist_v_vel(i_vel)
    end do
    do i_tem=1,n_tem
       hist_v_ltem(i_tem) = ltem_min + dltem*dble(i_tem-1)
    end do
    do i_sen=1,n_sen
       hist_v_lsen(i_sen) = lsen_min + dlsen*dble(i_sen-1)
    end do
    
    ! i_ye  = max(1,min(n_ye,int( (0.345d0-ye_min)/dye ))) + 1
    ! write(6,*) i_ye, hist_v_ye(i_ye), 0.345, hist_v_ye(i_ye+1)
    ! i_rho = max(1,min(n_rho,int( (log10(3.56d10)-lrho_min)/dlrho ))) + 1
    ! write(6,*) i_rho, hist_v_lrho(i_rho), log10(3.56d10), hist_v_lrho(i_rho+1)
    ! stop

    histogram_ye_outside_ah(:) = 0d0
    histogram_ye_ejecta(:) = 0d0
    histogram_vel_ejecta(:) = 0d0

    histogram_ye_ejectah(:) = 0d0
    histogram_vel_ejectah(:) = 0d0
    
    do lv=lv_min,lv_max

      !$omp parallel private(j,k,l,i_ye,i_rho,i_tem,i_sen, i_vel, dm, vel) &
      !$omp   reduction(+:histogram_rho_ye_total,histogram_rho_ye_ejecta,histogram_ye_ejecta,histogram_vel_ejecta)
      !$omp do
       do l=ld,lu
          do k=kd,ku
             do j=jd,ju

                dm = qb(j,k,l,lv)*vol3D(j,k,l,lv)
                
                i_ye  = max(1,min(n_ye,int( (ye(j,k,l,lv)-ye_min)/dye) + 1 ))
                i_rho = max(1,min(n_rho,int( (log10(qrho(j,k,l,lv))-lrho_min)/dlrho ) + 1 ))
                ! i_tem = max(1,min(n_tem,int( (log10( tem(j,k,l,lv))-ltem_min)/dltem ) + 1 ))
                ! i_sen = max(1,min(n_sen,int( (log10( sen(j,k,l,lv))-lsen_min)/dlsen ) + 1 ))
                !i_rho = max(1,min(n_rho,int( (log10(qb(j,k,l,lv))-lrho_min)/dlrho ))) + 1

                ! histogram_rho_ye_total(i_rho,i_ye) = histogram_rho_ye_total(i_rho,i_ye) + dm
                ! histogram_rho_tem_total(i_rho,i_tem) = histogram_rho_tem_total(i_rho,i_tem) + dm
                ! histogram_rho_sen_total(i_rho,i_sen) = histogram_rho_sen_total(i_rho,i_sen) + dm
                ! histogram_tem_ye_total(i_tem,i_ye) = histogram_tem_ye_total(i_tem,i_ye) + dm
                !if(flag_ah(j,k,l,lv) < 0.5d0)then
                   histogram_ye_outside_ah(i_ye) = histogram_ye_outside_ah(i_ye) + dm
                   histogram_rho_ye_outside_ah(i_rho,i_ye) = histogram_rho_ye_outside_ah(i_rho,i_ye) + dm
                   ! histogram_rho_tem_outside_ah(i_rho,i_tem) = histogram_rho_tem_outside_ah(i_rho,i_tem) + dm
                   ! histogram_rho_sen_outside_ah(i_rho,i_sen) = histogram_rho_sen_outside_ah(i_rho,i_sen) + dm
                   ! histogram_tem_ye_outside_ah(i_tem,i_ye) = histogram_tem_ye_outside_ah(i_tem,i_ye) + dm
                   ! if( ut(j,k,l,lv) > -1.d0 )then
                   
                   !    histogram_rho_tem_bound(i_rho,i_tem) = histogram_rho_tem_bound(i_rho,i_tem) + dm
                   !    histogram_rho_sen_bound(i_rho,i_sen) = histogram_rho_sen_bound(i_rho,i_sen) + dm
                   !    histogram_tem_ye_bound(i_tem,i_ye) = histogram_tem_ye_bound(i_tem,i_ye) + dm
                   ! endif
                   if( ut(j,k,l,lv) < -1.d0 )then
                      vel = sqrt(1d0-1d0/(-ut(j,k,l,lv))**2)
                      i_vel  = max(1,min(n_vel,int( (vel-vel_min)/dvel ) + 1 ))
                      histogram_vel_ejecta(i_vel) = histogram_vel_ejecta(i_vel) + dm
                      histogram_ye_ejecta(i_ye) = histogram_ye_ejecta(i_ye) + dm
                      histogram_rho_ye_ejecta(i_rho,i_ye) = histogram_rho_ye_ejecta(i_rho,i_ye) + dm
                      histogram_vel_ye_ejecta(i_vel,i_ye) = histogram_vel_ye_ejecta(i_vel,i_ye) + dm
                   endif
                   if( ut(j,k,l,lv)*hhh(j,k,l,lv)/hhh_min< -1.d0 )then
                      vel = sqrt(1d0-1d0/(-ut(j,k,l,lv)*hhh(j,k,l,lv)/hhh_min)**2)
                      i_vel  = max(1,min(n_vel,int( (vel-vel_min)/dvel ) + 1 ))
                      histogram_vel_ejectah(i_vel) = histogram_vel_ejectah(i_vel) + dm
                      histogram_ye_ejectah(i_ye) = histogram_ye_ejectah(i_ye) + dm
                   endif

                !endif
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel

    enddo

    ! write(6,*) sum(histogram_rho_ye_total(:,:)),sum(histogram_rho_ye_outside_ah(:,:)),sum(histogram_rho_ye_bound(:,:))
    !stop

    write(str1,'(i6.6)') it

    open(newunit=unit_num,file=trim(dir_out)//"/vel_ye_histogram.dat",status="replace")
    write(unit_num,'("# ",99es15.7)') t, sum(histogram_vel_ye_ejecta(:,:)), dvel,dye
    do i_ye=1,n_ye
       write(unit_num,*)
       do i_vel=1,n_vel
          write(unit_num,'(99es15.7)') hist_v_vel(i_vel),hist_v_ye(i_ye), histogram_vel_ye_ejecta(i_vel,i_ye)/dvel/dye
       enddo
    enddo
    close(unit_num)

    ! open(newunit=unit_num,file=trim(dir_out)//"/rho_ye_histogram_"//trim(str1)//".dat",status="replace")
    ! write(unit_num,'("# ",99es15.7)') t, sum(histogram_rho_ye_outside_ah(:,:)),sum(histogram_rho_ye_ejecta(:,:)), dlrho,dye
    ! do i_ye=1,n_ye
    !    write(unit_num,*)
    !    do i_rho=1,n_rho
    !       write(unit_num,'(99es15.7)') hist_v_lrho(i_rho),hist_v_ye(i_ye), histogram_rho_ye_outside_ah(i_rho,i_ye)/dlrho/dye, histogram_rho_ye_ejecta(i_rho,i_ye)/dlrho/dye
    !    enddo
    ! enddo
    ! close(unit_num)


    ! open(newunit=unit_num,file=trim(dir_out)//"/rho_sen_histogram_"//trim(str1)//".dat",status="replace")
    ! write(unit_num,'("# ",99es15.7)') t, sum(histogram_rho_sen_total(:,:)),sum(histogram_rho_tem_outside_ah(:,:)),sum(histogram_rho_sen_bound(:,:)), dlrho,dlsen
    ! do i_sen=1,n_sen
    !    write(unit_num,*)
    !    do i_rho=1,n_rho
    !       write(unit_num,'(99es15.7)') hist_v_lrho(i_rho),hist_v_lsen(i_sen), histogram_rho_sen_total(i_rho,i_sen)/dlrho/dlsen, histogram_rho_sen_outside_ah(i_rho,i_sen)/dlrho/dlsen, histogram_rho_sen_bound(i_rho,i_sen)/dlrho/dlsen
    !    enddo
    ! enddo
    ! close(unit_num)


    ! open(newunit=unit_num,file=trim(dir_out)//"/rho_tem_histogram_"//trim(str1)//".dat",status="replace")
    ! write(unit_num,'("# ",99es15.7)') t, sum(histogram_rho_tem_total(:,:)),sum(histogram_rho_tem_outside_ah(:,:)),sum(histogram_rho_tem_bound(:,:)), dlrho,dltem
    ! do i_tem=1,n_tem
    !    write(unit_num,*)
    !    do i_rho=1,n_rho
    !       write(unit_num,'(99es15.7)') hist_v_lrho(i_rho),hist_v_ltem(i_tem), histogram_rho_tem_total(i_rho,i_tem)/dlrho/dltem, histogram_rho_tem_outside_ah(i_rho,i_tem)/dlrho/dltem, histogram_rho_tem_bound(i_rho,i_tem)/dlrho/dltem
    !    enddo
    ! enddo
    ! close(unit_num)


    ! open(newunit=unit_num,file=trim(dir_out)//"/tem_ye_histogram_"//trim(str1)//".dat",status="replace")
    ! write(unit_num,'("# ",99es15.7)') t, sum(histogram_tem_ye_total(:,:)),sum(histogram_tem_ye_outside_ah(:,:)),sum(histogram_tem_ye_bound(:,:)), dltem,dye
    ! do i_ye=1,n_ye
    !    write(unit_num,*)
    !    do i_tem=1,n_tem
    !       write(unit_num,'(99es15.7)') hist_v_ltem(i_tem),hist_v_ye(i_ye), histogram_tem_ye_total(i_tem,i_ye)/dltem/dye, histogram_tem_ye_outside_ah(i_tem,i_ye)/dltem/dye, histogram_tem_ye_bound(i_tem,i_ye)/dltem/dye
    !    enddo
    ! enddo
    ! close(unit_num)


    open(newunit=unit_num,file=trim(dir_out)//"/ye_histogram.dat",status="replace")
    write(unit_num,'("# ",99es15.7)') t, sum(histogram_ye_ejecta(:)), sum(histogram_ye_ejectah(:)),sum(histogram_ye_outside_ah(:))
    do i_ye=1,n_ye
       write(unit_num,'(99es15.7)') hist_v_ye(i_ye), histogram_ye_ejecta(i_ye)/dye, histogram_ye_ejectah(i_ye)/dye, histogram_ye_outside_ah(i_ye)/dye
    enddo
    close(unit_num)

    open(newunit=unit_num,file=trim(dir_out)//"/vel_histogram.dat",status="replace")
    write(unit_num,'("# ",99es15.7)') t, sum(histogram_vel_ejecta(:)), sum(histogram_vel_ejectah(:))
    do i_vel=1,n_vel
       write(unit_num,'(99es15.7)') hist_v_vel(i_vel), histogram_vel_ejecta(i_vel)/dvel, histogram_vel_ejectah(i_vel)/dvel
    enddo
    close(unit_num)

    deallocate(histogram_rho_ye_total,histogram_rho_ye_outside_ah,histogram_rho_ye_ejecta,hist_v_lrho, hist_v_ye)
    deallocate(histogram_rho_sen_total,histogram_rho_sen_outside_ah,histogram_rho_sen_bound,hist_v_lsen)
    deallocate(histogram_rho_tem_total,histogram_rho_tem_outside_ah,histogram_rho_tem_bound)
    deallocate(histogram_tem_ye_total,histogram_tem_ye_outside_ah,histogram_tem_ye_bound,hist_v_ltem)
    
  end subroutine make_histogram_rho_ye
end module module_rho_ye
