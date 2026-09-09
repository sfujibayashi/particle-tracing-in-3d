program weak

  use unit
  use module_weak_interaction
  use module_etanu

  use, intrinsic :: ieee_arithmetic

  implicit none

  character(256) :: fn_eos, fn_eosb, dir_read, label
  character(256) :: fn_caprate, fn_caprate_neg, fn_etanu

  integer :: ntemp_in, nrho_in, nye_in
  integer :: itt_min, itt_max, np, ip, np_skip, np_start
  integer :: ntemp_list
  real(8) :: temp_min, temp_max
  logical :: mass_weighted

  real(8) :: buf(100)
  real(8), allocatable :: hut(:), ut1(:)

  block
    use inputparser

    call get_string_parameter("parameters", "dir_read", dir_read)
    call get_string_parameter("parameters", "label", label)
    call get_logical_parameter("parameters", "mass_weighted", mass_weighted)
    call get_string_parameter("parameters", "fn_eos", fn_eos)
    call get_string_parameter("parameters", "fn_eosb", fn_eosb)
    call get_integer_parameter("parameters", "nrho", nrho_in)
    call get_integer_parameter("parameters", "ntemp", ntemp_in)
    call get_integer_parameter("parameters", "nye", nye_in)

    write(6,'(a)') trim(fn_eos)
    write(6,*) nrho_in, ntemp_in, nye_in

    call get_integer_parameter("parameters", "ntemp_list", ntemp_list)
    call get_double_parameter("parameters", "temp_min", temp_min)
    call get_double_parameter("parameters", "temp_max", temp_max)
    write(6,*) ntemp_list, temp_min, temp_max

    itt_min = 1
    call get_integer_parameter("parameters", "itt_max", itt_max)
    call get_integer_parameter("parameters", "np", np)
    call get_integer_parameter("parameters", "np_start", np_start)
    call get_integer_parameter("parameters", "np_skip", np_skip)
    write(6,*) itt_min, itt_max
    write(6,*) np_start, np, np_skip

    call get_string_parameter("parameters", "fn_caprate", fn_caprate)
    write(6,'(a)') trim(fn_caprate)
    call get_string_parameter("parameters", "fn_caprate_neg", fn_caprate_neg)
    call get_string_parameter("parameters", "fn_etanu", fn_etanu)
  end block

  call init_weak_table(fn_caprate, fn_caprate_neg)
  call etanu_init(fn_etanu)

  block
    use module_eos
    call readeos(fn_eos, nrho_in, ntemp_in, nye_in)
  end block

  ! Preserve the original ana_traj.dat read, although these values are not
  ! currently used by calc_weak itself.
  allocate(ut1(np), hut(np))
  open(11, file=trim(dir_read)//"/ana_traj.dat", status="old", action="read")
  read(11,*)
  read(11,*)
  do ip = 1, np
    read(11,*) buf(1:30)
    ut1(ip) = buf(16)
    hut(ip) = buf(17)
  enddo
  close(11)

  ! Keep the same convention as the original code: reserve one extra slot
  ! so that the final EOF read cannot run beyond the allocated arrays.
  itt_max = itt_max + 1
  write(*,'("np, itt_min, itt_max=",3i7)') np, itt_min, itt_max
  write(*,'(a)') "Tracer loop uses OpenMP when compiled with OpenMP enabled."

  !$omp parallel do default(none) schedule(dynamic,4) &
  !$omp& shared(np_start,np,np_skip,dir_read,itt_min,itt_max) &
  !$omp& private(ip)
  do ip = np_start, np, np_skip
    call process_one_tracer(ip, dir_read, itt_min, itt_max)
  enddo
  !$omp end parallel do

contains

  subroutine process_one_tracer(ip, dir_read, itt_min, itt_max)

    implicit none

    integer, intent(in) :: ip, itt_min, itt_max
    character(*), intent(in) :: dir_read

    character(256) :: fn
    character(10) :: str1
    integer :: unit_traj, nunit_out
    integer :: it, nt, ios
    logical :: file_exists

    real(8) :: hoge, mass_dummy
    real(8) :: den_tmp, tem_tmp, ye_tmp
    real(8) :: ye_equil_cap, eta_dummy, xn_dummy, xp_dummy
    real(8) :: ecap_nrate, pcap_nrate
    real(8) :: r_tmp, vr_tmp, t_exp

    real(8), allocatable :: time(:)
    real(8), allocatable :: x_p(:), y_p(:), z_p(:)
    real(8), allocatable :: qrho_p(:), ye_p(:), tem_p(:)
    real(8), allocatable :: sen_p(:)
    real(8), allocatable :: vlx_p(:), vly_p(:), vlz_p(:)
    real(8), allocatable :: rne_p(:), rae_p(:)
    real(8), allocatable :: deptn_p(:), depta_p(:)
    real(8), allocatable :: yn_p(:), ya_p(:)

    allocate(time(itt_min:itt_max))
    allocate(x_p(itt_min:itt_max), y_p(itt_min:itt_max), z_p(itt_min:itt_max))
    allocate(qrho_p(itt_min:itt_max), ye_p(itt_min:itt_max), tem_p(itt_min:itt_max))
    allocate(sen_p(itt_min:itt_max))
    allocate(vlx_p(itt_min:itt_max), vly_p(itt_min:itt_max), vlz_p(itt_min:itt_max))
    allocate(rne_p(itt_min:itt_max), rae_p(itt_min:itt_max))
    allocate(deptn_p(itt_min:itt_max), depta_p(itt_min:itt_max))
    allocate(yn_p(itt_min:itt_max), ya_p(itt_min:itt_max))

    write(str1,'(i8.8)') ip
    fn = trim(dir_read)//"/traj_"//trim(str1)//".dat"

    !$omp critical(log_output)
    write(6,'(a)') trim(fn)
    !$omp end critical(log_output)

    open(newunit=unit_traj, file=fn, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: cannot open trajectory file:", trim(fn)
      !$omp end critical(log_output)
      return
    endif

    read(unit_traj,*,iostat=ios)
    if (ios == 0) read(unit_traj,*,iostat=ios)
    if (ios == 0) read(unit_traj,'(16x,es13.5)',iostat=ios) mass_dummy
    if (ios == 0) read(unit_traj,*,iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: incomplete trajectory header:", ip
      !$omp end critical(log_output)
      close(unit_traj)
      return
    endif

    nt = 0
    do
      if (nt + 1 > itt_max) then
        !$omp critical(log_output)
        write(6,*) "Error: trajectory exceeds itt_max:", ip, nt
        !$omp end critical(log_output)
        close(unit_traj)
        return
      endif

      it = nt + 1
      read(unit_traj,*,iostat=ios) &
           time(it), &
           x_p(it), &
           y_p(it), &
           z_p(it), &
           vlx_p(it), &
           vly_p(it), &
           vlz_p(it), &
           qrho_p(it), &
           tem_p(it), &
           ye_p(it), &
           sen_p(it), &
           rne_p(it), &
           rae_p(it), &
           deptn_p(it), &
           depta_p(it)

      if (ios < 0) exit
      if (ios > 0) then
        !$omp critical(log_output)
        write(6,*) "Error while reading trajectory:", ip, it
        !$omp end critical(log_output)
        close(unit_traj)
        return
      endif

      nt = it
    enddo
    close(unit_traj)

    ! Preserve the original misc-file read.  yn_p and ya_p are not currently
    ! used below, but keeping this makes the refactor behaviorally conservative.
    yn_p = 0d0
    ya_p = 0d0

    fn = trim(dir_read)//"/misc_"//trim(str1)//".dat"
    inquire(file=fn, exist=file_exists)
    if (file_exists) then
      open(newunit=unit_traj, file=fn, status="old", action="read", iostat=ios)
      if (ios /= 0) then
        !$omp critical(log_output)
        write(6,*) "Error: cannot open misc file:", trim(fn)
        !$omp end critical(log_output)
        return
      endif

      do it = 1, 4
        read(unit_traj,*,iostat=ios)
        if (ios /= 0) exit
      enddo

      if (ios == 0) then
        it = 0
        do
          if (it + 1 > itt_max) then
            !$omp critical(log_output)
            write(6,*) "Error: misc file exceeds itt_max:", ip, it
            !$omp end critical(log_output)
            close(unit_traj)
            return
          endif

          read(unit_traj,*,iostat=ios) hoge, yn_p(it+1), ya_p(it+1)
          if (ios < 0) exit
          if (ios > 0) then
            !$omp critical(log_output)
            write(6,*) "Error while reading misc file:", ip, it+1
            !$omp end critical(log_output)
            close(unit_traj)
            return
          endif
          it = it + 1
        enddo
      endif

      close(unit_traj)
    endif

    fn = "./weak/weak_"//trim(str1)//".dat"
    open(newunit=nunit_out, file=fn, status="replace", action="write", iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: cannot open output file:", trim(fn)
      !$omp end critical(log_output)
      return
    endif

    write(nunit_out, '("#",99a15)') &
         "t", "rho", "T", "Ye", "eta", "Xn", "Xp", &
         "lambda_ec*Xp", "lambda_pc*Xn", "r", "v^r", "t_exp", &
         "Ye(eq,cap)"

    do it = 1, nt
      den_tmp = qrho_p(it)
      tem_tmp = tem_p(it)/tem_uni
      ye_tmp  = ye_p(it)

      call ye_equilibrium_capture(den_tmp, tem_tmp, ye_equil_cap, &
           eta_dummy, xn_dummy, xp_dummy, ecap_nrate, pcap_nrate)
      call nrate_cap(den_tmp, tem_tmp, ye_tmp, eta_dummy, xn_dummy, &
           xp_dummy, ecap_nrate, pcap_nrate)

      r_tmp = sqrt(x_p(it)**2 + y_p(it)**2 + z_p(it)**2)
      vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + &
                z_p(it)*vlz_p(it))/r_tmp
      t_exp = r_tmp/abs(vr_tmp)

      write(nunit_out, '(" ",99es15.6e3)') &
           time(it), den_tmp, tem_tmp, ye_tmp, eta_dummy, xn_dummy, &
           xp_dummy, ecap_nrate, pcap_nrate, r_tmp, vr_tmp, t_exp, &
           ye_equil_cap
    enddo

    close(nunit_out)

  end subroutine process_one_tracer

end program weak
