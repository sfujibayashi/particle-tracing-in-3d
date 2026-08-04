program weak_freezeout

  use unit
  use, intrinsic :: ieee_arithmetic

  implicit none

  integer, parameter :: nresult = 18

  character(256) :: dir_read, fn, label, prefix, cumulative_dir
  integer :: itt_min, itt_max, np, ip, np_skip, np_start
  integer :: cumulative_stride, command_status
  logical :: write_cumulative
  integer :: nunit_out, nunit_out2, nunit_out3, nunit_out4

  real(8), allocatable :: mass_p(:)
  real(8), allocatable :: result_time(:,:), result_dye(:,:)
  real(8), allocatable :: result_dng(:,:), result_dnv(:,:)
  logical, allocatable :: valid(:)

  real(8) :: dye_fo, dng_fo, dnv_fo
  real(8) :: hoge(99)

  block
    use inputparser

    call get_string_parameter ("parameters", "dir_read", dir_read)
    call get_string_parameter ("parameters", "label",    label)

    itt_min = 1
    call get_integer_parameter("parameters", "itt_max",  itt_max)
    call get_integer_parameter("parameters", "np",       np)
    call get_integer_parameter("parameters", "np_start", np_start)
    call get_integer_parameter("parameters", "np_skip",  np_skip)

    call get_double_parameter("parameters", "dye_fo", dye_fo, 0.01d0)
    call get_double_parameter("parameters", "dng_fo", dng_fo, 0.01d0)
    call get_double_parameter("parameters", "dnv_fo", dnv_fo, 0.01d0)

    call get_string_parameter("parameters", "prefix", prefix, "")

    call get_logical_parameter("parameters", "write_cumulative", &
         write_cumulative,.false.)
    call get_integer_parameter("parameters", "cumulative_stride", &
         cumulative_stride,1)
    call get_string_parameter("parameters", "cumulative_dir", &
         cumulative_dir, "./cumulative")
  end block
  
  if (write_cumulative) then
     if (cumulative_stride < 1) then
        write(6,*) "Error: cumulative_stride must be positive:", &
             cumulative_stride
        stop 1
     endif
     
     call execute_command_line( &
          'mkdir -p "'//trim(cumulative_dir)//'"', &
          exitstat=command_status)
     if (command_status /= 0) then
        write(6,*) "Error: cannot create cumulative directory:", &
             trim(cumulative_dir)
        stop 1
     endif
  endif
  
  allocate(mass_p(np))
  
  open(11, file=trim(dir_read)//"/ana_traj.dat", &
       status="old", action="read")
  read(11,*)
  read(11,*)
  do ip = 1, np
     read(11,*) hoge(1:10)
     mass_p(ip) = hoge(3)
  enddo
  close(11)
  
  ! Keep the same convention as the original program: one extra array slot.
  itt_max = itt_max + 1
  write(*,'("np, itt_min, itt_max=",3i7)') np, itt_min, itt_max
  write(*,'(a)') "Tracer loop uses OpenMP when compiled with OpenMP enabled."

  allocate(result_time(nresult,np), result_dye(nresult,np))
  allocate(result_dng(nresult,np),  result_dnv(nresult,np))
  allocate(valid(np))

  result_time = 0d0
  result_dye  = 0d0
  result_dng  = 0d0
  result_dnv  = 0d0
  valid       = .false.

  if (len_trim(prefix) > 0) prefix = trim(prefix)//"_"

  fn = "./"//trim(prefix)//"weak_freezeout_timescale.dat"
  open(newunit=nunit_out, file=fn, status="replace", action="write")
  call write_header(nunit_out)

  fn = "./"//trim(prefix)//"weak_freezeout_dye.dat"
  open(newunit=nunit_out2, file=fn, status="replace", action="write")
  call write_header(nunit_out2)

  fn = "./"//trim(prefix)//"weak_freezeout_dng.dat"
  open(newunit=nunit_out3, file=fn, status="replace", action="write")
  call write_header(nunit_out3)

  fn = "./"//trim(prefix)//"weak_freezeout_dnv.dat"
  open(newunit=nunit_out4, file=fn, status="replace", action="write")
  call write_header(nunit_out4)

  ! Each tracer is independent.  The result arrays are shared, but each
  ! iteration writes only to its own ip column.  All work arrays inside
  ! process_one_tracer are local to the calling thread.
  !$omp parallel do default(none) schedule(dynamic,4) &
  !$omp& shared(np_start,np,np_skip,mass_p,dir_read,itt_min,itt_max) &
  !$omp& shared(dye_fo,dng_fo,dnv_fo) &
  !$omp& shared(write_cumulative,cumulative_stride,cumulative_dir) &
  !$omp& shared(result_time,result_dye,result_dng,result_dnv,valid) &
  !$omp& private(ip)
  do ip = np_start, np, np_skip
    call process_one_tracer( &
         ip, mass_p(ip), dir_read, itt_min, itt_max, &
         dye_fo, dng_fo, dnv_fo, &
         write_cumulative, cumulative_stride, cumulative_dir, &
         result_time(:,ip), result_dye(:,ip), &
         result_dng(:,ip), result_dnv(:,ip), valid(ip))
  enddo
  !$omp end parallel do

  ! Keep file output outside process_one_tracer.  This preserves particle
  ! ordering and avoids concurrent writes when the loop is parallelized.
  do ip = np_start, np, np_skip
    if (.not.valid(ip)) cycle

    write(nunit_out, '(" ",i15,99es15.6e3)') &
         ip, result_time(:,ip)
    write(nunit_out2,'(" ",i15,99es15.6e3)') &
         ip, result_dye(:,ip)
    write(nunit_out3,'(" ",i15,99es15.6e3)') &
         ip, result_dng(:,ip)
    write(nunit_out4,'(" ",i15,99es15.6e3)') &
         ip, result_dnv(:,ip)
  enddo

  close(nunit_out)
  close(nunit_out2)
  close(nunit_out3)
  close(nunit_out4)

contains

  subroutine write_header(nunit)
    implicit none

    integer, intent(in) :: nunit

    write(nunit, '("#",99a15)') &
         "ip", "mass", "t_FO", "x_FO", "y_FO", "z_FO", &
         "vx_FO", "vy_FO", "vz_FO", "rho_FO", "T_FO", &
         "Ye_FO", "Rec_FO", "Rpc_FO", "t_exp", "t_weak", &
         "dye(after FO)", "dng(after FO)", "dnv(after FO)"
  end subroutine write_header


  subroutine process_one_tracer( &
       ip, mass, dir_read, itt_min, itt_max, &
       dye_fo, dng_fo, dnv_fo, &
       write_cumulative, cumulative_stride, cumulative_dir, &
       result_time, result_dye, result_dng, result_dnv, valid)

    use, intrinsic :: ieee_arithmetic

    implicit none

    integer, intent(in) :: ip, itt_min, itt_max
    logical, intent(in) :: write_cumulative
    integer, intent(in) :: cumulative_stride
    real(8), intent(in) :: mass, dye_fo, dng_fo, dnv_fo
    character(*), intent(in) :: dir_read, cumulative_dir

    real(8), intent(out) :: result_time(nresult)
    real(8), intent(out) :: result_dye(nresult)
    real(8), intent(out) :: result_dng(nresult)
    real(8), intent(out) :: result_dnv(nresult)
    logical, intent(out) :: valid

    character(256) :: fn
    character(10) :: str1
    integer :: unit_traj, ios, it, nt
    integer :: it_fo_time, it_fo_dye, it_fo_dng, it_fo_dnv
    logical :: file_exists
    real(8) :: hoge(99), dt
    real(8) :: t_weak, t_exp, vr_tmp

    real(8), allocatable :: time(:)
    real(8), allocatable :: x_p(:), y_p(:), z_p(:)
    real(8), allocatable :: vlx_p(:), vly_p(:), vlz_p(:)
    real(8), allocatable :: qrho_p(:), tem_p(:), ye_p(:)
    real(8), allocatable :: rec_p(:), rpc_p(:)
    real(8), allocatable :: lambda_ec(:), lambda_pc(:)
    real(8), allocatable :: dye_rem(:), dng_rem(:), dnv_rem(:)

    result_time = 0d0
    result_dye  = 0d0
    result_dng  = 0d0
    result_dnv  = 0d0
    valid       = .false.

    allocate(time(itt_min:itt_max))
    allocate(x_p(itt_min:itt_max), y_p(itt_min:itt_max), &
             z_p(itt_min:itt_max))
    allocate(vlx_p(itt_min:itt_max), vly_p(itt_min:itt_max), &
             vlz_p(itt_min:itt_max))
    allocate(qrho_p(itt_min:itt_max), tem_p(itt_min:itt_max), &
             ye_p(itt_min:itt_max))
    allocate(rec_p(itt_min:itt_max), rpc_p(itt_min:itt_max))
    allocate(lambda_ec(itt_min:itt_max), lambda_pc(itt_min:itt_max))
    allocate(dye_rem(itt_min:itt_max), dng_rem(itt_min:itt_max), &
             dnv_rem(itt_min:itt_max))

    write(str1,'(i8.8)') ip

    ! Read the tracer trajectory.
    fn = trim(dir_read)//"/traj_"//trim(str1)//".dat"

    open(newunit=unit_traj, file=fn, status="old", &
         action="read", iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: cannot open trajectory file:", trim(fn)
      !$omp end critical(log_output)
      return
    endif

    do it = 1, 4
      read(unit_traj,*,iostat=ios)
      if (ios /= 0) then
        !$omp critical(log_output)
        write(6,*) "Error: incomplete trajectory header:", ip
        !$omp end critical(log_output)
        close(unit_traj)
        return
      endif
    enddo

    nt = 0
    do
      read(unit_traj,*,iostat=ios) hoge(1:10)
      if (ios < 0) exit
      if (ios > 0) then
        !$omp critical(log_output)
        write(6,*) "Error while reading trajectory:", ip, nt+1
        !$omp end critical(log_output)
        close(unit_traj)
        return
      endif

      nt = nt + 1
      if (nt > itt_max) then
        !$omp critical(log_output)
        write(6,*) "Error: trajectory exceeds itt_max:", ip, nt, itt_max
        !$omp end critical(log_output)
        close(unit_traj)
        return
      endif

      time(nt)   = hoge(1)
      x_p(nt)    = hoge(2)
      y_p(nt)    = hoge(3)
      z_p(nt)    = hoge(4)
      vlx_p(nt)  = hoge(5)
      vly_p(nt)  = hoge(6)
      vlz_p(nt)  = hoge(7)
      qrho_p(nt) = hoge(8)
      tem_p(nt)  = hoge(9)
      ye_p(nt)   = hoge(10)
    enddo
    close(unit_traj)

    if (nt < 1) then
      !$omp critical(log_output)
      write(6,*) "Error: empty trajectory:", ip
      !$omp end critical(log_output)
      return
    endif

    ! Read weak-interaction rates.
    rec_p     = 0d0
    rpc_p     = 0d0
    lambda_ec = 0d0
    lambda_pc = 0d0

    fn = "./weak/weak_"//trim(str1)//".dat"
    inquire(file=fn, exist=file_exists)
    if (.not.file_exists) then
      !$omp critical(log_output)
      write(6,*) "File does not exist:", trim(fn)
      !$omp end critical(log_output)
      return
    endif

    open(newunit=unit_traj, file=fn, status="old", &
         action="read", iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: cannot open weak-rate file:", trim(fn)
      !$omp end critical(log_output)
      return
    endif

    read(unit_traj,*,iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: incomplete weak-rate header:", ip
      !$omp end critical(log_output)
      close(unit_traj)
      return
    endif

    do it = 1, nt
      read(unit_traj,*,iostat=ios) hoge(1:9)
      if (ios /= 0) then
        !$omp critical(log_output)
        write(6,*) "Error: weak-rate row-count mismatch:", ip, it-1, nt
        !$omp end critical(log_output)
        close(unit_traj)
        return
      endif

      if (.not.ieee_is_finite(hoge(8))) then
        !$omp critical(log_output)
        write(6,*) "Warning: replacing invalid rec by zero:", &
             ip, it, hoge(8)
        !$omp end critical(log_output)
        hoge(8) = 0d0
      endif
      if (.not.ieee_is_finite(hoge(9))) then
        !$omp critical(log_output)
        write(6,*) "Warning: replacing invalid rpc by zero:", &
             ip, it, hoge(9)
        !$omp end critical(log_output)
        hoge(9) = 0d0
      endif

      rec_p(it) = hoge(8)
      rpc_p(it) = hoge(9)

      if (hoge(7) == 0d0) then
        lambda_ec(it) = 0d0
      else
        lambda_ec(it) = hoge(8)/hoge(7)
      endif

      if (hoge(6) == 0d0) then
        lambda_pc(it) = 0d0
      else
        lambda_pc(it) = hoge(9)/hoge(6)
      endif
    enddo
    close(unit_traj)

    ! Cumulative integrals from each time to the final time.
    dye_rem(nt) = 0d0
    dng_rem(nt) = 0d0
    dnv_rem(nt) = 0d0

    do it = nt-1, 1, -1
      dt = time(it+1) - time(it)

      dye_rem(it) = dye_rem(it+1) + 0.5d0*dt * &
           ((-rec_p(it)   + rpc_p(it)) + &
            (-rec_p(it+1) + rpc_p(it+1)))

      dng_rem(it) = dng_rem(it+1) + 0.5d0*dt * &
           ((rec_p(it)   + rpc_p(it)) + &
            (rec_p(it+1) + rpc_p(it+1)))

      dnv_rem(it) = dnv_rem(it+1) + 0.5d0*dt * &
           (abs(-rec_p(it)   + rpc_p(it)) + &
            abs(-rec_p(it+1) + rpc_p(it+1)))
    enddo

    if (write_cumulative) then
      call output_cumulative_history( &
           ip, nt, cumulative_dir, cumulative_stride, &
           time, tem_p, qrho_p, ye_p, &
           rec_p, rpc_p, lambda_ec, lambda_pc, &
           dye_rem, dng_rem, dnv_rem, &
           x_p, y_p, z_p, vlx_p, vly_p, vlz_p)
    endif

    ! Last point satisfying t_weak < t_exp during outward motion.
    it_fo_time = 0
    do it = nt, 1, -1
      call calculate_timescales( &
           it, x_p, y_p, z_p, vlx_p, vly_p, vlz_p, &
           lambda_ec, lambda_pc, t_exp, t_weak, vr_tmp)

      if (vr_tmp > 0d0 .and. t_weak < t_exp) then
        it_fo_time = it
        exit
      endif
    enddo

    if (it_fo_time == 0) then
      !$omp critical(log_output)
      write(6,*) "Warning: no timescale freeze-out point:", ip
      !$omp end critical(log_output)
      it_fo_time = nt
    endif

    ! First point after which |remaining net dYe| stays below dye_fo.
    it_fo_dye = 1
    do it = nt-1, 1, -1
      if (abs(dye_rem(it)) > dye_fo) then
        it_fo_dye = it + 1
        exit
      endif
    enddo

    ! First points where the monotonic gross/variation integrals are small.
    it_fo_dng = 1
    do it = 1, nt
      if (dng_rem(it) <= dng_fo) then
        it_fo_dng = it
        exit
      endif
    enddo

    it_fo_dnv = 1
    do it = 1, nt
      if (dnv_rem(it) <= dnv_fo) then
        it_fo_dnv = it
        exit
      endif
    enddo

    call make_result( &
         it_fo_time, mass, time, x_p, y_p, z_p, &
         vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, &
         rec_p, rpc_p, lambda_ec, lambda_pc, &
         dye_rem, dng_rem, dnv_rem, result_time)

    call make_result( &
         it_fo_dye, mass, time, x_p, y_p, z_p, &
         vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, &
         rec_p, rpc_p, lambda_ec, lambda_pc, &
         dye_rem, dng_rem, dnv_rem, result_dye)

    call make_result( &
         it_fo_dng, mass, time, x_p, y_p, z_p, &
         vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, &
         rec_p, rpc_p, lambda_ec, lambda_pc, &
         dye_rem, dng_rem, dnv_rem, result_dng)

    call make_result( &
         it_fo_dnv, mass, time, x_p, y_p, z_p, &
         vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, &
         rec_p, rpc_p, lambda_ec, lambda_pc, &
         dye_rem, dng_rem, dnv_rem, result_dnv)

    valid = .true.
  end subroutine process_one_tracer


  subroutine output_cumulative_history( &
       ip, nt, output_dir, output_stride, &
       time, tem_p, qrho_p, ye_p, &
       rec_p, rpc_p, lambda_ec, lambda_pc, &
       dye_rem, dng_rem, dnv_rem, &
       x_p, y_p, z_p, vlx_p, vly_p, vlz_p)

    implicit none

    integer, intent(in) :: ip, nt, output_stride
    character(*), intent(in) :: output_dir
    real(8), intent(in) :: time(:), tem_p(:), qrho_p(:), ye_p(:)
    real(8), intent(in) :: rec_p(:), rpc_p(:)
    real(8), intent(in) :: lambda_ec(:), lambda_pc(:)
    real(8), intent(in) :: dye_rem(:), dng_rem(:), dnv_rem(:)

    real(8), intent(in) :: x_p(:), y_p(:), z_p(:)
    real(8), intent(in) :: vlx_p(:), vly_p(:), vlz_p(:)

    character(256) :: fn
    character(10) :: str1
    integer :: it, nunit_cumulative, ios, stride
    real(8) :: t_exp, t_weak, vr_tmp, r_tmp
    integer :: i

    write(str1,'(i8.8)') ip
    fn = trim(output_dir)//"/cumulative_"//trim(str1)//".dat"

    open(newunit=nunit_cumulative, file=trim(fn), &
         status="replace", action="write", iostat=ios)
    if (ios /= 0) then
      !$omp critical(log_output)
      write(6,*) "Error: cannot open cumulative file:", trim(fn)
      !$omp end critical(log_output)
      return
    endif

    
    write(nunit_cumulative,'("#",i10,99i15)') (i,i=1,20)
    write(nunit_cumulative,'("#",a10,99a15)') &
         "it", "time", "T", "rho", "Ye", "Rec", "Rpc", "lambda_ec", "lambda_pc", "t_weak", &
         "rate_net", "rate_gross", "rate_variation", &
         "dye_rem", "abs_dye_rem", "dng_rem", "dnv_rem", &
         "t_exp", "v^r", "r"

    stride = max(1, output_stride)

    do it = 1, nt, stride
       call calculate_timescales( &
            it, x_p, y_p, z_p, vlx_p, vly_p, vlz_p, &
            lambda_ec, lambda_pc, t_exp, t_weak, vr_tmp)
       r_tmp = sqrt(x_p(it)**2 + y_p(it)**2 + z_p(it)**2)
       call write_cumulative_row( &
            nunit_cumulative, it, time, tem_p, qrho_p, ye_p, &
            rec_p, rpc_p, lambda_ec, lambda_pc, &
            dye_rem, dng_rem, dnv_rem, &
            t_exp, t_weak, vr_tmp, r_tmp)
    enddo
    
    ! Always include the final point, where all cumulative integrals are zero.
    if (mod(nt-1,stride) /= 0) then
       it=nt
       call calculate_timescales( &
            it, x_p, y_p, z_p, vlx_p, vly_p, vlz_p, &
            lambda_ec, lambda_pc, t_exp, t_weak, vr_tmp)
       r_tmp = sqrt(x_p(it)**2 + y_p(it)**2 + z_p(it)**2)
       call write_cumulative_row( &
            nunit_cumulative, it, time, tem_p, qrho_p, ye_p, &
            rec_p, rpc_p, lambda_ec, lambda_pc, &
            dye_rem, dng_rem, dnv_rem, &
            t_exp, t_weak, vr_tmp, r_tmp)
    endif

    close(nunit_cumulative)
  end subroutine output_cumulative_history


  subroutine write_cumulative_row( &
       nunit, it, time, tem_p, qrho_p, ye_p, &
       rec_p, rpc_p, lambda_ec, lambda_pc, &
       dye_rem, dng_rem, dnv_rem, &
       t_exp, t_weak, vr_tmp, r_tmp)

    implicit none

    integer, intent(in) :: nunit, it
    real(8), intent(in) :: time(:), tem_p(:), qrho_p(:), ye_p(:)
    real(8), intent(in) :: rec_p(:), rpc_p(:)
    real(8), intent(in) :: lambda_ec(:), lambda_pc(:)
    real(8), intent(in) :: dye_rem(:), dng_rem(:), dnv_rem(:)
    real(8), intent(in) :: t_exp, t_weak, vr_tmp, r_tmp

    real(8) :: rate_net, rate_gross, rate_variation

    rate_net       = -rec_p(it) + rpc_p(it)
    rate_gross     =  rec_p(it) + rpc_p(it)
    rate_variation = abs(rate_net)

    write(nunit,'(" ",i10,99es15.6e3)') &
         it, time(it), tem_p(it), qrho_p(it), ye_p(it), &
         rec_p(it), rpc_p(it), lambda_ec(it), lambda_pc(it), &
         t_weak, rate_net, rate_gross, rate_variation, &
         dye_rem(it), abs(dye_rem(it)), dng_rem(it), dnv_rem(it), &
         t_exp, vr_tmp, r_tmp
  end subroutine write_cumulative_row


  subroutine calculate_timescales( &
       it, x_p, y_p, z_p, vlx_p, vly_p, vlz_p, &
       lambda_ec, lambda_pc, t_exp, t_weak, vr_tmp)

    implicit none

    integer, intent(in) :: it
    real(8), intent(in) :: x_p(:), y_p(:), z_p(:)
    real(8), intent(in) :: vlx_p(:), vly_p(:), vlz_p(:)
    real(8), intent(in) :: lambda_ec(:), lambda_pc(:)
    real(8), intent(out) :: t_exp, t_weak, vr_tmp

    real(8) :: r_tmp, lambda_sum

    lambda_sum = lambda_ec(it) + lambda_pc(it)
    if (lambda_sum > 0d0) then
      t_weak = 1d0/lambda_sum
    else
      t_weak = huge(1d0)
    endif

    r_tmp = sqrt(x_p(it)**2 + y_p(it)**2 + z_p(it)**2)
    if (r_tmp > 0d0) then
      vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + &
                z_p(it)*vlz_p(it))/r_tmp
    else
      vr_tmp = 0d0
    endif

    if (abs(vr_tmp) > 0d0) then
      t_exp = r_tmp/abs(vr_tmp)
    else
      t_exp = huge(1d0)
    endif
  end subroutine calculate_timescales


  subroutine make_result( &
       it_fo, mass, time, x_p, y_p, z_p, &
       vlx_p, vly_p, vlz_p, qrho_p, tem_p, ye_p, &
       rec_p, rpc_p, lambda_ec, lambda_pc, &
       dye_rem, dng_rem, dnv_rem, result)

    implicit none

    integer, intent(in) :: it_fo
    real(8), intent(in) :: mass
    real(8), intent(in) :: time(:), x_p(:), y_p(:), z_p(:)
    real(8), intent(in) :: vlx_p(:), vly_p(:), vlz_p(:)
    real(8), intent(in) :: qrho_p(:), tem_p(:), ye_p(:)
    real(8), intent(in) :: rec_p(:), rpc_p(:)
    real(8), intent(in) :: lambda_ec(:), lambda_pc(:)
    real(8), intent(in) :: dye_rem(:), dng_rem(:), dnv_rem(:)
    real(8), intent(out) :: result(nresult)

    real(8) :: t_exp, t_weak, vr_tmp

    call calculate_timescales( &
         it_fo, x_p, y_p, z_p, vlx_p, vly_p, vlz_p, &
         lambda_ec, lambda_pc, t_exp, t_weak, vr_tmp)

    result(1)  = mass
    result(2)  = time(it_fo)
    result(3)  = x_p(it_fo)
    result(4)  = y_p(it_fo)
    result(5)  = z_p(it_fo)
    result(6)  = vlx_p(it_fo)
    result(7)  = vly_p(it_fo)
    result(8)  = vlz_p(it_fo)
    result(9)  = qrho_p(it_fo)
    result(10) = tem_p(it_fo)
    result(11) = ye_p(it_fo)
    result(12) = rec_p(it_fo)
    result(13) = rpc_p(it_fo)
    result(14) = t_exp
    result(15) = t_weak
    result(16) = dye_rem(it_fo)
    result(17) = dng_rem(it_fo)
    result(18) = dnv_rem(it_fo)
  end subroutine make_result

end program weak_freezeout
