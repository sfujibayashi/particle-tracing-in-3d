program weak_freezeout

  use unit

  use, intrinsic :: ieee_arithmetic

  implicit none

  real(8),parameter :: emev=0.51099996d0
  real(8),parameter :: mev2t9= 1.160445d1

  character(256) :: dir_read,fn,fn_out,label
  character(10) :: str1

  integer :: itt_min,itt_max,np,ip,it,nt,np_skip,np_start
  real(8),allocatable :: time(:),mass_p(:)
  real(8),allocatable :: x_p(:),y_p(:),z_p(:),&
       vlx_p (:),&
       vly_p (:),&
       vlz_p (:),&
       qrho_p(:),&
       tem_p (:),&
       ye_p  (:),&
       qb_p  (:),&
       rpc_p (:),&
       rec_p (:)

  real(8),allocatable :: lambda_ec(:), lambda_pc(:)

  integer :: unit_traj, nunit_out, nunit_out2, nunit_out3, nunit_out4
  logical :: file_exists
  integer :: idx
  real(8) :: hoge(99)

  real(8) :: dye_fo, dn_fo


  block
    use inputparser
    call get_string_parameter("parameters","dir_read",dir_read)
    call get_string_parameter("parameters","label",label)
    itt_min = 1
    call get_integer_parameter("parameters","itt_max",itt_max)
    call get_integer_parameter("parameters","np",np)
    call get_integer_parameter("parameters","np_start",np_start)
    call get_integer_parameter("parameters","np_skip",np_skip)

    call get_double_parameter("parameters","dye_fo",dye_fo,0.01d0)
    call get_double_parameter("parameters","dn_fo",dn_fo,0.01d0)
  end block

  allocate( mass_p(np) )
       

  open(11,file=trim(dir_read)//"/ana_traj.dat",status="old",action="read")
  read(11,*);read(11,*)
  do ip=1,np
     read(11,*) hoge(1:10)
     mass_p(ip) = hoge(3)
  enddo
  close(11)
  
  itt_max = itt_max + 1
  write(*,'("np, itt_min, itt_max=",3i7)') np,itt_min,itt_max

  allocate( &
       time     (itt_min:itt_max),&
       x_p      (itt_min:itt_max),&
       y_p      (itt_min:itt_max),&
       z_p      (itt_min:itt_max),&
       vlx_p    (itt_min:itt_max),&
       vly_p    (itt_min:itt_max),&
       vlz_p    (itt_min:itt_max),&
       qrho_p   (itt_min:itt_max),&
       tem_p   (itt_min:itt_max),&
       ye_p   (itt_min:itt_max),&
       rpc_p    (itt_min:itt_max),&
       rec_p    (itt_min:itt_max) )

  allocate(lambda_ec (itt_min:itt_max) )
  allocate(lambda_pc (itt_min:itt_max) )
  
  fn = "./weak_freezeout_timescale.dat"
  open(newunit=nunit_out, file=fn, status="replace", action="write")
  write(nunit_out, '("#",99a15)') "ip", "mass", "t_FO", "x_FO", "y_FO", "z_FO", "vx_FO", "vy_FO", "vz_FO", "rho_FO", "T_FO", "Ye_FO", "Rec_FO", "Rpc_FO", "t_exp", "t_weak", "dye(after FO)", "dng(after FO)", "dnv(after FO)"

  fn = "./weak_freezeout_dye.dat"
  open(newunit=nunit_out2, file=fn, status="replace", action="write")
  write(nunit_out2, '("#",99a15)') "ip", "mass", "t_FO", "x_FO", "y_FO", "z_FO", "vx_FO", "vy_FO", "vz_FO", "rho_FO", "T_FO", "Ye_FO", "Rec_FO", "Rpc_FO", "t_exp", "t_weak", "dye(after FO)", "dng(after FO)", "dnv(after FO)"

  fn = "./weak_freezeout_dng.dat"
  open(newunit=nunit_out3, file=fn, status="replace", action="write")
  write(nunit_out3, '("#",99a15)') "ip", "mass", "t_FO", "x_FO", "y_FO", "z_FO", "vx_FO", "vy_FO", "vz_FO", "rho_FO", "T_FO", "Ye_FO", "Rec_FO", "Rpc_FO", "t_exp", "t_weak", "dye(after FO)", "dng(after FO)", "dnv(after FO)"

  fn = "./weak_freezeout_dnv.dat"
  open(newunit=nunit_out4, file=fn, status="replace", action="write")
  write(nunit_out4, '("#",99a15)') "ip", "mass", "t_FO", "x_FO", "y_FO", "z_FO", "vx_FO", "vy_FO", "vz_FO", "rho_FO", "T_FO", "Ye_FO", "Rec_FO", "Rpc_FO", "t_exp", "t_weak", "dye(after FO)", "dng(after FO)", "dnv(after FO)"
  
  do ip=np_start,np,np_skip
     
     write(str1,'(i8.8)') ip
     fn = trim(dir_read)//"/traj_"//trim(str1)//".dat"
     write(6,'(a)') trim(fn)
     open(newunit=unit_traj,file=fn,status="old",action="read")
     read(unit_traj,*)
     read(unit_traj,*)
     read(unit_traj,*)
     read(unit_traj,*)
     it = 0
     do
        nt = it
        it = it + 1
        read(unit_traj,*,end=99) hoge(1:10)
        time(it) = hoge(1)
        x_p(it) = hoge(2)
        y_p(it) = hoge(3)
        z_p(it) = hoge(4)
        vlx_p(it) = hoge(5)
        vly_p(it) = hoge(6)
        vlz_p(it) = hoge(7)
        qrho_p(it) = hoge(8)
        tem_p(it) = hoge(9)
        ye_p(it) = hoge(10)
     enddo
99   continue
     close(unit_traj)
     
!!! misc. file
     rpc_p(:) = 0d0
     rec_p(:) = 0d0
     lambda_ec(:) = 0d0
     lambda_pc(:) = 0d0

     fn = "./weak/weak_"//trim(str1)//".dat"
     inquire(file=fn, exist=file_exists)
     if(.not.file_exists)then
        write(6,*) "file not exists"
        write(6,'(a)') trim(fn)
        cycle
     endif
        
     open(newunit=unit_traj,file=fn,status="old")
     read(unit_traj,*)

     do it=1,nt
        read(unit_traj,*,end=991) hoge(1:9)
        if (.not.ieee_is_finite(hoge(8))) then
           write(6,*) "Warning: replacing invalid rec by zero:",ip,it,hoge(8)
           hoge(8) = 0d0
        endif
        if (.not.ieee_is_finite(hoge(9))) then
           write(6,*) "Warning: replacing invalid rpc by zero:",ip,it,hoge(9)
           hoge(9) = 0d0
        endif
        
        rec_p(it) = hoge(8)
        rpc_p(it) = hoge(9)

        if(hoge(7) == 0d0)then
           lambda_ec(it) = 0d0
        else
           lambda_ec(it) = hoge(8)/hoge(7)
        endif

        if(hoge(6) == 0d0)then
           lambda_pc(it) = 0d0
        else
           lambda_pc(it) = hoge(9)/hoge(6)
        endif
     enddo
991  continue
     close(unit_traj)

     block
       real(8) :: t_weak, r_tmp, vr_tmp, t_exp
       integer :: it_fo, it_fo_dye, it_fo_dng, it_fo_dnv
       
       real(8) :: dt, dye, dng, dnv
       real(8) :: dye_after_fo , dng_after_fo , dnv_after_fo , t_weak_fo , t_exp_fo
       real(8) :: dye_after_fo2, dng_after_fo2, dnv_after_fo2, t_weak_dye, t_exp_dye
       real(8) :: dye_after_fo3, dng_after_fo3, dnv_after_fo3, t_weak_dng, t_exp_dng
       real(8) :: dye_after_fo4, dng_after_fo4, dnv_after_fo4, t_weak_dnv, t_exp_dnv

       it_fo     = 0
       find_fo:do it=nt,1,-1
          !write(6,*)it, ye_p(it),rec_p(it),rpc_p(it), x_p(it)
          ! t_weak = ye_p(it)/(rec_p(it) + rpc_p(it))
          t_weak = 1d0/(lambda_ec(it) + lambda_pc(it))
          
          r_tmp =  sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2)
          vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/r_tmp
          t_exp = r_tmp/abs(vr_tmp)
          
          if(vr_tmp > 0d0 .and. t_weak < t_exp)then
             it_fo = it
             exit find_fo
          endif
       enddo find_fo
       if (it_fo == 0) then
          write(6,*) "Warning: no timescale freeze-out point:", ip
          it_fo = nt
       endif

       it=it_fo
       ! t_weak_fo = ye_p(it)/(rec_p(it) + rpc_p(it))
       t_weak_fo = 1d0/(lambda_ec(it) + lambda_pc(it))
       r_tmp =  sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2)
       vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/r_tmp
       t_exp_fo = r_tmp/abs(vr_tmp)

       
       dye = 0d0
       dng = 0d0
       dnv = 0d0
       do it=it_fo,nt-1
          dt = time(it+1) - time(it)
          
          if (dt <= 0d0) then
             write(6,*) "Warning: non-positive dt:", ip, it, dt
             cycle
          endif
          
          dye = dye + 0.5d0*dt * &
               ( (-rec_p(it)   + rpc_p(it))  + &
                 (-rec_p(it+1) + rpc_p(it+1)) )
          dng  = dng + 0.5d0*dt * &
               ( (rec_p(it)   + rpc_p(it))  + &
                 (rec_p(it+1) + rpc_p(it+1)) )
          dnv = dnv + 0.5d0*dt * &
               ( abs(-rec_p(it)   + rpc_p(it))  + &
                 abs(-rec_p(it+1) + rpc_p(it+1)) )

       enddo
       dye_after_fo = dye
       dng_after_fo = dng
       dnv_after_fo = dnv

       ! FO condition 2 (dYe)
       it_fo_dye = 0
       dye = 0d0
       find_fo_dye:do it=nt-1,1,-1
          dt = time(it+1) - time(it)
          dye = dye + 0.5d0*dt * &
               ( (-rec_p(it)   + rpc_p(it)) + &
                 (-rec_p(it+1) + rpc_p(it+1)) )
          if( abs(dye) > dye_fo)then
             it_fo_dye = it
             exit find_fo_dye
          endif
       enddo find_fo_dye
       
       if (it_fo_dye == 0) then
          write(6,*) "Warning: total weak exposure is below dye_fo:", ip, dye
          it_fo_dye = 1
       endif

       dng = 0d0
       dnv = 0d0
       do it=it_fo_dye,nt-1
          dt = time(it+1) - time(it)
          dng  = dng + 0.5d0*dt * &
               ( (rec_p(it)   + rpc_p(it))  + &
                 (rec_p(it+1) + rpc_p(it+1)) )
          dnv = dnv + 0.5d0*dt * &
               ( abs(-rec_p(it)   + rpc_p(it))  + &
                 abs(-rec_p(it+1) + rpc_p(it+1)) )
       enddo
       dye_after_fo2 = dye
       dng_after_fo2 = dng
       dnv_after_fo2 = dnv


       it=it_fo_dye
       t_weak_dye = 1d0/(lambda_ec(it) + lambda_pc(it))
       r_tmp =  sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2)
       vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/r_tmp
       t_exp_dye = r_tmp/abs(vr_tmp)

       ! FO condition 3 (dNgross)
       it_fo_dng = 0
       dng = 0d0
       find_fo_dng:do it=nt-1,1,-1
          dt = time(it+1) - time(it)
          dng = dng + 0.5d0*dt * &
               ( ( rec_p(it)   + rpc_p(it)) + &
                 ( rec_p(it+1) + rpc_p(it+1)) )
          if( abs(dng) > dn_fo)then
             it_fo_dng = it
             exit find_fo_dng
          endif
       enddo find_fo_dng

       if (it_fo_dng == 0) then
          write(6,*) "Warning: total weak exposure is below dng_fo:", ip, dng
          it_fo_dng = 1
       endif

       dye = 0d0
       dnv = 0d0
       do it=it_fo_dng,nt-1
          dt = time(it+1) - time(it)
          dye = dye + 0.5d0*dt * &
               ( (-rec_p(it)   + rpc_p(it))  + &
                 (-rec_p(it+1) + rpc_p(it+1)) )
          dnv = dnv + 0.5d0*dt * &
               ( abs(-rec_p(it)   + rpc_p(it))  + &
                 abs(-rec_p(it+1) + rpc_p(it+1)) )

       enddo
       dye_after_fo3 = dye
       dng_after_fo3 = dng
       dnv_after_fo3 = dnv

       it=it_fo_dng
       ! t_weak_dng = ye_p(it)/(rec_p(it) + rpc_p(it))
       t_weak_dng = 1d0/(lambda_ec(it) + lambda_pc(it))
       r_tmp =  sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2)
       vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/r_tmp
       t_exp_dng = r_tmp/abs(vr_tmp)


       ! FO condition 4 (dNvariation)
       it_fo_dnv = 0
       dnv = 0d0
       find_fo_dnv:do it=nt-1,1,-1
          dt = time(it+1) - time(it)
          dnv = dnv + 0.5d0*dt * &
               ( abs(-rec_p(it)   + rpc_p(it)) + &
                 abs(-rec_p(it+1) + rpc_p(it+1)) )
          if( abs(dnv) > dn_fo)then
             it_fo_dnv = it
             exit find_fo_dnv
          endif
       enddo find_fo_dnv

       if (it_fo_dnv == 0) then
          write(6,*) "Warning: total weak exposure is below dnv_fo:", ip, dnv
          it_fo_dnv = 1
       endif

       dye = 0d0
       dng = 0d0
       do it=it_fo_dnv,nt-1
          dt = time(it+1) - time(it)
          dye = dye + 0.5d0*dt * &
               ( (-rec_p(it)   + rpc_p(it))  + &
                 (-rec_p(it+1) + rpc_p(it+1)) )
          dng = dng + 0.5d0*dt * &
               ( (rec_p(it)   + rpc_p(it))  + &
                 (rec_p(it+1) + rpc_p(it+1)) )
       enddo
       dye_after_fo4 = dye
       dng_after_fo4 = dng
       dnv_after_fo4 = dnv

       it=it_fo_dnv
       ! t_weak_dnv = ye_p(it)/(rec_p(it) + rpc_p(it))
       t_weak_dnv = 1d0/(lambda_ec(it) + lambda_pc(it))
       r_tmp =  sqrt(x_p(it)**2+y_p(it)**2+z_p(it)**2)
       vr_tmp = (x_p(it)*vlx_p(it) + y_p(it)*vly_p(it) + z_p(it)*vlz_p(it))/r_tmp
       t_exp_dnv = r_tmp/abs(vr_tmp)


       write(nunit_out, '(" ",i15,99es15.6e3)') ip, mass_p(ip), time(it_fo), x_p(it_fo), y_p(it_fo), z_p(it_fo), vlx_p(it_fo), vly_p(it_fo), vlz_p(it_fo), qrho_p(it_fo), tem_p(it_fo), ye_p(it_fo), rec_p(it_fo), rpc_p(it_fo), t_exp_fo, t_weak_fo, dye_after_fo, dng_after_fo, dnv_after_fo

       write(nunit_out2,'(" ",i15,99es15.6e3)') ip, mass_p(ip), time(it_fo_dye), x_p(it_fo_dye), y_p(it_fo_dye), z_p(it_fo_dye), vlx_p(it_fo_dye), vly_p(it_fo_dye), vlz_p(it_fo_dye), qrho_p(it_fo_dye), tem_p(it_fo_dye), ye_p(it_fo_dye), rec_p(it_fo_dye), rpc_p(it_fo_dye), t_exp_dye, t_weak_dye, dye_after_fo2, dng_after_fo2, dnv_after_fo2

       write(nunit_out3,'(" ",i15,99es15.6e3)') ip, mass_p(ip), time(it_fo_dng), x_p(it_fo_dng), y_p(it_fo_dng), z_p(it_fo_dng), vlx_p(it_fo_dng), vly_p(it_fo_dng), vlz_p(it_fo_dng), qrho_p(it_fo_dng), tem_p(it_fo_dng), ye_p(it_fo_dng), rec_p(it_fo_dng), rpc_p(it_fo_dng), t_exp_dng, t_weak_dng, dye_after_fo3, dng_after_fo3, dnv_after_fo3

       write(nunit_out4,'(" ",i15,99es15.6e3)') ip, mass_p(ip), time(it_fo_dnv), x_p(it_fo_dnv), y_p(it_fo_dnv), z_p(it_fo_dnv), vlx_p(it_fo_dnv), vly_p(it_fo_dnv), vlz_p(it_fo_dnv), qrho_p(it_fo_dnv), tem_p(it_fo_dnv), ye_p(it_fo_dnv), rec_p(it_fo_dnv), rpc_p(it_fo_dnv), t_exp_dnv, t_weak_dnv, dye_after_fo4, dng_after_fo4, dnv_after_fo4
       
       ! write(6, '(" ",i15,99es15.6e3)') ip, mass_p(ip), time(it_fo), x_p(it_fo), y_p(it_fo), z_p(it_fo), vlx_p(it_fo), vly_p(it_fo), vlz_p(it_fo), qrho_p(it_fo), tem_p(it_fo), ye_p(it_fo), rec_p(it_fo), rpc_p(it_fo), t_exp, t_weak, dye_after_fo

     end block
  enddo
  close(nunit_out)
  close(nunit_out2)
  close(nunit_out3)
  close(nunit_out4)
  
end program weak_freezeout
