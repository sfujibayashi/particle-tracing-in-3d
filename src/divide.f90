module divide

  type meshblock
     integer :: jd,ju,kd,ku,ld,lu
  end type meshblock

contains
  subroutine divide8(jd,ju,kd,ku,ld,lu,jd_divided,ju_divided,kd_divided,ku_divided,ld_divided,lu_divided)
    integer,intent(in) :: jd,ju,kd,ku,ld,lu
    integer,intent(out) :: jd_divided(8),ju_divided(8),kd_divided(8),ku_divided(8),ld_divided(8),lu_divided(8)
    integer :: j_half,k_half,l_half

    j_half = jd + (ju-jd)/2 +1
    k_half = kd + (ku-kd)/2 +1
    l_half = ld + (lu-ld)/2 +1

    jd_divided(1) = jd
    jd_divided(2) = j_half
    jd_divided(3) = jd
    jd_divided(4) = j_half
    jd_divided(5) = jd
    jd_divided(6) = j_half
    jd_divided(7) = jd
    jd_divided(8) = j_half

    kd_divided(1) = kd
    kd_divided(2) = kd
    kd_divided(3) = k_half
    kd_divided(4) = k_half
    kd_divided(5) = kd
    kd_divided(6) = kd
    kd_divided(7) = k_half
    kd_divided(8) = k_half

    ld_divided(1) = ld
    ld_divided(2) = ld
    ld_divided(3) = ld
    ld_divided(4) = ld
    ld_divided(5) = l_half
    ld_divided(6) = l_half
    ld_divided(7) = l_half
    ld_divided(8) = l_half

    ju_divided(1) = j_half - 1
    ju_divided(2) = ju
    ju_divided(3) = j_half - 1
    ju_divided(4) = ju
    ju_divided(5) = j_half - 1
    ju_divided(6) = ju
    ju_divided(7) = j_half - 1
    ju_divided(8) = ju

    ku_divided(1) = k_half - 1
    ku_divided(2) = k_half - 1
    ku_divided(3) = ku
    ku_divided(4) = ku
    ku_divided(5) = k_half - 1
    ku_divided(6) = k_half - 1
    ku_divided(7) = ku
    ku_divided(8) = ku

    lu_divided(1) = l_half - 1
    lu_divided(2) = l_half - 1
    lu_divided(3) = l_half - 1
    lu_divided(4) = l_half - 1
    lu_divided(5) = lu
    lu_divided(6) = lu
    lu_divided(7) = lu
    lu_divided(8) = lu

  end subroutine divide8

  subroutine divide4(jd,ju,kd,ku,ld,lu,jd_divided,ju_divided,kd_divided,ku_divided,ld_divided,lu_divided)
    integer,intent(in) :: jd,ju,kd,ku,ld,lu
    integer,intent(out) :: jd_divided(4),ju_divided(4),kd_divided(4),ku_divided(4),ld_divided(4),lu_divided(4)
    integer :: j_half,k_half,l_half

    j_half = jd + (ju-jd)/2 +1
    k_half = kd + (ku-kd)/2 +1
    l_half = ld + (lu-ld)/2 +1

    jd_divided(1) = jd
    jd_divided(2) = j_half
    jd_divided(3) = jd
    jd_divided(4) = j_half

    kd_divided(1) = kd
    kd_divided(2) = kd
    kd_divided(3) = k_half
    kd_divided(4) = k_half

    ld_divided(1) = ld
    ld_divided(2) = ld
    ld_divided(3) = ld
    ld_divided(4) = ld

    ju_divided(1) = j_half - 1
    ju_divided(2) = ju
    ju_divided(3) = j_half - 1
    ju_divided(4) = ju

    ku_divided(1) = k_half - 1
    ku_divided(2) = k_half - 1
    ku_divided(3) = ku
    ku_divided(4) = ku

    lu_divided(1) = lu
    lu_divided(2) = lu
    lu_divided(3) = lu
    lu_divided(4) = lu

  end subroutine divide4

end module divide
