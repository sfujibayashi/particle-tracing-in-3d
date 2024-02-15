subroutine print_data(t,job,it)
#include "macro.h"
  use io
  use simdata3D
  !$use omp_lib
  real(8),intent(in) :: t
  integer,intent(in) :: job,it

  integer :: j,k,l,lv
  character(256) :: str1,str2,str3

  real(8),parameter :: clight = 2.99792458d10, grav=6.674d-8, msun=1.989d33
  real(8),parameter :: r_unit_kkawa = grav*msun/clight**2,t_unit_kkawa = r_unit_kkawa/clight, rho_unit_kkawa=msun/r_unit_kkawa**3, e_unit_kkawa=rho_unit_kkawa*clight**2, b_unit_kkawa = sqrt(e_unit_kkawa)
  
  integer :: ld_print
  integer :: unit_num
 
ld_print=ld

! #ifdef STAGGERED
! #ifndef FULL
!   ld_print = ld + 1
! #endif
! #endif
 
  
  write(str2,'("job",i3.3)') job
  write(str3,'("it",i3.3)') it

  !$omp parallel &
  !$omp default(none) &
  !$omp shared(lv_min,lv_max,dir_out,jd,ju,kd,ku,lu,ld_print,str2,str3,t,x,y,z,qb,vlx,vly,vlz,pres,qrho) &
  !$omp private(str1,unit_num)
  !$omp do
  do lv=lv_min,lv_max
     write(str1,'(i10)') lv
     open(newunit=unit_num,file=trim(dir_out)//"/input3D_"//trim(str2)//"_"//trim(str3)//"_"//trim(adjustl(str1)),status="replace",action="write")
     
     !write(unit_num,'("#",99i5)') (ju-jd+1), (ku-kd+1), (lu-ld_print+1)
     !write(unit_num,'("# Normalization: c=1, G=",es13.5,"=1, Msun=",es13.5,"=1")') grav,msun
     !write(unit_num,'("# Extraction radius: ",es13.5,"cm, model: ",a)') r_ext*r_uni,trim(cmodel)
     !write(unit_num,'("# 1:time, 2:angle, 3:rho, 4:enthalpy, 5:pressure, 6:v^x, 7:v^y, 8:v^z, 9:Lorentz factor, 10:Lapse, 11:Conformal factor (Psi), 12:Ye (boundary), 13:S (boundary)")')
     write(unit_num,'(a)') "1: x, 2: y, 3: z, 4: rho_star, 5: v^x, 6: v^x, 7: v^x, 8: P, 9: rho_star/rho, 10: Bx, 11: By, 12:Bz"
     write(unit_num,*) t/t_unit_kkawa
     do l=ld_print,lu
        do k=kd,ku
           do j=jd,ju
              write(unit_num,'(99es15.7)') x(j,lv)/r_unit_kkawa,y(k,lv)/r_unit_kkawa,z(l,lv)/r_unit_kkawa, &
                   qb(j,k,l,lv)/rho_unit_kkawa,vlx(j,k,l,lv),vly(j,k,l,lv),vlz(j,k,l,lv),pres(j,k,l,lv)/e_unit_kkawa,qb(j,k,l,lv)/qrho(j,k,l,lv), &
                   0d0, 0d0, 0d0
              !bx(j,k,l,lv)/b_unit_kkawa,by(j,k,l,lv)/b_unit_kkawa,bz(j,k,l,lv)/b_unit_kkawa
           enddo
        enddo
     enddo
     close(unit_num)
  enddo
  !$omp end do
  !$omp end parallel
  
end subroutine print_data
