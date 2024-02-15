subroutine coorindex3D(xx,yy,zz,j1,k1,l1,lv)
  use simdata3D
  implicit none
  real(8),intent(in) :: xx,yy,zz
  integer,intent(out) :: lv,j1,k1,l1

  lv = lv_max
  do
     if((abs(xx)<x(ju,lv).and.abs(yy)<y(ku,lv).and.abs(zz)<z(lu,lv)).or.lv==lv_min)then
        exit
     else
        lv = lv - 1
     endif
  enddo

  j1 = jd + 1
  do while( ( xx - x(j1-1,lv) )*( xx - x(j1,lv) ) > 0.d0 .and.j1<ju)
   j1 = j1 + 1
  end do
  
  k1 = kd + 1
  do while( ( yy - y(k1-1,lv) )*( yy - y(k1,lv) ) > 0.d0 .and.k1<ku)
   k1 = k1 + 1
  end do
  
  l1 = ld + 1
  do while( ( zz - z(l1-1,lv) )*( zz - z(l1,lv) ) > 0.d0 .and.l1<lu)
     l1 = l1 + 1
  end do

  if(lv<lv_min .or. lv>lv_max .or. &
       j1 < jd + 1 .or. ju < j1 .or. &
       k1 < kd + 1 .or. ku < k1 .or. &
       l1 < ld + 1 .or. lu < l1 )then
     write(6,'(a,99i10)') "index out of range", j1,k1,l1,lv
  endif
  
  
  return
end subroutine coorindex3D
