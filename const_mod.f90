module const
  implicit real(8) (a-h,o-z) ! none
  real(8) f13, f23, f43, f53, f16, f56, f76, f17, f37, f18  & 
       ,f19, f111, f112, f114, f156, f715, f283, f203       &
       ,f316 &
       ,pi, pi2, pi3, pi4, pi8, pi16, pi32, pi332, pai      & 
       ,ppi, ppi2, ppi3, ppi4, ppi5, pi13, dln2 &
       ,pii, pihi
  parameter(f13  = 1.d0/3.d0             &
           ,f23  = 2.d0/3.d0             &
           ,f43  = 4.d0/3.d0             &
           ,f53  = 5.d0/3.d0             &
           ,f16  = 1.d0/6.d0             &
           ,f56  = 5.d0/6.d0             &
           ,f76  = 7.d0/6.d0             &
           ,f17  = 1.d0/7.d0             &
           ,f37  = 3.d0/7.d0             &
           ,f18  = 1.d0/8.d0             & 
           ,f19  = 1.d0/9.d0             &
           ,f111 = 1.d0/11.d0            & 
           ,f112 = 1.d0/12.d0            &
           ,f114 = 1.d0/14.d0            &
           ,f156 = 1.d0/56.d0            &
           ,f715 = 7.d0/15.d0            &
           ,f283 = 28.d0/3.d0            &
           ,f203 = 20.d0/3.d0            &
           ,f316 = 3.d0/16.d0            &
           ,pi   = 3.141592653589793d0   &
           ,pi2  = pi*2.d0               & 
           ,pi3  = pi*3.d0               & 
           ,pi4  = pi*4.d0               &
           ,pi8  = pi*8.d0               &
           ,pi16 = pi*16.d0              &
           ,pi32 = pi*32.d0              &
           ,pi332= pi*3.d0/32.d0         &
           ,pai  = pi                    & 
           ,ppi  = pi*pi                 &
           ,ppi2 = ppi *ppi              &
           ,ppi3 = ppi2*ppi              &
           ,ppi4 = ppi2*ppi2             &
           ,ppi5 = ppi2*ppi3             &
           ,pi13 = 1.d0/3.d0 *pi**2      &
           ,pii  = 1.d0/pi               &
           ,pihi = 1.d0/sqrt(pi)         &
           ,dln2 = 0.693147181d0         )
end module const
