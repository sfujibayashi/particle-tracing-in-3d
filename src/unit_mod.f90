module unit
!  K=c=G=1
  implicit none
  real(8) rho_uni, rhoi_uni,      t_uni,     ti_uni  &
        ,dedt_uni, dtde_uni, dmevdt_uni, dtdmev_uni  &
        ,dtdn_uni, dndt_uni,      r_uni,      v_uni  &
        , tem_uni,   chmerg,     chemev,      e_uni  &
        ,  dn_uni, dmev_uni
  parameter( rho_uni = 5.80783413d18     &
           ,rhoi_uni = 1.d0/rho_uni      &
           ,   t_uni = 1.60668219d-6     &
           ,  ti_uni = 1.d0/t_uni        &
           ,dedt_uni = 2.02775408d51     &
!           ,dedt_uni = 5.2052d39         &
           ,dtde_uni = 1.d0/dedt_uni     &
         ,dmevdt_uni = 5.2052d39         &
         ,dtdmev_uni = 1.d0/dmevdt_uni   &
!           ,dtdn_uni = 1.795486031d8     &
           ,dtdn_uni = dtde_uni          &
           ,dndt_uni = 1.d0/dtdn_uni     &
           ,   r_uni = 4.81671204d4      &  
           ,   v_uni = 2.99792458d10     &
           , tem_uni = 1.160445d10       &
           ,  chmerg = 1.60218d-6        &
           ,  chemev = 1.d0/1.60218d-6   &
           ,   e_uni = 5.829814713d53    &
           ,  dn_uni =   dndt_uni*t_uni  &
           ,dmev_uni =   dedt_uni*t_uni  )
end module unit
