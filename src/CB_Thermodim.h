!========================================================================
!     Common block Thermodim: dimensional thermodynamic coefficient
!     (defined on the C-grid) 
!========================================================================


      double precision Klat_ia, Klat_oa
      double precision Ksens_ai, Ksens_ao, Ksens_io
      double precision Kice, Kocn, Kadvo, Hocn
      double precision Kemis_i, Kemis_al, Kemis_o

      common/Thermodim/         &
                Klat_ia,        & ! LH transfer coefficient (ice/atm)
                Klat_oa,        & ! LH transfer coefficient (ocn/atm)
                Ksens_ai,       & ! SH transfer coefficient (ice/atm)
                Ksens_ao,       & ! SH transfer coefficient (ocn/atm)
                Ksens_io,       & ! SH transfer coefficient (ocn/ice)
                Kice,           & ! ice thermal conductivity
                Kocn,           & ! ocean diffusion coefficient [m2/s]
                Kadvo,          & ! adv heat transfer coefficient (ocn)
                Hocn,           & ! mixed layer depth
                Kemis_i,  &
                Kemis_al, &
                Kemis_o





