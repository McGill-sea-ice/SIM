!========================================================================
!     Common block for semi-lag: store alphamx and alphamy h,A thermo update
!========================================================================

      double precision                             &
                   alpmx         (1:nx,1:ny),      &
                   alpmy         (1:nx,1:ny)

      common/semilag/     &
                   alpmx                    ,      & ! alphamx
                   alpmy                             ! alphamy





