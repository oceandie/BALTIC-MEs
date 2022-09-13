MODULE zgrmes
   !!==============================================================================
   !!                       ***  MODULE zgrmes   ***
   !! Ocean initialization : Multiple Enveloped s coordinate (MES)
   !!==============================================================================
   !!  NEMO      4.0  ! 2021-07  (D. Bruciaferri)   
   !!----------------------------------------------------------------------
   !
   USE oce               ! ocean variables
   USE dom_oce           ! ocean domain
   USE depth_e3          ! depth <=> e3
   USE closea            ! closed seas
   USE mes               ! MEs-coordinates
   USE zgrloc            ! localised vert. coor. system 
   !
   USE in_out_manager    ! I/O manager
   USE iom               ! I/O library
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library
   USE wrk_nemo          ! Memory allocation
   USE timing            ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zgr_mes        ! called by domzgr.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"

CONTAINS

! =====================================================================================================

   SUBROUTINE zgr_mes
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE zgr_mes  ***
      !!
      !! ** Purpose : Wrap the steps to generate gloabal or localised 
      !!              MEs-coordinates systems
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk     ! dummy loop index
      !-----------------------------------------------------------------------
      !
      ! Initialise the mask to use MEs-coord. globally
      l2g_msk(:,:) = 2.0
     
      ! Generating a global MEs vertical grid
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL mes_build
      
      ! Local MEs
      ! ~~~~~~~~~
      IF( ln_loc_zgr )  CALL zgr_loc
      
      ! MEs with partial steps
      ! ~~~~~~~~~~~~~~~~~~~~~~
      IF(ln_pst_mes .OR. ln_pst_l2g)  CALL zgr_pst_mes
      !
      ! Computing gdep3w
      gde3w_0(:,:,1) = 0.5 * e3w_0(:,:,1)
      DO jk = 2, jpk
         gde3w_0(:,:,jk) = gde3w_0(:,:,jk-1) + e3w_0(:,:,jk)
      END DO
      !
      ! From here equal to sco code - domzgr.F90 line 2183
      !
      CALL lbc_lnk( e3t_0 , 'T', 1._wp )
      CALL lbc_lnk( e3u_0 , 'U', 1._wp )
      CALL lbc_lnk( e3v_0 , 'V', 1._wp )
      CALL lbc_lnk( e3f_0 , 'F', 1._wp )
      CALL lbc_lnk( e3w_0 , 'W', 1._wp )
      CALL lbc_lnk( e3uw_0, 'U', 1._wp )
      CALL lbc_lnk( e3vw_0, 'V', 1._wp )
      !
      WHERE (e3t_0   (:,:,:).eq.0.0)  e3t_0(:,:,:) = 1.0
      WHERE (e3u_0   (:,:,:).eq.0.0)  e3u_0(:,:,:) = 1.0
      WHERE (e3v_0   (:,:,:).eq.0.0)  e3v_0(:,:,:) = 1.0
      WHERE (e3f_0   (:,:,:).eq.0.0)  e3f_0(:,:,:) = 1.0
      WHERE (e3w_0   (:,:,:).eq.0.0)  e3w_0(:,:,:) = 1.0
      WHERE (e3uw_0  (:,:,:).eq.0.0)  e3uw_0(:,:,:) = 1.0
      WHERE (e3vw_0  (:,:,:).eq.0.0)  e3vw_0(:,:,:) = 1.0
      !
      IF ( ln_loc_zgr ) THEN
         IF ( .NOT. ln_pst_l2g ) THEN
            ! Only in the transition zone since MEs
            ! zone has been already taken care of
            IF ( lwp ) WRITE(numout,*) 'Refine mbathy'
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF ( l2g_msk(ji,jj) == 1._wp ) THEN
                     DO jk = 1, jpkm1
                        IF( scobot(ji,jj) >= gdept_0(ji,jj,jk) )   mbathy(ji,jj) = MAX( 2, jk )
                     END DO
                     IF( scobot(ji,jj) == 0._wp               )   mbathy(ji,jj) = 0
                  END IF
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE zgr_mes

! =====================================================================================================

   SUBROUTINE zgr_pst_mes
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE zgr_pst_mes  ***
     !!
     !! ** Purpose : Apply partial steps to MEs coordinates 
     !!              (and transition zone if using local MEs)
     !!
     !!----------------------------------------------------------------------
     INTEGER                             ::   ji, jj, jk, jp, jl   ! dummy loop index
     INTEGER                             ::   jbku, jbkv, jbkf     ! dummy loop index
     REAL(wp)                            ::   zdepth, e3psm, e3psr ! local variables
     REAL(wp)                            ::   zgdept, zgdepw, zgdepwp1
     REAL(wp)                            ::   ze3t, ze3w, e3u_j, e3u_jp1
     REAL(wp)                            ::   e3t_i, e3t_ip1, e3w_i, e3w_ip1
     REAL(wp)                            ::   e3t_j, e3t_jp1, e3w_j, e3w_jp1
     REAL(wp)                            ::   ps_i, ps_ip1, ps_j, ps_jp1
     REAL(wp)                            ::   ps_uj, ps_ujp1, ps_p, ps_s
     !
     REAL(wp), POINTER, DIMENSION(:,:)   :: pst_msk              !: mask for areas where to
                                                                 !: apply partial steps
     REAL(wp), POINTER, DIMENSION(:,:)   :: apst                 !: tracking points where
                                                                 !: we apply partial steps
     ! ARRAYS TO BACKUP ORIGINAL MEs COORDINATE SYSTEM 
     REAL(wp), POINTER, DIMENSION(:,:,:) :: e3t_mes, e3w_mes     !: scale factors [m]
     REAL(wp), POINTER, DIMENSION(:,:,:) :: gdept_mes, gdepw_mes !: depths [m]
     !!----------------------------------------------------------------------

     e3psm = rn_e3pst_min
     e3psr = rn_e3pst_rat

     CALL wrk_alloc( jpi, jpj, pst_msk, apst)
     CALL wrk_alloc( jpi, jpj, jpk, gdept_mes, gdepw_mes, e3t_mes, e3w_mes)

     ! Initialise the mask to identify where to apply partial steps
     pst_msk(:,:) = 0.0
     IF(ln_pst_mes)  WHERE (l2g_msk(:,:) == 2.0)  pst_msk(:,:) = 1.0
     IF(ln_pst_l2g)  WHERE (l2g_msk(:,:) == 1.0)  pst_msk(:,:) = 1.0


     ! Compute mbathy for ocean points (i.e. the number of ocean levels)
     ! find the number of ocean levels such that the last level thickness
     ! is larger than the minimum of e3zps_min and e3zps_rat * e3t_0
     DO jj = 1, jpj
        DO ji = 1, jpi
           IF ( pst_msk(ji,jj) == 1.0 .AND. bathy(ji,jj) > 0._wp) THEN
              DO jk = jpkm1, 1, -1
                 zdepth = gdepw_0(ji,jj,jk) + MIN( e3psm, e3t_0(ji,jj,jk)*e3psr )
                 IF ( bathy(ji,jj) <= zdepth )   mbathy(ji,jj) = jk-1
              END DO
            END IF
        END DO
     END DO

     ! ===================================================== 
     ! Scale factors and depth at T- and W- points
     ! =====================================================
     !
     ! Backup to the reference MEs-coordinate
     e3t_mes  (:,:,:) = e3t_0  (:,:,:)
     e3w_mes  (:,:,:) = e3w_0  (:,:,:)
     gdept_mes(:,:,:) = gdept_0(:,:,:)
     gdepw_mes(:,:,:) = gdepw_0(:,:,:)
     
     ! Mask identifying points where we apply partial steps:
     !  bottom cell is     a partial cell = 0 
     !  bottom cell is not a partial cell = 1
     apst(:,:) = 1._wp
     
     DO jj = 1, jpj
        DO ji = 1, jpi
           jk = mbathy(ji,jj)
           !
           ! We use partial steps only at ocean points where 
           ! envelopes differ from the bottom topography 
           IF ( pst_msk(ji,jj) == 1.0 .AND. jk > 0 ) THEN    
              !
              IF ( gdepw_0(ji,jj,jk+1) /= bathy(ji,jj) ) apst(ji,jj) = 0     
              !
              IF ( jk == jpkm1 ) THEN          ! max ocean level case
                 zgdepw = bathy(ji,jj)
                 ze3t   = bathy(ji,jj) - gdepw_mes(ji,jj,jk)
                 ze3w   = 0.5_wp * e3w_mes(ji,jj,jk) * ( 1._wp + ( ze3t/e3t_mes(ji,jj,jk) ) )
                 e3t_0(ji,jj,jk  ) = ze3t
                 e3w_0(ji,jj,jk  ) = ze3w
                 gdepw_0(ji,jj,jk+1) = zgdepw
                 gdept_0(ji,jj,jk  ) = gdept_mes(ji,jj,jk-1) + ze3w
                 gdept_0(ji,jj,jk+1) = gdept_0(ji,jj,jk) + ze3t
                 !
              ELSE                             ! standard case
                 zgdepw   = gdepw_mes(ji,jj,jk  )
                 zgdepwp1 = gdepw_mes(ji,jj,jk+1)
                 zgdept   = gdept_mes(ji,jj,jk  )
                 ze3t     = e3t_mes  (ji,jj,jk  )
                 ze3w     = e3w_mes  (ji,jj,jk  )
                 IF ( bathy(ji,jj) <= zgdepwp1 ) gdepw_0(ji,jj,jk+1) = bathy(ji,jj)
                 !gm Bug?  check the gdepw_1d
                 !       ... on jk
                 gdept_0(ji,jj,jk) = zgdepw + ( gdepw_0(ji,jj,jk+1) - zgdepw      ) &
                    &                       * ((zgdept-zgdepw) / (zgdepwp1-zgdepw))
                 e3t_0  (ji,jj,jk) = ze3t   * (gdepw_0(ji,jj,jk+1) - zgdepw) / (zgdepwp1 - zgdepw)
                 e3w_0  (ji,jj,jk) = 0.5_wp * (gdepw_0(ji,jj,jk+1) + zgdepwp1 - 2._wp * zgdepw) &
                    &                       * (ze3w / (zgdepwp1 - zgdepw))
                 !       ... on jk+1
                 gdept_0(ji,jj,jk+1) = gdept_0(ji,jj,jk) + e3t_0(ji,jj,jk)
              ENDIF
           ENDIF
        ENDDO
     ENDDO     
     !
     ! set value at mbathy+1
     DO jj = 1, jpj
        DO ji = 1, jpi
           jk = mbathy(ji,jj)
           IF ( pst_msk(ji,jj) == 1.0 .AND. jk > 0) THEN
               e3t_0  (ji,jj,jk+1) = e3t_0  (ji,jj,jk)
               e3w_0  (ji,jj,jk+1) = e3t_0  (ji,jj,jk)
           END IF
        END DO
     END DO
     !
     ! =========================================== 
     ! Scale factors at U-, V-, UW- and VW- points
     ! ===========================================
     DO jj = 1, jpjm1
        DO ji = 1, jpim1
           !
           DO jp = 0,1
              !  
              ! e3u at the bottom
              ! ~~~~~~~~~~~~~~~~~
              jbku = mbathy(ji+jp,jj)
              ! Only ocean points  
              IF (pst_msk(ji,jj) == 1.0 .AND. jbku > 0) THEN
                 ! work at bottom topography and level below
                 DO jk = jbku, jbku+1
                    ps_i   = apst(ji  ,jj  )   ! is bottom cell at T(ji  ,jj) point a partial cell or not
                    ps_ip1 = apst(ji+1,jj  )   ! is bottom cell at T(ji+1,jj) point a partial cell or not
                    ! correction for level jk
                    IF ( (jk /= mbathy(ji  ,jj)) .AND. (jk /= mbathy(ji  ,jj)+1) ) ps_i   = 1
                    IF ( (jk /= mbathy(ji+1,jj)) .AND. (jk /= mbathy(ji+1,jj)+1) ) ps_ip1 = 1
                    !
                    ps_p = ps_i * ps_ip1
                    ps_s = ps_i + ps_ip1
                    !
                    e3t_i   = e3t_0(ji  ,jj,jk)
                    e3t_ip1 = e3t_0(ji+1,jj,jk)
                    e3w_i   = e3w_0(ji  ,jj,jk)
                    e3w_ip1 = e3w_0(ji+1,jj,jk)
                    !
                    e3u_0(ji,jj,jk)  = MIN( e3u_0(ji,jj,jk) ,                              &
                     &                        ps_p    * e3u_0(ji,jj,jk)                    & ! no zps cell
                     &                     + (1-ps_p) * ( ps_ip1 * e3t_i + ps_i * e3t_ip1  & ! one of neighb is zps
                     &                                   + (1-ps_s) * MIN(e3t_i,e3t_ip1) ) & ! both neighb are zps
                     &                    )
                     
                    e3uw_0(ji,jj,jk) = MIN(e3uw_0(ji,jj,jk),                               &
                     &                        ps_p    * e3uw_0(ji,jj,jk)                   & ! no zps 
                     &                     + (1-ps_p) * ( ps_ip1 * e3w_i + ps_i * e3w_ip1  & ! one of neighb is zps
                     &                                   + (1-ps_s) * MIN(e3w_i,e3w_ip1) ) & ! both side are zps
                     &                    )
                 END DO
              END IF
              ! 
              ! e3v at the bottom
              ! ~~~~~~~~~~~~~~~~~
              jbkv = mbathy(ji,jj+jp)
              ! Only ocean points  
              IF (pst_msk(ji,jj) == 1.0 .AND. jbkv > 0) THEN
                 ! work at bottom topography and level below
                 DO jk = jbkv, jbkv+1
                    ps_j   = apst(ji, jj  )   ! is bottom cell at T(ji,jj  ) point a partial cell or not
                    ps_jp1 = apst(ji, jj+1)   ! is bottom cell at T(ji,jj+1) point a partial cell or not
                    ! correction for level jk
                    IF ( (jk /= mbathy(ji  ,jj)) .AND. (jk /= mbathy(ji  ,jj)+1) ) ps_j   = 1
                    IF ( (jk /= mbathy(ji+1,jj)) .AND. (jk /= mbathy(ji+1,jj)+1) ) ps_jp1 = 1
                    !
                    ps_p = ps_j * ps_jp1
                    ps_s = ps_j + ps_jp1
                    !
                    e3t_j   = e3t_0(ji,jj  ,jk)
                    e3t_jp1 = e3t_0(ji,jj+1,jk)
                    e3w_j   = e3w_0(ji,jj  ,jk)
                    e3w_jp1 = e3w_0(ji,jj+1,jk)
                    !
                    e3v_0(ji,jj,jk)  = MIN( e3v_0(ji,jj,jk) ,                              &
                     &                        ps_p    * e3v_0(ji,jj,jk)                    & ! no zps cell
                     &                     + (1-ps_p) * ( ps_jp1 * e3t_j + ps_j * e3t_jp1  & ! one of neighb is zps
                     &                                   + (1-ps_s) * MIN(e3t_j,e3t_jp1) ) & ! both neighb are zps
                     &                    )

                    e3vw_0(ji,jj,jk) = MIN(e3vw_0(ji,jj,jk),                               &
                     &                        ps_p    * e3vw_0(ji,jj,jk)                   & ! no zps 
                     &                     + (1-ps_p) * ( ps_jp1 * e3w_j + ps_j * e3w_jp1  & ! one of neighb is zps
                     &                                   + (1-ps_s) * MIN(e3w_j,e3w_jp1) ) & ! both side are zps
                     &                    )
                 END DO
              END IF
           END DO
        END DO
     END DO

     ! ========================= 
     ! Scale factors at F-points
     ! =========================   
     DO jj = 1, jpjm1
        DO ji = 1, jpim1
           DO jp = 0,1
              DO jl = 0,1
                 jbkf = mbathy(ji+jl,jj+jp)
                 ! Only ocean points  
                 IF (pst_msk(ji,jj) == 1.0 .AND. jbkf > 0) THEN
                    ! work at bottom topography and level below
                    DO jk = jbkf, jbkf+1
                       ps_uj   = MIN(apst(ji,jj  ),apst(ji+1,jj  )) ! is bot cell at U(ji,jj  ) pnt a part. cell or not
                       ps_ujp1 = MIN(apst(ji,jj+1),apst(ji+1,jj+1)) ! is bot cell at U(ji,jj+1) pnt a part. cell or not
                       ! correction for level jk
                       IF (      (jk /= mbathy(ji  ,jj  )) .AND. (jk /= mbathy(ji  ,jj  )+1) &
                         & .AND. (jk /= mbathy(ji+1,jj  )) .AND. (jk /= mbathy(ji+1,jj  )+1) ) ps_uj   = 1
                       IF (      (jk /= mbathy(ji  ,jj+1)) .AND. (jk /= mbathy(ji  ,jj+1)+1) &
                         & .AND. (jk /= mbathy(ji+1,jj+1)) .AND. (jk /= mbathy(ji+1,jj+1)+1) ) ps_ujp1 = 1
                       !
                       ps_p = ps_uj * ps_ujp1
                       ps_s = ps_uj + ps_ujp1
                       !
                       e3u_j   = e3u_0(ji,jj  ,jk)
                       e3u_jp1 = e3u_0(ji,jj+1,jk)
                       !
                       e3f_0(ji,jj,jk) = MIN( e3f_0(ji,jj,jk) ,                                &
                        &                       ps_p    * e3f_0(ji,jj,jk)                      & ! no zps cell
                        &                    + (1-ps_p) * ( ps_ujp1 * e3u_j + ps_uj * e3u_jp1  & ! one of neighb is zps
                        &                                  + (1-ps_s) * MIN(e3u_j,e3u_jp1)   ) & ! both neighb are zps
                        &                   )
                    END DO
                 END IF
              END DO
           END DO
        END DO
     END DO
     !
     ! set to z-scale factor if zero (i.e. along closed boundaries)
     DO jk = 1, jpk
        WHERE( e3f_0(:,:,jk) == 0._wp ) e3f_0(:,:,jk) = e3t_0(ji,jj,jk)
     END DO
     !
     e3t_0(:,mj0(1),:) = e3t_0(:,mj0(2),:)     ! we duplicate factor scales for jj = 1 and jj = 2
     e3w_0(:,mj0(1),:) = e3w_0(:,mj0(2),:)
     e3u_0(:,mj0(1),:) = e3u_0(:,mj0(2),:)
     e3v_0(:,mj0(1),:) = e3v_0(:,mj0(2),:)
     e3f_0(:,mj0(1),:) = e3f_0(:,mj0(2),:)

     ! Control of the sign
     IF( MINVAL( e3t_0  (:,:,:) ) <= 0._wp )   CALL ctl_stop( '  zgr_pst_mes : e r r o r   e3t_0 <= 0' )
     IF( MINVAL( e3w_0  (:,:,:) ) <= 0._wp )   CALL ctl_stop( '  zgr_pst_mes : e r r o r   e3w_0 <= 0' )
     IF( MINVAL( gdept_0(:,:,:) ) <  0._wp )   CALL ctl_stop( '  zgr_pst_mes : e r r o r   gdept_0 < 0' )
     IF( MINVAL( gdepw_0(:,:,:) ) <  0._wp )   CALL ctl_stop( '  zgr_pst_mes : e r r o r   gdepw_0 < 0' )

     CALL wrk_dealloc( jpi, jpj, pst_msk, apst)
     CALL wrk_dealloc( jpi, jpj, jpk, gdept_mes, gdepw_mes, e3t_mes, e3w_mes)

   END SUBROUTINE zgr_pst_mes

END MODULE zgrmes
