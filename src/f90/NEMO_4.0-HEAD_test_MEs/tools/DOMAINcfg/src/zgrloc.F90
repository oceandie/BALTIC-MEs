MODULE zgrloc
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
   !
   USE in_out_manager    ! I/O manager
   USE iom               ! I/O library
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library
   USE wrk_nemo          ! Memory allocation
   USE timing            ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zgr_loc        ! called by zgrmes.F90
   !
   ! INPUT GLOBAL VERTICAL COORDINATE SYSTEM
   !
   REAL(wp), POINTER, DIMENSION(:)     :: e3t_1d_in  , e3w_1d_in   !: ref. scale factors (m)
   REAL(wp), POINTER, DIMENSION(:)     :: gdept_1d_in, gdepw_1d_in !: ref. depth (m)
   !
   REAL(wp), POINTER, DIMENSION(:,:)   :: mbathy_in
   !
   REAL(wp), POINTER, DIMENSION(:,:,:) :: e3t_in, e3u_in , e3v_in , e3f_in !: scale factors [m]
   REAL(wp), POINTER, DIMENSION(:,:,:) :: e3w_in, e3uw_in, e3vw_in         !:   -      -
   ! 
   REAL(wp), POINTER, DIMENSION(:,:,:) :: gdept_in, gdepw_in               !: depths [m]
   !
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"

CONTAINS

! =====================================================================================================

   SUBROUTINE zgr_loc
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE zgr_loc  ***
      !!
      !! ** Purpose : Wrap the steps to generate a localised vertical 
      !!              coordinate systems
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk     ! dummy loop index
      !-----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpk, gdept_1d_in, gdepw_1d_in, e3t_1d_in  , e3w_1d_in)
      CALL wrk_alloc( jpi, jpj, mbathy_in)
      CALL wrk_alloc( jpi, jpj, jpk, gdept_in, gdepw_in, e3t_in, e3w_in)
      CALL wrk_alloc( jpi, jpj, jpk, e3u_in  , e3v_in  , e3f_in, e3uw_in, e3vw_in)
      !
      ! Reading the input vertical grid that will be used globally
      CALL zgr_dom_read
      ! Creating the local cooridnate system within the global input vertical grid
      IF ( ln_mes ) CALL zgr_mes_local
      !
      CALL wrk_dealloc( jpk, gdept_1d_in, gdepw_1d_in, e3t_1d_in  , e3w_1d_in)
      CALL wrk_dealloc( jpi, jpj, mbathy_in)
      CALL wrk_dealloc( jpi, jpj, jpk, gdept_in, gdepw_in, e3t_in, e3w_in)
      CALL wrk_dealloc( jpi, jpj, jpk, e3u_in  , e3v_in  , e3f_in, e3uw_in, e3vw_in)

   END SUBROUTINE zgr_loc

! =====================================================================================================

   SUBROUTINE zgr_dom_read
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE zgr_dom_read  ***
     !!
     !! ** Purpose :   Read the vertical information in the domain configuration file
     !!
     !!----------------------------------------------------------------------
     INTEGER  ::   jk     ! dummy loop index
     INTEGER  ::   inum   ! local logical unit
     REAL(WP) ::   z_cav
     REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
     !!----------------------------------------------------------------------
     !
     IF(lwp) THEN
       WRITE(numout,*)
       WRITE(numout,*) '   zgr_dom_read : read the vertical coordinates in domain_cfg_global.nc file'
       WRITE(numout,*) '   ~~~~~~~~'
     ENDIF
     !
     CALL iom_open( 'domain_cfg_global.nc', inum )
     !
     !                          !* ocean cavities under iceshelves
     CALL iom_get( inum, 'ln_isfcav', z_cav )
     IF( z_cav == 0._wp ) THEN   ;   ln_isfcav = .false.   ;   ELSE   ;   ln_isfcav = .true.   ;   ENDIF
     !
     !                          !* vertical scale factors
     CALL iom_get( inum, jpdom_unknown, 'e3t_1d'  , e3t_1d_in ) ! 1D reference coordinate
     CALL iom_get( inum, jpdom_unknown, 'e3w_1d'  , e3w_1d_in )
     !
     CALL iom_get( inum, jpdom_global, 'bottom_level' , mbathy_in )     ! 2D mbathy
     !
     CALL iom_get( inum, jpdom_global, 'e3t_0'  , e3t_in  )     ! 3D coordinate
     CALL iom_get( inum, jpdom_global, 'e3u_0'  , e3u_in  )
     CALL iom_get( inum, jpdom_global, 'e3v_0'  , e3v_in  )
     CALL iom_get( inum, jpdom_global, 'e3f_0'  , e3f_in  )
     CALL iom_get( inum, jpdom_global, 'e3w_0'  , e3w_in  )
     CALL iom_get( inum, jpdom_global, 'e3uw_0' , e3uw_in )
     CALL iom_get( inum, jpdom_global, 'e3vw_0' , e3vw_in )
     !
     !                          !* depths
     !                                   !- old depth definition (obsolescent feature)
     IF(  iom_varid( inum, 'gdept_1d', ldstop = .FALSE. ) > 0  .AND.  &
        & iom_varid( inum, 'gdepw_1d', ldstop = .FALSE. ) > 0  .AND.  &
        & iom_varid( inum, 'gdept_0' , ldstop = .FALSE. ) > 0  .AND.  &
        & iom_varid( inum, 'gdepw_0' , ldstop = .FALSE. ) > 0    ) THEN
        CALL ctl_warn( 'zgr_dom_read : old definition of depths and scale factors used ', & 
           &           '           depths at t- and w-points read in the domain configuration file')
        CALL iom_get( inum, jpdom_unknown, 'gdept_1d', gdept_1d_in )   
        CALL iom_get( inum, jpdom_unknown, 'gdepw_1d', gdepw_1d_in )
        CALL iom_get( inum, jpdom_global , 'gdept_0' , gdept_in )
        CALL iom_get( inum, jpdom_global , 'gdepw_0' , gdepw_in )
        !
     ELSE                                !- depths computed from e3. scale factors
        CALL e3_to_depth( e3t_1d_in, e3w_1d_in, gdept_1d_in, gdepw_1d_in )    ! 1D reference depth
        CALL e3_to_depth( e3t_in   , e3w_in   , gdept_in   , gdepw_in    )    ! 3D depths
        IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) '              GLOBAL reference 1D z-coordinate depth and scale factors:'
          WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
          WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, gdept_1d_in(jk), gdepw_1d_in(jk), e3t_1d_in(jk), e3w_1d_in(jk), jk = 1, jpk )
        ENDIF
     ENDIF
     !
     CALL iom_close( inum )
     !
   END SUBROUTINE zgr_dom_read

! =====================================================================================================

   SUBROUTINE zgr_mes_local
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE zgr_dom_merge  ***
     !!
     !! ** Purpose :   Create a local MEs grid within a gloabal grid 
     !!                using different vertical coordinates.
     !!
     !!----------------------------------------------------------------------
     INTEGER                      ::   ji, jj, jk ! dummy loop index
     INTEGER                      ::   inum       ! local logical unit
     REAL(wp), DIMENSION(jpi,jpj) ::   l2g_wgt    ! weigths for computing model
                                                  ! levels in transition area
     !!----------------------------------------------------------------------

     IF(lwp) THEN
       WRITE(numout,*)
       WRITE(numout,*) '   zgr_mes_merge : read the vertical coordinates in domain_cfg_global.nc file'
       WRITE(numout,*) '   ~~~~~~~~'
     ENDIF
     !
     CALL iom_open( 'bathy_meter.nc', inum )
     CALL iom_get( inum, jpdom_data, 's2z_msk', l2g_msk)
     CALL iom_get( inum, jpdom_data, 's2z_wgt', l2g_wgt)

     DO jj = 1,jpj
        DO ji = 1,jpi
           SELECT CASE (INT(l2g_msk(ji,jj)))
             CASE (0) ! global zps area
                  mbathy (ji,jj  ) = mbathy_in(ji,jj)
                  gdept_0(ji,jj,:) = gdept_in(ji,jj,:)
                  gdepw_0(ji,jj,:) = gdepw_in(ji,jj,:)
                  e3t_0  (ji,jj,:) = e3t_in (ji,jj,:)
                  e3w_0  (ji,jj,:) = e3w_in (ji,jj,:)
                  e3u_0  (ji,jj,:) = e3u_in (ji,jj,:)
                  e3v_0  (ji,jj,:) = e3v_in (ji,jj,:)
                  e3f_0  (ji,jj,:) = e3f_in (ji,jj,:)
                  e3uw_0 (ji,jj,:) = e3uw_in (ji,jj,:)
                  e3vw_0 (ji,jj,:) = e3vw_in (ji,jj,:)
             CASE (1) ! MEs to zps transition area 
                  gdept_0(ji,jj,:) =           l2g_wgt(ji,jj)   * gdept_0(ji,jj,:) + &
                    &                ( 1._wp - l2g_wgt(ji,jj) ) * gdept_in(ji,jj,:)
                  gdepw_0(ji,jj,:) =           l2g_wgt(ji,jj)   * gdepw_0(ji,jj,:) + &
                    &                ( 1._wp - l2g_wgt(ji,jj) ) * gdepw_in(ji,jj,:)
             CASE (2) ! MEs area
                  CYCLE     
           END SELECT
        END DO
     END DO
     !
     ! e3t, e3w for transition zone
     ! as finite differences
     DO jj = 1, jpj
        DO ji = 1, jpi
           IF ( l2g_msk(ji,jj) == 1._wp ) THEN
              DO jk = 1,jpkm1
                 e3t_0(ji,jj,jk)   = gdepw_0(ji,jj,jk+1) - gdepw_0(ji,jj,jk)
                 e3w_0(ji,jj,jk+1) = gdept_0(ji,jj,jk+1) - gdept_0(ji,jj,jk)
              ENDDO
              ! Surface
              jk = 1
              e3w_0(ji,jj,jk) = 2.0_wp * (gdept_0(ji,jj,1) - gdepw_0(ji,jj,1))
              !
              ! Bottom
              jk = jpk
              e3t_0(ji,jj,jk) = 2.0_wp * (gdept_0(ji,jj,jk) - gdepw_0(ji,jj,jk))
           END IF
        END DO
     END DO
     !
     ! MBATHY transition zone
     DO jj = 1, jpj
        DO ji = 1, jpi
           IF ( l2g_msk(ji,jj) == 1._wp ) THEN
              DO jk = 1, jpkm1
                 IF( scobot(ji,jj) >= gdept_0(ji,jj,jk) ) mbathy(ji,jj) = MAX( 2, jk )
                 IF( scobot(ji,jj) == 0.e0              ) mbathy(ji,jj) = 0
                 IF( scobot(ji,jj) < 0.e0               ) mbathy(ji,jj) = INT( scobot(ji,jj)) ! do we need it?
              END DO
           END IF
        END DO
     END DO
     !
     ! Computing e3u_0, e3v_0, e3f_0, e3uw_0, e3vw_0 
     ! for transition zone
     !
     DO jj = 1, jpjm1
        DO ji = 1, jpim1
           IF ( l2g_msk(ji,jj) == 1._wp ) THEN
              DO jk = 1, jpk
                 e3u_0(ji,jj,jk)=(MIN(1,mbathy(ji,jj))* e3t_0(ji,jj,jk)      +         &
                                  MIN(1,mbathy(ji+1,jj))*e3t_0(ji+1,jj,jk) ) /         &
                                  MAX( 1, MIN(1,mbathy(ji,jj))+MIN(1,mbathy(ji+1,jj))  )

                 e3v_0(ji,jj,jk)=(MIN(1,mbathy(ji,jj))* e3t_0(ji,jj,jk)      +         &
                                  MIN(1,mbathy(ji,jj+1))*e3t_0(ji,jj+1,jk) ) /         &
                                  MAX( 1, MIN(1,mbathy(ji,jj))+MIN(1,mbathy(ji,jj+1))  )

                 e3uw_0(ji,jj,jk)=(MIN(1,mbathy(ji,jj))* e3w_0(ji,jj,jk)      +         &
                                   MIN(1,mbathy(ji+1,jj))*e3w_0(ji+1,jj,jk) ) /         &
                                   MAX( 1, MIN(1,mbathy(ji,jj))+MIN(1,mbathy(ji+1,jj))  )

                 e3vw_0(ji,jj,jk)=(MIN(1,mbathy(ji,jj))* e3w_0(ji,jj,jk)      +         &
                                   MIN(1,mbathy(ji,jj+1))*e3w_0(ji,jj+1,jk) ) /         &
                                   MAX( 1, MIN(1,mbathy(ji,jj))+MIN(1,mbathy(ji,jj+1))  )

                 e3f_0(ji,jj,jk)=(MIN(1,mbathy(ji,jj))* e3t_0(ji,jj,jk)         +       &
                                  MIN(1,mbathy(ji+1,jj))*e3t_0(ji+1,jj,jk)      +       &
                                  MIN(1,mbathy(ji+1,jj+1))* e3t_0(ji+1,jj+1,jk) +       &
                                  MIN(1,mbathy(ji,jj+1))*e3t_0(ji,jj+1,jk) )    /       &
                                  MAX(  1, MIN(1,mbathy(ji,jj))+MIN(1,mbathy(ji,jj+1))  &
                                +   MIN(1,mbathy(ji+1,jj))+MIN(1,mbathy(ji+1,jj+1))  )
              END DO
           END IF
        END DO
     END DO

   END SUBROUTINE zgr_mes_local

END MODULE zgrloc
