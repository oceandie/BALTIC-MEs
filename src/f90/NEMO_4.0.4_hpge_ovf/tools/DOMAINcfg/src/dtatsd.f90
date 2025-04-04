MODULE dtatsd
   !!======================================================================
   !!                     ***  MODULE  dtatsd  ***
   !! Ocean data  :  read ocean Temperature & Salinity Data from gridded data
   !!======================================================================
   !! History :  OPA  ! 1991-03  ()  Original code
   !!             -   ! 1992-07  (M. Imbard)
   !!            8.0  ! 1999-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module 
   !!            3.3  ! 2010-10  (C. Bricaud, S. Masson)  use of fldread
   !!            3.4  ! 2010-11  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dta_tsd      : read and time interpolated ocean Temperature & Salinity Data
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_tsd_init   ! called by opa.F90
   PUBLIC   dta_tsd        ! called by istate.F90 and tradmp.90

   LOGICAL , PUBLIC ::   ln_tsd_init      !: T & S data flag
   LOGICAL , PUBLIC ::   ln_tsd_tradmp    !: internal damping toward input data flag

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tsd   ! structure of input SST (file informations, fields read)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dtatsd.F90 6140 2015-12-21 11:35:23Z timgraham $ 
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_tsd_init( ld_tradmp )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tsd_init  ***
      !!                    
      !! ** Purpose :   initialisation of T & S input data 
      !! 
      !! ** Method  : - Read namtsd namelist
      !!              - allocates T & S data structure 
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in), OPTIONAL ::   ld_tradmp   ! force the initialization when tradp is used
      !
      INTEGER ::   ierr0, ierr1, ierr2, ierr3   ! temporary integers
      !
      CHARACTER(len=100)            ::   cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION( jpts) ::   slf_i           ! array of namelist informations on the fields to read
      TYPE(FLD_N)                   ::   sn_tem, sn_sal
      !!
      NAMELIST/namtsd/   ln_tsd_init, ln_tsd_tradmp, cn_dir, sn_tem, sn_sal
      INTEGER  ::   ios
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_tsd_init')
      !
      !  Initialisation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0
      !
      REWIND( numnam_ref )              ! Namelist namtsd in reference namelist : 
      READ  ( numnam_ref, namtsd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtsd in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namtsd in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namtsd, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtsd in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtsd )

      IF( PRESENT( ld_tradmp ) )   ln_tsd_tradmp = .TRUE.     ! forces the initialization when tradmp is used
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dta_tsd_init : Temperature & Salinity data '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtsd'
         WRITE(numout,*) '      Initialisation of ocean T & S with T &S input data   ln_tsd_init   = ', ln_tsd_init
         WRITE(numout,*) '      damping of ocean T & S toward T &S input data        ln_tsd_tradmp = ', ln_tsd_tradmp
         WRITE(numout,*)
         IF( .NOT.ln_tsd_init .AND. .NOT.ln_tsd_tradmp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   T & S data not used'
         ENDIF
      ENDIF
      !
      IF( ln_rstart .AND. ln_tsd_init ) THEN
         CALL ctl_warn( 'dta_tsd_init: ocean restart and T & S data intialisation, ',   &
            &           'we keep the restart T & S values and set ln_tsd_init to FALSE' )
         ln_tsd_init = .FALSE.
      ENDIF
      !
      !                             ! allocate the arrays (if necessary)
      IF(  ln_tsd_init .OR. ln_tsd_tradmp  ) THEN
         !
         ALLOCATE( sf_tsd(jpts), STAT=ierr0 )
         IF( ierr0 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd_init: unable to allocate sf_tsd structure' )   ;   RETURN
         ENDIF
         !
                                ALLOCATE( sf_tsd(jp_tem)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
         IF( sn_tem%ln_tint )   ALLOCATE( sf_tsd(jp_tem)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
                                ALLOCATE( sf_tsd(jp_sal)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
         IF( sn_sal%ln_tint )   ALLOCATE( sf_tsd(jp_sal)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
         !
         IF( ierr0 + ierr1 + ierr2 + ierr3 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd : unable to allocate T & S data arrays' )   ;   RETURN
         ENDIF
         !                         ! fill sf_tsd with sn_tem & sn_sal and control print
         slf_i(jp_tem) = sn_tem   ;   slf_i(jp_sal) = sn_sal
         CALL fld_fill( sf_tsd, slf_i, cn_dir, 'dta_tsd', 'Temperature & Salinity data', 'namtsd' )
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dta_tsd_init')
      !
   END SUBROUTINE dta_tsd_init


   SUBROUTINE dta_tsd( kt, ptsd )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tsd  ***
      !!                    
      !! ** Purpose :   provides T and S data at kt
      !! 
      !! ** Method  : - call fldread routine
      !!              - ORCA_R2: add some hand made alteration to read data  
      !!              - 'key_orca_lev10' interpolates on 10 times more levels
      !!              - s- or mixed z-s coordinate: vertical interpolation on model mesh
      !!              - ln_tsd_tradmp=F: deallocates the T-S data structure
      !!                as T-S data are no are used
      !!
      !! ** Action  :   ptsd   T-S data on medl mesh and interpolated at time-step kt
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt     ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   ptsd   ! T & S data
      !
      INTEGER ::   ji, jj, jk, jl, jkk   ! dummy loop indicies
      INTEGER ::   ik, il0, il1, ii0, ii1, ij0, ij1   ! local integers
      REAL(wp)::   zl, zi
      REAL(wp), POINTER, DIMENSION(:) ::  ztp, zsp   ! 1D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_tsd')
      !
      CALL fld_read( kt, 1, sf_tsd )      !==   read T & S data at kt time step   ==!
      !
      !
      !                                   !==   ORCA_R2 configuration and T & S damping   ==! 
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 .AND. ln_tsd_tradmp ) THEN    ! some hand made alterations
         !
         ij0 = 101   ;   ij1 = 109                       ! Reduced T & S in the Alboran Sea
         ii0 = 141   ;   ii1 = 155
         DO jj = mj0(ij0), mj1(ij1)
            DO ji = mi0(ii0), mi1(ii1)
               sf_tsd(jp_tem)%fnow(ji,jj,13:13) = sf_tsd(jp_tem)%fnow(ji,jj,13:13) - 0.20_wp
               sf_tsd(jp_tem)%fnow(ji,jj,14:15) = sf_tsd(jp_tem)%fnow(ji,jj,14:15) - 0.35_wp
               sf_tsd(jp_tem)%fnow(ji,jj,16:25) = sf_tsd(jp_tem)%fnow(ji,jj,16:25) - 0.40_wp
               !
               sf_tsd(jp_sal)%fnow(ji,jj,13:13) = sf_tsd(jp_sal)%fnow(ji,jj,13:13) - 0.15_wp
               sf_tsd(jp_sal)%fnow(ji,jj,14:15) = sf_tsd(jp_sal)%fnow(ji,jj,14:15) - 0.25_wp
               sf_tsd(jp_sal)%fnow(ji,jj,16:17) = sf_tsd(jp_sal)%fnow(ji,jj,16:17) - 0.30_wp
               sf_tsd(jp_sal)%fnow(ji,jj,18:25) = sf_tsd(jp_sal)%fnow(ji,jj,18:25) - 0.35_wp
            END DO
         END DO
         ij0 =  87   ;   ij1 =  96                          ! Reduced temperature in Red Sea
         ii0 = 148   ;   ii1 = 160
         sf_tsd(jp_tem)%fnow( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ,  4:10 ) = 7.0_wp
         sf_tsd(jp_tem)%fnow( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 11:13 ) = 6.5_wp
         sf_tsd(jp_tem)%fnow( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 14:20 ) = 6.0_wp
      ENDIF
      !
      ptsd(:,:,:,jp_tem) = sf_tsd(jp_tem)%fnow(:,:,:)    ! NO mask
      ptsd(:,:,:,jp_sal) = sf_tsd(jp_sal)%fnow(:,:,:) 
      !
      IF( ln_sco .OR. ln_mes ) THEN                   !==   s- or mixed s-zps-coordinate   ==!
         !
         CALL wrk_alloc( jpk, ztp, zsp )
         !
         IF( kt == nit000 .AND. lwp )THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dta_tsd: interpolates T & S data onto the s- or mixed s-z-coordinate mesh'
         ENDIF
         !
         DO jj = 1, jpj                         ! vertical interpolation of T & S
            DO ji = 1, jpi
               DO jk = 1, jpk                        ! determines the intepolated T-S profiles at each (i,j) points
                  zl = gdept_0(ji,jj,jk)
                  IF(     zl < gdept_1d(1  ) ) THEN          ! above the first level of data
                     ztp(jk) =  ptsd(ji,jj,1    ,jp_tem)
                     zsp(jk) =  ptsd(ji,jj,1    ,jp_sal)
                  ELSEIF( zl > gdept_1d(jpk) ) THEN          ! below the last level of data
                     ztp(jk) =  ptsd(ji,jj,jpkm1,jp_tem)
                     zsp(jk) =  ptsd(ji,jj,jpkm1,jp_sal)
                  ELSE                                      ! inbetween : vertical interpolation between jkk & jkk+1
                     DO jkk = 1, jpkm1                                  ! when  gdept(jkk) < zl < gdept(jkk+1)
                        IF( (zl-gdept_1d(jkk)) * (zl-gdept_1d(jkk+1)) <= 0._wp ) THEN
                           zi = ( zl - gdept_1d(jkk) ) / (gdept_1d(jkk+1)-gdept_1d(jkk))
                           ztp(jk) = ptsd(ji,jj,jkk,jp_tem) + ( ptsd(ji,jj,jkk+1,jp_tem) - ptsd(ji,jj,jkk,jp_tem) ) * zi 
                           zsp(jk) = ptsd(ji,jj,jkk,jp_sal) + ( ptsd(ji,jj,jkk+1,jp_sal) - ptsd(ji,jj,jkk,jp_sal) ) * zi
                        ENDIF
                     END DO
                  ENDIF
               END DO
               DO jk = 1, jpkm1
                  ptsd(ji,jj,jk,jp_tem) = ztp(jk) * tmask(ji,jj,jk)     ! mask required for mixed zps-s-coord
                  ptsd(ji,jj,jk,jp_sal) = zsp(jk) * tmask(ji,jj,jk)
               END DO
               ptsd(ji,jj,jpk,jp_tem) = 0._wp
               ptsd(ji,jj,jpk,jp_sal) = 0._wp
            END DO
         END DO
         ! 
         CALL wrk_dealloc( jpk, ztp, zsp )
         ! 
      ELSE                                !==   z- or zps- coordinate   ==!
         !                             
         ptsd(:,:,:,jp_tem) = ptsd(:,:,:,jp_tem) * tmask(:,:,:)    ! Mask
         ptsd(:,:,:,jp_sal) = ptsd(:,:,:,jp_sal) * tmask(:,:,:)
         !
         IF( ln_zps ) THEN                      ! zps-coordinate (partial steps) interpolation at the last ocean level
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ik = mbkt(ji,jj) 
                  IF( ik > 1 ) THEN
                     zl = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
                     ptsd(ji,jj,ik,jp_tem) = (1.-zl) * ptsd(ji,jj,ik,jp_tem) + zl * ptsd(ji,jj,ik-1,jp_tem)
                     ptsd(ji,jj,ik,jp_sal) = (1.-zl) * ptsd(ji,jj,ik,jp_sal) + zl * ptsd(ji,jj,ik-1,jp_sal)
                  ENDIF
                  ik = mikt(ji,jj)
                  IF( ik > 1 ) THEN
                     zl = ( gdept_0(ji,jj,ik) - gdept_1d(ik) ) / ( gdept_1d(ik+1) - gdept_1d(ik) ) 
                     ptsd(ji,jj,ik,jp_tem) = (1.-zl) * ptsd(ji,jj,ik,jp_tem) + zl * ptsd(ji,jj,ik+1,jp_tem)
                     ptsd(ji,jj,ik,jp_sal) = (1.-zl) * ptsd(ji,jj,ik,jp_sal) + zl * ptsd(ji,jj,ik+1,jp_sal)
                  END IF
               END DO
            END DO
         ENDIF
         !
      ENDIF
      !
      IF( .NOT.ln_tsd_tradmp ) THEN                   !==   deallocate T & S structure   ==! 
         !                                              (data used only for initialisation)
         IF(lwp) WRITE(numout,*) 'dta_tsd: deallocte T & S arrays as they are only use to initialize the run'
                                        DEALLOCATE( sf_tsd(jp_tem)%fnow )     ! T arrays in the structure
         IF( sf_tsd(jp_tem)%ln_tint )   DEALLOCATE( sf_tsd(jp_tem)%fdta )
                                        DEALLOCATE( sf_tsd(jp_sal)%fnow )     ! S arrays in the structure
         IF( sf_tsd(jp_sal)%ln_tint )   DEALLOCATE( sf_tsd(jp_sal)%fdta )
                                        DEALLOCATE( sf_tsd              )     ! the structure itself
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dta_tsd')
      !
   END SUBROUTINE dta_tsd

   !!======================================================================
END MODULE dtatsd
