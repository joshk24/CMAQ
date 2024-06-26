!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!

!------------------------------------------------------------------------!
! WORK IN PROGRESS:                                                      !
! AJ Lenzen started this as a combination of commands in driver.F,       !
!    opGSI.F and rdGSI.F                                                 !
! JD East updated it and merged it into v532+                            !
! BH Henderson made it a single module and simplified the interface      !
!   status: conceptually works, but needs to be compiled and tested      !
!           assimilation scripts will need to be updated substantially   !
!------------------------------------------------------------------------!
      MODULE ASSIMUTIL

#ifdef parallel
      USE SE_MODULES            ! stenex (using SE_UTIL_MODULE)
#else
      USE NOOP_MODULES          ! stenex (using NOOP_UTIL_MODULE)
#endif



      ! ASSIMRUN - run assimilation
      ! opASSIM - output state for assimilation
      ! rdASSIM - read in state from assimilation
c     IMPLICIT NONE
      integer :: ASSIMTSTEP = 10000
      character*1024 :: ASSIMSCRIPT
      logical :: LDOASSIM
     
      CONTAINS

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ASSIMINIT()
          ! ASSIMINIT declares retrieves all global variables from the environment.
          use get_env_module
          ! DOASSIM, ASSIMSCRIPT and ASSIMTSTEP from environment 

          IMPLICIT NONE


          INTEGER          :: VARDEV = -1
          CALL GET_ENV( LDOASSIM, 'DOASSIM', LDOASSIM, VARDEV )
          CALL GET_ENV( ASSIMSCRIPT, 'ASSIMSCRIPT', ASSIMSCRIPT, VARDEV )
          CALL GET_ENV( ASSIMTSTEP, 'ASSIMTSTEP', ASSIMTSTEP, VARDEV )
      END SUBROUTINE

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ASSIMRUN (CGRID, JDATE, JTIME, TSTEP )
          ! ASSIMRUN checks if assimilation should run and executes:
          ! 1) opASSIM
          ! 2) ASSIMSCRIPT
          ! 3) rdASSIM
          USE UTILIO_DEFN
          USE RUNTIME_VARS
          USE CENTRALIZED_IO_MODULE

          IMPLICIT NONE

          REAL, POINTER :: CGRID( :,:,:,: )  ! for initial CONC
          INTEGER      :: JDATE        ! starting date (YYYYDDD)
          INTEGER      :: JTIME        ! starting time (HHMMSS)
          INTEGER      :: TSTEP        ! output timestep (HHMMSS)
          ! Variables to control ASSIM output and script
          LOGICAL ::  exists, STATUS
          INTEGER :: checkcount
          LOGICAL, SAVE :: FIRST = .true.
          character*1024 :: XMSG ! JTK
          character*1024 :: PNAME ! JTK

          IF (FIRST) THEN
              ! make sure LDOASSIM, ASSIMTSTEP and ASSIMSCRIPT are defined.
              call ASSIMINIT()
              FIRST = .false.
          ENDIF

          IF ( LDOASSIM ) then ! switch for ASSIM on or off
              ! check if CMAQ is at an ASSSIMILATION time
              IF ( MOD( TIME2SEC( JTIME ), TIME2SEC( ASSIMTSTEP ) ) .EQ. 0 ) THEN
                  ! Output CTM_ASSIM_1
                  call opASSIM(cgrid,jdate,jtime,tstep(1))
                  IF (mype.eq.0) THEN ! IF PRIMARY PROCESSOR
                    ! ASSIMSCRIPT should assume the existence of the CTM_ASSIM_1
                    ! environmental variable that points to the current state.
                    ! from that file, it can identify the JDATE, JTIME and available
                    ! variables.
                    call system(ASSIMSCRIPT, STATUS)
                    IF (STATUS.ne.0) THEN
                        XMSG='***Failure, ASSIM script failed'
                        CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                    ENDIF
                  ENDIF ! if primary processor
                  ! WAIT for all processors
                  CALL SUBST_BARRIER
                  ! Update CGRID from CTM_ASSIM_1
                  call rdASSIM(jdate,jtime,cgrid)
              ENDIF ! IF ISASSIMTIME
          ENDIF ! END ASSIM SWITCH
      END SUBROUTINE ASSIMRUN

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPASSIM (CGRID, JDATE, JTIME, TSTEP )
          ! OPASSIM outputs CGRID to a netcdf file for assimilation
          ! script to read in
          USE GRID_CONF             ! horizontal & vertical domain specifications
          USE WVEL_DEFN             ! derived vertical velocity component -- BHH: not using this... right?
          USE STD_CONC              ! standard CONC
          USE UTILIO_DEFN           ! M3WARN etc.
          USE SE_MODULES            ! stenex

          IMPLICIT NONE

          INCLUDE SUBST_FILES_ID    ! file name parameters

    ! Arguments:
          REAL, POINTER :: CGRID( :,:,:,: )  ! for initial CONC
          INTEGER      JDATE        ! starting date (YYYYDDD)
          INTEGER      JTIME        ! starting time (HHMMSS)
          INTEGER      TSTEP        ! output timestep (HHMMSS)

    ! Local Variables:

          REAL, ALLOCATABLE :: DBUFF( :,:,: )
          INTEGER      ALLOCSTAT
          INTEGER      STATUS               !  ENVINT status

          CHARACTER( 16 ) :: PNAME = 'OPASSIM'
          CHARACTER( 96 ) :: XMSG = ' '
          CHARACTER( 28 ) :: SSTR = ' species saved to ASSIM file:'

    ! environment variable for no. of layers from bottom to save on CONC file
          CHARACTER( 16 ) :: NLAYS_CONC = 'NLAYS_CONC'
          character(16) CTMASSIM

    !  environment variable description
          CHARACTER( 80 ) :: VARDESC

    !JDE      INTEGER      LOGDEV

          INTEGER      K
          INTEGER :: KD=5
          INTEGER :: L=1
          INTEGER :: SPC=0   ! loop counters
          INTEGER      INDX
          INTEGER ::     INO2=-1
          INTEGER ::     VNO2=-1 ! pointer to NO2 information

    !-----------------------------------------------------------------------

          CALL SUBST_BARRIER

    ! Set output file characteristics based on COORD.EXT and open it
    ! need to replace output with a revolving time file
          CTMASSIM='CTM_ASSIM_1'

          FTYPE3D = GRDDED3
          SDATE3D = JDATE
          STIME3D = JTIME
          TSTEP3D = TSTEP
          NVARS3D = 1
          NCOLS3D = GL_NCOLS
          NROWS3D = GL_NROWS
          NLAYS3D = NLAYS
          NTHIK3D = 1
          GDTYP3D = GDTYP_GD
          P_ALP3D = P_ALP_GD
          P_BET3D = P_BET_GD 
          P_GAM3D = P_GAM_GD
          XORIG3D = XORIG_GD
          YORIG3D = YORIG_GD
          XCENT3D = XCENT_GD
          YCENT3D = YCENT_GD
          XCELL3D = XCELL_GD
          YCELL3D = YCELL_GD
          VGTYP3D = VGTYP_GD
          VGTOP3D = VGTOP_GD
          DO L = 1, NLAYS3D + 1
             VGLVS3D( L ) = VGLVS_GD( L )
          END DO
          GDNAM3D = GRID_NAME  ! from HGRD_DEFN

          FDESC3D( 1 ) = 'Concentration file output'
          FDESC3D( 2 ) = 'From CMAQ model dyn alloc version CTM'
          FDESC3D( 3 ) = 'Set of variables just NO2'
          FDESC3D( 4 ) = 'For next scenario continuation runs,'
          FDESC3D( 5 ) = 'use the "one-step" CGRID file'
          KD = 5
          DO L = 1, MIN ( NLAYS3D, MXDESC3 - KD )
             WRITE( FDESC3D( L + KD ),'( "Layer", I3, " to", I3, " " )' )
     &   L, L
          END DO
          IF ( ( KD + 1 + L ) .LT. MXDESC3 ) THEN
             DO K = KD + 1 + L, MXDESC3
                FDESC3D( K ) = ' '
             END DO
          END IF

          ino2=-1
          vno2=-1
          DO SPC = 1, N_C_GC_SPC
             if( C_GC_SPC( SPC ).eq.'NO2')then
               INO2=SPC
               VNO2=SPC
               ! FOUND NO2, so exit loop
               exit
             endif
          END DO
          if (INO2<1)then
                XMSG = 'Could not find NO2'
                CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          endif
          VTYPE3D( 1 ) = M3REAL
          VNAME3D( 1 ) = C_GC_SPC( VNO2 )
          UNITS3D( 1 ) = 'ppmV'
          VDESC3D( 1 ) = 'Variable ' // VNAME3D( 1 )

    ! create header
          IF ( MYPE .EQ. 0 ) THEN   ! open new
             IF ( .NOT. OPEN3( CTMASSIM, FSNEW3, PNAME ) ) THEN
                XMSG = 'opASSIM Could not open ' // CTMASSIM // '. File may already exist.' 
                CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
             END IF
          END IF
          CALL SUBST_BARRIER
          if ( .not. io_pe_inclusive ) then !JDE
            IF ( .NOT. OPEN3( CTMASSIM, FSREAD3, PNAME ) ) THEN
               XMSG = 'opASSIM Could not open ' // CTMASSIM
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
          end if !JDE if PE0
          CALL SUBST_BARRIER

    ! write the initial concentrations as step 0 on the conc file
    ! (inital data assumed to be in correct output units)

          ALLOCATE ( DBUFF( NCOLS,NROWS,NLAYS ), STAT = ALLOCSTAT )
          IF ( ALLOCSTAT .NE. 0 ) THEN
             XMSG = 'Failure allocating DBUFF'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

          WRITE( LOGDEV,* ) ' '
          INDX = CONC_MAP( VNO2 )

          DBUFF = CGRID( 1:NCOLS,1:NROWS,1:NLAYS,INDX )

          IF ( .NOT. WRITE3( CTMASSIM, C_GC_SPC( VNO2 ),
     &                    JDATE, JTIME, DBUFF ) ) THEN
              XMSG = 'opASSIM Could not write ' //
     &           TRIM( C_GC_SPC( VNO2 ) ) //
     &           ' to ' // CTMASSIM
              CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

          DEALLOCATE ( DBUFF )

          WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &  'Timestep written to', CTMASSIM,
     &  'for date and time', JDATE, JTIME
          WRITE( LOGDEV, '(  5X,  A, 1X, I8, ":", I6.6 )' )
     &  'from timestep on initial data files for date and time',
     &   JDATE, JTIME
          IF ( .NOT. CLOSE3( CTMASSIM ) ) THEN 
             XMSG = 'opassim Could not close file "' // TRIM( CTMASSIM ) // '"'
             CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT2 )
          END IF

          RETURN
      END SUBROUTINE

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE rdASSIM(JDATE,JTIME,CGRID)
          ! OPASSIM outputs CGRID to a netcdf file for assimilation
          ! script to read in
          USE GRID_CONF             ! horizontal & vertical domain specifications
          USE CGRID_SPCS            ! CGRID mechanism species
          USE UTILIO_DEFN
          USE STD_CONC

          IMPLICIT NONE

          INCLUDE SUBST_CONST       ! constants
          INCLUDE SUBST_FILES_ID    ! file name parameters

          CHARACTER( 16 ) :: FNAME = 'CTM_ASSIM_1'
          INTEGER      JDATE
          INTEGER      JTIME,l,r,c
          REAL, POINTER :: CGRID( :,:,:,: )  ! for initial CONC
          REAL      :: DBUFF( NCOLS,NROWS,NLAYS )
          INTEGER       :: STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI
          INTEGER      GXOFF, GYOFF               ! global origin offset from file
          CHARACTER( 96 ) :: XMSG = ' '
          CHARACTER( 16 ) :: PNAME = 'RDASSIM'
          INTEGER :: INDX
          INTEGER :: INO2=-1,VNO2=-1 ! pointer to NO2 information
          INTEGER :: SPC=0

          IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
             XMSG = 'Could not open ' // FNAME // ' file'
             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF
          CALL SUBHFILE ( fname, GXOFF, GYOFF,
     &                STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI )
          DO SPC = 1, N_C_GC_SPC
             if( C_GC_SPC( SPC ).eq.'NO2')then
               INO2=SPC
               VNO2=SPC
               exit
             endif
          END DO
          if (INO2<1)then
                XMSG = 'Could not find NO2'
                CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          endif
          IF ( .NOT. XTRACT3( FNAME, VNAME3D ( 1 ),
     &           1,NLAYS, STRTROWINI,ENDROWINI, STRTCOLINI,ENDCOLINI,
     &           JDATE, JTIME, DBUFF ) ) THEN
                   XMSG = 'Could not read ' // TRIM( VNAME3D( 1 ) )
     &              // ' from ' // FNAME
                   CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF
          INDX = CONC_MAP( VNO2 )
          DO L = 1, NLAYS
            DO R = 1, NROWS
              DO C = 1, NCOLS
                CGRID( C,R,L,INDX ) = DBUFF( C,R,L )
              END DO
            END DO
          END DO
          IF ( .NOT. CLOSE3( fname ) ) THEN 
            XMSG = 'rdassim Could not close file "' // TRIM( fname ) // '"'
            CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT2 )
          END IF
          RETURN
      END SUBROUTINE rdASSIM


      END MODULE ASSIMUTIL


