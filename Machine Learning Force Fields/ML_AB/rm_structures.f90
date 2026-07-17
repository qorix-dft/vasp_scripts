! Read in ML_ABN, eliminate all structures in the section(s) of configuration
! numbers entered by the user, and renumber the remaining structures/basis
! sets accordingly. The result is written to ML_ABN_purge.
! The total number of configurations is printed first; based on it, enter the
! first and last configuration of the section to eliminate (e.g. "30 100"
! removes structures 30 to 100). Sections can be entered repeatedly; enter
! "0 0" to finish.
! NB: -Basis sets belonging to eliminated structures are removed as well, and
!      the number of configurations (line 5) and the numbers of basis sets
!      per atom type in ML_ABN_purge are updated automatically, so no manual
!      editing is needed
!     -Copy/rename the ML_ABN* files as needed
!
! Compile: gfortran -o purge_structures purge_structures.f90
!
  PROGRAM PURGE_STRUCTURES
    IMPLICIT NONE
    INTEGER                         :: NCONFM   ! total no. of configurations
    INTEGER                         :: MTYP     ! no. of species
    INTEGER, ALLOCATABLE            :: NB(:)    ! no. of basis sets before purge
    INTEGER, ALLOCATABLE            :: NBNEW(:) ! no. of basis sets after purge
    INTEGER                         :: ICOUNT   ! no. of configurations kept
    INTEGER                         :: IWRITTEN ! no. of configurations written
    INTEGER, ALLOCATABLE            :: ISTR(:,:), IB(:,:), ICONF(:)
    INTEGER                         :: IOS, I, ITYP, IPOS
    INTEGER                         :: I1, I2   ! section to eliminate
    LOGICAL                         :: LWRITE, LFOUND, LPEND
    LOGICAL, ALLOCATABLE            :: LCONF(:)
    CHARACTER(LEN=200)              :: LINE, LADJ
    CHARACTER(LEN=200), ALLOCATABLE :: ALINE(:)
    !
    ! read and write old/new ML_AB file
    OPEN(10,FILE='ML_ABN',STATUS='OLD',IOSTAT=IOS)
    IF (IOS /= 0) ERROR STOP 'Error: could not open ML_ABN'
    OPEN(11,FILE='ML_ABN_purge',STATUS='REPLACE')
    DO I=1,4
      READ(10,'(a)',IOSTAT=IOS) LINE
      IF (IOS /= 0) ERROR STOP 'Error: unexpected end of ML_ABN in header'
      WRITE(11,'(a)') TRIM(LINE)
    ENDDO
    READ(10,*,IOSTAT=IOS) NCONFM
    IF (IOS /= 0 .OR. NCONFM < 1) ERROR STOP 'Error: could not read the number of configurations'
    ALLOCATE(ICONF(NCONFM), LCONF(NCONFM))
    !
    ! choose the configurations to eliminate
    LCONF = .TRUE.
    WRITE(*,*) 'Total number of configurations = ', NCONFM
    WRITE(*,*) 'Enter the first and last configuration of the section to eliminate'
    WRITE(*,*) '(e.g. "30 100" eliminates structures 30 to 100).'
    WRITE(*,*) 'Sections can be entered repeatedly; enter "0 0" when done:'
    DO
      READ(*,*,IOSTAT=IOS) I1, I2
      IF (IOS > 0) THEN
        WRITE(*,*) 'Invalid input: enter two integers, e.g. "30 100" ("0 0" to finish)'
        CYCLE
      ELSEIF (IOS < 0) THEN
        ERROR STOP 'Error: unexpected end of input'
      ENDIF
      IF (I1 == 0 .AND. I2 == 0) EXIT
      IF (I1 > I2) THEN
        ! swap so that I1 <= I2
        I  = I1
        I1 = I2
        I2 = I
      ENDIF
      IF (I1 < 1 .OR. I2 > NCONFM) THEN
        WRITE(*,*) 'Invalid section: configurations must lie between 1 and ', NCONFM
        CYCLE
      ENDIF
      LCONF(I1:I2) = .FALSE.
      WRITE(*,*) 'Eliminating configurations ', I1, ' to ', I2
    ENDDO
    !
    ! create new config. ind. for the kept configurations
    ICOUNT = 0
    ICONF  = 0
    DO I=1,NCONFM
      IF (LCONF(I)) THEN
        ICOUNT = ICOUNT+1
        ICONF(I) = ICOUNT
      ENDIF
    ENDDO
    IF (ICOUNT == 0) ERROR STOP 'Error: all configurations would be eliminated'
    IF (ICOUNT == NCONFM) WRITE(*,*) 'Warning: nothing eliminated, ML_ABN_purge will equal ML_ABN'
    ! the header directly gets the reduced number of configurations
    WRITE(11,*) ICOUNT
    DO I=1,3
      READ(10,'(a)',IOSTAT=IOS) LINE
      IF (IOS /= 0) ERROR STOP 'Error: unexpected end of ML_ABN in header'
      WRITE(11,'(a)') TRIM(LINE)
    ENDDO
    READ(10,*,IOSTAT=IOS) MTYP
    IF (IOS /= 0 .OR. MTYP < 1) ERROR STOP 'Error: could not read the number of atom types'
    WRITE(11,*) MTYP
    ALLOCATE(NB(MTYP), NBNEW(MTYP), ALINE(MTYP))
    !
    ! copy everything up to the basis sets, then purge/renumber them
    LFOUND = .FALSE.
    DO
      READ(10,'(a)',IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      WRITE(11,'(a)') TRIM(LINE)
      IF (INDEX(LINE,'The numbers of basis') > 0) THEN
        LFOUND = .TRUE.
        READ(10,'(a)',IOSTAT=IOS) LINE
        IF (IOS /= 0) ERROR STOP 'Error: unexpected end of ML_ABN in basis sets'
        WRITE(11,'(a)') TRIM(LINE)
        ! read old no. of basis sets
        READ(10,*,IOSTAT=IOS) (NB(I),I=1,MTYP)
        IF (IOS /= 0) ERROR STOP 'Error: could not read the numbers of basis sets'
        ALLOCATE( ISTR(MTYP,MAXVAL(NB)), IB(MTYP,MAXVAL(NB)) )
        DO ITYP=1,MTYP
          READ(10,'(a)',IOSTAT=IOS) LINE         ! **** separator
          READ(10,'(a)',IOSTAT=IOS) ALINE(ITYP)  ! 'Basis set for ...'
          READ(10,'(a)',IOSTAT=IOS) LINE         ! ---- separator
          IF (IOS /= 0) ERROR STOP 'Error: unexpected end of ML_ABN in basis sets'
          ! read config and basis sets
          DO I=1,NB(ITYP)
            READ(10,*,IOSTAT=IOS) ISTR(ITYP,I), IB(ITYP,I)
            IF (IOS /= 0) ERROR STOP 'Error: could not read a basis set'
            IF (ISTR(ITYP,I) < 1 .OR. ISTR(ITYP,I) > NCONFM) &
              ERROR STOP 'Error: basis set points to a configuration out of range'
          ENDDO
        ENDDO
        ! keep only the basis sets of configurations that survive the purge
        DO ITYP=1,MTYP
          NBNEW(ITYP) = COUNT( LCONF(ISTR(ITYP,1:NB(ITYP))) )
          IF (NBNEW(ITYP) == 0) WRITE(*,*) 'Warning: no basis sets left for atom type ', ITYP
        ENDDO
        WRITE(11,*) (NBNEW(I),I=1,MTYP)
        ! write new config ind and basis sets
        DO ITYP=1,MTYP
          WRITE(11,'(a)') REPEAT('*', 50)
          WRITE(11,'(a)') TRIM(ALINE(ITYP))
          WRITE(11,'(a)') REPEAT('-', 50)
          DO I=1,NB(ITYP)
            IF (LCONF(ISTR(ITYP,I))) WRITE(11,*) ICONF(ISTR(ITYP,I)), IB(ITYP,I)
          ENDDO
        ENDDO
        EXIT
      ENDIF ! done with basis sets
    ENDDO
    IF (.NOT.LFOUND) ERROR STOP 'Error: basis set section not found in ML_ABN'
    !
    ! read config, energies etc. and write the selected ones renumbered
    LWRITE   = .FALSE.
    LPEND    = .FALSE. ! a '****' separator line is pending
    IWRITTEN = 0
    DO
      READ(10,'(a)',IOSTAT=IOS) LINE
      IF (IOS /= 0) EXIT
      IPOS = INDEX(LINE,'Configuration num.')
      IF (IPOS > 0) THEN
        ! config. ind.
        READ(LINE(IPOS+18:),*,IOSTAT=IOS) I
        IF (IOS /= 0 .OR. I < 1 .OR. I > NCONFM) ERROR STOP 'Error: could not read a configuration number'
        IF (LCONF(I)) THEN
          ! write only selected config. with their new config. ind.
          IF (LPEND) WRITE(11,'(a)') REPEAT('*', 50)
          WRITE(11,'(a,i7)') '     Configuration num.', ICONF(I)
          IWRITTEN = IWRITTEN+1
          LWRITE = .TRUE.
        ELSE
          LWRITE = .FALSE.
        ENDIF
        LPEND = .FALSE.
      ELSE
        ! '****' separators are held back: written in front of the next kept
        ! config., or as part of the current config. if more content follows
        IF (LPEND .AND. LWRITE) WRITE(11,'(a)') REPEAT('*', 50)
        LPEND = .FALSE.
        LADJ = ADJUSTL(LINE)
        IF (LADJ(1:4) == '****') THEN
          LPEND = .TRUE.
        ELSEIF (LWRITE) THEN
          WRITE(11,'(a)') TRIM(LINE)
        ENDIF
      ENDIF
    ENDDO ! done with config.
    !
    WRITE(*,*) 'Total number of configurations before = ', NCONFM
    WRITE(*,*) 'Configurations eliminated             = ', NCONFM-ICOUNT
    WRITE(*,*) 'Total number of configurations after  = ', ICOUNT
    WRITE(*,*) 'Basis sets per atom type before       = ', NB
    WRITE(*,*) 'Basis sets per atom type after        = ', NBNEW
    IF (IWRITTEN /= ICOUNT) THEN
      WRITE(*,*) 'WARNING: number of configurations actually written = ', IWRITTEN
      WRITE(*,*) '!CHANGE the number of configurations in ML_ABN_purge to this number! (line 5)'
    ELSE
      WRITE(*,*) 'Number of configurations in ML_ABN_purge (line 5) already updated'
    ENDIF
    !
    CLOSE(10)
    CLOSE(11)
    !
  END PROGRAM PURGE_STRUCTURES
