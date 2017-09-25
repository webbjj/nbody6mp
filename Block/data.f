      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'common6.h'
      REAL*8  RAN2
*
*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      RN1 = RAN2(KDUM)
*       Skip the first random numbers (IDUM1 specified at input).
      DO 1 K = 1,IDUM1
          RN1 = RAN2(KDUM)
    1 CONTINUE
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
*       Read IMF parameters, # primordials, Z-abundance, epoch & HR interval.
      READ (5,*) ALPHAS, BODY1, BODYN, NBIN0, NHI0, ZMET, EPOCH0, DTPLOT
      IF (N + 2*NBIN0 + 2*NHI0.GE.NMAX - 2) THEN
          WRITE (6,2)  N, NBIN0, NHI0
    2     FORMAT (' FATAL ERROR!    BAD INPUT    N NBIN0 NHI0 ',I7,2I6)
          STOP
      END IF
      IF (ZMET.LE.0.0D0) ZMET = 1.0D-04
      IF (ZMET.GT.0.03) ZMET = 0.03
      IF (KZ(12).GT.0) DTPLOT = MAX(DTPLOT,DELTAT)
*
*     BEGIN NBODY6MP - READ IN SUB-POPULATION PARAMETERS
      IF (KZ(50).GT.1) THEN
        IF (KZ(22).NE.5)THEN
          PRINT *,'MULTIPLE POPULATION ERROR'
          RETURN
        END IF

        DO I=1,KZ(50)
          READ(5,*) NPOPS(I),YPOP(I),ZPOP(I)
          IF (ZPOP(I).LE.0.0D0) ZPOP(I) = 1.0D-04
          IF (ZPOP(I).GT.0.03) ZPOP(I) = 0.03
        END DO
*       DEBUG: TEMP SET ZMET EQUAL TO ZPOP(1)
        ZMET=ZPOP(1)
      END IF

      IF (KZ(19).GE.3) THEN
        IF(KZ(50).LE.1) THEN
            CALL zcnsts(ZMET,ZPARS)
            WRITE (6,4)  ZPARS(11), ZPARS(12), ZMET
    4       FORMAT (//,12X,'ABUNDANCES:  X =',F6.3,'  Y =',F6.3,
     &                                           '  Z =',F7.4)
        ELSE
*           Assign sub-population type its own ZPARS
            DO J=1,KZ(50)
              CALL zcnsts(ZPOP(J),ZPARS)
              ZPARSP(1:20,J)=ZPARS
              ZPARS=0.0
*             If Y is not assigned then take value from ZPARS
              IF(YPOP(J).EQ.0) YPOP(J)=ZPARS(12)

              WRITE (6,*) 'POPULATION #',J,': ',NPOPS(J),' STARS'
              WRITE (6,4)  ZPARSP(11,J), ZPARSP(12,J),
     &                     ZPOP(J)
            END DO
          END IF

      END IF
*     END NBODY6MP
*
*       Check options for reading initial conditions from input file.
*       NBODY6MP - ADD KZ(22) != 5
      IF (KZ(22).GE.2.AND.KZ(22).NE.5.OR.KZ(22).EQ.-1) THEN
          ZMASS = 0.0
          DO 5 I = 1,N
              READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
              ZMASS = ZMASS + BODY(I)
    5     CONTINUE
*       Include possibility of a new IMF via option #20.
          IF (KZ(22).GT.2.OR.KZ(22).EQ.-1) GO TO 50
      END IF

*       BEGIN NBODY6MP - MPTYPE =1 UNLESS OTHERWISE SPECIFIED
      IF (KZ(22).EQ.5)THEN
          DO 6 I = 1,N
              READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3),
     &MPTYPE(I)
              ZMASS = ZMASS + BODY(I)
    6     CONTINUE
          GO TO 50
      END IF
*       END NBODY6MP
*
*       Include optional initial conditions on #10 in astrophysical units.
*     IF (KZ(22).EQ.-1) THEN
*         CALL SETUP2
*         GO TO 30
*     END IF
*
*       Include the case of equal masses (ALPHAS = 1 or BODY1 = BODYN).
      IF (ALPHAS.EQ.1.0.OR.BODY1.EQ.BODYN) THEN
          DO 10 I = 1,N
              BODY(I) = 1.0
   10     CONTINUE
*       Set provisional total mass (rescaled in routine SCALE).
          ZMASS = FLOAT(N)
          GO TO 40
      END IF
*
*       Choose between two realistic IMF's and standard Salpeter function.
      IF (KZ(20).EQ.1) THEN
          CALL IMF(BODY1,BODYN)
          GO TO 30
      ELSE IF (KZ(20).GE.2) THEN
          CALL IMF2(BODY1,BODYN)
          GO TO 30
      END IF
*       Assign #20 = 0 for no new IMF and no scaling with input on fort.10.
      IF (KZ(22).EQ.2.AND.KZ(20).EQ.0) GO TO 50
*
      WRITE (6,15)  ALPHAS, BODY1, BODYN
   15 FORMAT (/,12X,'STANDARD IMF    ALPHAS =',F5.2,
     &              '  BODY1 =',F5.1,'  BODYN =',F5.2)
*
*       Generate a power-law mass function with exponent ALPHAS.
      ALPHA1 = ALPHAS - 1.0
      FM1 = 1.0/BODY1**ALPHA1
      FMN = (FM1 - 1.0/BODYN**ALPHA1)/(FLOAT(N) - 1.0)
      ZMASS = 0.0D0
      CONST = 1.0/ALPHA1
*
*       Assign individual masses sequentially.
      DO 20 I = 1,N
          FMI = FM1 - FLOAT(I - 1)*FMN
          BODY(I) = 1.0/FMI**CONST
          ZMASS = ZMASS + BODY(I)
   20 CONTINUE
*
*       Scale the masses to <M> = 1 for now and set consistent total mass.
   30 ZMBAR1 = ZMASS/FLOAT(N)
      DO 35 I = 1,N
          BODY(I) = BODY(I)/ZMBAR1
   35 CONTINUE
      ZMASS = FLOAT(N)
*
*       Set up initial coordinates & velocities (uniform or Plummer model).
   40 IF (KZ(22).EQ.0.OR.KZ(22).EQ.1) THEN
          CALL SETUP
      END IF
*
   50 RETURN
*
      END
