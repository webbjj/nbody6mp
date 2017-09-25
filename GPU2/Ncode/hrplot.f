      SUBROUTINE HRPLOT
*
*
*       HR diagnostics of evolving stars.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  M0,M1,M2,LUM,LUM2,MC,ME,OSPIN,OSPIN2,EP1,EP2
      REAL*8  MENV,RENV,RCC,K2,K3
      PARAMETER(K3=0.21d0)
      REAL*8  XS(3),VS(3)
      REAL*8  RDD(NMAX),XDD(NMAX),MTOT,RDMAX,RHALF,RHOAVE
      INTEGER II,MLIST(NMAX),NCNT
      INTEGER NGHOST,JGHOST(MMAX+10),IGHOST(NMAX),IX
*
*
      TPHYS = (TIME+TOFF)*TSTAR + EPOCH0
*      WRITE(82,100)NPAIRS,TPHYS
      WRITE(82,*)NPAIRS,TPHYS
 100  FORMAT(I8,F14.5)
      NS = N - 2*NPAIRS
      IMERGE = 0
      NCNT = 0
      II = 0
      MTOT = 0.D0
      RDMAX = 0.D0
      NGHOST = 0
*JW - sometimes need full version of parameters
*      WRITE(83,101)NS,TPHYS,TCR,TSCALE
      WRITE(83,*)NS,TPHYS,TCR,TSCALE
*      WRITE(83,102)NC,RC,RBAR,RTIDE,(RDENS(K),K=1,3)
      WRITE(83,*) NC,RC,RBAR,RTIDE,(RDENS(K),K=1,3)   
*      WRITE(83,103)ZMBAR,MIN(TURN,99.999),RSCALE
      WRITE(83,*)ZMBAR,MIN(TURN,99.999),RSCALE

 101  FORMAT(I8,F14.5,2F6.2)
 102  FORMAT(I8,3F6.2,3F10.5)
 103  FORMAT(F12.4,2F6.2)
*
      write(6,*)' splits ',ifirst,n,ntot,nzero
      DO 20 I = 1,N
          M0 = BODY0(I)*ZMBAR
          M1 = BODY(I)*ZMBAR
*       Replace ghost mass of single star with original value from merged KS.
          IX = I
          IF (M1.EQ.0.0.AND.I.GE.IFIRST) THEN
              IM = 0
*       Search merger table for ghost to identify corresponding index.
              DO 2 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(I)) THEN
                      IM = K
                  END IF
    2         CONTINUE
*       Skip any ghosts associated with chain regularization.
              IF (IM.EQ.0) THEN
                  WRITE (6,3)  I, NCH
    3             FORMAT (' WARNING!    HRPLOT   I NCH ',I5,I4)
                  GO TO 20
              END IF
              M1 = CM(2,IM)*ZMBAR
              IM = 0
              DO 22 K = 1,NGHOST
                  IF(JGHOST(K).EQ.I) IM = K
   22         CONTINUE
              if(i.eq.7121) write(6,*)' im ',im
              IF (IM.GT.0) THEN
                 IX = IGHOST(I)
                 M1 = BODY(IX)*ZMBAR
              ENDIF
          END IF
*
*       Obtain stellar parameters at current epoch.
          IF(KZ(19).GT.2)THEN
             KW = KSTAR(I)
             AGE = MAX(TPLOT,TEV0(I))*TSTAR - EPOCH(I)
             MC = 0.D0
*            BEGIN NBODY6MP
             ZPARS=ZPARSP(1:20,MPTYPE(I))
*            END NBODY6MP
             CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
             CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                   RM,LUM,KW,MC,RCC,ME,RE,K2)
             OSPIN = MAX(SPIN(I)*SPNFAC,1.0D-10)/
     &                   (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
             EP1 = MIN(EPOCH(I)+EPOCH0+TOFF*TSTAR,999999.999)
             EP1 = MAX(EP1,-999999.999)
          ELSE
             KW = 0
             RM = 1.0D-10
             LUM = 1.0D-10
             OSPIN = 0.D0
             EP1 = 0.D0
          ENDIF
          IF (I.LT.IFIRST) THEN
              JPAIR = KVEC(I)
              J2 = 2*JPAIR
              IF (I.EQ.J2) GO TO 20
              J1 = 2*JPAIR - 1
              ICM = N + JPAIR
              M2 = BODY(J2)*ZMBAR
              RI = (X(1,ICM) - RDENS(1))**2 + (X(2,ICM) - RDENS(2))**2 +
     &                                        (X(3,ICM) - RDENS(3))**2
              DO K = 1,3
*                XS(K) = X(K,ICM) - RDENS(K)
                 XS(K) = X(K,ICM)
                 VS(K) = XDOT(K,ICM)
              ENDDO
              if(i.eq.4195) write(6,*)' jpair ',jpair,name(n+jpair),
     &                      j1,j2,m1,m2
*       Check for ghost binary.
              IF (M1.EQ.0.0) THEN
                  IM = 0
*       Search merger table to identify corresponding index of c.m. name.
                  DO 4 K = 1,NMERGE
                      IF (NAMEM(K).EQ.NAME(ICM)) THEN
                          IM = K
                      END IF
    4             CONTINUE
                  IF (IM.EQ.0) GO TO 20
*       Copy masses and obtain evolution parameters for first component.
                  M1 = CM(3,IM)*ZMBAR
                  M2 = CM(4,IM)*ZMBAR
                  IF(KZ(19).GT.2)THEN
                     KW = KSTAR(J1)
                     AGE = MAX(TPLOT,TEV0(J1))*TSTAR - EPOCH(J1)
                     MC = 0.D0
*                    BEGIN NBODY6MP
                     ZPARS=ZPARSP(1:20,MPTYPE(J1))
*                    END NBODY6MP
                     CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
                     CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                           RM,LUM,KW,MC,RCC,ME,RE,K2)
                     OSPIN = MAX(SPIN(J1)*SPNFAC,1.0D-10)/
     &                           (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
                     EP1 = MIN(EPOCH(J1)+EPOCH0+TOFF*TSTAR,999999.999)
                     EP1 = MAX(EP1,-999999.999)
                  ELSE
                     KW = 0
                     RM = 1.0D-10
                     LUM = 1.0D-10
                     OSPIN = 0.D0
                     EP1 = 0.D0
                  ENDIF
              END IF
              RJ = R(JPAIR)
              HJ = H(JPAIR)
*       Determine merger & ghost index for negative c.m. name.
              IF (NAME(N+JPAIR).LT.0) THEN
                  CALL FINDJ(J1,J,IMERGE)
                  NGHOST = NGHOST + 1
                  JGHOST(NGHOST) = J
                  IGHOST(J) = J2
                  write(6,*)' new ghost ',nghost,i,j,j2,name(j)
*       Skip second binary of quadruple.
                  IF (NAME(J).GT.NZERO) GO TO 20
                  M1 = CM(1,IMERGE)*ZMBAR
                  M2 = CM(2,IMERGE)*ZMBAR
                  HJ = HM(IMERGE)
                  RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
     &                                          XREL(3,IMERGE)**2)
                  DO K = 1,3
                     XS(K) = X(K,J1)
                     VS(K) = XDOT(K,J1)
                  ENDDO
*       Re-define index of second component and obtain parameters of M1.
                  J2 = J
                  IF(KZ(19).GT.2)THEN
                     AGE = MAX(TPLOT,TEV0(J1))*TSTAR - EPOCH(J1)
                     KW = KSTAR(J1)
                     MC = 0.D0
*                    BEGIN NBODY6MP
                     ZPARS=ZPARSP(1:20,MPTYPE(J1))
*                    END NBODY6MP
                     CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
                     CALL HRDIAG(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                           RM,LUM,KW,MC,RCC,ME,RE,K2)
                     OSPIN = MAX(SPIN(J1)*SPNFAC,1.0D-10)/
     &                           (K2*RM*RM*(M1-MC)+K3*RCC*RCC*MC)
                     EP1 = MIN(EPOCH(J1)+EPOCH0+TOFF*TSTAR,999999.999)
                     EP1 = MAX(EP1,-999999.999)
                  ELSE
                     KW = 0
                     RM = 1.0D-10
                     LUM = 1.0D-10
                     OSPIN = 0.D0
                     EP1 = 0.D0
                  ENDIF
              END IF
              M0 = BODY0(J2)*ZMBAR
              IF(KZ(19).GT.2)THEN
                 KW2 = KSTAR(J2)
                 AGE = MAX(TPLOT,TEV0(J2))*TSTAR - EPOCH(J2)
                 MC = 0.D0
*                BEGIN NBODY6MP
                 ZPARS=ZPARSP(1:20,MPTYPE(J1))
*                END NBODY6MP
                 CALL STAR(KW2,M0,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
                 CALL HRDIAG(M0,AGE,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                       RM2,LUM2,KW2,MC,RCC,ME,RE,K2)
                 OSPIN2 = MAX(SPIN(J2)*SPNFAC,1.0D-10)/
     &                        (K2*RM2*RM2*(M2-MC)+K3*RCC*RCC*MC)
                 EP2 = MIN(EPOCH(J2)+EPOCH0+TOFF*TSTAR,999999.999)
                 EP2 = MAX(EP2,-999999.999)
              ELSE
                 KW2 = 0
                 RM2 = 1.0D-10
                 LUM2 = 1.0D-10
                 OSPIN2 = 0.D0
                 EP2 = 0.D0
              ENDIF
              RI = SQRT(RI)/RC
*       Specify relevant binary mass.
              IF (BODY(J1).GT.0.0D0) THEN
                  BODYI = (M1 + M2)/ZMBAR
              ELSE
                  BODYI = CM(3,IMERGE) + CM(4,IMERGE)
              END IF 
              SEMI = -0.5*BODYI/HJ 
              ECC2 = (1.0 - RJ/SEMI)**2
              ECC = SQRT(ECC2)
              PB = DAYS*SEMI*SQRT(ABS(SEMI)/BODYI)
              PB = MIN(PB,999999.9D0)
              PB = LOG10(ABS(PB))
              SEMI = LOG10(ABS(SEMI*SU))
              R1 = LOG10(RM)
              R2 = LOG10(RM2)
              ZL1 = LOG10(LUM)
              ZL2 = LOG10(LUM2)
              TE1 = 0.25*(ZL1 - 2.0*R1) + 3.7
              TE2 = 0.25*(ZL2 - 2.0*R2) + 3.7
*             WRITE (82,5)  NAME(J1), NAME(J2), KW, KW2, KSTAR(ICM),
*    &            RI, ECC, PB, SEMI, M1, M2, ZL1, ZL2, R1, R2, TE1, TE2
*   5         FORMAT (2I6,2I3,I4,F6.1,F6.3,10F7.3)
* JW - Change Write(82,123) to 82,* for all decimal places 
              WRITE(82,*)NAME(J1),NAME(J2),KW,KW2,KSTAR(ICM),
     &            ECC,PB,SEMI,M1,M2,ZL1,ZL2,R1,R2,
     &            (XS(K),K=1,3),(VS(K),K=1,3), 
     &            EP1,EP2,OSPIN,OSPIN2,
*    BEGIN NBODY6MP
     &MPTYPE(J1),MPTYPE(J2)
*    END NBODY6MP  
              NCNT = NCNT + 2
              if(name(j1).eq.13.or.name(j2).eq.13)then
                 write(6,*)' bin ',name(j1),name(j2),j1,j2,i
              endif
              if(name(j1).eq.14.or.name(j2).eq.14)then
                 write(6,*)' bin ',name(j1),name(j2),j1,j2,i
              endif
              if(name(j1).eq.4173.or.name(j2).eq.4173)then
                 write(6,*)' bin ',name(j1),name(j2),j1,j2,i
              endif
              IF(KZ(12).GE.3)THEN
                 II = II + 1
                 RI = XS(1)**2 + XS(2)**2 + XS(3)**2
                 IF(RI.GT.0.0)THEN
                    RI = SQRT(RI)
                    RI = MIN(RI,2.D0*RTIDE)
                    RDMAX = MAX(RDMAX,RI)
                 ENDIF
                 RDD(II) = RI
                 XDD(II) = M1 + M2
                 MLIST(II) = II
                 MTOT = MTOT + XDD(II)
              ENDIF
          ELSE
*       Create output file for single stars (skip rare chain subsystem).
              IF (NAME(I).EQ.0) GO TO 20
              IF (M1.LE.0.0)THEN
                 write(6,*)' zero ',i,name(ix),kw,m1
*                GO TO 20
              ENDIF
              DO K = 1,3
*                XS(K) = X(K,I) - RDENS(K)
                 XS(K) = X(K,IX)
                 VS(K) = XDOT(K,IX)
              ENDDO
              RI = (XS(1) - RDENS(1))**2 + (XS(2) - RDENS(2))**2 +
     &                                     (XS(3) - RDENS(3))**2
              RI = SQRT(RI)/RC
              RI = MIN(RI,99.0D0)
              R1 = LOG10(RM)
              ZL1 = LOG10(LUM)
*       Form LOG(Te) using L = 4*pi*R**2*\sigma*T**4 and solar value 3.7.
              TE = 0.25*(ZL1 - 2.0*R1) + 3.7
*             WRITE (83,10)  NAME(I), KW, RI, M1, ZL1, R1, TE
*  10         FORMAT (I6,I3,F6.1,4F7.3)
* JW - Change WRITE(83,121) to (83,*) if necessary
              WRITE(83,*)NAME(IX),KW,M1,ZL1,R1, 
     &                     (XS(K),K=1,3),(VS(K),K=1,3),EP1,OSPIN,
* BEGIN NBODY6MP
     &                     MPTYPE(IX) 
* END NBODY6MP
              if(xs(1).gt.9999.9)then
                 write(6,*)' BIGX ',i,name(i),m1,xs(1)
*                STOP
              endif
              if(name(ix).eq.13.or.name(ix).eq.4173.or.
     &           name(ix).eq.14)then
                 write(6,*)' sgl ',name(ix),ix,i
              endif
              NCNT = NCNT + 1
              IF(KZ(12).GE.3)THEN
                 II = II + 1
                 RI = XS(1)**2 + XS(2)**2 + XS(3)**2
                 IF(RI.GT.0.0)THEN
                    RI = SQRT(RI)
                    RI = MIN(RI,2.D0*RTIDE)
                    RDMAX = MAX(RDMAX,RI)
                 ENDIF
                 RDD(II) = RI
                 XDD(II) = M1
                 MLIST(II) = II
                 MTOT = MTOT + XDD(II)
              ENDIF
          END IF
   20 CONTINUE
*
      KW = -1000
      KW2 = 1
      SEMI = 0.0
      WRITE(82,*)KW,KW,KW2,KW2,KW2,
     &             SEMI,SEMI,SEMI,M1,M1,ZL1,ZL1,R1,R1,
     &             (XS(K),K=1,3),(VS(K),K=1,3),SEMI,SEMI,OSPIN,OSPIN,
*    BEGIN NBODY6MP
     &             MPTYPE(1),MPTYPE(2)
*    END NBODY6MP
      WRITE(83,*)KW,KW2,M1,ZL1,R1,(XS(K),K=1,3),(VS(K),K=1,3),
     &             SEMI,OSPIN,
*    BEGIN NBODY6MP
     &             MPTYPE(1)
*    END NBODY6MP
 121  FORMAT(I7,I3,3F8.3,1P,6E14.6,0P,F12.3,1P,E12.4)
 123  FORMAT(2I7,2I3,I4,F6.3,F10.2,7F8.3,1P,6E14.6,0P,2F12.3,1P,2E12.4)
*
      WRITE(6,*)' HRPLOT CHECK: N NCNT ',N,NCNT
      call flush(6)
      IF(KZ(12).GE.3)THEN
*
*       Sort radii into increasing order.
         CALL SORT1(II,RDD,MLIST)
*
*       Determine the half-mass radius (pc) and relaxation time (Myrs).
         MC = 0.D0
         DO 80 , I = 1,II
            MC = MC + XDD(MLIST(I))
            IF(MC.GE.0.5*MTOT) GOTO 85
 80      CONTINUE
 85      CONTINUE
         RHALF = RDD(I)*RBAR
         TRH = FLOAT(II)/LOG10(0.4D0*FLOAT(II))
         TRH = 0.858D0*TRH*SQRT(RHALF**3/MTOT)
         RHOAVE = 0.23873*ZMBAR/(RDMAX**3)
         WRITE(6,125)RHALF,TRH
 125     FORMAT('  HALF-MASS RH TRH ',F6.2,F12.4)
         WRITE(6,*)' RT RMAX RHO ',RTIDE*RBAR,RDMAX*RBAR,RHOAVE
      ENDIF
*
*       Update next plotting time.
      TPLOT = TPLOT + DTPLOT
      CALL FLUSH(6)
      CALL FLUSH(82)
      CALL FLUSH(83)
*
      RETURN
*
      END
