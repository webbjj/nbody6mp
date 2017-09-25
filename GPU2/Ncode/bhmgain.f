      SUBROUTINE BHMGAIN
*
*
*       Mass gain for SMBHs
*       ------------------------------
*
*       Original scheme of Stone, Kupper, and Ostriker 2016.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 TBH

      DO I=1,N

          TBH=(-3.0)*(BODY(I)**(-1.0/3.0)-BODY0(I)**(-1.0/3.0))/
     &(TIME-TOFF)

          BODY(I)=-1.0*(TIME-TOFF)/(3.0*TBH)+BODY0(I)**(-1.0/3.0)
          BODY(I)=BODY(I)**(-3.0)

      END DO 


      END 
