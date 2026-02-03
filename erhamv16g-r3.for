      PROGRAM ERHAM
      IMPLICIT REAL(8) (A-H,O-Z)
C
C    Version V16g-R3 (20 may 2013)
C
C    Revision R1: Tunneling parameter matrix elements corrections
C       Coriolis type operators (odd powers of Px, Pz) [1/7/09]
C       High powers of Pz as used in octic CD constants [10/12/09]
C
C    Revision R3: Added option (using IFPR) to print predictions in JPL
C       catalog file format as separate file, corrected formula for
c       intensities, threshold now applies to calculated intensity
c       [previously to S*mu^2], changed dimensions in ORDER to 140000
c       transitions, added calculation of sum of states, increased
c       number of tunneling parameters per state to 37, introduced
c       parameter scaling for least-squares fit [5/20/13]
c
c: 5/20/13
c:      DIMENSION A(50,6),INPAR(34,6),IVR(50,6),JMIN(6),JMAX(6),FMIN(6),
c:     * FMAX(6),THRES(6),NTUP(6),NTE(6),NSIG(6),ISIG(4,10,6),FRQ(8191),
c:     * WT(8191),BL(8191),ITRA(8,8191),ILEV(2,16383),DIP(3)
      DIMENSION A(60,6),INPAR(38,6),IVR(60,6),JMIN(6),JMAX(6),FMIN(6),
     * FMAX(6),THRES(6),NTUP(6),NTE(6),NSIG(6),ISIG(4,10,6),FRQ(8191),
     * WT(8191),BL(8191),ITRA(8,8191),ILEV(2,16383),DIP(3),SCP(80)
ce 5/20/13
      DIMENSION IVL(8)
      CHARACTER*10 IB
      CHARACTER*4 MONTH(12)
      DATA MONTH /' jan',' feb',' mar',' apr',' may',' jun',
     *            ' jul',' aug',' sep',' oct',' nov',' dec'/
      CHARACTER*20 FILEIN,FILEOUT
      WRITE (*,*) 'Enter input file name !'
      READ (*,*) FILEIN
      OPEN(5,FILE=FILEIN)
      WRITE (*,*) 'Enter output file name !'
      READ (*,*) FILEOUT
      OPEN(6,FILE=FILEOUT)
      WRITE (*,*) 'program ERHAM V16g-R3 20 may 2013'
      CALL DATE_AND_TIME (IB,IB,IB,IVL)
      WRITE (6,4001) IVL(3),MONTH(IVL(2)),IVL(1),IVL(5),IVL(6),IVL(7)
 4001 FORMAT(' program ERHAM V16g-R3 20 may 2013  ***  date and time:'
     *  ,i3,a4,2i5,2(':',i2)/)
      CALL INPUT(A,MQ,N1,N2,NC,INPAR,IVR,NTUP,JMIN,JMAX,NTRA,FRQ,BL,WT,
     *           ITRA,ILEV,NIT,FMIN,FMAX,NSIG,ISIG,NTE,THRES,LP,NIV,
c: 5/20/13
c:     *           IFPR,TEMP,NUNC,DIP,ISCD)
c:      IF (NIT.GT.0) CALL ITER(A,MQ,N1,N2,NC,INPAR,IVR,NTUP,NTE,NTRA,
c:     *                        FRQ,BL,WT,ITRA,ILEV,NIT,LP,NIV,IFPR,ISCD)
     *           IFPR,TEMP,NUNC,DIP,ISCD,SCP)
      IF (NIT.GT.0) CALL ITER(A,MQ,N1,N2,NC,INPAR,IVR,NTUP,NTE,NTRA,
     *                     FRQ,BL,WT,ITRA,ILEV,NIT,LP,NIV,IFPR,ISCD,SCP)
ce 5/20/13
      CALL FTEST(NIT,MQ,N1,N2,NC,NTE,A,INPAR,IVR,NIV,LP)
      CALL PREDIC(A,MQ,N1,N2,NC,INPAR,NTUP,NTE,JMIN,JMAX,FMIN,
     *     FMAX,NSIG,ISIG,THRES,NIV,TEMP,IVR,LP,NUNC,DIP,IFPR,ISCD)
      CALL DATE_AND_TIME (IB,IB,IB,IVL)
      WRITE (6,4001) IVL(3),MONTH(IVL(2)),IVL(1),IVL(5),IVL(6),IVL(7)
      END
      SUBROUTINE FTEST(NIT,MQ,N1,N2,NC,NTE,A,INPAR,IVR,NIV,LP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION RHO(6,2),B(3),IQ1(10),IQ2(10),IVX(6,9),EPS(6,9),RES(30),
c: 5/20/13
c:     *          STER(27),A(50,6),INPAR(34,6),NTE(6),IVR(50,6)
     *          STER(27),A(60,6),INPAR(38,6),NTE(6),IVR(60,6)
ce 5/20/13
      COMMON D(2,243,243),D1(243,243),D2(243,243),PHI1(970),C(243,243),
     *       U(2,243,243),E(243),EJ(243),H(2,243,243),EW(243,243)
      DATA PI,C0,ZERO/3.141592653589793D0,29979.2458D0,0D0/
      DEGRAD=180D0/PI
c     MQ:   equivalence parameter (=1 for equiv, =2 for non-equiv)
c     N1,N2:periodicities of 1st and 2nd rotor
c     NC:   direction cosine parameter
c     NTE:  number of tunneling energy parameters
c     RHO:  rho-vectors in the order rho1,rho2,beta1,beta2,alpha1,alpha2
c     B:    rotational constants
c     IQ1,IQ2: q's of the epsilon parameter
c     EPS:  epsilon parameter
c     EW:   variance-covariance matrix (= inverted normal matrix) times
c           standard deviation (transmitted in COMMON block)
      NTEE=0
      I2=0
      DO 3 I=4,6
         DO 2 K=1,2
            I2=I2+1
    2       RHO(I,K)=A(I2,1)
    3    B(I-3)=A(I+3,1)
c: 5/20/13
      IF (NC.EQ.1) RHO(6,2)=RHO(6,2)+PI
ce 5/20/13
      DO 5 IV=1,NIV
         NTEE=NTEE+NTE(IV)
         DO 5 I=1,NTE(IV)
            K=INPAR(I,IV)/256
            IQ2(I)=MOD(K/64,16)-8
            IQ1(I)=K/1024-8
            IVX(IV,I)=IVR(21+I,IV)
    5       EPS(IV,I)=A(21+I,IV)
      IF (NC.EQ.-1) RHO(5,2)=PI-RHO(5,2)
      DO 10 I=1,(MQ+1)*3+NTEE
         DO 10 K=1,(MQ+1)*3+NTEE
   10       C(K,I)=ZERO
      IF (NIT.EQ.0) GOTO 50
      L=0
      I2=9
      DO 15 I=1,I2
         IF (MQ.EQ.1.AND.(I.EQ.2.OR.I.EQ.4.OR.I.EQ.6)) GOTO 15
         N=IVR(I,1)
         L=L+1
         IF (N.EQ.0) THEN
            DO 11 K=1,LP
   11          D1(L,K)=ZERO
         ELSE
            DO 12 K=1,LP
   12          D1(L,K)=EW(N,K)
         ENDIF
   15    CONTINUE
      DO 20 IV=1,NIV
         DO 20 IP=1,NTE(IV)
            N=IVR(21+IP,IV)
            IF (N.GT.0) THEN
               L=L+1
               DO 19 K=1,LP
   19          D1(L,K)=EW(N,K)
            ENDIF
   20    CONTINUE
      I1=0
      DO 35 I=1,I2
         IF (MQ.EQ.1.AND.(I.EQ.2.OR.I.EQ.4.OR.I.EQ.6)) GOTO 35
         N=IVR(I,1)
         I1=I1+1
         IF (N.NE.0) THEN
            DO 32 K=1,L
   32          C(K,I1)=D1(K,N)
         ENDIF
   35    CONTINUE
      DO 40 IV=1,NIV
         DO 40 IP=1,NTE(IV)
            N=IVR(21+IP,IV)
            IF (N.GT.0) THEN
               I1=I1+1
               DO 39 K=1,L
   39          C(K,I1)=D1(K,N)
            ENDIF
   40       CONTINUE
   50 CALL DERPAR(NTE,NTEE,N1,N2,MQ,NC,B,RHO,IQ1,IQ2,IVX,EPS,C,RES,STER,
     *            L,NIV)
      WRITE (6,951) (RES(K),STER(K),K=1,5)
  951 FORMAT(/'DERIVED PARAMETERS',20X,'VALUE',7X,'STD ERROR'//
     *       'INTERNAL MOMENTS OF INERTIA (u*A**2)',2(/35X,2F12.7)/
     *    'REDUCED INTERNAL MOMENTS OF INERTIA (u*A**2)',3(/35X,2F12.7))
      WRITE (6,953) (RES(K)*DEGRAD,STER(K)*DEGRAD,K=6,14)
  953 FORMAT('ANGLES  (A,1)',22X,2F12.7/8X,'(B,1)',22X,2F12.7/8X,
     *  '(C,1)',22X,2F12.7/8X,'(A,2)',22X,2F12.7/8X,'(B,2)',22X,2F12.7/
     *  8X,'(C,2)',22X,2F12.7/'ALPHAQ1',28X,2F12.7/'ALPHAQ2',28X,2F12.7/
     *  'OMEGA',30X,2F12.7)
      WRITE (6,952) (RES(K+14),STER(K+14),K=1,3),(RES(K+14)/C0,
     *               STER(K+14)/C0,K=1,3)
  952 FORMAT('F-NUMBERS (F1,F2,F'') (MHz)',3(/35X,2F12.2)/
     *  'F-NUMBERS (F1,F2,F'') (CM-1)',3(/35X,2F12.7))
      WRITE (6,954) (RES(K),STER(K),K=18,L)
  954 FORMAT('TORSIONAL ENERGY DIFFERENCES (MHz)'/(F47.4,F12.4))
      END
      SUBROUTINE DERPAR(NTE,NTEE,N1,N2,MQ,NC,B,RHO,IQ1,IQ2,IVX,EPS,C,
     *                  RES,STER,L,NIV)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION DD(27,30),BIN(3),DIN(3),B(3),RHO(6,2),BITAU(2),RES(30),
     *          RITAU(3),DIRC(7,2),F(3),U(3,3),IQ1(10),IQ2(10),IVX(6,9),
     *          EPS(6,9),C(243,243),AB(27),AA(27,27),STER(27),NTE(6)
c                 1:RHO(a), 2:RHO(b), 3:RHO(c), 4:RHO, 5:BETA, 6:ALPHA
c                 DIRC: direction cosines (1-3) and angles (4-6), alphaQ(7)
c                 BIN: principal moments of inertia
c                 DIN: derivative of BIN w.r.t. rot constants B
c                 BITAU: moment of inertia of internal rotor
c                 RITAU: reduced moments of inertia
c                 F: internal rotation constants
c                 OMEG: angle between internal rotation axes
c                 C: correlation matrix
      DATA CON,PI,ZERO,ONE,TWO/505379.05D0,3.141592653589793D0,0D0,
     *           1D0,2D0/
      DEGRAD=180D0/PI
      DO 3 K=1,9
         DO 3 I=1,17+NTEE
    3       DD(I,K)=ZERO
      DO 10 I=1,3
         BIN(I)=CON/B(I)
   10    DIN(I)=-BIN(I)/B(I)
      LL=5
      DO 20 K=1,2
         S=SIN(RHO(5,K))
         RHO(2,K)=COS(RHO(6,K))*S
         RHO(3,K)=SIN(RHO(6,K))*S
         RHO(1,K)=COS(RHO(5,K))
         SI=ZERO
         SR=ZERO
         DO 15 I=1,3
            RHO(I,K)=RHO(I,K)*RHO(4,K)
            SQ=RHO(I,K)*RHO(I,K)
            SI=SI+BIN(I)*BIN(I)*SQ
   15       SR=SR+BIN(I)*SQ
         BITAU(K)=DSQRT(SI)
         RES(K)=BITAU(K)
         IF (RHO(4,K).EQ.ZERO) BITAU(K)=ONE
         RITAU(K)=BITAU(K)-SR
         RES(K+2)=RITAU(K)
         DO 19 I=1,3
            DIRC(I,K)=BIN(I)*RHO(I,K)/BITAU(K)
            DIRC(I+3,K)=DACOS(DIRC(I,K))
            LL=LL+1
            RES(LL)=DIRC(I+3,K)
            DD(K,I+6)=DIRC(I,K)*RHO(I,K)
            DD(K+2,I+6)=DD(K,I+6)-RHO(I,K)*RHO(I,K)
            DD(K,I+(K-1)*3)=BIN(I)*DIRC(I,K)
   19       DD(K+2,I+(K-1)*3)=DD(K,I+(K-1)*3)-2*BIN(I)*RHO(I,K)
            IF (RHO(4,K).EQ.ZERO) THEN
               DIRC(7,K)=ONE
            ELSE
               DIRC(7,K)=DATAN(DIRC(3,K)/DIRC(2,K))
            ENDIF
   20       RES(11+K)=DIRC(7,K)
      SR=ZERO
      OMEG=ZERO
      DO 25 I=1,3
         OMEG=OMEG+DIRC(I,1)*DIRC(I,2)
         S=RHO(I,1)*RHO(I,2)
         SR=SR+BIN(I)*S
         DD(5,I+6)=-S
         DD(5,I)=-RHO(I,2)*BIN(I)
   25    DD(5,I+3)=-RHO(I,1)*BIN(I)
      RITAU(3)=-SR
      RES(5)=RITAU(3)
      RES(14)=DACOS(OMEG)
      DO 110 K=1,2
         DO 105 I=1,3
            DD(2+K*3+I,I+6)=RHO(I,K)/BITAU(K)
            DD(2+K*3+I,I+(K-1)*3)=BIN(I)/BITAU(K)
            DO 105 J=1,9
  105          DD(2+K*3+I,J)=DD(2+K*3+I,J)-DIRC(I,K)*DD(K,J)/BITAU(K)
         SR=ONE-DIRC(1,K)*DIRC(1,K)
         DO 110 J=1,9
  110       DD(11+K,J)=(DIRC(2,K)*DD(5+K*3,J)-DIRC(3,K)*DD(4+K*3,J))/SR
      SR=DSIN(RES(14))
      DO 120 J=1,9
         DO 115 I=1,3
  115       DD(14,J)=DD(14,J)+DIRC(I,2)*DD(5+I,J)+DIRC(I,1)*DD(8+I,J)
  120    DD(14,J)=-DD(14,J)/SR
      DO 125 K=1,2
         DO 125 I=1,3
            SR=DSQRT(DABS(ONE-DIRC(I,K)*DIRC(I,K)))
            DO 125 J=1,9
  125          DD(2+K*3+I,J)=-DD(2+K*3+I,J)/SR
      DO 130 K=1,17
         DO 130 I=1,3
  130       DD(K,I+6)=DD(K,I+6)*DIN(I)
      D=RITAU(1)*RITAU(2)-RITAU(3)*RITAU(3)
      F(1)=CON*RITAU(2)/D
      F(2)=CON*RITAU(1)/D
      F(3)=-CON*RITAU(3)/D
      DO 131 I=1,3
  131    RES(14+I)=F(I)
      U(1,1)=-F(1)*F(1)
      U(1,2)=-F(3)*F(3)
      U(2,1)=U(1,2)
      U(2,2)=-F(2)*F(2)
      U(3,1)=-F(1)*F(3)
      U(1,3)=TWO*U(3,1)
      U(3,2)=-F(2)*F(3)
      U(2,3)=TWO*U(3,2)
      U(3,3)=-F(1)*F(2)-F(3)*F(3)
      DO 30 I=1,3
         DO 30 K=1,3
   30       U(K,I)=U(K,I)/CON
      DO 40 I=1,3
         DO 40 J=1,9
            S=ZERO
            DO 35 K=1,3
   35          S=S+U(I,K)*DD(K+2,J)
   40       DD(I+14,J)=S
      DO 50 K=1,2
         IF (RHO(4,K).EQ.ZERO) THEN
            RES(14+K)=ZERO
            RES(2+K)=ZERO
         ENDIF
         DO 45 I=1,3
            IF (RHO(4,K).EQ.ZERO) THEN
               U(I,1)=ZERO
            ELSE
               U(I,1)=RHO(I,K)/RHO(4,K)
            ENDIF
   45    CONTINUE
         SR=ONE
         IF (K.EQ.2) SR=DBLE(NC)
         U(2,2)=RHO(1,K)*COS(RHO(6,K))*SR
         U(3,2)=RHO(1,K)*SIN(RHO(6,K))*SR
         U(1,2)=-RHO(4,K)*SIN(RHO(5,K))*SR
         U(2,3)=-RHO(3,K)
         U(3,3)=RHO(2,K)
         U(1,3)=ZERO
         DO 48 J=1,17
            DO 47 I=1,3
               S=ZERO
               DO 46 L=1,3
   46             S=S+DD(J,L+3*(K-1))*U(L,I)
   47          DIN(I)=S
            DO 48 I=1,3
   48          DD(J,I+3*(K-1))=DIN(I)
   50    CONTINUE
      DO 54 J=1,17
         S=DD(J,2)
         DD(J,2)=DD(J,4)
         DD(J,4)=DD(J,5)
         DD(J,5)=DD(J,3)
   54    DD(J,3)=S
      PIN1=TWO*PI/DBLE(N1)
      PIN2=TWO*PI/DBLE(N2)
      L=17
      K2=9
      DO 60 IV=1,NIV
      DO 59 I=1,N1/2+1
         PIN1I=PIN1*DBLE(I-1)
         J2=N2
         IF (I.EQ.1.OR.2*(I-1).EQ.N1) J2=N2/2+1
         J1=1
         IF (MQ.EQ.1) J1=I
         DO 59 J=J1,J2
         IF (I*J.EQ.1) GOTO 59
            PIN2J=PIN2*DBLE(J-1)
            L=L+1
            RES(L)=ZERO
            K3=K2
            DO 55 K=1,NTE(IV)
               S1=TWO*(DCOS(PIN1I*IQ1(K)+PIN2J*IQ2(K))-ONE)
               DD(L,K+K3)=S1
               S2=TWO*(DCOS(PIN1I*IQ2(K)+PIN2J*IQ1(K))-ONE)
               IF (MQ.EQ.1.AND.IABS(IQ2(K)).NE.IQ1(K)) THEN
                  S1=S1+S2
                  DD(L,K+K3)=DD(L,K+K3)+S2
               ENDIF
               IF (IVX(IV,K).EQ.0) K3=K3-1
   55          RES(L)=RES(L)+EPS(IV,K)*S1
   59       CONTINUE
   60    K2=K3+NTE(IV)
      I2=(MQ+1)*3+NTEE
      IF (MQ.EQ.1) THEN
         DO 70 I=1,L
         DD(I,1)=DD(I,1)+DD(I,2)
         DD(I,2)=DD(I,3)+DD(I,4)
         DD(I,3)=DD(I,5)+DD(I,6)*SIGN(ONE,RHO(6,1)*RHO(6,2))
         DO 65 K=4,I2
   65       DD(I,K)=DD(I,K+3)
   70    CONTINUE
      ENDIF
      DO 100 I=1,L
         DO 95 K=1,I2
            S=ZERO
            DO 90 J=1,I2
   90          S=S+DD(I,J)*C(J,K)
   95       AB(K)=S
         DO 98 K=1,L
            S=ZERO
            DO 97 J=1,I2
   97          S=S+AB(J)*DD(K,J)
   98       AA(I,K)=S
C                   WRITE (6,*) (AA(I,K),K=1,I)
         IF (AA(I,I).LT.ZERO) WRITE (6,*) 'negative AA(I,I)   I =',I
  100    STER(I)=DSQRT(DABS(AA(I,I)))
      RETURN
      END
      SUBROUTINE INDMAT(NU,DC,JDM,IS1,IS2,PIN1,PIN2,NTE,INPAR,A,MQ)
      IMPLICIT REAL(8) (A-H,O-Z)
c     initialization of D matrix, phase angles and energy matrix EV
      DIMENSION INPAR(1),A(1),DC(1)
      COMMON B(2,243,243),D1(243,243),D2(243,243),PHI1(485),PHI2(485),
     *       EV(243,243)
      DATA ZERO,ONE,TWO,SQRT2/0D0,1D0,2D0,1.414213562373095D0/
      NL=(NU+1)/2
      DO 10 K=1,NU
         DO 10 L=1,NU
            EV(L,K)=ZERO
            D1(L,K)=ZERO
   10       D2(L,K)=ZERO
      D1(NL,NL)=ONE
      D2(NL,NL)=ONE
      DC(1)=DCOS(A(3))
      DC(2)=(ONE+DC(1))/TWO
      DC(3)=(ONE-DC(1))/TWO
      DC(4)=DSIN(A(3))/SQRT2
      DC(5)=DCOS(A(4))
      DC(6)=(ONE+DC(5))/TWO
      DC(7)=(ONE-DC(5))/TWO
      DC(8)=DSIN(A(4))/SQRT2
      JDM=0
      SS1=IS1*2
      SS2=IS2*2
      DO 20 K=-NU+1,NU-1
         PHI1(K+NU)=(SS1-A(1)*DBLE(K))*PIN1
   20    PHI2(K+NU)=(SS2-A(2)*DBLE(K))*PIN2
      DO 30 I=1,NTE
         IQ=INPAR(I)/256
         Q2=MOD(IQ/64,16)-8
         Q1=IQ/1024-8
         DO 30 K2=-NL+1,NL-1
            DO 30 K1=-NL+1,NL-1
   30          EV(K1+NL,K2+NL)=EV(K1+NL,K2+NL)+A(21+I)*
     *                   TRIG(2*K1,2*K2,+1,Q1,Q2,1D0,PHI1,PHI2,-1,NU,MQ)
      RETURN
      END
      SUBROUTINE ITER(A,MQ,N1,N2,NC,INPAR,IVR,NTUP,NTE,NTRA,FRQ,BL,WT,
c: 5/20/13
c:     *                ITRA,ILEV,NIT,LP,NIV,IFPR,ISCD)
     *                ITRA,ILEV,NIT,LP,NIV,IFPR,ISCD,SCP)
ce 5/20/13
c  solves the non-linear least-squares problem of determining spectroscopic
c  parameters from the observed transition frequencies by iteration
      IMPLICIT REAL(8) (A-H,O-Z)
c: 5/20/13
c:      DIMENSION A(50,1),FRQ(1),ITRA(8,1),WT(1),INPAR(34,1),IVR(50,1),
c:     *    ILEV(2,1),NTUP(1),NTE(1),DC(8),BL(1),MSYM(10,3)
      DIMENSION A(60,1),FRQ(1),ITRA(8,1),WT(1),INPAR(38,1),IVR(60,1),
     *    ILEV(2,1),NTUP(1),NTE(1),DC(8),BL(1),MSYM(10,3),SCP(1)
ce 5/20/13
      CHARACTER*8 TT(21)
      CHARACTER*1 LBL(241)
      COMMON B(2,243,243),D1(243,243),D2(243,243),PHI1(485),PHI2(485),
     *     EV(243,243),U(2,243,243),E(243),EJ(243),H(2,243,243),
     *     EW(243,243),EX(243,243),EZ(972),DER(8191,80),CALC(8191),
     *     AUX(8191),STER(80),CORR(160),JCR(80)
      DATA NTRX,MS1,MS2,MS4/8191,8,128,16384/
      DATA NU,ZERO,ONE,PI/243,0D0,1D0,3.141592653589793D0/
      DATA TT/'RHO1','RHO2','BETA1','BETA2','ALPHA1','ALPHA2','A','B',
     *'C','DELTA J','DELTA JK','DELTA K','DDELTA J','DDELTA K','PHI J',
     *'PHI JK','PHI KJ','PHI K','PPHI J','PPHI JK','PPHI K'/
      DATA MSYM/
     *    3,2,4,4,2,3,0,3,1,0 ,3,2,1,3,0,0,0,3,0,0, 0,0,2,1,0,0,0,0,1,0/
  901 FORMAT(/20X,14('*')/20X,'*  CYCLE',I3,'  *'/20X,14('*')/)
  902 FORMAT(I3,6D14.6/(3X,6D14.6))
  903 FORMAT(I4,I2,I1,2(I3,2I4),F13.4,F6.2,F8.2,F14.4,2F10.4)
  904 FORMAT(I5,' NORMAL EQUATIONS'//' STATE',7X,'OLD PARAMETER',16X,
c: 5/20/13
c:     *       'STANDARD ERROR',6X,'CHANGE',9X,'PREC'/)
c:  905 FORMAT(I8,A12,D20.12,3D15.5)
c:  906 FORMAT(8X,A12,D20.12)
c:  907 FORMAT(2I4,6I2,D20.12,3D15.5)
c:  908 FORMAT(I4,4X,6I2,D20.12)
c:  909 FORMAT(2I4,A12,D20.12,3D15.5)
     *       'STANDARD ERROR',6X,'CHANGE',9X,'PREC',7X,'SCALE FAC'/)
  905 FORMAT(I8,A12,D20.12,3D15.5,D12.3)
  906 FORMAT(8X,A12,D20.12)
  907 FORMAT(2I4,6I2,D20.12,3D15.5,D12.3)
  908 FORMAT(I4,4X,6I2,D20.12)
  909 FORMAT(2I4,A12,D20.12,3D15.5,D12.3)
ce 5/20/13
  910 FORMAT(I4,4X,A12,D20.12)
  911 FORMAT(/' CORRELATION MATRIX')
  912 FORMAT(/' PAR',12I7)
  913 FORMAT(I4,2X,12I7)
      PIN1=PI/DBLE(N1)
      PIN2=PI/DBLE(N2)
      MS3=MS1*MS1*MS2
      NTR=NTRA+NTRA+1
      DO 200 IT=1,NIT
         WRITE (*,*) 'CYCLE',IT
         WRITE (6,901) IT
         DO 10 I=1,NTRA
            CALC(I)=ZERO
            DO 10 K=1,80
   10          DER(I,K)=ZERO
         M1=1
         IQ=ILEV(1,1)
         IX1=-1
         DO 50 I=1,NTR
            IN=ILEV(1,I)
            IF (IN.EQ.IQ) GOTO 50
            IX=IQ/MS2
            J=MOD(IQ,MS2)
            IV=MOD(IX,MS1)
            IS2=MOD(IX/MS1,MS1)
            IS1=IQ/MS3
            IF (IS1.EQ.IS2) THEN
               ISZ=2
               IF (IS1.EQ.0) ISZ=1
            ELSE
               ISZ=4
               IF (MOD(IS1+IS2,N1).EQ.0) ISZ=3
            ENDIF
            JSYM=0
            IF (ISZ.LT.4) JSYM=MSYM(ISCD+6,ISZ)
            IF (NC.EQ.+1.OR.JSYM.EQ.0.OR.IABS(ISCD).LT.2.OR.IABS(ISCD)
     *          .GT.3.OR.ISCD.EQ.-3.OR.JSYM.EQ.4) GOTO 14
            JSYM=4-JSYM
   14       IF (IX.EQ.IX1) GOTO 15
            CALL INDMAT(NU,DC,JDM,IS1,IS2,PIN1,PIN2,NTE(IV),INPAR(1,IV),
     *                                                       A(1,IV),MQ)
            IX1=IX
            WRITE (*,*) ' IS1',IS1,' IS2',IS2,' IV',IV
   15       CALL HAMILT(A(1,IV),NC,J,J1,J21,RJJ1,NTUP(IV),NTE(IV),
     *                              INPAR(1,IV),DC,JDM,JSYM,LBL,MQ,ISCD)
            DO 20 M=M1,M2
               L=ILEV(2,M)
               NQ=L/MS4
               L=MOD(L,MS4)-MS4/2
               IL=IABS(L)
   20          CALC(IL)=CALC(IL)+E(NQ)*DBLE(L/IL)
            KP1=-1
            DO 30 IP=1,21+NTUP(IV)
               IF (IVR(IP,IV).EQ.0) GOTO 30
               KP=IVR(IP,IV)
c: 5/20/13
c:               S=ONE
c:               IF (IP.EQ.6.AND.KP.EQ.KP1) S=SIGN(ONE,A(5,1)*A(6,1))
               S=SCP(KP)
               IF (IP.EQ.6.AND.KP.EQ.KP1) S=SIGN(S,A(5,1)*A(6,1))
ce 5/20/13
               CALL DERIV(A(1,IV),J,J1,J21,RJJ1,NTUP(IV),INPAR(1,IV),IP,
     *                    PIN1,PIN2,NC,NTE(IV),MQ,ISCD)
               DO 25 M=M1,M2
                  L=ILEV(2,M)
                  NQ=L/MS4
                  L=MOD(L,MS4)-MS4/2
                  IL=IABS(L)
   25             DER(IL,KP)=DER(IL,KP)+E(NQ)*DBLE(L/IL)*S
               KP1=KP
   30          CONTINUE
            M1=I
            IQ=IN
   50       M2=I
         IF (MOD(IFPR,2).EQ.0) GOTO 60
         DO 55 I=1,NTRA
   55       WRITE (6,902) I,(DER(I,K),K=1,LP)
   60    NT=0
         IBL=1
         DO 90 I=1,NTRA
            CALC(I)=FRQ(I)-CALC(I)
            AUX(I)=CALC(I)*WT(I)
            IF (WT(I)*BL(I).GE.ZERO) THEN
		     WRITE (6,903) I,(ITRA(K,I),K=1,8),FRQ(I),BL(I),WT(I),
     *                       CALC(I),AUX(I)
               IF (WT(I).GT.ZERO.AND.BL(I).EQ.ZERO) THEN
                  NT=NT+1
                  CALC(NT)=AUX(I)
                  DO 80 K=1,LP
   80                DER(NT,K)=DER(I,K)*WT(I)
               ELSE
                  IF (WT(I)*BL(I).GT.ZERO) GOTO 90
               ENDIF
            ELSE
               NT=NT+1
               S=ZERO
               DO 85 IB=IBL,I
   85             S=S+AUX(IB)*DABS(BL(IB))
               CALC(NT)=S
               WRITE (6,903) I,(ITRA(K,I),K=1,8),FRQ(I),BL(I),WT(I),
     *                       CALC(I),AUX(I),S
               DO 87 K=1,LP
                  S=ZERO
                  DO 86 IB=IBL,I
   86                S=S+DER(IB,K)*WT(IB)*DABS(BL(IB))
   87             DER(NT,K)=S
            ENDIF
            IBL=I+1
   90       CONTINUE
         CALL LEASQU1(0,DER,NTRX,CALC,AUX,LP,NT,CORR,STER,S,SSQ,AUX,DER,
     *                NTRX)
c: 5/20/13
         DO 110 K=1,LP
            CORR(K)=CORR(K)*SCP(K)
            CORR(K+LP)=CORR(K+LP)*SCP(K)
  110       STER(K)=STER(K)*SCP(K)
ce 5/20/13
         WRITE (6,904) NT
         WRITE (*,*) 'STANDARD DEVIATION', DSQRT(S)
         S1=SIGN(ONE,A(5,1)*A(6,1))
         L1=-1
         DO 140 I=1,6
            S=ONE
            IF (I.GT.2) S=180D0/PI
            L=IVR(I,1)
            IF (L.EQ.0) THEN
               WRITE (6,906) TT(I),A(I,1)*S
            ELSE
               IF (L1.EQ.L.AND.I.EQ.6) CORR(L)=CORR(L)*S1
c: 5/20/13
c:         WRITE (6,905) L,TT(I),A(I,1)*S,STER(L)*S,CORR(L)*S,CORR(L+LP)*S
           WRITE (6,905) L,TT(I),A(I,1)*S,STER(L)*S,CORR(L)*S,
     1		 CORR(L+LP)*S,SCP(L)
ce 5/20/13
               A(I,1)=A(I,1)+CORR(L)
            ENDIF
  140       L1=L
         IF (MQ.EQ.1) A(2,1)=A(1,1)
         IF (MQ.EQ.1) A(4,1)=A(3,1)
         IF (MQ.EQ.1) A(6,1)=A(5,1)*S1
         DO 180 IV=1,NIV
            WRITE (6,*)
            A(1,IV)=A(1,1)
            A(2,IV)=A(2,1)
            A(3,IV)=A(3,1)
            A(4,IV)=A(4,1)
            A(5,IV)=A(5,1)
            A(6,IV)=A(6,1)
            DO 170 I=7,21+NTUP(IV)
               L=IVR(I,IV)
               IF (I.LE.21) GOTO 150
               M=I-21
               IQ=INPAR(M,IV)/256
               KAP=MOD(IQ,16)
               MEG=MOD(IQ/16,4)
               IMG=2-MOD(MEG,2)
               IMGS=2*MOD(MEG,2)-1
               MEG=MEG+IMG-3
               IQ2=MOD(IQ/64,16)-8
               IQ1=MOD(IQ/1024,16)-8
               JP=MOD(INPAR(M,IV)/16,16)
               KP=MOD(INPAR(M,IV),16)
               IF (L.GT.0) THEN
                  WRITE (6,907) IV,L,IQ1,IQ2,MEG,IMGS*KAP,JP,KP,A(I,IV),
c: 5/20/13
c:     *                          STER(L),CORR(L),CORR(L+LP)
     *                          STER(L),CORR(L),CORR(L+LP),SCP(L)
ce 5/20/13
                  A(I,IV)=A(I,IV)+CORR(L)
               ELSE
                  WRITE (6,908) IV,IQ1,IQ2,MEG,IMGS*KAP,JP,KP,A(I,IV)
               ENDIF
               GOTO 170
  150          IF (L.GT.0) THEN
                  WRITE (6,909) IV,L,TT(I),A(I,IV),STER(L),CORR(L),
c: 5/20/13
c:     *                          CORR(L+LP)
     *                          CORR(L+LP),SCP(L)
ce 5/20/13
                  A(I,IV)=A(I,IV)+CORR(L)
               ELSE
                  WRITE (6,910) IV,TT(I),A(I,IV)
               ENDIF
  170          CONTINUE
  180       CONTINUE
         WRITE (6,911)
         DO 190 K=1,LP,12
            J1=MIN0(LP,K+11)
            WRITE (6,912) (J,J=K,J1)
            WRITE (6,913)
            DO 190 I=K,LP
               J2=MIN0(J1,I)
               DO 185 J=K,J2
  185             JCR(J)=DER(I,J)*1D5+.5D0
  190          WRITE (6,913) I,(JCR(J),J=K,J2)
  200    CONTINUE
      DO 230 K=1,LP
         DO 220 L=K,LP
c: 5/20/13
c:            EW(K,L)=DER(K,L)
c:  220       EW(L,K)=DER(K,L)
            EW(K,L)=DER(K,L)*SCP(K)*SCP(L)
  220       EW(L,K)=EW(K,L)
ce 5/20/13
  230    EW(K,K)=STER(K)*STER(K)
      RETURN
      END
      SUBROUTINE HAMILT(A,NC,J,J1,J21,RJJ1,NTUP,NTE,INPAR,DC,JDM,JSYM,
     *                                                      LBL,MQ,ISCD)
C  sets up effective rotational Hamiltonian for a single (vibrational) state of
C  of a molecule with two equivalent internal rotors (periodicity N)
C     array A contains the regular asymmetric rigid rotor constants,the internal
C             rotation parameters (RHO, BETA, ALPHA) and the tunneling parameters
C     J = rotational quantum number
C     NTUP = number of tunneling parameters and higher order distortion constants
C     array INPAR contains the information about the tunneling parameters
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(1),INPAR(1),DC(8),MCODE(10,2)
      CHARACTER*1 LBL(1),LK
c new
      LOGICAL LJ,LA,LB,LC
c new
      COMMON B(2,243,243),D1(243,243),D2(243,243),PHI1(485),PHI2(485),
     *      EV(243,243),U(2,243,243),E(243),EJ(243),H(2,243,243),
     *      EW(243,243),EC(243),AUX(2673),SK(121)
      COMPLEX(8) SC,QC,HC(243,243),UC(243,243),DAL1,DAL2,EC,BC(243,243)
      EQUIVALENCE (H,HC),(U,UC),(B,BC)
cccc      DATA NU,ZERO,ONE,TWO/243,0D0,1D0,2D0/
      DATA NU,ZERO,ONE,TWO,PI/243,0D0,1D0,2D0,3.141592653589793D0/
c new
      DATA SQRT2/1.414213562373095D0/
      DATA MCODE/2,2,3,1,1,3,3,3,1,2, 0,0,4,1,1,4,4,1,3,0/
c new
      LA=NC.EQ.1
      LB=NC.EQ.-1.AND.(A(5).EQ.ZERO.OR.A(5).EQ.PI)
      LC=NC.EQ.-1.AND.(A(5).NE.ZERO.AND.A(5).NE.PI)
c new
      NL=(NU+1)/2
      J1=J+1
      J21=J+J1
      RJJ1=J*J1
      DO 5 KQ=1,J21
         DO 5 K=1,KQ
    5       HC(KQ,K)=(0D0,0D0)
      DO 10 K=-J,J
   10    EJ(K+NL)=DSQRT(DBLE((J-K)*(J+K+1)))
C     set up Asymmetric Rigid Rotor part of Hamiltonian
      CALL ASYMRO(A(7),J,J1,J21,RJJ1,H,NU,EJ,NL)
C     add the tunneling contributions to the matrix elements
C                               first the pure vibrational contributions
      MSC=MCODE(ISCD+6,(NC+3)/2)
      MSCM=MOD(MSC,2)
      N1=NL-J
      N2=NL+J
      IF (JDM.EQ.J) GOTO 20
      DO 15 JJ=JDM,J-1
         CALL DMAT(JJ,JJ+1,NL,DC(1),D1,D1,0)
   15    CALL DMAT(JJ,JJ+1,NL,DC(5),D2,D1,MQ)
      JDM=J
   20 CALL BMAT(J,NL,A(6)-A(5),MSC)
      DAL1=A(5)*(0D0,1D0)
      DAL2=-A(6)*(0D0,1D0)
      N11=N1-1
      DO 50 K=N1,N2
         DO 40 K1=N1,N2
            SC=(0D0,0D0)
            DO 30 K2=N1,N2
   30          SC=SC+BC(K1,K2)*EV(K1,NC*(K2-NL)+NL)*D2(K2,K)
   40       EC(K1)=SC
         F2=1-2*MOD(IABS(NL-K),2)*MSCM
         DO 50 KQ=N1,N2
            SC=(0D0,0D0)
            DO 45 K1=N1,N2
   45          SC=SC+D1(K1,KQ)*EC(K1)
            SC=EXP(DAL1*DBLE(KQ-NL))*SC*F2*EXP(DAL2*DBLE(K-NL))
            IF (KQ.EQ.K) HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+SC
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
   50    CONTINUE
C                    then the contributions to the rotational parameters
      M1=NTE+1
      IQ=INPAR(NTE+1)/256
C     loop over tunneling parameters and higher order distortion constants
      DO 100 I=NTE+1,NTUP+1
         INP1=INPAR(I)/256
         IF (INP1.EQ.IQ) GOTO 100
C     get Q1, Q2, KAP, MEG for each group of parameters
         KAP=MOD(IQ,16)
         MEG=MOD(IQ/16,4)
         IMG=2-MOD(MEG,2)
         MEG=MEG+IMG-3
         Q2=MOD(IQ/64,16)-8
         Q1=MOD(IQ/1024,16)-8
         S2=MEG*(1-2*MOD(KAP,2))
         IF (IABS(ISCD).EQ.4) S2=MEG
         IF (ISCD.EQ.2.AND.IMG.EQ.2.AND.NC.EQ.-1) S2=-S2
         IF (ISCD.EQ.3.AND.IMG.EQ.2.AND.NC.EQ.+1) S2=-S2
         DO 90 K=-J,J-KAP
            KQ=K+KAP
            KK=K+KQ
c: 10/12/09
            SKK=K
            SKQ=KQ
ce 10/12/09
C     M1,M2 are the boundaries of the group of tunneling parameters with
C     the same Q1, Q2, KAP, MEG
            S=ZERO
            DO 60 M=M1,M2
               JP=MOD(INPAR(M)/16,16)
               KP=MOD(INPAR(M),16)
c: 10/12/09
c:   60          S=S+A(21+M)*RJJ1**(JP/2)*(K**KP+KQ**KP)
   60          S=S+A(21+M)*RJJ1**(JP/2)*(SKK**KP+SKQ**KP)
ce 10/12/09
            DO 70 K1=0,KAP-1
   70          S=S*EJ(K+K1+NL)
   90       H(IMG,KQ+J1,K+J1)=H(IMG,KQ+J1,K+J1)+
     1              S*TRIG(KK,NC*KK,MEG,Q1,Q2,S2,PHI1,PHI2,-1,NU,MQ)/TWO
         M1=I
         IQ=INP1
  100    M2=I
      DO 101 I=1,J21
  101    LBL(I)=' '
C     symmetrize blocks (for N=3, blocks 00 and 11 (NC=-1) or 12 (NC=+1)
      IF (JSYM.EQ.0) THEN
c     diagonalize full matrix
         S=ZERO
         CALL SEIGCX(H,J21,NU,E,1,J21,S,J21,U,AUX)
      ELSE
c: 5/20/13
         LL1=(J+MOD(J,2))/2
         LL2=J1*(1-MOD(J,2))
         LL3=LL2+LL1
         IF (JSYM.EQ.4) THEN
            LL2=(LL1+1)*(1-MOD(J,2))
            LL3=LL1
            LJ=MOD(J,2).EQ.0
            IF ((LB.AND.LJ).OR.(LC.AND..NOT.LJ)) LL2=LL2+LL1
            IF ((LC.AND.LJ).OR.(LA.AND..NOT.LJ)) LL3=2*LL3
            LL3=LL2+LL3
         ENDIF
ce 5/20/13
      SELECT CASE(JSYM)
      CASE (1:2,4)
c     symmetrize: 1,4 -> (e+,o+)(e-,o-); 2 -> (e+,o-)(e-,o+)
         J2=2*J1
         DO 108 K=1,J
  108       SK(K)=ONE
         IF (JSYM.EQ.2) THEN
            FJ=1-2*MOD(J,2)
            DO 109 K=1,J
               SK(K)=FJ
  109          FJ=-FJ
cc: 5/20/13
c         ELSE
c            IF (MOD(J,2).EQ.1) THEN
c               LL1=J1
c               LL2=0
c            ENDIF
cce 5/20/13
	   ENDIF
         DO 110 KQ=1,J
            DO 105 K=1,KQ
               SC=HC(KQ,K)+CONJG(HC(J2-K,J2-KQ))*SK(K)*SK(KQ)
               QC=HC(J2-KQ,K)*SK(KQ)+CONJG(HC(J2-K,KQ))*SK(K)
               HC(KQ,K)=(SC+QC)/TWO
  105          HC(J2-K,J2-KQ)=CONJG(SC-QC)/TWO
            SC=HC(J1,KQ)+CONJG(HC(J2-KQ,J1))*SK(KQ)
  110       HC(J1,KQ)=SC/SQRT2
         IF (JSYM.NE.4) THEN
            S=ZERO
            CALL SEIGCX(H,J1,NU,E,1,J1,S,J1,U,AUX)
         ELSE
            IF (J1.GT.0) CALL EVENODD(HC,NU,J1,UC,E,AUX)
         ENDIF
         DO 115 I=1,J1
            DO 115 K=1,J
               UC(K,I)=UC(K,I)/SQRT2
  115          UC(J2-K,I)=UC(K,I)*SK(K)
         DO 120 KQ=1,J
            DO 120 K=1,KQ
  120          HC(KQ,K)=HC(J1+KQ,J1+K)
         IF (JSYM.NE.4) THEN
            S=ZERO
            CALL SEIGCX(H,J,NU,E(J1+1),1,J,S,J,U(1,1,J1+1),AUX)
         ELSE
            IF (J.GT.0) CALL EVENODD(HC,NU,J,UC(1,J1+1),E(J1+1),AUX)
         ENDIF
         DO 130 I=J1+1,J21
            DO 125 K=1,J
  125          UC(J1+K,I)=-UC(K,I)/SQRT2*SK(J1-K)
            DO 126 K=1,J
  126          UC(K,I)=-UC(J2-K,I)*SK(K)
  130       UC(J1,I)=ZERO
      CASE (3)
c     symmetrize even/odd Ka
         CALL EVENODD(HC,NU,J21,UC,E,AUX)
      END SELECT
c: 5/20/13
c:         DO 135 I=1,J
c:  135       LBL(I+J1)='-'
         DO 135 I=1,LL1
            LBL(I+LL3)='-'
  135       LBL(I+LL2)='-'
ce 5/20/13
         DO 150 I=1,J21
            KQ=I
            S=E(KQ)
            DO 145 K=I,J21
               IF (E(K).LT.S) KQ=K
  145          S=E(KQ)
            E(KQ)=E(I)
            E(I)=S
            LK=LBL(KQ)
            LBL(KQ)=LBL(I)
            LBL(I)=LK
            DO 150 K=1,J21
               SC=UC(K,KQ)
               UC(K,KQ)=UC(K,I)
  150          UC(K,I)=SC
      ENDIF
      RETURN
      END
      SUBROUTINE EVENODD(HC,NH,M,UC,E,AUX)
      IMPLICIT REAL(8) (A-H,O-Z)
      COMPLEX(8) HC(NH,1),UC(NH,1)
      DIMENSION E(1),AUX(1)
      DATA ZERO/0D0/
      M1=(M+1)/2
      M2=MOD(M,2)*M/2
      M3=1-MOD(M,2)
      M5=M-M1
      M4=M/2-M2
      DO 160 K=1,M
         KK=(K-M3)/2
         DO 160 KQ=K,M,2
            KK=KK+1
  160       HC(KK,K)=HC(KQ,K)
      DO 170 K=2-M3,M,2
         KK=(K-M3)/2+1
         DO 170 KQ=KK,M1
            HC(KQ+M2,KK+M2)=HC(KQ,K)
  170       HC(KQ+M4,KK+M4)=HC(KQ,K+1)
      S=ZERO
      CALL SEIGCX(HC,M1,NH,E,1,M1,S,M1,UC,AUX)
      DO 175 I=1,M1
         L=M+M3
         IF (MOD(M,2).EQ.1) UC(L,I)=UC(M1,I)
         DO 175 K=M5,1,-1
            L=L-1
            UC(L,I)=ZERO
            L=L-1
  175       UC(L,I)=UC(K,I)
      DO 180 KQ=1,M5
          DO 180 K=1,KQ
  180        HC(KQ,K)=HC(M1+KQ,M1+K)
      S=ZERO
      CALL SEIGCX(HC,M5,NH,E(M1+1),1,M5,S,M5,UC(1,M1+1),AUX)
      DO 185 I=M1+1,M
         L=M+M3
         IF (MOD(M,2).EQ.1) UC(L,I)=ZERO
         DO 185 K=M5,1,-1
            L=L-1
            UC(L,I)=UC(K,I)
            L=L-1
  185       UC(L,I)=ZERO
      RETURN
      END
      SUBROUTINE DMAT(J,J1,NL,DC,D,DD,M)
      IMPLICIT REAL(8) (A-H,O-Z)
C     calculation of matrix D(J+1) from D(J)
C                    note:   J1 = J+1
C                            where    J = ROT. Q.N. OF LOWER J
      COMMON B(243,243)
      DIMENSION D(243,243),DD(243,243),DC(4)
      F1(L,K)=DSQRT(DBLE((L+K)*(L+K+1)))
      F2(L,K)=DSQRT(DBLE(2*(L-K+1)*(L+K+1)))
      IF (M.NE.1) THEN
         BJ=(J+J1)*(J1+J1)
         DO 30 K=0,J1
            FA=F1(J,-K)
            FB=F2(J,K)
            FC=F1(J,K)
            L=NL+K
            DO 30 KQ=-K,K
               LQ=NL+KQ
               BR=FA*D(LQ+1,L+1)*DC(2)-FB*D(LQ+1,L)*DC(4)+FC*D(LQ+1,L-1)
     *                                                            *DC(3)
               BQ=FA*D(LQ,L+1)*DC(4)  +FB*D(LQ,L)*DC(1)  -FC*D(LQ,L-1)
     *                                                            *DC(4)
               BP=FA*D(LQ-1,L+1)*DC(3)+FB*D(LQ-1,L)*DC(4)+FC*D(LQ-1,L-1)
     *                                                            *DC(2)
   30          B(LQ,L)=(F1(J,-KQ)*BR+F2(J,KQ)*BQ+F1(J,KQ)*BP)
         J2=NL+NL
         DO 40 L=NL,NL+J1
            DO 40 LQ=J2-L,L
               D(LQ,L)=B(LQ,L)/BJ
               D(J2-L,J2-LQ)=D(LQ,L)
               D(L,LQ)=B(LQ,L)*DBLE(1-2*MOD(L-LQ,2))/BJ
   40          D(J2-LQ,J2-L)=D(L,LQ)
      ELSE
   50    DO 55 K=NL-J1,NL+J1
            DO 55 KQ=NL-J1,NL+J1
   55          D(KQ,K)=DD(KQ,K)
      ENDIF
      RETURN
      END
      SUBROUTINE BMAT(J,NL,DELALP,M)
      IMPLICIT REAL(8) (A-H,O-Z)
      COMMON B(2,243,243),D1(243,243),D2(243,243)
      DIMENSION F(243)
      COMPLEX(8) DALP,SC,BC(243,243)
      EQUIVALENCE (B,BC)
      DATA ZERO,TWO/0D0,2D0/
      DALP=(0D0,1D0)*DELALP
      J1=NL-J
      J2=NL+J
      SELECT CASE (M)
      CASE (1,3)
      DO 10 K=J1,J2
   10    F(K)=1-2*MOD(IABS(K-NL),2)
      CASE (2,4)
      DO 15 K=J1,J2
   15    F(K)=1D0
      END SELECT
      SELECT CASE (M)
      CASE (1:2)
c             most common case, B is real  (# 1,2)
      DO 30 K1=J1,J2
         DO 30 K2=J1,J2
            S=ZERO
            DO 25 K=J1,J2
   25          S=S+F(K)*D1(K1,K)*D2(K2,K)
            B(1,K1,K2)=S/TWO
   30       B(2,K1,K2)=ZERO
      CASE (3:4)
c             B is complex  (# 3,4)
      DO 40 K1=J1,J2
         DO 40 K2=J1,J2
            SC=(0D0,0D0)
            DO 35 K=J1,J2
   35          SC=SC+F(K)*D1(K1,K)*D2(K2,K)*EXP(DBLE(K-NL)*DALP)
   40       BC(K1,K2)=SC/TWO
      END SELECT
      RETURN
      END
      FUNCTION TRIG(KP,KM,MEG,Q1,Q2,S2,PHI1,PHI2,NOD,NU,MQ)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION PHI1(1),PHI2(1)
      DATA ZERO,ONE,TWO/0D0,1D0,2D0/
      S=ONE
      IF (Q1.EQ.ZERO.AND.Q2.EQ.ZERO) GOTO 50
      PHI1Q=PHI1(KP+NU)*Q1+PHI2(KM+NU)*Q2
      PHI2Q=PHI1(KP+NU)*Q2+PHI2(KM+NU)*Q1
c: 1/7/09
c:      DER1=KP*Q1
c:      DER2=KP*Q2
      DER1=-KP*Q1
      DER2=-KP*Q2
ce 1/7/09
      FAC=TWO 
c: 1/7/09
c:      IF (MEG.EQ.-1) FAC=DBLE(-NOD)*FAC
      IF (MEG.EQ.-1) FAC=DBLE(NOD)*FAC
ce 1/7/09
      S=ZERO
      IX=2
      IF (MQ.EQ.1) IX=MIN0(1,IDINT((Q1+Q2)*(Q1-Q2)))+2
      IF (NOD.EQ.1) GOTO 5
      DER1=ONE
      DER2=ONE
    5 GOTO (20,10,40,30) IX+MEG*NOD
   10 S=S2*DCOS(PHI2Q)*DER2
   20 S=FAC*(S+DCOS(PHI1Q)*DER1)
      GOTO 50
   30 S=S2*DSIN(PHI2Q)*DER2
c: 1/7/09
c:   40 S=-FAC*(S-DSIN(PHI1Q)*DER1)
   40 S=-FAC*(S+DSIN(PHI1Q)*DER1)
ce 1/7/09
   50 TRIG=S
      RETURN
      END
      SUBROUTINE DERIV(A,J,J1,J21,RJJ1,NTUP,INPAR,IP,PIN1,PIN2,NC,NTE,
     *                  MQ,ISCD)
C     calculates the derivatives of the energy levels with respect to the
C                    spectroscopic parameters
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(1),INPAR(1),MCODE(10,2),EW(243,243)
      COMMON B(2,243,243),D1(243,243),D2(243,243),PHI1(485),PHI2(485),
     *     EV(243,243),U(2,243,243),E(243),EJ(243),H(2,243,243),
     *     EX(243,243),EC(243),B1C(243)
      COMPLEX(8) UC(243,243),HC(243,243),SC,BC(243,243),EC,B1C,S1C,DAL1,
     *     DAL2,EX,DALP
      EQUIVALENCE (H,HC),(U,UC),(B,BC),(EX,EW)
      DATA MCODE/2,2,3,1,1,3,3,3,1,2,0,0,4,1,1,4,4,1,3,0/
      DATA NU,ZERO,ONE,TWO/243,0D0,1D0,2D0/
      NL=(NU+1)/2
      MSC=MCODE(ISCD+6,(NC+3)/2)
      MSCM=MOD(MSC,2)
      DO 10 KQ=1,J21
         DO 10 K=1,KQ
   10       HC(KQ,K)=(0D0,0D0)
      N1=NL-J
      N2=NL+J
      N11=N1-1
      DAL1=A(5)*(0D0,1D0)
      DAL2=-A(6)*(0D0,1D0)
      SELECT CASE (IP)
c
c     derivatives with respect to the rotational and distortion constants
c
      CASE (7:21)
      DO 30 I=1,15
   30    E(I)=ZERO
      E(IP-6)=ONE
      CALL ASYMRO(E,J,J1,J21,RJJ1,H,NU,EJ,NL)
c
c     derivative with respect to RHO1 or RHO2
c
      CASE (1:2)                                   
         DO 41 K2=-J,J
            DO 41 K1=-J,J
   41          EW(K1+NL,K2+NL)=ZERO
      DO 42 I=1,NTE
         IQ=INPAR(I)/256
         Q2=MOD(IQ/64,16)-8
         Q1=IQ/1024-8
         QA=Q1+(Q1-Q2)*DBLE(1-MQ)
         QB=Q1+Q2-QA
         DO 42 K2=-J,J
            DO 42 K1=-J,J
               IF (IP.EQ.1) EW(K1+NL,K2+NL)=EW(K1+NL,K2+NL)+A(21+I)*
     *               TRIG(2*K1,2*K2,+1,Q1,Q2,1D0,PHI1,PHI2,1,NU,MQ)*PIN1
   42          IF (IP.EQ.MQ) EW(K1+NL,K2+NL)=EW(K1+NL,K2+NL)+A(21+I)*
     *               TRIG(2*K2,2*K1,+1,QA,QB,1D0,PHI2,PHI1,1,NU,MQ)*PIN2
      DO 60 K=N1,N2
         DO 50 K1=N1,N2
            SC=(0D0,0D0)
            DO 45 K2=N1,N2
   45          SC=SC+BC(K1,K2)*EW(K1,NC*(K2-NL)+NL)*D2(K2,K)
   50       EC(K1)=SC
         F2=1-2*MOD(IABS(NL-K),2)*MSCM
         DO 60 KQ=N1,N2
            SC=(0D0,0D0)
            DO 55 K1=N1,N2
   55          SC=SC+D1(K1,KQ)*EC(K1)
            SC=EXP(DAL1*DBLE(KQ-NL))*SC*F2*EXP(DAL2*DBLE(K-NL))
            IF (KQ.EQ.K) HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+SC
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
   60       CONTINUE
      M1=NTE+1
      IQ=INPAR(NTE+1)/256
      DO 100 I=NTE+1,NTUP+1
         INP1=INPAR(I)/256
         IF (INP1.EQ.IQ) GOTO 100
         KAP=MOD(IQ,16)
         MEG=MOD(IQ/16,4)
         IMG=2-MOD(MEG,2)
         MEG=MEG+IMG-3
         Q2=MOD(IQ/64,16)-8
         Q1=MOD(IQ/1024,16)-8
         QA=Q1+(Q1-Q2)*DBLE(1-MQ)
         QB=Q1+Q2-QA
         S2=MEG*(1-2*MOD(KAP,2))
         IF (IABS(ISCD).EQ.4) S2=MEG
         IF (ISCD.EQ.2.AND.IMG.EQ.2.AND.NC.EQ.-1) S2=-S2
         IF (ISCD.EQ.3.AND.IMG.EQ.2.AND.NC.EQ.+1) S2=-S2
         DO 90 K=-J,J-KAP
            KQ=K+KAP
            KK=K+KQ
c: 10/12/09
            SKQ=KQ
            SKK=K
ce 10/12/09
            S=ZERO
            DO 70 M=M1,M2
               JP=MOD(INPAR(M)/16,16)
               KP=MOD(INPAR(M),16)
c: 10/12/09
c:   70          S=S+A(21+M)*RJJ1**(JP/2)*(K**KP+KQ**KP)
   70          S=S+A(21+M)*RJJ1**(JP/2)*(SKK**KP+SKQ**KP)
ce 10/12/09
            DO 80 K1=0,KAP-1
   80          S=S*EJ(K+K1+NL)
   89         IF (IP.EQ.1) H(IMG,KQ+J1,K+J1)=H(IMG,KQ+J1,K+J1)+S*
     1            TRIG(KK,NC*KK,MEG,Q1,Q2,S2,PHI1,PHI2,1,NU,MQ)*PIN1/TWO
   90         IF (IP.EQ.MQ) H(IMG,KQ+J1,K+J1)=H(IMG,KQ+J1,K+J1)+S*
     1            TRIG(NC*KK,KK,MEG,QA,QB,S2,PHI2,PHI1,1,NU,MQ)*PIN2/TWO
         M1=I
         IQ=INP1
  100    M2=I
c
c     derivatives with respect to the tunneling parameters
c
      CASE (22:)                                   
  120 IF (IP.LE.21+NTE) THEN
      IQ=INPAR(IP-21)/256
      Q2=MOD(IQ/64,16)-8
      Q1=IQ/1024-8
      DO 130 K2=-J,J
         DO 130 K1=-J,J
  130    EW(K1+NL,K2+NL)=TRIG(2*K1,2*K2,+1,Q1,Q2,1D0,PHI1,PHI2,-1,NU,MQ)
      DO 160 K=N1,N2
         DO 145 K1=N1,N2
            SC=(0D0,0D0)
            DO 140 K2=N1,N2
  140          SC=SC+BC(K1,K2)*EW(K1,NC*(K2-NL)+NL)*D2(K2,K)
  145       EC(K1)=SC
         F2=1-2*MOD(IABS(NL-K),2)*MSCM
         DO 160 KQ=N1,N2
            SC=(0D0,0D0)
            DO 150 K1=N1,N2
  150          SC=SC+D1(K1,KQ)*EC(K1)
            SC=EXP(DAL1*DBLE(KQ-NL))*SC*F2*EXP(DAL2*DBLE(K-NL))
            IF (KQ.EQ.K) HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+SC
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
  160       CONTINUE
      ELSE
  170 IQ=INPAR(IP-21)/256
      KAP=MOD(IQ,16)
         MEG=MOD(IQ/16,4)
         IMG=2-MOD(MEG,2)
         MEG=MEG+IMG-3
      Q2=MOD(IQ/64,16)-8
      Q1=MOD(IQ/1024,16)-8
      S2=MEG*(1-2*MOD(KAP,2))
         IF (IABS(ISCD).EQ.4) S2=MEG
         IF (ISCD.EQ.2.AND.IMG.EQ.2.AND.NC.EQ.-1) S2=-S2
         IF (ISCD.EQ.3.AND.IMG.EQ.2.AND.NC.EQ.+1) S2=-S2
      DO 190 K=-J,J-KAP
         KQ=K+KAP
         KK=K+KQ
c: 10/12/09
            SKQ=KQ
            SKK=K
ce 10/12/09
         JP=MOD(INPAR(IP-21)/16,16)
         KP=MOD(INPAR(IP-21),16)
c: 10/12/09
c:         S=RJJ1**(JP/2)*(K**KP+KQ**KP)
         S=RJJ1**(JP/2)*(SKK**KP+SKQ**KP)
ce 10/12/09
         DO 180 K1=0,KAP-1
  180          S=S*EJ(K+K1+NL)
  190    H(IMG,KQ+J1,K+J1)=S*TRIG(KK,NC*KK,MEG,Q1,Q2,S2,PHI1,PHI2,-1,NU,
     *                        MQ)/TWO
      ENDIF
c
c     derivative with respect to BETA1 or BETA2
c
      CASE (3:4)
  200 JN=NL+NL
      IF (IP.EQ.4) GOTO 255
      DO 250 K=N1,N2
         DO 230 K1=N1,N2
            SC=(0D0,0D0)
            S1C=SC
            DO 220 K2=N1,N2
               L2=NC*(K2-NL)+NL
               FF=EV(K1,L2)*D2(K2,K)
               S1C=S1C+(EJ(JN-K1)*BC(K1-1,K2)-EJ(K1)*BC(K1+1,K2))*FF
  220          SC=SC+BC(K1,K2)*FF
            B1C(K1)=S1C
  230       EC(K1)=SC
         F2=1-2*MOD(IABS(NL-K),2)*MSCM
         DO 250 KQ=N1,N2
            SC=(0D0,0D0)
            DO 240 K1=N1,N2
  240          SC=SC-(EJ(K1)*D1(K1+1,KQ)-EJ(JN-K1)*D1(K1-1,KQ))*EC(K1)
     *            +D1(K1,KQ)*B1C(K1)
            SC=EXP(DAL1*DBLE(KQ-NL))*SC*F2*EXP(DAL2*DBLE(K-NL))/TWO
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
  250       IF (KQ.EQ.K) HC(K-N11,K-N11)=HC(K-N11,K-N11)+SC
      IF (MQ.EQ.2) GOTO 400
  255 DO 290 K=N1,N2
         DO 270 K2=N1,N2
            SC=(0D0,0D0)
            S1C=SC
            DO 260 K1=N1,N2
               L2=NC*(K2-NL)+NL
               FF=EV(K1,L2)*D1(K1,K)
               S1C=S1C+(EJ(JN-K2)*BC(K1,K2-1)-EJ(K2)*BC(K1,K2+1))*FF
  260          SC=SC+BC(K1,K2)*FF
            B1C(K2)=S1C
  270       EC(K2)=SC
         DO 290 KQ=N1,N2
            F2=1-2*MOD(IABS(NL-KQ),2)*MSCM
            SC=(0D0,0D0)
            DO 280 K2=N1,N2
  280          SC=SC-(EJ(K2)*D2(K2+1,KQ)-EJ(JN-K2)*D2(K2-1,KQ))*EC(K2)
     *            +D2(K2,KQ)*B1C(K2)
            SC=EXP(DAL2*DBLE(KQ-NL))*SC*F2*EXP(DAL1*DBLE(K-NL))/TWO
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
  290       IF (KQ.EQ.K) HC(K-N11,K-N11)=HC(K-N11,K-N11)+SC
c
c     derivative with respect to ALPHA1 or ALPHA2
c
      CASE (5:6)
         DALP=(0D0,1D0)*(A(6)-A(5))
         IF (MSC.EQ.3) THEN
         DO 305 K=N1,N2
  305       E(K)=(K-NL)*(1-2*MOD(IABS(K-NL),2))
         ELSE
         DO 308 K=N1,N2
  308       E(K)=K-NL
         ENDIF
         DO 315 K1=N1,N2
            DO 315 K2=N1,N2
               SC=(0D0,0D0)
               IF (MSC.GE.3) THEN
                  DO 310 K=N1,N2
  310                SC=SC+E(K)*D1(K1,K)*D2(K2,K)*EXP(DBLE(K-NL)*DALP)
               ENDIF
  315          EX(K1,K2)=(0D0,1D0)*SC/TWO
      S=ONE
      IF (IP.EQ.6) GOTO 355
      DO 350 K=N1,N2
         DO 330 K1=N1,N2
            SC=(0D0,0D0)
            S1C=SC
            DO 320 K2=N1,N2
               L2=NC*(K2-NL)+NL
               FF=EV(K1,L2)*D2(K2,K)
               S1C=S1C+BC(K1,K2)*FF
  320          SC=SC+EX(K1,K2)*FF
            B1C(K1)=S1C
  330       EC(K1)=SC
         F2=1-2*MOD(IABS(NL-K),2)*MSCM
         DO 350 KQ=N1,N2
            SC=(0D0,0D0)
            S1C=SC
            DO 340 K1=N1,N2
               SC=SC+D1(K1,KQ)*EC(K1)
  340          S1C=S1C+D1(K1,KQ)*B1C(K1)
            S1C=S1C*(0D0,1D0)*DBLE(KQ-NL)
            SC=EXP(DAL1*DBLE(KQ-NL))*(S1C-SC)*F2*EXP(DAL2*DBLE(K-NL))
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
  350       IF (KQ.EQ.K) HC(K-N11,K-N11)=HC(K-N11,K-N11)+SC
      IF (MQ.EQ.2) GOTO 400
      S=SIGN(ONE,A(6)*A(5))
  355 DO 390 K=N1,N2
         DO 370 K1=N1,N2
            SC=(0D0,0D0)
            S1C=SC
            DO 360 K2=N1,N2
               L2=NC*(K2-NL)+NL
               FF=EV(K1,L2)*D2(K2,K)
               S1C=S1C+BC(K1,K2)*FF
  360          SC=SC+EX(K1,K2)*FF
            B1C(K1)=S1C
  370       EC(K1)=SC
         F2=(1-2*MOD(IABS(NL-K),2)*MSCM)*S
         DO 390 KQ=N1,N2
            SC=(0D0,0D0)
            S1C=SC
            DO 380 K1=N1,N2
               SC=SC+D1(K1,KQ)*EC(K1)
  380          S1C=S1C+D1(K1,KQ)*B1C(K1)
            S1C=S1C*(0D0,-1D0)*DBLE(K-NL)
            SC=EXP(DAL1*DBLE(KQ-NL))*(S1C+SC)*F2*EXP(DAL2*DBLE(K-NL))
            IF (KQ.LE.K) THEN
               HC(K-N11,KQ-N11)=HC(K-N11,KQ-N11)+CONJG(SC)
            ELSE
               HC(KQ-N11,K-N11)=HC(KQ-N11,K-N11)+SC
            ENDIF
  390       IF (KQ.EQ.K) HC(K-N11,K-N11)=HC(K-N11,K-N11)+SC
      END SELECT                                    
C     transform derivatives to eigenbasis
  400 DO 420 I=1,J21
         S=ZERO
         DO 415 K=1,J21
            SC=-CONJG(UC(K,I))*HC(K,K)/TWO
            DO 410 KQ=K,J21
  410          SC=SC+CONJG(UC(KQ,I))*HC(KQ,K)
  415       S=S+REAL(SC*UC(K,I))
  420    E(I)=S*TWO
      RETURN
      END
      SUBROUTINE ASYMRO(B,J,J1,J21,RJJ1,H,NU,EJ,NL)
C     sets up the matrix of the asymmetric rotor in the symmetric rotor basis
C     (a-reduction)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION B(1),H(2,NU,1),EJ(1)
      DATA ONE,TWO,SIX,EM3/1D0,2D0,6D0,1D-3/
      QTM=RJJ1*EM3
      B1=(B(2)+B(3))/TWO
      A1=     RJJ1*(B1+QTM*(-B(4)+QTM*B(9)))
      A2=      B(1)-B1+QTM*(-B(5)+QTM*B(10))
      A3=              EM3*(-B(6)+QTM*B(11))
      A4=              EM3*       EM3*B(12)
      A5=(B1-B(3))/TWO+QTM*(-B(7)+QTM*B(13))
      A6=              EM3*(-B(8)+QTM*B(14))
      A7=              EM3*       EM3*B(15)
      DO 15 K=1,J21
         RK2=(K-J1)*(K-J1)
   15    H(1,K,K)=A1+RK2*(A2+RK2*(A3+RK2*A4))
      DO 20 K=1,J21-2
         RK2=(K-J)*(K-J)
   20    H(1,K+2,K)=EJ(K-J1+NL)*EJ(K-J+NL)*
     *                          (A5+(ONE+RK2)*A6+(ONE+RK2*(SIX+RK2))*A7)
      RETURN
      END
      SUBROUTINE PREDIC(A,MQ,N1,N2,NC,INPAR,NTUP,NTE,JMIN,JMAX,FMIN,
     *     FMAX,NSIG,ISIG,THRES,NIV,TEMP,IVR,LP,NUNC,DIP,IFPR,ISCD)
C     predicts rotational energy levels, transition frequencies and moments
      IMPLICIT REAL(8) (A-H,O-Z)
c: 5/20/13
c:      DIMENSION A(50,1),INPAR(34,1),ISIG(4,10,1),NTUP(1),NTE(1),JMIN(1),
c:     *   JMAX(1),FMIN(1),FMAX(1),NSIG(1),THRES(1),IVR(50,1),DIP(1),DC(8)
c:      DIMENSION MSYM(10,3)
      DIMENSION A(60,1),INPAR(38,1),ISIG(4,10,1),NTUP(1),NTE(1),JMIN(1),
     *   JMAX(1),FMIN(1),FMAX(1),NSIG(1),THRES(1),IVR(60,1),DIP(1),DC(8)
      DIMENSION MSYM(10,3),SST(121,10,2)
ce 5/20/13
      CHARACTER*1 LBL(241)
      COMMON B(2,243,243),D1(243,243),D2(243,243),PHI1(485),PHI2(485),
     *     EV(243,243),U(2,243,243),E(243),EJ(243),H(2,243,243),
     *     EW(243,243),EX(243,243),EZ(972),EE(243),DER(243,80),
     *     DERJ(243,80),VCM(80,80),LBJ(241),LB(241)
      COMPLEX(8) UC(243,243),HC(243,243),S1C,S2C,S3C
      EQUIVALENCE (H,HC),(U,UC)
c: 5/20/13
c:      DATA NU,ZERO,ONE,TWO,PI,PLANCK,BOLTZ,C0/243,0D0,1D0,2D0,
c:     *        3.141592653589793D0,6.6260755D-34,1.380658D-23,29979.2458/
      DATA NU,ZERO,ONE,TWO,PI,PLANCK,BOLTZ,C0,XINTF/243,0D0,1D0,2D0,
     *        3.141592653589793D0,6.6260755D-34,1.380658D-23,29979.2458,
     *        4.16231D-5/
ce 5/20/13
      DATA MSYM/
     *    3,2,4,4,2,3,0,3,1,0 ,3,2,1,3,0,0,0,3,0,0, 0,0,2,1,0,0,0,0,1,0/
  980 FORMAT(' J =',I3/(8(F14.3,A1)))
  981 FORMAT(/' VIBRATIONAL STATE',I3,6X,'IS1',I2,4X,'IS2',I2/)
  982 FORMAT(4I4,F15.4,I4,4F10.5,D15.4)
  983 FORMAT(10F9.5)
      F1(J,K)=DSQRT(DBLE((J-K)*(J+K+1)))/TWO
      F2(J,K)=DSQRT(DBLE((J+K)*(J+K+1)))
      PIN1=PI/DBLE(N1)
      PIN2=PI/DBLE(N2)
      FAC=-PLANCK*1D6/(TEMP*BOLTZ)
      OPEN(1,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(2,STATUS='SCRATCH',FORM='UNFORMATTED')
      IF (NUNC.EQ.0) GOTO 10
      DO 5 K=1,LP
         DO 5 L=1,LP 
    5       VCM(L,K)=EW(L,K)
   10 WRITE (*,*) 'PREDICTIONS'
      WRITE (6,984)
  984 FORMAT(/' PREDICTIONS')
      DO 100 IV=1,NIV
      FIV=FMIN(IV)
      FAV=FMAX(IV)
c: 5/20/13
      EMIN=ZERO
ce 5/20/13
      DO 90 IS=1,NSIG(IV)
         IS1=ISIG(1,IS,IV)
         IS2=ISIG(2,IS,IV)
c: 5/20/13
         ISG1=ISIG(3,IS,IV)
         ISG2=ISIG(4,IS,IV)
c     Initialize sum of states
         SST1=ZERO
         SST2=ZERO
ce:5/20/13
         WRITE (*,*) ' IS1',IS1,' IS2',IS2,' IV',IV
         WRITE (6,981) IV,IS1,IS2
         IF (IS1.EQ.IS2) THEN
            ISZ=2
            IF (IS1.EQ.0) ISZ=1
         ELSE
            ISZ=4
            IF (MOD(IS1+IS2,N1).EQ.0) ISZ=3
         ENDIF
         JSYM=0
         IF (ISZ.LT.4) JSYM=MSYM(ISCD+6,ISZ)
         IF (NC.EQ.+1.OR.JSYM.EQ.0.OR.IABS(ISCD).LT.2.OR.IABS(ISCD).GT.3
     *       .OR.ISCD.EQ.-3.OR.JSYM.EQ.4) GOTO 14
         JSYM=4-JSYM
   14    CALL INDMAT(NU,DC,JDM,IS1,IS2,PIN1,PIN2,NTE(IV),INPAR(1,IV),
     *                                                       A(1,IV),MQ)
c: 5/20/13
c:         DO 90 J=JMIN(IV),JMAX(IV)
         DO 80 J=JMIN(IV),JMAX(IV)
ce 5/20/13
            CALL HAMILT(A(1,IV),NC,J,J1,J21,RJJ1,NTUP(IV),NTE(IV),
     *                              INPAR(1,IV),DC,JDM,JSYM,LBL,MQ,ISCD)
            WRITE (6,980) J,(E(K),LBL(K),K=1,J21)
c: 5/20/13
c:            IF (IFPR.LT.2) GOTO 12
            IF (IFPR.LT.2.OR.IFPR.EQ.4) GOTO 12
ce 5/20/13
c     Printing complex eigenvectors !!!!!!!!!!!!!!!!!!!!!!!!!!
            DO 11 I=1,J21
   11          WRITE (6,983) (U(1,K,I),K=1,J21)
   12       DO 13 I=1,J21
               LB(I)=I
               IF (LBL(I).EQ.'-') LB(I)=-I
   13          CONTINUE
            DO 15 M=1,J21
   15          EE(M)=E(M)
            IF (NUNC.EQ.0) GOTO 19
            DO 18 IP=1,21+NTUP(IV)
               IF (IVR(IP,IV).EQ.0) GOTO 18
               KP=IVR(IP,IV)
               CALL DERIV(A(1,IV),J,J1,J21,RJJ1,NTUP(IV),INPAR(1,IV),IP,
     *                    PIN1,PIN2,NC,NTE(IV),MQ,ISCD)
               DO 16 M=1,J21
   16             DER(M,KP)=E(M)
   18          CONTINUE
c     calculate transition frequencies and transition moments : R-transitions
   19       SJ=4*J
            J22=J21+1
            TH=THRES(IV)
            IF (J.EQ.JMIN(IV)) GOTO 35
            WRITE (6,*) '  R-TRANSITIONS',J,J-1
            READ (1) EJ,LBJ,H,DERJ
            DO 30 I=1,J21-2
c: 5/20/13
c:               ISP=MOD(I/2+(J21-1-I)/2,2)+3
c:               ISP=ISIG(ISP,IS,IV)
c:               SS=ISP
c:               SS=DEXP(FAC*EJ(I))*SS
               ISP=ISG1
c:               IF (LBJ(I).EQ.'-') ISP=ISG2
               IF (LBJ(I).LT.0) ISP=ISG2
               SS=ISP
               SS=DEXP(FAC*EJ(I))*SS*XINTF
ce 5/20/13
               DO 30 M=1,J21
                  S1C=(0D0,0D0)
                  S2C=S1C
                  S3C=S1C
                  DEL=EE(M)-EJ(I)
                  IF (DABS(DEL).LT.FIV.OR.DABS(DEL).GT.FAV) GOTO 30
                  DO 20 KQ=1,J21-2
                     K=KQ-J
                     S1C=S1C+CONJG(UC(KQ,M))*F2(J,-K)*HC(KQ,I)
                     S2C=S2C+CONJG(UC(KQ+2,M))*F2(J,K)*HC(KQ,I)
   20                S3C=S3C+CONJG(UC(KQ+1,M))*DSQRT(DBLE((J-K)*(J+K)))
     *                                                         *HC(KQ,I)
                  S3C=S3C*TWO
                  SM=(REAL(CONJG(S3C)*(S1C-S2C))*DIP(2)
     *               -AIMAG(CONJG(S3C)*(S1C+S2C))*DIP(3))*DIP(1)
     *               +AIMAG(CONJG(S2C)*S1C)*DIP(2)*DIP(3)
                  SC=CONJG(S1C+S2C)*(S1C+S2C)/SJ
                  SB=CONJG(S1C-S2C)*(S1C-S2C)/SJ
                  SA=CONJG(S3C)*S3C/SJ
                  SM=SM/SJ+SA*DIP(1)**2+SB*DIP(2)**2+SC*DIP(3)**2
c: 5/20/13
c:                  IF (SM.GT.TH) THEN
c:                     SBB=SM*SS*(ONE-DEXP(FAC*DEL))
                  SBB=SM*SS*(ONE-DEXP(FAC*DEL))*DEL
                  IF (SBB.GT.TH) THEN
ce 5/20/13
                     UP=DMAX1(EJ(I),EE(M))/C0
                     S=ZERO
                     IF (NUNC.NE.0) THEN
                        DO 25 KP=1,LP
                           E(KP)=DER(M,KP)-DERJ(I,KP)
                           S2=ZERO
                           DO 24 L=1,KP
   24                         S2=S2+VCM(L,KP)*E(L)
   25                      S=S+E(KP)*(S2*TWO-VCM(KP,KP)*E(KP))
                        S=DSQRT(S) 
c         DEL = transition frequency
c         ISP = spin multiplicity
c         S1,S2,S3 = line strengths for mu(a),mu(b),mu(c)
c         SBB = line strength * spin weight * Boltzmann factor at temperature TEMP
c         S  = uncertainty of transition frequency
c         UP = upper state energy (in cm-1)
                        WRITE (6,982) J,LB(M),J-1,LBJ(I),DEL,ISP,SA,SB,
     *                                SC,SBB,S
                     ELSE
                        WRITE (6,982) J,LB(M),J-1,LBJ(I),DEL,ISP,SA,SB,
     *                                SC,SBB
                     ENDIF
                        WRITE (2) 8*IS1+IS2,IV,J,LB(M),J-1,LBJ(I),DEL,
     *                            ISP,SA,SB,SC,SBB,S,UP
                  ENDIF
   30          CONTINUE
            REWIND 1
c: 5/20/13
   35       SS=DEXP(FAC*EE(J21))*J21
            EMIN=MIN(EMIN,EE(1))
            IF (LB(J21).GT.0) THEN
               SST1=SST1+SS
            ELSE
               SST2=SST2+SS
            ENDIF
c:   35       IF (J.EQ.0) GOTO 60
            IF (J.EQ.0) GOTO 60
ce 5/20/13
c         Q-transitions
            WRITE (6,*) '  Q-TRANSITIONS',J
            RJJ1=DBLE(J21)/RJJ1
            DO 50 I=1,J21-1
c: 5/20/13
c:               ISP=MOD(I/2+(J22-I)/2,2)+3
c:               ISP=ISIG(ISP,IS,IV)
c:               SS=ISP
c:               SS=DEXP(FAC*EE(I))*SS
               SS=DEXP(FAC*EE(I))
               ISP=ISG1
               IF (LB(I).GT.0) THEN
                  SST1=SST1+SS*J21
               ELSE
                  SST2=SST2+SS*J21
                  ISP=ISG2
               ENDIF
               SS=SS*XINTF*ISP
ce 5/20/13
               DO 50 M=I+1,J21
                  S1C=(0D0,0D0)
                  S2C=S1C
                  S3C=CONJG(UC(1,I))*DBLE(-J)*UC(1,M)
                  DEL=EE(M)-EE(I)
                  IF (DEL.LT.FIV.OR.DEL.GT.FAV) GOTO 50
                  DO 40 KQ=2,J21
                     K=KQ-J1
                     S1C=S1C+CONJG(UC(KQ-1,I))*F1(J,-K)*UC(KQ,M)
                     S2C=S2C+CONJG(UC(KQ,I))*F1(J,-K)*UC(KQ-1,M)
   40                S3C=S3C+CONJG(UC(KQ,I))*DBLE(K)*UC(KQ,M)
                  SM=(REAL(CONJG(S3C)*(S1C+S2C))*DIP(2)
     *               -AIMAG(CONJG(S3C)*(S1C-S2C))*DIP(3))*DIP(1)
     *               +AIMAG(CONJG(S2C)*S1C)*DIP(2)*DIP(3)
                  SC=CONJG(S1C-S2C)*(S1C-S2C)*RJJ1
                  SB=CONJG(S1C+S2C)*(S1C+S2C)*RJJ1
                  SA=CONJG(S3C)*S3C*RJJ1
                  SM=SM/SJ+SA*DIP(1)**2+SB*DIP(2)**2+SC*DIP(3)**2
c: 5/20/13
c:                  IF (SM.GT.TH) THEN
c:                     SBB=SM*SS*(ONE-DEXP(FAC*DEL))
                  SBB=SM*SS*(ONE-DEXP(FAC*DEL))*DEL
                  IF (SBB.GT.TH) THEN
ce 5/20/13
                     UP=EE(M)/C0
                     S=ZERO 
                     IF (NUNC.NE.0) THEN
                        DO 45 KP=1,LP
                           E(KP)=DER(M,KP)-DER(I,KP)
                           S2=ZERO 
                           DO 44 L=1,KP
   44                         S2=S2+VCM(L,KP)*E(L)
   45                      S=S+E(KP)*(S2*TWO-VCM(KP,KP)*E(KP))
                        S=DSQRT(S)
                        WRITE (6,982) J,LB(M),J,LB(I),DEL,ISP,SA,SB,SC,
     *                                SBB,S
                     ELSE
                        WRITE (6,982) J,LB(M),J,LB(I),DEL,ISP,SA,SB,SC,
     *                                SBB
                     ENDIF
                     WRITE (2) 8*IS1+IS2,IV,J,LB(M),J,LB(I),DEL,ISP,SA,
     *                         SB,SC,SBB,S,UP
                  ENDIF
   50          CONTINUE
   60       WRITE (1) EE,LB,U,DER
            REWIND 1
c: 5/20/13
c:   90    CONTINUE
            SST(J+1,IS,1)=SST1*ISG1
            SST(J+1,IS,2)=SST2*ISG2
   80    CONTINUE
   90    CONTINUE
      IF (JMIN(IV).EQ.0) THEN
         WRITE (6,985) ((ISIG(1,IS,IV),ISIG(2,IS,IV),I=1,2),
     1                     IS=1,NSIG(IV))
	   DO 95 J=JMIN(IV)+1,JMAX(IV)+1
   95       WRITE (6,986) (SST(J,IS,1),SST(J,IS,2),IS=1,NSIG(IV))
	   SS=DEXP(-FAC*EMIN)
         J=JMAX(IV)+1
         WRITE (6,986)
         WRITE (6,986) (SST(J,IS,1)*SS,SST(J,IS,2)*SS,IS=1,NSIG(IV))
      ENDIF
  985 FORMAT(/' SUM OF STATES'//'IS1,IS2',8(I3,I1,7X))
  986 FORMAT(8F12.3)
ce 5/20/13
      J=0
      WRITE (2) J,J,J,J,J,J,ZERO,J,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO
c: 5/20/13
c:  100 CALL ORDER
  100 CALL ORDER(IFPR)
ce 5/20/13
      RETURN
      END
c: 5/20/13
c:      SUBROUTINE ORDER
c:      IMPLICIT REAL(8) (A-H,O-Z)
      SUBROUTINE ORDER(IFPR)
      IMPLICIT REAL(8) (A-H,O-Z)
c     Subroutine contains new code to print predictions in the format of a
c     JPL catalog file (SPCAT). Suggested by Brian Drouin (JPL).
c     Also larger dimensions to store information for 140000 transitions.
      REAL(8) LOGSTR0,LOGSTR1
      CHARACTER*20 FILECAT
c:      COMMON FR(50000),STR(3,50000),BL(50000),IS(50000),J(50000),IV(5000
c:     *0),N(50000),J1(50000),N1(50000),UN(50000),ISP(50000),UP(50000)
      COMMON FR(140000),STR(3,140000),BL(140000),IS(140000),J(140000),
     1 IV(140000),N(140000),J1(140000),N1(140000),UN(140000),
     2 ISP(140000),UP(140000)
      DATA ZERO,C0/0D0,29979.2458/
      IOUT=6
      IF (IFPR.EQ.4) THEN
         IOUT2=7
         WRITE (*,*)  'Enter catalog filename !'
         READ (*,*)    FILECAT
         WRITE (*,*)  'Enter catalog ID !'
         READ (*,*)    ICATID
         WRITE (*,*)  'Enter Partition function !'
         READ (*,*)    Q
         WRITE (*,*)  'Enter first LOG Intensity Cutoff (LOGSTR0) !'
         READ (*,*)    LOGSTR0
         WRITE (*,*)  'Enter second LOG Intensity Cutoff (LOGSTR1) !'
         READ (*,*)    LOGSTR1
         OPEN(IOUT2,FILE=FILECAT)
         WRITE (IOUT,996) ICATID,Q,LOGSTR0,LOGSTR1
      ENDIF
  996 FORMAT(//80('*')/' CATALOG ENTRY INFORMATION'/' Catalog ID',I21/
     1 ' PARTITION FUNCTION', D22.8/ ' 1st INT CUTOFF (LOGSTR0)',D16.8/
     2 ' 2nd INT CUTOFF (LOGSTR1)',D16.8/80('*')/)
c:      IOUT=6
ce 5/20/13
      REWIND 2
c: 5/20/13
c:      DO 10 I=1,50000
      DO 10 I=1,140000
ce 5/20/13
         READ (2) IS(I),IV(I),J(I),N(I),J1(I),N1(I),FR(I),ISP(I),
     1      (STR(K,I),K=1,3),BL(I),UN(I),UP(I)
c         FR = transition frequency
c         STR = line strengths
c         BL = line strength * spin weight * Boltzmann factor at temperature TEMP
c         ISP = spin multiplicity
c         UN  = uncertainty of transition frequency
c         UP = upper state energy (in cm-1)
         IF (N(I)*N1(I).EQ.0) GOTO 20
         IF (FR(I).GT.0D0) GOTO 10
         K=J(I)
         J(I)=J1(I)
         J1(I)=K
         K=N(I)
         N(I)=N1(I)
         N1(I)=K
         FR(I)=-FR(I) 
c: 5/20/13
c:         BL(I)=-BL(I)
c:   10 CONTINUE
c:      I=50001
   10 CONTINUE
      I=140001
ce 5/20/13
   20 I=I-1
      IF (I.EQ.0) RETURN 
      WRITE (6,997)
  997 FORMAT(/' TRANSITIONS ORDERED BY FREQUENCY'/)
      DO 50 M=1,I
         L=M
         S=FR(L)
         DO 40 K=M,I
            IF (FR(K).LT.S) L=K
   40       S=FR(L)
         FR(L)=FR(M)
         FR(M)=S
         DO 41 K=1,3
            S=STR(K,L)
            STR(K,L)=STR(K,M)
   41       STR(K,M)=S
         S=BL(L)
         BL(L)=BL(M)
         BL(M)=S
         S=UN(L)
         UN(L)=UN(M)
         UN(M)=S
         S=UP(L)
         UP(L)=UP(M)
         UP(M)=S
         K=IS(L)
         IS(L)=IS(M)
         IS(M)=K
         K=J(L)
         J(L)=J(M)
         J(M)=K
         K=IV(L)
         IV(L)=IV(M)
         IV(M)=K
         K=J1(L)
         J1(L)=J1(M)
         J1(M)=K
         K=N(L)
         N(L)=N(M)
         N(M)=K
         K=N1(L)
         N1(L)=N1(M)
         N1(M)=K
         K=ISP(L)
         ISP(L)=ISP(M)
         ISP(M)=K
   50    CONTINUE
      DO 65 M=1,I
         KA1=IABS(N(M))/2
         KC1=(2*J(M)+2-IABS(N(M)))/2
         KA2=IABS(N1(M))/2
         KC2=(2*J1(M)+2-IABS(N1(M)))/2
         K=IV(M)
c: 5/20/13
         IF (IFPR.EQ.4) THEN
            EL=UP(M)-FR(M)/C0
            IGUP=ISP(M)*(2*J(M)+1)
            IF (BL(M).EQ.ZERO) THEN
               BLM=-99.
            ELSE
               BLM=LOG10(BL(M)/Q)
            ENDIF
            THRLOG=LOG10(10**LOGSTR0+(FR(M)/300000)**2*10**LOGSTR1)
            IF (BLM.LT.THRLOG) GOTO 65
            IGUP=MIN(IGUP,999)
            UNM=DMIN1(UN(M),999.9999D0)
            WRITE (IOUT2,998) FR(M),UNM,BLM,3,EL,IGUP,ICATID,1404,J(M),
     1         KA1,KC1,IS(M)/8,J1(M),KA2,KC2,MOD(IS(M),8)
  998 FORMAT(F13.4,2F8.4,I2,F10.4,I3,I7,I4,2(4I2,4X))
         ENDIF
ce 5/20/13
   65    WRITE (IOUT,999) IS(M)/8,MOD(IS(M),8),K,J(M),N(M),KA1,KC1,K,
     1      J1(M),N1(M),KA2,KC2,FR(M),UN(M),ISP(M),(STR(K,M),K=1,3),
     2      BL(M),UP(M)
  999 FORMAT(I2,I1,2(I3,4I4),F12.3,F9.3,I4,4F8.4,F10.3)
      REWIND 2
      RETURN
      END
      SUBROUTINE INPUT(A,MQ,N1,N2,NC,INPAR,IVR,NTUP,JMIN,JMAX,LM,FRQ,BL,
     1                 WT,ITRA,ILEV,NIT,FMN,FMX,NSIG,ISIG,NTE,THRES,LP,
c: 5/20/13
c:     2                 NIV,IFPR,TEMP,NUNC,DIP,ISCD)
     2                 NIV,IFPR,TEMP,NUNC,DIP,ISCD,SCP)
ce 5/20/13
      IMPLICIT REAL(8) (A-H,O-Z)
      CHARACTER*4 WTWT
c: 5/20/13
c:      DIMENSION A(50,1),INPAR(34,1),IVR(50,1),JMIN(1),JMAX(1),FMN(1),
c:     1   FMX(1),THRES(1),NTUP(1),NTE(1),FRQ(1),WT(1),BL(1),ITRA(8,1),
c:     2   ILEV(2,1),NSIG(1),ISIG(4,10,1),SPAR(33),IGR(33),DIP(1)
      DIMENSION A(60,1),INPAR(38,1),IVR(60,1),JMIN(1),JMAX(1),FMN(1),
     1   FMX(1),THRES(1),NTUP(1),NTE(1),FRQ(1),WT(1),BL(1),ITRA(8,1),
     2   ILEV(2,1),NSIG(1),ISIG(4,10,1),SPAR(37),IGR(37),DIP(1),SCP(80),
     3   SCC(60,6)
ce 5/20/13
      LOGICAL LISC
      DATA NTRX,MS1,MS2,MS4/8191,8,128,16384/
      DATA JMX,ZERO,ONE,TWO,PI/120,0D0,1D0,2D0,3.141592653589793D0/
  900 FORMAT(1X,10A8)
  901 FORMAT(/' **** CALCULATION FOR ',A4,'EQUIVALENT MOTIONS ****'//
     1' SYMMETRY PARAMETER (ISCD)',I12/' DIRECTION COSINE PARAMETER (NC
     2)',I6/' NUMBER OF VIBRATIONAL STATES (NIV)',I3/' NUMBER OF ITERATI
     3ONS (NIT)',I11/' PRINT OPTION PARAMETER (IFPR)',I8/' UNCERTAINTY P
     4ARAMETER (NUNC)',I9/' DIPOLE MOMENT COMPONENTS (DIP) ',3F6.3/' TEM
     5PERATURE/KELVIN (TEMP)',F12.2/' VIB STATE READ (KVQ,1=YES,ELSE=NO)
     6',I3)
c: 5/20/13
c:  902 FORMAT(6I4,F15.8,I4)
c:  903 FORMAT(' INPUT ERROR',6I4,F15.8,I4,'  INPUT DELETED')
  902 FORMAT(6I4,F15.8,I4,D12.3)
  903 FORMAT(' INPUT ERROR',6I4,F15.8,I4,D12.3'  INPUT DELETED')
ce 5/20/13
  904 FORMAT(' INPUT ERROR',6I4,F12.4,2F8.3,'    TRANSITION DELETED')
  905 FORMAT(8I4,F15.4,2F8.3)
  906 FORMAT(/I4,' TRANSITIONS,',I4,' WITH NON-ZERO WEIGHT')
  907 FORMAT(/'VIBRATIONAL STATE',I4/' JMIN',I4,'  JMAX',I4,6X,'FMIN',
     1   F10.1,6X,'FMAX',F10.1,'  INTENSITY CUTOFF',F6.2/3F15.8/
     2   2(5F15.8/),2F15.8/15I4)
  908 FORMAT(' ILLEGAL INPUT OF VARIABLE PARAMETERS'/17I4)
  909 FORMAT('INPUT ERROR',2I4,6X,'INPUT DELETED')
  910 FORMAT(8X,'SYMMETRY BLOCK',3X,2I1,'  SPIN WEIGHTS',2I4)
  911 FORMAT(/' PERIODICITY OF INTERNAL ROTOR (N)',I13,I14/' RHO PARAMET
     1ER',F33.8,F14.8/' VARIATION PARAMETER FOR RHO, IVRHO',I12,I14/
     2' RHO AXIS ANGLE, BETA',F26.8,F14.8/' VARIATION PARAMETER FOR BETA
     3, IVBET',I11,I14/' RHO AXIS ANGLE, ALPHA',F25.8,F14.8/' VARIATION 
     4PARAMETER FOR ALPHA, IVALP',I10,I14)
      READ (5,900) (WT(I),I=1,10)
c     problem information
      WRITE (6,900) (WT(I),I=1,10) 
      WRITE (*,900) (WT(I),I=1,10)
      KVQ=0
      IIN=5
      READ (5,*) ISCD,NC,NIV,NIT,IFPR,NUNC,(DIP(I),I=1,3),TEMP
c     ISCD = symmetry parameter
c            ISCD = 1 => C'1 (non-equivalent) = -1 => C's (non-equivalent)
c                   2 => C2 (equivalent)      = -2 => C2v (equivalent)
c                   3 => Cs (equivalent)      = -3 => C2v (equivalent)
c                   4 => Ci (equivalent)      = -4 => C2h (equivalent)
c                               0 => C's (non-equivalent, in bc plane)
c                              -3 => C2v (equivalent, in bc plane)
c                              -5 => C2h (equivalent, in bc plane)
c     NC  = direction cosine parameter (with respect to principal a axis)
c           NC = +1 ==> direction cosines of internal rotors are equal
c           NC = -1 ==> direction cosines of internal rotors are opposite
c     NIV = number of vibrational states (min=1, max=6)
c     NIT = number of iterations (NIT.LE.0 => no iterations)
c     IFPR = print option parameter
c           IFPR = 0  no additional printout
c                = 1  prints derivatives,   a LOT of printout !!!!
c                = 2  prints eigenvectors of rotational energy levels
c                = 3  prints eigenvectors and derivatives 
c                = 4  prints predictions also in JPL catalogue file
c                     format as separate file (no derivatives or eigenvectors)
c     NUNC = uncertainties parameter
c           NUNC.NE.0 : calculate uncertainties of predicted frequencies
c           NUNC = 0  : no uncertainties of predicted frequencies 
c     DIP = components of electric dipole moment (mu(a),mu(b),mu(c)
c     TEMP = temperature in K(elvin) for Boltzmann population factor
c     KVQ = switch for vibrational label
c           KVQ = 0 : don't read vibrational states label
c           KVQ = 1 : read vibrational states label
c
c: 5/20/13
c      IF (MIN(5+ISCD,4-ISCD,NIV-1,6-NIV,NIT,IFPR,3-IFPR).LT.0.OR.
      IF (MIN(5+ISCD,4-ISCD,NIV-1,6-NIV,NIT,IFPR,4-IFPR).LT.0.OR.
ce 5/20/13
     1     IABS(NC).NE.1) GOTO 220
      IF (IABS(ISCD).GT.3) NC=-1
      MQ=3-MIN(MAX(IABS(ISCD),1),2)
      LISC=(ISCD*NC.EQ.3.OR.ISCD*NC.EQ.-2).AND.ISCD.GT.0
      WTWT='NON-'
      IF (MQ.EQ.1) WTWT='    '
      IF (NIT.LE.0) NUNC=0
      WRITE (6,901) WTWT,ISCD,NC,NIV,NIT,IFPR,NUNC,(DIP(I),I=1,3),TEMP,
     1 KVQ
      WRITE (*,901) WTWT,ISCD,NC,NIV,NIT,IFPR,NUNC,(DIP(I),I=1,3),TEMP,
     1 KVQ
      LXY=0
      IF (ISCD.EQ.0.OR.ISCD.EQ.-3.OR.ISCD.EQ.-5) LXY=1
      READ (5,*) N1,RHO1,IVRHO1,BETA1,IVBET1,ALPHA1,IVALP1
      IF (IVRHO1.NE.1) IVRHO1=0
      IF (IVBET1.NE.1) IVBET1=0
      IF (IVALP1.NE.1) IVALP1=0
      IF (ISCD.LT.0.AND.LXY.EQ.0) THEN
         ALPHA1=ZERO
         IVALP1=0
      ENDIF
      IF (LXY.EQ.1) THEN
         BETA1=PI/TWO
         IVBET1=0
      ENDIF      
      IF (MQ.EQ.2) THEN
         READ (5,*) N2,RHO2,IVRHO2,BETA2,IVBET2,ALPHA2,IVALP2
         IF (IVRHO2.LT.1.OR.IVRHO2.GT.2) IVRHO2=0
         IF (IVBET2.LT.1.OR.IVBET2.GT.2) IVBET2=0
         IF (IVALP2.LT.1.OR.IVALP2.GT.2) IVALP2=0
         IF (ISCD.EQ.-1) THEN
            ALPHA2=ALPHA1
            IVALP2=IVALP1
         ELSE IF (ISCD.EQ.0) THEN
               BETA2=BETA1
               IVBET2=IVBET1
         ENDIF
      ELSE
         N2=N1
         RHO2=RHO1
         IVRHO2=0
         BETA2=BETA1
         IVBET2=0
         ALPHA2=ALPHA1
         IF (LISC.OR.ISCD.EQ.-3) ALPHA2=-ALPHA2
         IVALP2=0
      ENDIF
c     N1,N2 =   periodicities of internal rotors (0 < N1, 0 < N2)
c     RHO1, RHO2 = internal rotation parameters
c     IVRHO1, IVRHO2 = variation parameter for RHO
c           IVRHO  = 0 : keep RHO constant
c           IVRHO  = 1 : RHO is a variable in least-squares fit
c           IVRHO2 = 2 : RHO2 is equal to RHO1 (non-equivalent rotors)
c     BETA1, BETA2 = RHO axis polar angle
c           between RHO axis and "a" principal axis
c     IVBET1, IVBET2 = variation parameter for BETA
c           IVBET = 0 : keep BETA constant
c           IVBET = 1 : BETA is a variable in least-squares fit
c           IVBET2 = 2 : BETA2 is equal to BETA1 (non-equivalent rotors)
c     ALPHA1, ALPHA2 = RHO axis angle
c           of RHO axis w.r.t. "ab" principal plane
c     IVALP1, IVALP2 = variation parameter for ALPHA
c           IVALP = 0 : keep ALPHA constant
c           IVALP = 1 : ALPHA is a variable in least-squares fit
c           IVALP2 = 2 : ALPHA2 is equal to +/- ALPHA1 (non-equivalent rotors)
c     Note: N2, RHO2, IVRHO2, BETA2, IVBET2, ALPHA2, IVALP2 are needed only
c                                                                    if MQ = 2
      IF (MIN(N1,N2).LT.1) GOTO 230
      WRITE (*,911) N1,N2,RHO1,RHO2,IVRHO1,IVRHO2,BETA1,BETA2,IVBET1,
     1       IVBET2,ALPHA1,ALPHA2,IVALP1,IVALP2
      WRITE (6,911) N1,N2,RHO1,RHO2,IVRHO1,IVRHO2,BETA1,BETA2,IVBET1,
     1       IVBET2,ALPHA1,ALPHA2,IVALP1,IVALP2
      BETA1=BETA1*PI/180D0
      BETA2=BETA2*PI/180D0
      ALPHA1=ALPHA1*PI/180D0
      ALPHA2=ALPHA2*PI/180D0
      LM=0
      LL=0
      DO 100 IV=1,NIV
         DO 5 KP=1,6
    5       IVR(KP,IV)=-1
         READ (5,*) JMIN(IV),JMAX(IV),FMN(IV),FMX(IV),THRES(IV),
     1              (A(I,IV),I=7,21),(IVR(I,IV),I=7,21)
c
c     JMIN, JMAX = range of J quantum number for predictions ( < 121 )
c               No predictions result from JMIN > JMAX
c     FMN, FMX = lower and upper limit of frequencies for predictions
c     THRES = intensity threshold for predictions
c     array A : 3 rotational constants,
c               5 quartic centrifugal distortion constants (A-reduction)
c               7 sextic centrifugal distortion constants (A-reduction)
c     array IVR : variation parameters for rotational and distortion constants
c               = 0 : keep constant at input value
c               > 0 : variable during least-squares fit
c               < 0 : constant is kept identical to corresponding parameter of
c                     the vibrational state for IV = ABS(IVR)
c
         JMAX(IV)=MIN(JMAX(IV),JMX)
         WRITE (6,907) IV,JMIN(IV),JMAX(IV),FMN(IV),FMX(IV),THRES(IV),
     1              (A(I,IV),I=7,21),(IVR(I,IV),I=7,21)
c     input of tunneling parameters             (less than 32 per state)
         I=0
         II=0
c: 5/20/13
c:   10    READ (5,*) IQ1,IQ2,MEG,KAP,JP,KP,PAR,IVAR
   10    READ (5,*) IQ1,IQ2,MEG,KAP,JP,KP,PAR,IVAR,SCPP
ce 5/20/13
c     IQ1,IQ2 = indicate to which localized state a tunneling parameter is
c               related to
c               max(IQ1,IQ2) = 7, min(IQ1) = 0
c               min(IQ2) = 0 for equiv, min(IQ2) = -7 for nonequiv motions
c     MEG     = omega (omega = +1 for even order of angular momenta,
c                              -1     odd                          ,
c                              0  terminates input of tunneling parameters)
c     KAP = distance of matrix element from main diagonal (KAPPA = KQ - K)
c           KAP>0 real, imaginary part of matrix element for MEG = +1, -1, resp.
c           KAP<0 imaginary, real part of matrix element for MEG = +1, -1, resp.
c
c     JP  = exponent of P operator (total angular momentum)
c     KP  = exponent of Pz operator (component of P along z)
c     PAR = value of tunneling parameter
c     IVAR = variation of tunneling parameter
c          =+1 : yes
c          =-1 : yes, same value as tunneling parameter with IQ1 and IQ2
c                interchanged
c          else: no 
c     SCPP = scaling factor for parameter
c
c     comment: the vibrational energy tunneling parameters have
c              KAP = 0, JP = 0, KP = 0, MEG = +1
c
c  ISCD =  1 => C'1 (non-equivalent) effective
c               KAP >=< 0  => IMG = 0,1
c               MEG = +1 => SumC1
c               MEG = -1 => SumS2
c  ISCD = -1 => C's (non-equivalent) effective
c               KAP >= 0 => IMG = 0
c               MEG = +1 => SumC1
c               MEG = -1 => SumS2
c  ISCD =  2 => C2 (equivalent) effective NC = +1
c               KAP >=< 0  => IMG = 0,1
c               MEG = +1, KAP = even => SumC3
c                         KAP = odd  => SumC4
c               MEG = -1, KAP = even => SumS5
c                         KAP = odd  => SumS6
c  ISCD =  2 => C2 (equivalent) effective NC = -1
c               KAP >=< 0  => IMG = 0,1
c               MEG = +1, KAP = even, IMG = 0 => SumC3
c                                     IMG = 1 => SumC4
c                         KAP = odd,  IMG = 0 => SumC4
c                                     IMG = 1 => SumC3
c               MEG = -1, KAP = even, IMG = 0 => SumS6
c                                     IMG = 1 => SumS5
c                         KAP = odd,  IMG = 0 => SumS5
c                                     IMG = 1 => SumS6
c  ISCD =  3 => C3 (equivalent) effective NC = +1
c               KAP >=< 0  => IMG = 0,1
c               MEG = +1, KAP = even, IMG = 0 => SumC3
c                                     IMG = 1 => SumC4
c                         KAP = odd,  IMG = 0 => SumC4
c                                     IMG = 1 => SumC3
c               MEG = -1, KAP = even, IMG = 0 => SumS5
c                                     IMG = 1 => SumS6
c                         KAP = odd,  IMG = 0 => SumS6
c                                     IMG = 1 => SumS5
c  ISCD =  3 => Cs (equivalent) effective NC = -1
c               KAP >=< 0  => IMG = 0,1
c               MEG = +1, KAP = even => SumC3
c                         KAP = odd  => SumC4
c               MEG = -1, KAP = even => SumS6
c                         KAP = odd  => SumS5
c  ISCD = -2 => C2v (equivalent) effective NC = +1
c               KAP >= 0 => IMG = 0
c               MEG = +1, KAP = even => SumC3
c                         KAP = odd  => SumC4
c               MEG = -1, KAP = even => SumS5
c                         KAP = odd  => SumS6
c  ISCD = -2 => C2v (equivalent) effective NC = -1
c               KAP >= 0 => IMG = 0
c               MEG = +1, KAP = even => SumC3
c                         KAP = odd  => SumC4
c               MEG = -1, KAP = even => SumS6
c                         KAP = odd  => SumS5
c  SumC1: IQ1+IQ2 > 0 or
c         IQ1+IQ2 = 0 .AND. IQ1-IQ2 >= 0
c  SumS2: IQ1+IQ2 > 0 or
c         IQ1+IQ2 = 0 .AND. IQ1-IQ2 > 0
c  SumC3: IQ1+IQ2 >= 0 .AND. IQ1-IQ2 >= 0
c  SumC4: IQ1+IQ2 >  0 .AND. IQ1-IQ2 >  0
c  SumS5: IQ1+IQ2 >  0 .AND. IQ1-IQ2 >= 0
c  SumS6: IQ1+IQ2 >= 0 .AND. IQ1-IQ2 >  0
c
         IF (MEG.EQ.0) GOTO 20
c   Original version: for MEG = +1 and MEG = -1, only real and imaginary
c     parts of tunneling integrals were allowed, respectively -> IMG = 0
c   Newer versions: imaginary and real parts, respectively, were made
c     available by allowing MEG to be of opposite sign -> IMG = 1
c   Modification Version16: Use KAP < 0 to indicate IMG = 1
c   NOTE: Elsewhere in the program, IMG is automatically increased by 1;
c         that means IMG values 0 & 1 used here become 1 & 2, resp.
         MEG=SIGN(1,MEG)
         IMGS=SIGN(1,KAP)
         IMG=0
         IF (KAP.NE.0.AND.ISCD.GT.0) IMG=(1-IMGS)/2
      KAP=IABS(KAP)
      L=MOD(KAP,2)
      IRO=JP+KP+KAP
         IQ=IQ1+IQ2
         IQQ=IQ1-IQ2
c   Allow for S-reduction Hamiltonians
         IF (MIN(7-IQ1,7-IABS(IQ2),IQ,JP,KP,14-KAP,ISCD*IMG).LT.0)
     *       GOTO 14
         IF (ISCD*IMG.LT.0.OR.KAP+1-IMG.EQ.0) GOTO 14
         IF (ISCD*NC*IMG.EQ.3.OR.ISCD*NC*IMG.EQ.-2) L=1-L
         IF (IABS(ISCD).EQ.4) L=0
         IF (MQ.EQ.1) THEN
            IF (MIN(IQQ,N1-IQ1+NC*IQ2).GE.0.AND.
     *       (1+IQQ)*(L*(IQ-1)+1)+MEG*(IQQ-1)*(L*(IQ+1)-1).GT.0) GOTO 15
         ELSE
            IF (IQ.GT.0.OR.(IQ.EQ.0.AND.IQQ.GE.1-MEG)) GOTO 15
         ENDIF
c: 5/20/13
c:   14    WRITE (6,903) IQ1,IQ2,MEG,IMGS*KAP,JP,KP,PAR,IVAR
c:         WRITE (*,903) IQ1,IQ2,MEG,IMGS*KAP,JP,KP,PAR,IVAR
c:   15    IF (IABS(IVAR).NE.1) IVAR=0
c:         IF (I.LT.31) GOTO 16
   14    WRITE (6,903) IQ1,IQ2,MEG,IMGS*KAP,JP,KP,PAR,IVAR,SCPP
         WRITE (*,903) IQ1,IQ2,MEG,IMGS*KAP,JP,KP,PAR,IVAR,SCPP
   15    IF (IABS(IVAR).NE.1) IVAR=0
         IF (I.LT.37) GOTO 16
ce 5/20/13
         WRITE (6,*) ' TOO MANY TUNNELING PARAMETERS FOR STATE IV =',IV
         STOP
   16    I=I+1
c   Generate parameter information
c   IRO = 0: tunneling energy parameter
c   IRO = 1: other tunneling parameter
         IRO=MIN(IRO,1)
         II=II+1-IRO
         INPAR(I,IV)=(((((IRO*16+8+IQ1)*16+8+IQ2)*4+MEG+2-IMG)*16+KAP)*
     1               16+JP)*16+KP
         SPAR(I)=PAR
         IGR(I)=IVAR
c: 5/20/13
c:         WRITE (6,902) IQ1,IQ2,MEG,IMGS*KAP,JP,KP,PAR,IVAR
         SCP(I)=SCPP
         WRITE (6,902) IQ1,IQ2,MEG,IMGS*KAP,JP,KP,PAR,IVAR,SCPP
ce 5/20/13
         GOTO 10
   20    DO 50 KP=1,I 
            JP=KP
            IQ=INPAR(JP,IV)
            DO 45 K=KP,I
               IF (INPAR(K,IV).LT.IQ) JP=K
   45          IQ=INPAR(JP,IV)
            S=SPAR(JP)
            SPAR(JP)=SPAR(KP)
            SPAR(KP)=S
c: 5/20/13
            S=SCP(JP)
            SCP(JP)=SCP(KP)
            SCP(KP)=S
ce 5/20/13
            L=IGR(JP)
            IGR(JP)=IGR(KP)
            IGR(KP)=L
            INPAR(JP,IV)=INPAR(KP,IV)
   50       INPAR(KP,IV)=IQ
         INPAR(I+1,IV)=0 
         DO 60 KP=1,I
            A(KP+21,IV)=SPAR(KP)
c: 5/20/13
            SCC(KP+21,IV)=SCP(KP)
ce 5/20/13
   60       IVR(KP+21,IV)=IGR(KP)
         NTUP(IV)=I
         NTE(IV)=II
C     input of symmetry blocks for predictions
         I=0
         IVQ=IV
         IV1=IV
   70    READ (5,*) IS1,IS2,ISP1,ISP2
c     IS1,IS2 = symmetry numbers of symmetry block
c     ISP1,ISP2 = spin multiplicity for (sym) and (antisym), resp.
         IF (IS1.LT.0) GOTO 79
         IF (MQ.EQ.1) THEN
         IF (MIN(IS2-IS1,N1/2-IS1,N1-IS1-IS2,N1/2+IS1*(1+IS2)-IS2).GE.0)
     *       GOTO 75
         ELSE 
            IF (MIN(N1/2-IS1,IS2,N2/2-IS2).GE.0) GOTO 75
            IF ((IS1.EQ.0.OR.IS1.EQ.N1/2).AND.IS2.LT.N2) GOTO 75
         ENDIF
         WRITE (6,909) IS1,IS2
         WRITE (*,909) IS1,IS2
         GOTO 70
   75    I=I+1
         ISIG(1,I,IV)=IS1
         ISIG(2,I,IV)=IS2
         ISIG(3,I,IV)=ISP1
         ISIG(4,I,IV)=ISP2
         WRITE (6,910) IS1,IS2,ISP1,ISP2
         NSIG(IV)=I
         GOTO 70
c     everything is in but the transitions
c     input of transitions
   79 BLS=ZERO
      IBL=1
      FREQQ=-ONE
   80 IF (KVQ.NE.1) THEN
         READ (IIN,*) IS1,IS2,JQ,NQ,J,NN,FREQ,BLE,ER
      ELSE
         READ (IIN,*) IS1,IS2,IVQ,JQ,NQ,IV1,J,NN,FREQ,ER,BLE
      ENDIF
c
c     IS1, IS2  =  symmetry label of transition
c     JQ, NQ    = J quantum number and label of energy level of upper state
c     J, NN     = J quantum number and label of energy level of lower state
c     FREQ      = frequency of transitions (MHz)
c     BLE       = coefficient in linear combinations for blends
c                 BLE = 0. ==> no blend
c                 BLE > 0. ==> transition(s) of blended line
c                 BLE < 0. ==> last transition of blended line
c                 only one frequency (the weighted average with weight equal
c                 to BLE/sum(BLE) ) is fit
c     ER        = estimated uncertainty of transition (MHz) (its inverse square
c                 becomes the weight of the transition in the least-square fit)
c                 ER = 0. ==> zero weight in fit for this transition
c
         IF (IS1.LT.0) GOTO 98
         IF (MAX(NN-2*J-1,NQ-2*JQ-1,IABS(JQ-J)-1,JQ-JMX,J-JMX,-JQ,-NQ,
     *           -J,-NN,-IS2,IVQ-NIV,IV1-NIV).GT.0) GOTO 84
         IF (MQ.EQ.1) THEN
            IF (MAX(IS1-N1/2,IS1+IS2-N1,IS1-IS2,IS2-N1/2-IS1*(1+IS2))
     *          .LE.0) GOTO 85
         ELSE 
            IF (MIN(IS1-N1/2,IS2-N2/2).GE.0) GOTO 85
            IF ((IS1.EQ.0.OR.IS1.EQ.N1/2).AND.IS2.LT.N2) GOTO 85
         ENDIF
         ER=DABS(ER)
   84    WRITE (6,904) IS1,IS2,JQ,NQ,J,NN,FREQ,BLE,ER
         WRITE (*,904) IS1,IS2,JQ,NQ,J,NN,FREQ,BLE,ER
         GOTO 80
   85    LM=LM+1
         ITRA(1,LM)=IS1
         ITRA(2,LM)=IS2
         ITRA(3,LM)=IVQ
         ITRA(4,LM)=JQ
         ITRA(5,LM)=NQ
         ITRA(6,LM)=IV1
         ITRA(7,LM)=J
         ITRA(8,LM)=NN
         FRQ(LM)=FREQ
         WT(LM)=ZERO
         IF (ER.EQ.ZERO) GOTO 90
         LL=LL+1
         WT(LM)=ONE/ER
   90    BL(LM)=BLE
         IF (IIN.EQ.5) THEN
            IF (BLE.EQ.ZERO) THEN
               IF (BLS.GT.ZERO) THEN
                  DO 91 I=IBL,LM
   91                BL(I)=ZERO
                  BLS=ZERO
               ENDIF
               IBL=LM+1
            ELSE
c   Generate information for blended transitions
               IF (IBL.EQ.LM) THEN
                  FREQQ=FREQ
                  WTT=WT(LM)
               ENDIF
               FRQ(LM)=FREQQ
               WT(LM)=WTT
               BLS=BLS+DABS(BLE)
               IF (BLE.LT.ZERO) THEN
                  DO 92 I=IBL,LM
   92                BL(I)=BL(I)/BLS
                  BLS=ZERO
                  IBL=LM+1
               ENDIF
            ENDIF
         ELSE
            IF (FREQ.NE.FREQQ) THEN
               IF (LM.GT.1+IBL) THEN
                  IF (BLS.EQ.ZERO) BLS=ONE
                  K=-1
                  DO 93 I=IBL,LM-1
                     IF (BL(I).NE.ZERO) K=I
   93                BL(I)=DABS(BL(I))/BLS
                  IF (K.NE.-1) BL(K)=-BL(K)
               ENDIF
               BLS=DABS(BL(LM))
               WTT=WT(LM)
               FREQQ=FREQ
               IBL=LM
            ELSE
               BLS=BLS+DABS(BL(LM))
               WT(LM)=WTT
            ENDIF
         ENDIF
c   Generate energy level information
         K=(IS1*MS1+IS2)*MS1
         ILEV(1,LM+LM-1)=(K+IVQ)*MS2+JQ
         ILEV(2,LM+LM-1)=NQ*MS4+MS4/2+LM
         ILEV(1,LM+LM)  =(K+IV1)*MS2+J
         ILEV(2,LM+LM)  =NN*MS4+MS4/2-LM
         WRITE (6,905) IS1,IS2,IVQ,JQ,NQ,IV1,J,NN,FREQ,BLE,ER
         IF (LM.LT.NTRX) GOTO 80
   98    IF (IIN.EQ.4.AND.LM.GT.IBL) THEN
            IF (BLS.EQ.ZERO) BLS=ONE
            K=-1
            DO 99 I=IBL,LM
               IF (BL(I).NE.ZERO) K=I
   99          BL(I)=DABS(BL(I))/BLS
            IF (K.NE.-1) BL(K)=-BL(K)
         ENDIF
  100    CONTINUE
c   Print summary of transitions and sort energy level information
      WRITE (6,906) LM,LL
      DO 130 KP=1,2*LM
         JP=KP
         IQ=ILEV(1,JP)
         NN=ILEV(2,JP)
         DO 120 K=KP,2*LM
      IF (ILEV(1,K).LT.IQ.OR.(ILEV(1,K).EQ.IQ.AND.ILEV(2,K).LT.NN)) JP=K
            NN=ILEV(2,JP)
  120       IQ=ILEV(1,JP)
         ILEV(1,JP)=ILEV(1,KP)
         ILEV(2,JP)=ILEV(2,KP)
         ILEV(1,KP)=IQ
  130    ILEV(2,KP)=NN
      ILEV(1,2*LM+1)=0
      A(1,1)=RHO1
      A(2,1)=RHO2
      A(3,1)=BETA1
      A(4,1)=BETA2
      A(5,1)=ALPHA1
      A(6,1)=ALPHA2
      IVR(1,1)=IVRHO1
      IVR(2,1)=IVRHO2
      IVR(3,1)=IVBET1
      IVR(4,1)=IVBET2
      IVR(5,1)=IVALP1
      IVR(6,1)=IVALP2
c: 5/20/13
      SCC(1,1)=1D-3
      SCC(2,1)=1D-3
      DO 135 IV=1,NIV
         DO 135 KP=3,21
  135       SCC(KP,IV)=ONE
ce 5/20/13
      LP=0
      DO 150 KP=7,21
         DO 140 IV=1,NIV
  140       IVR(KP,IV)=MIN(IVR(KP,IV),1)
  150    LP=MIN(LP,IVR(KP,1))
      IF (LP.GE.0) GOTO 160
      WRITE (6,908) (IVR(KP,1),KP=1,21)
      STOP
  160 LP=0
      DO 200 IV=1,NIV
         DO 190 KP=1,21+NTUP(IV)
  170       IF (IVR(KP,IV).GT.0) THEN
               IF (IVR(KP,IV).EQ.1) LP=LP+1
               IVR(KP,IV)=LP
c: 5/20/13
               SCP(LP)=SCC(KP,IV)
ce 5/20/13
            ELSE IF (IVR(KP,IV).LT.0) THEN
c
c   RE: Rotational and centrifugal distortion constants in different states
c
               IF (KP.LE.21) THEN
                  I=-IVR(KP,IV)
                  IVR(KP,IV)=IVR(KP,I)
                  A(KP,IV)=A(KP,I)
               ELSE
c
c   RE: Tunneling parameters (non-equivalent rotors)
c
                  IQ=MOD(INPAR(KP-21,IV)/256,16384)
                  IQ=IABS(MOD(IQ/1024,16)-MOD(IQ/64,16))*245760
                  DO 180 JP=22,21+NTUP(IV)
                   IF (IABS(INPAR(KP-21,IV)-INPAR(JP-21,IV)).EQ.IQ) THEN
                        IF (JP.LT.KP) THEN
                           IVR(KP,IV)=IVR(JP,IV)
                           A(KP,IV)=A(JP,IV)
                        ELSE IF (JP.GT.KP) THEN
                           IVR(JP,IV)=-1
                           IVR(KP,IV)=1
                           GOTO 170
                        ENDIF
                     ENDIF
  180             CONTINUE
               ENDIF
            END IF
  190       CONTINUE
  200    CONTINUE
      IF (LP.GT.0) GOTO 210
      NIT=0
      NUNC=0
  210 IF (LP.LE.80) RETURN
      WRITE (6,*) ' TOO MANY VARIABLE PARAMETERS',LP
      GOTO 240
  220 WRITE (6,*) ' INPUT ERROR'
      WRITE (6,*) ' ISCD=',ISCD,' NC=',NC,' NIV=',NIV,' NIT=',NIT,
     1            ' IFPR=',IFPR
  230 WRITE (6,*) ' INPUT ERROR','   N1=',N1,'   N2=',N2
  240 END
      SUBROUTINE LEASQU1(LSMD,A,NA,RES,DELP,N,M,P,SP,SIGMA2,SS,B,AC,NC)
C  FILE LEASQU1                                 VERSION FROM 25 SEP 1990
C                        CHANGES MADE IN INITIAL COMMENTS ON 13 FEB 1991
C  LSMD    I : CONTROL PARAMETER = -1 LINEAR LEAST SQUARES
C                                =  0 REGULAR LEAST SQUARES
C                                =  1 DAMPED LEAST SQUARES
C                                =  2 MODIFIED DAMPED LEAST SQUARES
C                                =  3 DIAGNOSTIC LEAST SQUARES
C  A(M,N)  I : JACOBIAN MATRIX (M.GE.N)
C   (N,N)  O : CORRELATION MATRIX
C  NA      I : LEADING DIMENSION OF A IN DIMENSION STATEMENT (NA.GE.M)
C  RES(M)  I : VECTOR OF RESIDUALS
C          O : FINAL VECTOR OF RESIDUALS IF LSMD = -1
C  DELP(M) I : VECTOR OF RESIDUALS OF ESTIMATED PARAMETERS IF LSMD = 2,
C              DUMMY ARGUMENT UNLESS LSMD = 2
C  P(2*N)  O : FIRST N COMPONENTS = VECTOR OF CHANGES TO PARAMETERS
C          O : LAST N COMPONENTS = VECTOR FOR PRECISION ESTIMATE
C  SP(N)   O : VECTOR OF STANDARD ERRORS
C  SIGMA2  O : VARIANCE
C  SS      O : SUM OF SQUARES
C  B(M)      : WORK AREA (BECOMES RESIDUAL VECTOR IN EIGENSPACE)
C  AC(M,N)   : WORK AREA (TEMPORARY STORAGE FOR A)
C  NC      I : LEADING DIMENSION OF AC IN DIMENSION STATEMENT (NC.GE.M)
C              IF LSMD.NE.-1 AND A IS NO LONGER NEEDED, A AND AC MAY
C              BE IDENTICAL WITH NC=NA
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NA,1),RES(1),DELP(1),P(1),SP(1),B(1),AC(NC,1)
      DATA ZERO,NPRIN/0D0,6/
  910 FORMAT(D12.4,12F10.5/(12X,12F10.5))
  911 FORMAT(/' SINGULAR VALUES AND RIGHT SINGULAR VECTORS FOR JACOBIAN
     *MATRIX'/)
  912 FORMAT(/I3,' SINGULAR VALUES USED OUT OF',I3,', STANDARD DEVIATION
     * =',D16.8)
      SS=ZERO
      DO 10 L=1,M
         B(L)=RES(L)
         SS=SS+RES(L)*RES(L)
         DO 10 I=1,N
   10       AC(L,I)=A(L,I)
      CALL LSVDF(A,NA,M,N,B,M,1,SP,P,IER)
      WRITE (NPRIN,911)
      DO 20 I=1,N
         WRITE (NPRIN,910) SP(I),(A(K,I),K=1,N)
         P(I+N)=ZERO
         DO 20 K=1,N
   20       P(I+N)=P(I+N)+(SP(K)*A(I,K))**2
      IF (LSMD.EQ.3) GOTO 70
      NN=N
      S1=M
      IF (LSMD.LE.0) S1=M-N
      IF (LSMD.NE.2) GOTO 40
      DO 30 I=1,N
   30    S1=S1-DELP(I)*DELP(I)
   40 SIGMA2=SS/S1
      S1=ZERO
      IF (LSMD.GT.0) S1=SIGMA2
      DO 50 I=1,N
   50    SP(I)=DSQRT(SP(I)*SP(I)+S1)
      S1=DSQRT(SIGMA2)
      GOTO 100
   70 NN=MIN(M-1,N)
   75 S1=DSQRT(SS/DFLOAT(M-NN))
      IF (SP(NN).GT.S1.OR.NN.EQ.1) GOTO 80
      NN=NN-1
      GOTO 75
   80 SIGMA2=S1*S1
      DO 90 I=1,N
   90    SP(I)=DMAX1(SP(I),S1)
  100 WRITE (NPRIN,912) NN,N,S1
      DO 120 K=1,N
         P(K)=ZERO
         DO 110 I=1,N
  110       A(K,I)=A(K,I)/SP(I)
         DO 120 I=1,NN
  120       P(K)=P(K)+A(K,I)*B(I)
      DO 140 K=1,N
         P(K+N)=DSQRT(SIGMA2*DFLOAT(M)/P(K+N))/DFLOAT(N)
         DO 130 L=K,N
            B(L)=ZERO
            DO 130 I=1,N
  130          B(L)=B(L)+A(K,I)*A(L,I)
         DO 140  L=K,N
  140       A(K,L)=B(L)*SIGMA2
      DO 160 K=1,N
         SP(K)=DSQRT(A(K,K))
         DO 150 L=1,K
  150       A(K,L)=A(L,K)
         IF (LSMD.NE.2) GOTO 156
         DO 155 L=1,N
  155       P(K)=P(K)+A(K,L)*DELP(L)
  156    DO 160 L=1,K
  160       A(K,L)=A(K,L)/(SP(K)*SP(L))
      IF (LSMD.NE.-1) GOTO 190
      SS=ZERO
      DO 180 L=1,M
         DO 170 I=1,N
  170       RES(L)=RES(L)-AC(L,I)*P(I)
  180    SS=SS+RES(L)*RES(L)
  190 RETURN
      END
C   IMSL ROUTINE NAME - LSVDB
C-----------------------------------------------------------------------
C   COMPUTER          - IBM77/DOUBLE
C   LATEST REVISION   - NOVEMBER 1, 1984
C   PURPOSE           - SINGULAR VALUE DECOMPOSITION OF A BIDIAGONAL
C                         MATRIX.
C   USAGE             - CALL LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C   ARGUMENTS    D    - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                       ON INPUT, D CONTAINS THE DIAGONAL ELEMENTS OF
C                         THE BIDIAGONAL MATRIX B. D(I)=B(I,I),
C                         I=1,...,N.
C                       ON OUTPUT, D CONTAINS THE N (NONNEGATIVE)
C                         SINGULAR VALUES OF B IN NONINCREASING ORDER.
C                E    - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                       ON INPUT, E CONTAINS THE SUPERDIAGONALELEMENTS
C                         OF B. E(1) IS ARBITRARY,
C                         E(I)=B(I-1,I), I=2,...,N.
C                       ON OUTPUT, THE CONTENTS OF E ARE MODIFIED BY THE
C                         SUBROUTINE.
C                N    - ORDER OF THE MATRIX B. (INPUT)
C                V    - NRV BY N MATRIX. (INPUT/OUTPUT)
C                         IF NRV.LE.0, V IS NOT USED. OTHERWISE, V IS
C                         REPLACED BY THE NRV BY N PRODUCT MATRIX V*VB.
C                         SEE REMARKS.
C                IV   - ROW DIMENSION OF MATRIX V EXACTLY AS SPECIFIED
C                         IN THE DIMENSION STATEMENT IN THE CALLING
C                         PROGRAM. (INPUT)
C                NRV  - NUMBER OF ROWS OF V. (INPUT)
C                C    - N BY NCC MATRIX. (INPUT/OUTPUT)
C                         IF NCC.LE.0 C IS NOT USED. OTHERWISE, C IS
C                         REPLACED BY THE N BY NCC PRODUCT
C                         MATRIX UB**(T) * C. SEE REMARKS.
C                IC   - ROW DIMENSION OF MATRIX C EXACTLY AS SPECIFIED
C                         IN THE DIMENSION STATEMENT IN THE CALLING
C                         PROGRAM. (INPUT)
C                NCC  - NUMBER OF COLUMNS IN C. (INPUT)
C                IER  - ERROR PARAMETER. (INPUT)
C                       WARNING ERROR
C                         IER=33 INDICATES THAT MATRIX B IS NOT FULL
C                           RANK OR VERY ILL-CONDITIONED. SMALL SINGULAR
C                           VALUES MAY NOT BE VERY ACCURATE.
C                       TERMINAL ERROR
C                         IER=129 INDICATES THAT CONVERGENCE WAS NOT
C                           ATTAINED AFTER 10*N QR SWEEPS. (CONVERGENCE
C                           USUALLY OCCURS IN ABOUT 2*N SWEEPS).
C   REQD. IMSL ROUTINES - SINGLE/LSVG2,VHS12,UERTST,UGETIO,VBLA=SROTG
C                       - DOUBLE/LSVG2,VHS12,UERTST,UGETIO,VBLA=DROTG
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C   REMARKS      LSVDB COMPUTES THE SINGULAR VALUE DECOMPOSITION OF
C                AN N BY N BIDIAGONAL MATRIX
C                     B = UB * S * VB**(T)       WHERE UB AND VB ARE
C                N BY N ORTHOGONAL MATRICES AND S IS DIAGONAL.
C                IF ARGUMENTS V AND C ARE N BY N IDENTITY MATRICES, ON
C                EXIT THEY ARE REPLACED BY VB AND UB**T, RESPECTIVELY.
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C-----------------------------------------------------------------------
      SUBROUTINE LSVDB (D,E,N,V,IV,NRV,C,IC,NCC,IER)
C                             SPECIFICATIONS FOR ARGUMENTS
      INTEGER          N,IV,NRV,IC,NCC,IER
      DOUBLE PRECISION D(N),E(N),V(IV,1),C(IC,1)
C                             SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER          I,II,J,K,KK,L,LL,LP1,NQRS,N10
      LOGICAL          WNTV,HAVERS,FAIL
      DOUBLE PRECISION DNORM,ZERO,ONE,TWO,CS,F,SQINF,FTEMP,G,H,HTEMP,SN,
     *                 T,X,Y,Z
      DATA             SQINF/0.8507059173023461D+38/
      DATA             ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/
C                             FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (N.LE.0) GO TO 9005
      N10 = 10*N
      WNTV = NRV.GT.0
      HAVERS = NCC.GT.0
      FAIL = .FALSE.
      NQRS = 0
      E(1) = ZERO
      DNORM = ZERO
      DO 5 J=1,N
    5 DNORM = DMAX1(DABS(D(J))+DABS(E(J)),DNORM)
      DO 100 KK=1,N
         K = N+1-KK
C                             TEST FOR SPLITTING OR RANK DEFICIENCIES
C                               FIRST MAKE TEST FOR LAST DIAGONAL TERM,
C                               D(K), BEING SMALL.
   10    IF (K.EQ.1) GO TO 25
         T = DNORM+D(K)
         IF (T.NE.DNORM) GO TO 25
C                             SINCE D(K) IS SMALL WE WILL MAKE A SPECIAL
C                               PASS TO TRANSFORM E(K) TO ZERO.
         CS = ZERO
         SN = -ONE
         DO 20 II=2,K
            I = K+1-II
            F = -SN*E(I+1)
            E(I+1) = CS*E(I+1)
            T = D(I)
            FTEMP = F
            CALL DROTG (D(I),FTEMP,CS,SN)
C                             TRANSFORMATION CONSTRUCTED TO ZERO
C                               POSITION (I,K).
            IF (.NOT.WNTV) GO TO 20
            DO 15 J=1,NRV
   15       CALL LSVG2 (CS,SN,V(J,I),V(J,K))
C                             ACCUMULATE RT. TRANSFORMATIONS IN V.
   20    CONTINUE
C                             THE MATRIX IS NOW BIDIAGONAL, AND OF
C                               LOWER ORDER SINCE E(K) .EQ. ZERO
   25    DO 30 LL=1,K
            L = K+1-LL
            T = DNORM+E(L)
            IF (T.EQ.DNORM) GO TO 50
            T = DNORM+D(L-1)
            IF (T.EQ.DNORM) GO TO 35
   30    CONTINUE
C                             THIS LOOP CANT COMPLETE SINCE E(1) = ZERO.
         GO TO 50
C                             CANCELLATION OF E(L), L.GT.1.
   35    CS = ZERO
         SN = -ONE
         DO 45 I=L,K
            F = -SN*E(I)
            E(I) = CS*E(I)
            T = DNORM+F
            IF (T.EQ.DNORM) GO TO 50
            T = D(I)
            FTEMP = F
            CALL DROTG (D(I),FTEMP,CS,SN)
            IF (.NOT.HAVERS) GO TO 45
            DO 40 J=1,NCC
   40       CALL LSVG2 (CS,SN,C(I,J),C(L-1,J))
   45    CONTINUE
C                             TEST FOR CONVERGENCE
   50    Z = D(K)
         IF (L.EQ.K) GO TO 85
C                             SHIFT FROM BOTTOM 2 BY 2 MINOR OF B**(T)*B
         X = D(L)
         Y = D(K-1)
         G = E(K-1)
         H = E(K)
         F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
         G = DABS(F)
         IF (DABS(F) .LT. SQINF) G = DSQRT(ONE+F**2)
         IF (F.LT.ZERO) GO TO 55
         T = F+G
         GO TO 60
   55    T = F-G
   60    F = ((X-Z)*(X+Z)+H*(Y/T-H))/X
C                             NEXT QR SWEEP
         CS = ONE
         SN = ONE
         LP1 = L+1
         DO 80 I=LP1,K
            G = E(I)
            Y = D(I)
            H = SN*G
            G = CS*G
            HTEMP = H
            CALL DROTG (F,HTEMP,CS,SN)
            E(I-1) = F
            F = X*CS+G*SN
            G = -X*SN+G*CS
            H = Y*SN
            Y = Y*CS
            IF (.NOT.WNTV) GO TO 70
C                             ACCUMULATE ROTATIONS (FROM THE RIGHT) IN V
            DO 65 J=1,NRV
   65       CALL LSVG2 (CS,SN,V(J,I-1),V(J,I))
            HTEMP = H
   70       CALL DROTG (F,HTEMP,CS,SN)
            D(I-1) = F
            F = CS*G+SN*Y
            X = -SN*G+CS*Y
            IF (.NOT.HAVERS) GO TO 80
            DO 75 J=1,NCC
   75       CALL LSVG2 (CS,SN,C(I-1,J),C(I,J))
C                             APPLY ROTATIONS FROM THE LEFT TO RIGHT
C                               HAND SIDES IN C
   80    CONTINUE
         E(L) = ZERO
         E(K) = F
         D(K) = X
         NQRS = NQRS+1
         IF (NQRS.LE.N10) GO TO 10
C                             RETURN TO TEST FOR SPLITTING.
         FAIL = .TRUE.
C                             CUTOFF FOR CONVERGENCE FAILURE. NQRS WILL
C                               BE 2*N USUALLY.
   85    IF (Z.GE.ZERO) GO TO 95
         D(K) = -Z
         IF (.NOT.WNTV) GO TO 95
         DO 90 J=1,NRV
   90    V(J,K) = -V(J,K)
   95    CONTINUE
C                             CONVERGENCE. D(K) IS MADE NONNEGATIVE
  100 CONTINUE
      IF (N.EQ.1) GO TO 140
      DO 105 I=2,N
         IF (D(I).GT.D(I-1)) GO TO 110
  105 CONTINUE
      GO TO 140
C                             EVERY SINGULAR VALUE IS IN ORDER
  110 DO 135 I=2,N
         T = D(I-1)
         K = I-1
         DO 115 J=I,N
            IF (T.GE.D(J)) GO TO 115
            T = D(J)
            K = J
  115    CONTINUE
         IF (K.EQ.I-1) GO TO 135
         D(K) = D(I-1)
         D(I-1) = T
         IF (.NOT.HAVERS) GO TO 125
         DO 120 J=1,NCC
            T = C(I-1,J)
            C(I-1,J) = C(K,J)
  120    C(K,J) = T
  125    IF (.NOT.WNTV) GO TO 135
         DO 130 J=1,NRV
            T = V(J,I-1)
            V(J,I-1) = V(J,K)
  130    V(J,K) = T
  135 CONTINUE
C                             END OF ORDERING ALGORITHM.
  140 IER = 129
      IF (FAIL) GO TO 9000
C                             CHECK FOR POSSIBLE RANK DEFICIENCY
      IER = 33
      T = 0.0D0
      IF (D(1).NE.ZERO) T=D(N)/D(1)
      F=100.0D0+T
      IF (F.EQ.100.0D0) GO TO 9000
      IER = 0
      GO TO 9005
 9000 CONTINUE
cc      CALL UERTST (IER,'LSVDB ')
      write (6,9999) ier
 9999 format(/' error code =',i4,' in subroutine LSVDB')     
 9005 RETURN
      END
C   IMSL ROUTINE NAME - LSVDF
C-----------------------------------------------------------------------
C   COMPUTER          - IBM77/DOUBLE
C   LATEST REVISION   - NOVEMBER 1, 1984
C   PURPOSE           - SINGULAR VALUE DECOMPOSITION OF A REAL
C                           MATRIX
C   USAGE             - CALL LSVDF (A,IA,M,N,B,IB,NB,S,WK,IER)
C   ARGUMENTS    A    - REAL M BY N MATRIX. (INPUT/OUTPUT)
C                       ON INPUT, A CONTAINS THE MATRIX TO BE
C                       ON OUTPUT, A CONTAINS THE N BY N MATRIX V IN ITS
C                         ITS FIRST N ROWS. SEE REMARKS. EITHER M.GE.N
C                         OR M.LT.N IS PERMITTED.
C                IA   - ROW DIMENSION OF MATRIX A EXACTLY AS SPECIFIED
C                         IN THE DIMENSION STATEMENT IN THE CALLING
C                         PROGRAM. (INPUT)
C                         A IS USED BY LSVDF AS WORK STORAGE FOR AN
C                         N BY N MATRIX. THEREFORE, IA MUST BE GREATER
C                         THAN OR EQUAL TO MAX(M,N).
C                M    - NUMBER OF ROWS IN A. (INPUT)
C                N    - NUMBER OF COLUMNS IN A. (INPUT)
C                B    - M BY NB MATRIX. (INPUT/OUTPUT)
C                         B IS NOT USED IF NB.LE.0. OTHERWISE, B IS
C                         REPLACED BY THE MATRIX PRODUCT U**(T) * B.
C                         SEE REMARKS.
C                IB   - ROW DIMENSION OF MATRIX B EXACTLY AS SPECIFIED
C                         IN THE DIMENSION STATEMENT IN THE CALLING
C                         PROGRAM. (INPUT)
C                NB   - NUMBER OF COLUMNS IN B. (INPUT)
C                         IF NB.LE.0, B IS NOT USED.
C                S    - VECTOR OF LENGTH N. (OUTPUT)
C                       ON OUTPUT, S CONTAINS THE ORDERED SINGULAR
C                         VALUES OF A.  S(1) .GE. S(2),..., .GE. S(N)
C                         .GE. 0.
C                WK   - WORK VECTOR OF LENGTH 2N.
C                IER  - ERROR PARAMETER. (OUTPUT)
C                       WARNING ERROR
C                         IER=33 INDICATES THAT MATRIX A IS NOT FULL
C                           RANK OR VERY ILL-CONDITIONED. SMALL SINGULAR
C                           VALUES MAY NOT BE VERY ACCURATE.
C                         IER=34 INDICATES THAT EITHER N.LE.0 OR M.LE.0.
C                       TERMINAL ERROR
C                         IER=129 INDICATES THAT CONVERGENCE WAS NOT
C                           OBTAINED BY LSVDB AND COMPUTATION WAS
C                           DISCONTINUED.
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C   REQD. IMSL ROUTINES - SINGLE/LSVDB,LSVG2,VHS12,UERSET,UERTST,UGETIO
C                           VBLA=SROTG
C                       - DOUBLE/LSVDB,LSVG2,VHS12,UERSET,UERTST,UGETIO
C                           VBLA=DROTG
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C   REMARKS  1. LSVDF COMPUTES THE SINGULAR VALUE DECOMPOSITION OF
C               A REAL M BY N MATRIX
C                    A = U * Q * V**(T)  WHERE
C               U IS AN M BY M ORTHOGONAL MATRIX,
C               V IS AN N BY N ORTHOGONAL MATRIX, AND
C               Q IS AN M BY N MATRIX WITH ALL ELEMENTS ZERO EXCEPT
C                    Q(I,I) = S(I) I=1,...,MIN(M,N).
C               V IS RETURNED IN THE FIRST N ROWS OF A.
C               U IS OBTAINED BY SETTING B TO THE M BY M IDENTITY
C               MATRIX, ON INPUT, AND SETTING NB=M. ON OUTPUT, B IS
C               REPLACED BY U**(T).
C            2. THE NOTATION U**(T) AND V**(T) REPRESENTS U TRANSPOSE
C               AND V TRANSPOSE, RESPECTIVELY. Q**(+) DENOTES THE
C               GENERALIZED INVERSE OF Q.
C            3. LSVDF IS USEFUL IN ANALYZING AND SOLVING THE LEAST
C               SQUARES PROBLEM A*X.APPR.B (WHERE .APPR. MEANS
C               APPROXIMATELY EQUALS). IN THIS CASE B IS A VECTOR OF
C               LENGTH M AND LSVDF IS CALLED WITH IB=M, NB=1. THE
C               SOLUTION IS X=V*Q**(+)*U**(T)*B. U**(T)*B REPLACES B ON
C               OUTPUT. THE SOLUTION X IS OBTAINED AS FOLLOWS... (THE
C               USER MAY WISH TO SET SMALL SINGULAR VALUES SUCH AS S(I)
C               TO ZERO JUST PRIOR TO COMPUTING X. SEE REFERENCE FOR
C               DETAILS.)
C               C                  COMPUTE Q**(+) * U**(T) * B
C                     L=MIN0(M,N)
C                     DO 10 I=1,L
C                        T=0.0
C                        IF (S(I).NE.0.0) T=B(I)/S(I)
C                        B(I)=T
C                  10 CONTINUE
C               C                  COMPUTE V * Q**(+) * U**(T) * B
C                  15 DO 25 I=1,N
C                        X(I)=0.0
C                        DO 20 J=1,L
C                  20    X(I)=X(I)+A(I,J)*B(J)
C                  25 CONTINUE
C               IF B IS SET TO THE J-TH COLUMN OF THE M BY M IDENTITY
C               MATRIX ON INPUT, X IS THE J=TH COLUMN OF A**(+) (THE
C               GENERALIZED INVERSE OF A).
C            4. THE USER SHOULD BE AWARE OF SEVERAL PRACTICAL ASPECTS OF
C               THE SINGULAR VALUE ANALYSIS. ONE OF THESE IS THE EFFECT
C               OF THE UNCERTAINTY OF THE DATA AND THE NEED FOR SCALING.
C               SEE THE LAWSON-HANSON REFERENCE PAGES 180-198 FOR
C               DETAILS.
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C-----------------------------------------------------------------------
      SUBROUTINE LSVDF  (A,IA,M,N,B,IB,NB,S,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER          IA,M,N,IB,NB,IER
      DOUBLE PRECISION A(IA,N),B(IB,1),S(N),WK(N,2)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER          I,J,JP1,K,L,MM,NN,NNP1,NS,NSP1
      DOUBLE PRECISION ZERO,ONE,T
      DATA             ZERO/0.0D0/,ONE/1.0D0/
C                             FIRST EXECUTABLE STATEMENT
      IER = 0
C                             BEGIN SPECIAL FOR ZERO ROWS AND COLS. PACK
C                               THE NONZERO COLS TO THE LEFT
      NN = N
      IER = 34
      IF (NN.LE.0.OR.M.LE.0) GO TO 9000
      IER = 0
      J = NN
    5 CONTINUE
      DO 10 I=1,M
         IF (A(I,J).NE.ZERO) GO TO 25
   10 CONTINUE
C                             COL J IS ZERO. EXCHANGE IT WITH COL N
      IF (J.EQ.NN) GO TO 20
      DO 15 I=1,M
   15 A(I,J) = A(I,NN)
   20 CONTINUE
      A(1,NN) = J
      NN = NN-1
   25 CONTINUE
      J = J-1
      IF (J.GE.1) GO TO 5
C                             IF N=0 THEN A IS ENTIRELY ZERO AND SVD
C                               COMPUTATION CAN BE SKIPPED
      NS = 0
      IF (NN.EQ.0) GO TO 120
C                             PACK NONZERO ROWS TO THE TOP. QUIT PACKING
C                               IF FIND N NONZERO ROWS
      I = 1
      MM = M
   30 IF (I.GT.N.OR.I.GE.MM) GO TO 75
      IF (A(I,I).NE.ZERO) GO TO 40
      DO 35 J=1,NN
         IF (A(I,J).NE.ZERO) GO TO 40
   35 CONTINUE
      GO TO 45
   40 I = I+1
      GO TO 30
C                             ROW I IS ZERO EXCHANGE ROWS I AND M
   45 IF (NB.LE.0) GO TO 55
      DO 50 J=1,NB
         T = B(I,J)
         B(I,J) = B(MM,J)
         B(MM,J) = T
   50 CONTINUE
   55 DO 60 J=1,NN
   60 A(I,J) = A(MM,J)
      IF (MM.GT.NN) GO TO 70
      DO 65 J=1,NN
   65 A(MM,J) = ZERO
   70 CONTINUE
C                             EXCHANGE IS FINISHED
      MM = MM-1
      GO TO 30
   75 CONTINUE
C                             END SPECIAL FOR ZERO ROWS AND COLUMNS
C                             BEGIN SVD ALGORITHM..
C                             (1) REDUCE THE MATRIX TO UPPER BIDIAGONAL
C                               FORM WITH HOUSEHOLDER TRANSFORMATIONS
C                               H(N)...H(1)AQ(1)...Q(N-2) = (D**T,0)**T
C                               WHERE D IS UPPER BIDIAGONAL.
C                             (2) APPLY H(N)...H(1) TO B. HERE
C                               H(N)...H(1)*B REPLACES B IN STORAGE.
C                             (3) THE MATRIX PRODUCT W=Q(1)...Q(N-2)
C                               OVERWRITES THE FIRST N ROWS OF A IN
C                               STORAGE.
C                             (4) AN SVD FOR D IS COMPUTED. HERE K
C                               ROTATIONS RI AND PI ARE COMPUTED SO THAT
C                               RK...R1*D*P1**(T)...PK**(T) =
C                               DIAG(S1,...,SM) TO WORKING ACCURACY. THE
C                               SI ARE NONNEGATIVE AND NONINCREASING.
C                               HERE RK...R1*B OVERWRITES B IN STORAGE
C                               WHILE A*P1**(T)...PK**(T) OVERWRITES A
C                               IN STORAGE.
C                             (5) IT FOLLOWS THAT,WITH THE PROPER
C                               DEFINITIONS, U**(T)*B OVERWRITES B,
C                               WHILE V OVERWRITES THE FIRST N ROW
C                               AND COLUMNS OF A.
      L = MIN0(MM,NN)
C                             THE FOLLOWING LOOP REDUCES A TO UPPER
C                               BIDIAGONAL AND ALSO APPLIES THE
C                               PREMULTIPLYING TRANSFORMATIONS TO B.
      DO 85 J=1,L
         IF (J.GE.MM) GO TO 80
         JP1 = MIN0(J+1,NN)
         CALL VHS12 (1,J,J+1,MM,A(1,J),1,T,A(1,JP1),1,IA,NN-J)
         CALL VHS12 (2,J,J+1,MM,A(1,J),1,T,B,1,IB,NB)
   80    IF (J.GE.NN-1) GO TO 85
         CALL VHS12 (1,J+1,J+2,NN,A(J,1),IA,WK(J,2),A(J+1,1),IA,1,MM-J)
   85 CONTINUE
C                             COPY THE BIDIAGONAL MATRIX INTO THE ARRAY
C                               S FOR LSVDB
      IF (L.EQ.1) GO TO 95
      DO 90 J=2,L
         S(J) = A(J,J)
         WK(J,1) = A(J-1,J)
   90 CONTINUE
   95 S(1) = A(1,1)
      NS = NN
      IF (MM.GE.NN) GO TO 100
      NS = MM+1
      S(NS) = ZERO
      WK(NS,1) = A(MM,MM+1)
  100 CONTINUE
C                             CONSTRUCT THE EXPLICIT N BY N PRODUCT
C                               MATRIX, W=Q1*Q2*...*QL*I IN THE ARRAY A
      DO 115 K=1,NN
         I = NN+1-K
         IF (I.GT.MIN0(MM,NN-2)) GO TO 105
         CALL VHS12 (2,I+1,I+2,NN,A(I,1),IA,WK(I,2),A(1,I+1),1,IA,NN-I)
  105    DO 110 J=1,NN
  110    A(I,J) = ZERO
         A(I,I) = ONE
  115 CONTINUE
C                             COMPUTE THE SVD OF THE BIDIAGONALMATRIX
cc      LEVEL=1
cc      CALL UERSET(LEVEL,LEVOLD)
      CALL LSVDB (S(1),WK(1,1),NS,A,IA,NN,B,IB,NB,IER)
C                             TEST FOR IER=33
      IF (IER.GT.128) GO TO 9000
cc      CALL UERSET(LEVOLD,LEVOLD)
      IF (IER.NE.33) GO TO 120
      T=0.0D0
      NM=MIN0(M,N)
      IF (S(1).NE.ZERO) T=S(NM)/S(1)
      F=100.0D0+T
      IF (F.EQ.100.0D0) GO TO 120
      IER=0
  120 CONTINUE
      IF (NS .LT. MIN0(M,N)) IER = 33
      IF (NS.GE.NN) GO TO 130
      NSP1 = NS+1
      DO 125 J=NSP1,NN
  125 S(J) = ZERO
  130 CONTINUE
      IF (NN.EQ.N) GO TO 155
      NNP1 = NN+1
C                             MOVE RECORD OF PERMUTATIONS AND STORE
C                               ZEROS
      DO 140 J=NNP1,N
         S(J) = A(1,J)
         IF (NN.LT.1) GO TO 140
         DO 135 I=1,NN
  135    A(I,J) = ZERO
  140 CONTINUE
C                             PERMUTE ROWS AND SET ZERO SINGULAR VALUES
      DO 150 K=NNP1,N
         I = S(K)
         S(K) = ZERO
         DO 145 J=1,N
            A(K,J) = A(I,J)
  145    A(I,J) = ZERO
         A(I,K) = ONE
  150 CONTINUE
C                             END SPECIAL FOR ZERO ROWS AND COLUMNS
  155 IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
cc      CALL UERTST (IER,'LSVDF ')
      write (6,9999) ier
 9999 format(/' error code =',i4,' in subroutine LSVDF')     
 9005 RETURN
      END
C   IMSL ROUTINE NAME   - LSVG2
C-----------------------------------------------------------------------
C   COMPUTER            - IBM77/DOUBLE
C   LATEST REVISION     - JULY 1, 1983
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE LSVDB
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C   REQD. IMSL ROUTINES - NONE REQUIRED
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C-----------------------------------------------------------------------
      SUBROUTINE LSVG2  (CS,SN,X,Y)
C                             SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION CS,SN,X,Y
C                             SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION XR
C                             FIRST EXECUTABLE STATEMENT
      XR=CS*X+SN*Y
      Y=-SN*X+CS*Y
      X=XR
      RETURN
      END
C   IMSL ROUTINE NAME - VBLA=DROTG
C-----------------------------------------------------------------------
C   COMPUTER          - IBM77/DOUBLE
C   LATEST REVISION   - NOVEMBER 1, 1984
C   PURPOSE           - CONSTRUCT GIVENS PLANE ROTATION
C                         (DOUBLE PRECISION)
C   USAGE             - CALL DROTG (DA,DB,DC,DS)
C   ARGUMENTS    DA   - FIRST ELEMENT OF DOUBLE PRECISION VECTOR.
C                         (INPUT/OUTPUT)
C                       ON OUTPUT, R=(+/-)DSQRT(DA**2 + DB**2)
C                         OVERWRITES DA.
C                DB   - SECOND ELEMENT OF DOUBLE PRECISION VECTOR.
C                         (INPUT/OUTPUT)
C                       ON OUTPUT, Z OVERWRITES DB.
C                         Z IS DEFINED TO BE..
C                         DS,     IF DABS(DA).GT.DABS(DB)
C                         1.0D0/DC, IF DABS(DB).GE.DABS(DA) AND
C                                   DC.NE.0.0D0
C                         1.0D0,    IF DC.EQ.0.0D0.
C                DC   - DOUBLE PRECISION ELEMENT OF OUTPUT
C                         TRANSFORMATION MATRIX. SEE REMARKS.
C                DS   - DOUBLE PRECISION ELEMENT OF OUTPUT
C                         TRANSFORMATION MATRIX. SEE REMARKS.
C   PRECISION/HARDWARE  - DOUBLE/ALL
C   REQD. IMSL ROUTINES - NONE REQUIRED
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C   REMARKS      DROTG CONSTRUCTS THE GIVENS TRANSFORMATION
C                    ( DC  DS )
C                G = (        ) ,  DC**2 + DS**2 = 1 ,
C                    (-DS  DC )
C                WHICH ZEROS THE SECOND ELEMENT OF (DA,DB)**T.
C   COPYRIGHT         - 1984 BY IMSL, INC. ALL RIGHTS RESERVED.
C   WARRANTY          - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                         APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                         EXPRESSED OR IMPLIED, IS APPLICABLE.
C-----------------------------------------------------------------------
      SUBROUTINE DROTG (DA,DB,DC,DS)
C                             SPECIFICATIONS FOR ARGUMENTS
      DOUBLE PRECISION DA,DB,DC,DS
C                             SPECIFICATIONS FOR LOCAL VARIABLES
      DOUBLE PRECISION U,V,R
C                             FIRST EXECUTABLE STATEMENT
      IF (DABS(DA).LE.DABS(DB)) GO TO 5
C                             HERE DABS(DA) .GT. DABS(DB)
      U = DA+DA
      V = DB/U
C                             NOTE THAT U AND R HAVE THE SIGN OF DA
      R = DSQRT(.25D0+V**2)*U
C                             NOTE THAT DC IS POSITIVE
      DC = DA/R
      DS = V*(DC+DC)
      DB = DS
      DA = R
      RETURN
C                             HERE DABS(DA) .LE. DABS(DB)
    5 IF (DB.EQ.0.D0) GO TO 15
      U = DB+DB
      V = DA/U
C                             NOTE THAT U AND R HAVE THE SIGN OF DB
C                               (R IS IMMEDIATELY STORED IN DA)
      DA = DSQRT(.25D0+V**2)*U
C                             NOTE THAT DS IS POSITIVE
      DS = DB/DA
      DC = V*(DS+DS)
      IF (DC.EQ.0.D0) GO TO 10
      DB = 1.D0/DC
      RETURN
   10 DB = 1.D0
      RETURN
C                             HERE DA = DB = 0.D0
   15 DC = 1.D0
      DS = 0.D0
      DA = 0.D0
      DB = 0.D0
      RETURN
      END
C   IMSL ROUTINE NAME - VHS12
C-----------------------------------------------------------------------
C   COMPUTER          - IBM77/DOUBLE
C   LATEST REVISION   - JANUARY 1, 1978
C   PURPOSE           - REAL HOUSEHOLDER TRANSFORMATION -
C                         COMPUTATION AND APPLICATIONS.
C   USAGE             - CALL VHS12 (MODE,LP,L1,M,U,INCU,UP,C,INCC,
C                         ICV,NCV)
C   ARGUMENTS    MODE - OPTION PARAMETER. (INPUT)
C                       IF MODE=1, THE SUBROUTINE COMPUTES A HOUSEHOLDER
C                         TRANSFORMATION AND IF NCV.GT.0, MULTIPLIES IT
C                         BY THE SET OF NCV VECTORS (EACH OF LENGTH M)
C                         STORED IN C. FOR A GIVEN VECTOR V OF LENGTH M
C                         AND TWO INTEGER INDICES LP AND L1 THAT SATISFY
C                         1 .LE. LP .LT. L1 .LE. M, THE SUBROUTINE
C                         DEFINES AN M BY M HOUSEHOLDER TRANSFORMATION Q
C                         WHICH SATIFIES QV=W WHERE
C                         W(I)=V(I) FOR I.LT.LP
C                         W(LP)=-SIG*SQRT(V(LP)**2+V(L1)**2+...
C                           +V(M)**2)
C                           SIG=1  IF V(LP).GE.0
C                           SIG=-1 IF V(LP).LT.0
C                         W(I)=V(I) FOR LP.LT.I.LT.L1
C                         W(I)=0    FOR I.GE.L1.
C                       IF MODE=2, THE SUBROUTINE ASSUMES THAT A
C                         HOUSEHOLDER TRANSFORMATION HAS ALREADY BEEN
C                         DEFINED BY A PREVIOUS CALL WITH MODE=1, AND IF
C                         NCV.GT.0, MULTIPLIES IT BY THE SET OF NCV
C                         VECTORS (EACH OF LENGTH M) STORED IN C.
C                LP   - PARAMETERS THAT DEFINE THE DESIRED
C                L1       HOUSEHOLDER TRANSFORMATION. (INPUT)
C                M        IF THE CONDITION 1.LE.LP.LT.L1.LE.M IS NOT
C                         SATISFIED, THE SUBROUTINE RETURNS TO THE
C                         CALLING PROGRAM WITHOUT PERFORMING ANY
C                         COMPUTATIONS.
C                U    - VECTOR OF M ELEMENTS. (INPUT, AND OUTPUT IF
C                         MODE=1)
C                         THE STORAGE INCREMENT BETWEEN ELEMENTS OF U IS
C                         INCU. (I.E., U(1+(J-1)*INCU), J=1,...,M). IF
C                         MODE=1, THE ARRAY V IS DEFINED AS
C                         V(J)=U(1+(J-1)*INCU), J=1,...,M.
C                       ON OUTPUT, U(1+(LP-1)*INCU) IS SET TO W(LP) (AS
C                         DEFINED ABOVE IN THE DESCRIPTION OF MODE=1).
C                INCU - INCREMENT BETWEEN ELEMENTS OF U. (INPUT)
C                UP   - SCALAR SET TO V(LP)-W(LP) TO DEFINE THE
C                         HOUSEHOLDER TRANSFORMATION Q. (INPUT IF
C                         MODE=2, OUTPUT IF MODE=1)
C                C    - VECTOR OF NCV*M ELEMENTS. (INPUT/OUTPUT)
C                         IF NCV.LE.0, C IS NOT USED.
C                         IF NCV.GT.0, C CONTAINS NCV VECTORS OF LENGTH
C                         M WITH INCREMENT INCC BETWEEN ELEMENTS OF
C                         VECTORS AND INCREMENT ICV BETWEEN VECTORS.
C                         ELEMENT I OF VECTOR J IS DEFINED AS
C                         C(1+(I-1)*INCC+(J-1)*ICV), I=1,...,M AND
C                         J=1,...,NCV.
C                       ON OUTPUT, C CONTAINS THE SET OF NCV VECTORS
C                         RESULTING FROM MULTIPLYING THE GIVEN VECTORS
C                         BY Q.
C                INCC - INCREMENT BETWEEN ELEMENTS OF VECTORS IN C.
C                         (INPUT)
C                ICV  - INCREMENT BETWEEN VECTORS IN C. (INPUT)
C                NCV  - NUMBER OF VECTORS STORED IN C. (INPUT)
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C   REQD. IMSL ROUTINES - NONE REQUIRED
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C   REMARKS  1. IF U IS A SINGLE SUBSCRIPTED ARRAY OR THE J-TH COLUMN OF
C               A MATRIX, THEN INCU=1. IF U IS THE I-TH ROW OF A MATRIX
C               THEN INCU IS THE ROW DIMENSION OF THE MATRIX EXACTLY AS
C               SPECIFIED IN THE CALLING PROGRAM.
C            2. IF C IS A DOUBLE SUBSCRIPTED MATRIX AND THE VECTORS ARE
C               THE FIRST NCV COLUMNS OF C, THEN INCC=1 AND ICV IS THE
C               ROW DIMENSION OF C EXACTLY AS SPECIFIED IN THE CALLING
C               PROGRAM. IN THIS CASE C IS REPLACED BY QC. IF THE
C               VECTORS ARE SUCCESSIVE ROWS OF C THEN INCC IS THE ROW
C               DIMENSION OF C EXACTLY AS SPECIFIED IN THE CALLING
C               PROGRAM AND ICV=1. IN THIS CASE C IS REPLACED BY CQ.
C   COPYRIGHT         - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C   WARRANTY          - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                         APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                         EXPRESSED OR IMPLIED, IS APPLICABLE.
C-----------------------------------------------------------------------
      SUBROUTINE VHS12  (MODE,LP,L1,M,U,INCU,UP,C,INCC,ICV,NCV)
C                             SPECIFICATIONS FOR ARGUMENTS
      INTEGER          MODE,LP,L1,M,INCU,INCC,ICV,NCV
      DOUBLE PRECISION U(1),UP,C(1)
C                             SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER          IJ,ILP,IL1,IM,INCR,I2,I3,I4,J
      DOUBLE PRECISION SM,B
      DOUBLE PRECISION ONE,CL,CLINV,SM1
C                             FIRST EXECUTABLE STATEMENT
      ONE = 1.D0
      IF (0.GE.LP.OR.LP.GE.L1.OR.L1.GT.M) GO TO 9005
      ILP = (LP-1)*INCU+1
      IL1 = (L1-1)*INCU+1
      IM = (M-1)*INCU+1
      CL = DABS(U(ILP))
      IF (MODE.EQ.2) GO TO 15
C                             CONSTRUCT THE TRANSFORMATION.
      DO 5 IJ=IL1,IM,INCU
    5 CL = DMAX1(DABS(U(IJ)),CL)
      IF (CL.LE.0.0D0) GO TO 9005
      CLINV = ONE/CL
      SM = (U(ILP)*CLINV)**2
      DO 10 IJ=IL1,IM,INCU
   10 SM = SM+(U(IJ)*CLINV)**2
C                             CONVERT DBLE. PREC. SM TO SNGL. PREC. SM1
      SM1 = SM
      CL = CL*DSQRT(SM1)
      IF (U(ILP).GT.0.0D0) CL = -CL
      UP = U(ILP)-CL
      U(ILP) = CL
      GO TO 20
C                             APPLY THE TRANSFORMATION I+U*(U**T)/B TO C
   15 IF (CL.LE.0.0D0) GO TO 9005
   20 IF (NCV.LE.0) GO TO 9005
      B = UP*U(ILP)
C                             B MUST BE NONPOSITIVE HERE. IF B = 0.,
C                               RETURN.
      IF (B.GE.0.0D0) GO TO 9005
      B = ONE/B
      I2 = 1-ICV+INCC*(LP-1)
      INCR = INCC*(L1-LP)
      DO 35 J=1,NCV
         I2 = I2+ICV
         I3 = I2+INCR
         I4 = I3
         SM = C(I2)*UP
         DO 25 IJ=IL1,IM,INCU
            SM = SM+C(I3)*U(IJ)
            I3 = I3+INCC
   25    CONTINUE
         IF (SM.EQ.0.0D0) GO TO 35
         SM = SM*B
         C(I2) = C(I2)+SM*UP
         DO 30 IJ=IL1,IM,INCU
            C(I4) = C(I4)+SM*U(IJ)
            I4 = I4+INCC
   30    CONTINUE
   35 CONTINUE
 9005 RETURN
      END
      SUBROUTINE SEIGCX(A,N,ND,VAL,K1,K2,TOL,NVEC,VEC,C)
c     modified 03/28/01
C  DIAGONALIZES FULL HERMITIAN MATRIX A, USES ONLY LOWER TRIANGLE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(ND),VAL(1),A(2,ND,1),C(1)
      IF (N.GT.1) GOTO 10
      VAL(1)=A(1,1,1)
      IF (NVEC.GT.0) VEC(1)=1D0
      IF (NVEC.GT.0) VEC(2)=0D0
      GOTO 20
   10 CALL HOUSCX(A,N,ND,C,C(ND+1),C(2*ND+1),C(3*ND+1),C(4*ND+1))
      CALL BISECT(C,C(ND+1),C(2*ND+1),VAL,N,K1,K2,TOL,ERRBND,C(6*ND+1))
      TOL=ERRBND
      IF(NVEC.LE.0) GOTO 20
      CALL INVITR(C,C(ND+1),N,VAL,NVEC,VEC,ND,C(2*ND+1),C(7*ND+1),
     *  C(8*ND+1),C(9*ND+1),C(6*ND+1),C(10*ND+1))
      K=NVEC*ND+1
      NDMN=ND-N
      DO 1 J=1,NVEC
      K=K-NDMN
      DO 1 I=1,N
      K=K-1
      VEC(2*K-1)=VEC(K)
    1 VEC(2*K)=0D0
      CALL REVCX(A,N,ND,VEC,NVEC,C(4*ND+1),C(3*ND+1))
   20 RETURN
      END
      SUBROUTINE HOUSCX(A,N,ND,C,B,B2,H,TAU)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(2,ND,1),TAU(2,1),C(1),B(1),B2(1),H(1)
c TOL = smallest floating point number = .2225073858507201-307 = 16**(-255)/4 = 2**(-1022)              
c (at usc) TOL=/z0d10000000000000/=2**(-208)=16**(-52)=.24308653429145084800d-62
      DATA TOL,ZERO,ONE/1d-100,0D0,1D0/
      DO 11 K=2,N
      KM1=K-1
      Q1=A(1,K,KM1)**2+A(2,K,KM1)**2
      Q2=DSQRT(Q1)
      IF(Q1.GT.0.D0)GO TO 1
      TAU(1,KM1)=-ONE
      TAU(2,KM1)=ZERO
      GO TO 2
    1 TAU(1,KM1)=-A(1,K,KM1)/Q2
      TAU(2,KM1)=-A(2,K,KM1)/Q2
    2 IF(K.EQ.N)GO TO 4
      V=ZERO
      DO 3 I=K+1,N
    3 V=V+A(1,I,KM1)**2+A(2,I,KM1)**2
      IF(V.GT.TOL)GO TO 5
      DELTA=ZERO
    4 H(KM1)=ZERO
      B2(KM1)=Q1
      B(KM1)=Q2
      TAU(1,KM1)=-TAU(1,KM1)
      TAU(2,KM1)=-TAU(2,KM1)
      GO TO 11
    5 V=V+Q1
      Q1=DSQRT(V)
      DELTA=V+Q1*Q2
      H(KM1)=DELTA
      B2(KM1)=V
      B(KM1)=Q1
      A(1,K,KM1)=A(1,K,KM1)-Q1*TAU(1,KM1)
      A(2,K,KM1)=A(2,K,KM1)-Q1*TAU(2,KM1)
      RHO=ZERO
      DO 9 I=K,N
      Q1=ZERO
      Q2=ZERO
      DO 6 J=K,I
      Q1=Q1+A(1,I,J)*A(1,J,KM1)-A(2,I,J)*A(2,J,KM1)
    6 Q2=Q2+A(1,I,J)*A(2,J,KM1)+A(2,I,J)*A(1,J,KM1)
      IF(I.EQ.N)GO TO 8
      DO 7 J=I+1,N
      Q1=Q1+A(1,J,I)*A(1,J,KM1)+A(2,J,I)*A(2,J,KM1)
    7 Q2=Q2+A(1,J,I)*A(2,J,KM1)-A(2,J,I)*A(1,J,KM1)
    8 TAU(1,I)=Q1/DELTA
      TAU(2,I)=Q2/DELTA
    9 RHO=RHO+Q1*A(1,I,KM1)+Q2*A(2,I,KM1)
      RHO=RHO/DELTA**2
      DO 10 I=K,N
      T1=A(1,I,KM1)
      T2=A(2,I,KM1)
      Q1=TAU(1,I)-RHO*T1
      Q2=TAU(2,I)-RHO*T2
      DO 10 J=K,I
      A(1,I,J)=A(1,I,J)-T1*TAU(1,J)-T2*TAU(2,J)-Q1*A(1,J,KM1)-Q2*
     * A(2,J,KM1)
   10 A(2,I,J)=A(2,I,J)-T2*TAU(1,J)+T1*TAU(2,J)-Q2*A(1,J,KM1)+Q1*
     * A(2,J,KM1)
   11 CONTINUE
      B2(N)=ZERO
      B(N)=ZERO
      DO 12 I=1,N
   12 C(I)=A(1,I,I)
      RETURN
      END
      SUBROUTINE REVCX(A,N,ND,X,NVEC,TAU,DELTA)
c     modified to avoid problems with the complex arithmetic 
      IMPLICIT REAL*8(A-H,O-Z)
c      COMPLEX*16 A,X,TAU,T,DCONJG
c      DIMENSION X(ND,1),A(ND,1),TAU(1),DELTA(1)
      DIMENSION A(2,ND,1),X(2,ND,1),TAU(2,1),DELTA(1)
      DATA ZERO,ONE/0D0,1D0/
      T1=ONE
      T2=ZERO
      DO 1 I=2,N
      X1=T1
      T1=TAU(1,I-1)*X1-TAU(2,I-1)*T2
      T2=TAU(1,I-1)*T2+TAU(2,I-1)*X1 
      DO 1 J=1,NVEC
      X1=X(1,I,J)
      X(1,I,J)=T1*X1-T2*X(2,I,J)
    1 X(2,I,J)=T1*X(2,I,J)+T2*X1
      IF(N.EQ.2)GOTO 5
      NM2=N-2
      DO 4 M=1,NM2
      K=N-M
      IF(DELTA(K-1).EQ.ZERO)GO TO 4
      DO 3 J=1,NVEC
      T1=ZERO
      T2=ZERO
      DO 2 I=K,N
      T1=T1+A(1,I,K-1)*X(1,I,J)+A(2,I,K-1)*X(2,I,J)
    2 T2=T2+A(1,I,K-1)*X(2,I,J)-A(2,I,K-1)*X(1,I,J)
      T1=T1/DELTA(K-1)
      T2=T2/DELTA(K-1)
      DO 3 I=K,N
      X(1,I,J)=X(1,I,J)-A(1,I,K-1)*T1+A(2,I,K-1)*T2
    3 X(2,I,J)=X(2,I,J)-A(1,I,K-1)*T2-A(2,I,K-1)*T1
    4 CONTINUE
    5 RETURN
      END
      SUBROUTINE INVITR(C,B,N,VAL,NVEC,VEC,ND,P,Q,R,U,W,CIN)
c  modified 03/28/01
C  MODIFIED TO YIELD ORTHOGONAL EIGENVECTORS FOR DEGENERATE EIGENVALUES
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 R4
      DIMENSION VEC(ND,1),VAL(1),C(1),B(1),P(1),Q(1),R(1),U(1),W(1),
     * CIN(1)
c   EPS = machine precision =   .2220446049250313D-15 = 16**(-13) =  2**(-52)
c  (at usc) EPS=/z3410000000000000/=2**(-52)=16**(-13)=.22204460492503130800d-15
      DATA EPS,ZERO,ONE,TWO,THREE/.2221d-15,0D0,1D0,2D0,3D0/
   11 NM1=N-1
      H=DMAX1(DABS(C(1))+DABS(B(1)),DABS(B(NM1))+DABS(C(N)))
      DO 12 I=2,NM1
   12 H=DMAX1(H,DABS(C(I))+DABS(B(I-1))+DABS(B(I)))
      H=H*EPS
      HT3=THREE*H
      B(N)=ZERO
      NP1=N+1
      DO 32 I=1,NVEC
      IF(I.EQ.1)GO TO 13
      IF(DABS(VAL(I)-XMULT).GT.HT3)GO TO 13
      MULT=MULT+1
      T=(MULT+1)/2
      IF(2*(MULT/2).EQ.MULT)T=-T
      T=XMULT+T*HT3
      GO TO 14
   13 MULT=1
      XMULT=VAL(I)
      T=XMULT
   14 DO 15 J=1,N
      P(J)=ZERO
      Q(J)=B(J)
      R(J)=C(J)-T
   15 U(J)=ONE
      IF(MULT.EQ.1)GO TO 34
      DO 33 J=1,N
      CALL RANDOM(R4)
   33 U(J)=TWO*R4-ONE
   34 DO 18 J=1,NM1
      IF(DABS(R(J)).LT.DABS(B(J)))GO TO 16
      IF(R(J).EQ.ZERO)R(J)=H
      XR=B(J)/R(J)
      CIN(J)=-ONE
      GO TO 17
   16 XR=R(J)/B(J)
      CIN(J)=ONE
      R(J)=B(J)
      XM=R(J+1)
      R(J+1)=Q(J)
      Q(J)=XM
      P(J)=Q(J+1)
      Q(J+1)=ZERO
   17 W(J)=XR
      Q(J+1)=Q(J+1)-XR*P(J)
   18 R(J+1)=R(J+1)-XR*Q(J)
      IF(R(N).EQ.ZERO)R(N)=H
      IT=0
      GO TO 23
   20 DO 22 J=1,NM1
      IF (CIN(J).EQ.ONE) GOTO 21
      U(J+1)=U(J+1)-W(J)*U(J)
      GO TO 22
   21 XR=U(J)
      U(J)=U(J+1)
      U(J+1)=XR-W(J)*U(J+1)
   22 CONTINUE
   23 U(N)=U(N)/R(N)
      U(NM1)=(U(NM1)-U(N)*Q(NM1))/R(NM1)
      IF(N.EQ.2)GO TO 25
      DO 24 J=3,N
      K=NP1-J
   24 U(K)=(U(K)-U(K+1)*Q(K)-U(K+2)*P(K))/R(K)
   25 IF(IT.NE.0)XM=ONE/U(M)
      XR=DABS(U(1))
      M=1
      DO 26 J=2,N
      XL=DABS(U(J))
      IF(XR.GE.XL)GO TO 26
      XR=XL
      M=J
   26 CONTINUE
      XR=ONE/U(M)
      IF(IT.EQ.0)GO TO 28
      T=ZERO
      DO 27 J=1,N
      XL=DABS(VEC(J,I)-XM*U(J))
   27 IF(XL.GT.T)T=XL
      IF(T*DABS(XR).LE.H)GO TO 300
   28 DO 29 J=1,N
      U(J)=XR*U(J)
   29 VEC(J,I)=U(J)
      IT=IT+1
      IF(IT.LT.10)GO TO 20
      WRITE(6,100)I
  100 FORMAT('0NO CONVERGENCE TOOK PLACE IN COMPUTING EIGENVECTOR',I3)
  300 IF (MULT.EQ.1) GOTO 30
      DO 302 K=I-MULT+1,I-1
      T=ZERO
      DO 301 J=1,N
  301 T=T+U(J)*VEC(J,K)
      DO 302 J=1,N
  302 U(J)=U(J)-VEC(J,K)*T
   30 T=ZERO
      DO 31 J=1,N
   31 T=T+U(J)**2
      T=ONE/DSQRT(T)
      DO 32 J=1,N
   32 VEC(J,I)=T*U(J)
   40 RETURN
      END
      SUBROUTINE BISECT(C,B,B2,VAL,N,K1,K2,TOL,ERRBND,W)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(1),B(1),B2(1),W(1),VAL(1)
c   EPS = machine precision =   .2220446049250313D-15 = 16**(-13) =  2**(-52)
c  (at usc) EPS=/z3410000000000000/=2**(-52)=16**(-13)=.22204460492503130800d-15
c      DATA EPS,ZERO,TWO/Z3410000000000000,0D0,2D0/
      DATA EPS,ZERO,TWO/.2221d-15,0D0,2D0/
      NM1=N-1
      H=DABS(B(1))
      XM=DABS(B(NM1))
      XMAX=DMAX1(C(1)+H,C(N)+XM)
      XMIN=DMIN1(C(1)-H,C(N)-XM)
      DO 1 I=2,NM1
      H=DABS(B(I-1))+DABS(B(I))
      XMAX=DMAX1(XMAX,C(I)+H)
    1 XMIN=DMIN1(XMIN,C(I)-H)
      H=EPS*DMAX1(XMAX,-XMIN)
      ERRBND=DMAX1(7.5*H,TOL)
      EPS1=TWO*(ERRBND-7.*H)
      K2P1=K2+1
      NVAL=K2P1-K1
      DO 2 I=1,NVAL
      VAL(I)=XMAX
    2 W(I)=XMIN
      XR=XMAX
      DO 11 K=K1,K2
      XL=XMIN
      IND=K2P1-K
      DO 3 I=1,IND
      IF(XL.GE.W(I))GO TO 3
      XL=W(I)
      GO TO 4
    3 CONTINUE
    4 IF(XR.GT.VAL(IND))XR=VAL(IND)
    5 XM=(XL+XR)/TWO
      IF(XR-XL.LE.(EPS+EPS)*(DABS(XR)+DABS(XL))+EPS1)GO TO 11
      M=0
      T=C(1)-XM
      IF(T.GE.ZERO)M=1
      DO 8 I=2,N
      IF(T.NE.ZERO)GO TO 6
      H=DABS(B(I-1))/EPS
      GO TO 7
    6 H=B2(I-1)/T
    7 T=C(I)-XM-H
    8 IF(T.GE.ZERO)M=M+1
      IF(M.GE.K)GO TO 9
      XR=XM
      GO TO 5
    9 XL=XM
      IF(M.LT.K2)GO TO 10
      W(1)=XM
      GO TO 5
   10 I=K2P1-M
      W(I)=XM
      IF(VAL(I-1).GT.XM)VAL(I-1)=XM
      GO TO 5
   11 VAL(IND)=XM
      RETURN
      END
