C   Last modified 08/06/2005

C   This code computes Fourier components of the reflection
C   function for a semi-infinite particulate slab characterized by the single
C   scattering albedo ALB and coefficients appearing in the expansion
C   of the phase function in Legendre polynomials (file 10).
C   The computation is based on the rigorous solution of the radiative
C   transfer equation as described in the paper
C                                                                  
C   M. I. MISHCHENKO, J. M. DLUGACH, E. G. YANOVITSKIJ,          
C   AND N. T. ZAKHAROVA, BIDIRECTIONAL REFLECTANCE OF FLAT,       
C   OPTICALLY THICK PARTICULATE LAYERS:  AN EFFICIENT RADIATIVE      
C   TRANSFER SOLUTION AND APPLICATIONS TO SNOW AND SOIL SURFACES,        
C   J. QUANT. SPECTROSC. RADIAT. TRANSFER, VOL. 63, 409-432 (1999).
C
C   The paper is available upon reuest from Michael Mishchenko
C   (send a mesage to michael.i.mishchenko@nasa.gov) and 
C   is also available in the .pdf format at
C   http://www.giss.nasa.gov/staff/mmishchenko/publications/index.html                                                                       
C   The code was developed by Michael Mishchenko and 
C   Nadia Zakharova at the NASA      
C   Goddard Institute for Space Studies, New York.  The development    
C   of the code was supported by the NASA Radiation Science Program.  
                                                                       
C   The code can be used without limitations in any not-for-       
C   profit scientific research.  The only request is that in any       
C   publication using the code the source of the code be acknowledged  
C   and relevant references be made.                                   
                                                                       
C   We would highly appreciate informing us of any problems encountered 
C   with this code.  Please send your message to Michael Mishchenko         
C   at michael.i.mishchenko@nasa.gov.                              
                                                                   
C   MMAX1 Fourier components of the reflection function are computed 
C   sequentially and each of them is given by the array R(I,J), where
C   I numbers the values of the cosine of the reflection zenith
C   angle (mu) and J numbers the values of the cosine of the illumination
C   zenith angle (mu_0).  Both mu and mu_0 take the values of the 
C   division points of a quadrature formula of order NG
C   on the interval (0,1].  The Fourier components can be used to compute  
C   the reflection function for any azimuth angle phi via Eq. (2), in
C   which m_max=MMAX1-1.  To compute the reflection function for
C   mu and/or mu_0 values different from the quadrature division points
C   on the interval (0,1] used to solve the radiative transfer equation,
C   an interpolation procedure must be used.  The interpolation code is
C   available at the same web site.
C   The Fourier components of the reflection function are writen in file 11.  
 
C   By default, MMAX1 is equal to the number of numerically
C   significant coefficients in the Legendre expansion of the
C   phase function [i.e., MMAX1=s_max+1, where s_max is the upper
C   summation limit in Eq. (6)].  However, if the specified value of
C   the parameter NSEP is 1, then the first-order-scattering
C   contribution to the Fourier components of the reflection function is 
C   separated (Section 3.7 of the above paper) and only the first MMAX2
C   modified Fourier components of the reflection function are computed
C   [MMAX2=m_1+1, where m_1 is the upper summation limit in Eq. (47)].
 
C   The code also computes the spherical (Bond) albedo and the values
C   of the plane albedo for MU_0 given by the quadrature division points
C   [Eqs. (12) and (11) of the above paper].
C   If only these albedos are required (not the complete reflection
C   function), then the line
C
C     DO M1=1,MMAX1
C
C   in the main program should be modified as follows:
C
C     DO M1=1,1
C
C   This will save a major fraction of CPU time.
 
C   The parameter NG should be increased simultaneously in all
C   PARAMETER statements until the reflection function values converge
C   within the desired accuracy.  Messages 
C
C     WARNING:  I(MU) HAS NOT CONVERGED.
C
C   and
C
C     WARNING:  NUMERICAL CONERGENCE HAS NOT BEEN REACHED.
C     TRY A LARGER NG.
C
C   may be indicators that the NG value chosen is too small and should
C   be increased.

C   WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
C   IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
C   INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
C   VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
C   THE AUTHORS AND THEIR ORGANIZATION DISCLAIM ALL LIABILITY FOR
C   ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM. 

      SUBROUTINE BRF(ALB, MMAX1, AL1, NSEP,
     *   IERR, X, GLOBAL, ALBEDO, RZERO_, R_)
!f2py intent(in) AL1, MMAX1, AL1, NSEP
!f2py intent(out) IERR, X, GLOBAL, ALBEDO, R_, RZERO_

      IMPLICIT REAL*8 (A-H,O-Z)                                         
      PARAMETER (NG=100, LMAX1=700)
      REAL*8 X(NG),W(NG),WW(NG),AL1(LMAX1),X1(NG),W1(NG)
      REAL*4 R(NG,NG),ALBEDO(NG),RZERO(NG,NG)
      REAL*4 R_(LMAX1,NG,NG),RZERO_(NG,NG)
      COMMON /RMATR/ R,RZERO   

      IERR = 0                                               

c      open (6,file='refl.print')
C      open (10,file='spher.write')
c      open (11,file='refl.write')

C***  SEPARATE FIRST-ORDER SCATTERING? (1 FOR YES)

C      NSEP=1
C      IF (NSEP.EQ.1) PRINT 1020
C 1020 FORMAT ('FIRST-ORDER SCATTERING IS SEPARATED')
C      WRITE (11,*) NSEP

C***  SPECIFY NUMERICAL ACCURACY

      EP=1D-5                                                        

C***  COMPUTATION OF THE QUADRATURE DIVISION POINTS ON THE INTERVAL
C***  (0,1) (GAUSS) OR (0,1] (MARK) (Section 3.5)

      NMU=NG
c     CALL MARK (NMU,X,W,1,0)    
c     CALL GAUSS (NMU,1,0,X,W)        

C***  SPECIAL QUADRATURE FORMULA GIVEN BY EQ. (43)

      CALL GAUSS (NMU,0,0,X,W)        
      PI=DACOS(-1D0)
      DO I=1,NG
         X1(I)=PI*X(I)/4D0+PI/4D0
         W1(I)=PI*W(I)/4D0
         W1(I)=W1(I)*DSIN(X1(I))
         X1(I)=DCOS(X1(I))
      ENDDO
      DO I=1,NG
         X(I)=X1(NG-I+1)
         W(I)=W1(NG-I+1)
      ENDDO

C      WRITE (11,1501) NMU
C      WRITE (11,1502) X
C      PRINT 7446,NG                                                    
C 7446 FORMAT('NUMBER OF MU- AND MU_0-VALUES = ',I3)                       
C      PRINT 7447
C 7447 FORMAT ('MU- AND MU_0-VALUES:')
C      WRITE (6,7448) X
C 7448 FORMAT (5D15.8)
      DO I=1,NG
         WW(I)=DSQRT(W(I))                                             
      ENDDO                                                             

C***  READ IN SINGLE-SCATTERING ALBEDO AND LEGENDRE 
C***  EXPANSION COEFFICIENTS  

C      READ (10,*) ALB,MMAX1
      IF (ALB.GT.1D0) PRINT 1010
 1010 FORMAT ('SINGLE-SCATTERING ALBEDO EXCEEDS 1. ',
     $        'EXECUTION TERMINATED.')
      IF (ALB.GT.1D0) IERR = 1010
      IF (ALB.GT.1D0) RETURN
      IF (MMAX1.GT.LMAX1) IERR = 1000
      IF (MMAX1.GT.LMAX1) RETURN
C      AL1(2)=0D0
C      DO L=1,MMAX1
C         READ (10,*) AL1(L)
C      ENDDO   
 1000 FORMAT ('INCREASE LMAX1 IN ALL PARAMETER STATEMENTS!')
      G=AL1(2)/3D0
      PRINT 1200, ALB,G                                                   
 1200 FORMAT('SINGLE-SCATTERING ALBEDO = ',D18.12,
     +       '  G = ', D12.6)      
C      WRITE (11,1501) MMAX1,ALB
C      WRITE (11,1502) (AL1(M1), M1=1,MMAX1)
C 1501 FORMAT (I5,D15.8)
C 1502 FORMAT (6D15.8)

C***  COMPUTATION OF FOURIER COMPONENTS OF THE REFLECTION FUNCTION 
C***  (SECTION 3.3)

      IFLAG=1
      DO M1=1,MMAX1
         print * ,m1
         M=M1-1
         CALL PHAS (M,MMAX1,AL1,X,WW)
         IF (M.EQ.0) CALL F_I (ALB,X,W,WW,EP)
         IF (M.NE.0) CALL RENRM2                   
         CALL RMINF (M,X,WW,EP,ALB,MMAX1,AL1,IFLAG,NSEP)                       
         IF (M.NE.0) GO TO 80
         GLOBAL=0D0
         DO J=1,NG                                                     
            ALBEDO(J)=0D0                                     
            DO I=1,NG                                                     
               ALBEDO(J)=ALBEDO(J)+W(I)*X(I)*RZERO(I,J)*2D0                
            ENDDO                                                           
            GLOBAL=GLOBAL+ALBEDO(J)*W(J)*X(J)*2D0
         ENDDO   
         PRINT 6100, GLOBAL  
C         WRITE (6,*) RZERO(NG,NG)
         DO I=1,NG
            DO J=1,NG     
               RZERO_(I, J) = RZERO(I, J)
            ENDDO
         ENDDO

 6100    FORMAT('SPHERICAL (BOND) ALBEDO = ',D12.6)                  
   80    CONTINUE
C         WRITE (11,*) IFLAG
C         DO I=1,NG
C            WRITE (11,1502) (R(I,J), J=1,I)
C         ENDDO
         DO I=1,NG
            DO J=1,NG     
               R_(M1, I, J) = R(I, J)
            ENDDO
         ENDDO
         MMAX2=M1
         IF (IFLAG.EQ.0) GO TO 450
      ENDDO   
  450 CONTINUE
C      PRINT 6200
C 6200 FORMAT ('PLANE ALBEDO VALUES:')
C      WRITE (6,7448) ALBEDO
C      itime=mclock()
C      time=dfloat(itime)/6000d0
C      print 1001,time
C 1001 format (' time =',f8.2,' min')
C     copy results in the returned variables


      RETURN
      END          

C*******************************************************************            

      SUBROUTINE F_I(ALB,X,W,WW,EPSL)
      PARAMETER (NG=100)
      IMPLICIT REAL*8 (A-H,O-Z)                          
      REAL*8 PH1(NG,NG),PH2(NG,NG),RNRM(NG),
     &       X(NG),DI(NG),DNI(NG),DI1(NG),DNI1(NG),  
     &       XMN(NG),XMNN(NG),W(NG),WW(NG),
     &       PR(NG,NG),PT(NG,NG)                           
      COMMON /PH/ PH1,PH2,RNRM
      COMMON /IPARA/ DI,DNI

C__ RENORMALIZATION OF THE ZEROTH FOURIER COMPONENT OF THE PHASE 
C__ FUNCTION (SECTION 3.6)

      DO I=1,NG
         WI=1D0/WW(I)
         DO J=1,NG
            WJ=WI/WW(J)
            PR(I,J)=PH1(I,J)*WJ
            PT(I,J)=PH2(I,J)*WJ
         ENDDO
      ENDDO
      DO K=1,NG
         S=0D0
         DO I=1,NG
            S=S+W(I)*(PR(I,K)+PT(I,K))
         ENDDO
         EP=(2D0-S)/(PT(K,K)*W(K)) + 1D0
         RNRM(K)=EP
         PH2(K,K)=PH2(K,K)*EP
      ENDDO

C__ CALCULATION OF THE ASYMPTOTIC INTERNAL FIELD (SECTION 3.4)

      IF (ALB.LT.0.7999D0) GO TO 2
      IF (DABS(ALB-1D0).GT.1D-10) GO TO 1
      DO I=1,NG 
         DI(I) = WW(I)
         DNI(I) = WW(I)
      END DO  
      GO TO 2

    1 DO I=1,NG 
         DI(I) = 2D0*WW(I)
         DNI(I) = 0.5D0*WW(I)
      END DO  
      EPS = 0.1D0*EPSL 
      DDK = DSQRT(1D0-ALB)
      DDK1 = DDK 
      LMAX=400
      DO L = 1,LMAX

C_______________ Eq. (38) 
         DDK=0.5D0*(DDK+DDK1)
         DDK1=DDK
         DO I=1,NG 
            DI1(I) = DI(I)
            DNI1(I) = DNI(I)
            XMN(I) = ALB*0.5D0/(1D0-DDK*X(I))
            XMNN(I) = ALB*0.5D0/(1D0+DDK*X(I))
         END DO  
         DO I=1,NG
            SUM = 0D0
            SUMN= 0D0
            DO J=1,NG
               SUM  = SUM + DI1(J)*PH2(I,J) + DNI1(J)*PH1(I,J)    
               SUMN = SUMN + DI1(J)*PH1(I,J) + DNI1(J)*PH2(I,J)    
            END DO
            DI(I) = SUM*XMN(I)
            DNI(I)= SUMN*XMNN(I)
         END DO

C_______________ Eq. (39) 
         SINT= 0D0
         DO I=1,NG
            SINT = SINT + (DI(I)+DNI(I))*WW(I)    
         END DO
         SINT=SINT*ALB*0.5D0
         DO I=1,NG
            DI(I) = DI(I)/SINT
            DNI(I)= DNI(I)/SINT
         END DO

C________________ Eq.(40) 
         SUMK = 0D0
         DO I=1, NG
            SUMK = SUMK+WW(I)*X(I)*(DI(I)-DNI(I))
         END DO
         ZN = ALB*SUMK
         IF (ABS(ZN).GT.1D-8) DDK=2D0*(1D0-ALB)/ZN
C        WRITE (6,*) L,SUM,DDK
         G = 0D0
         DO I=1,NG
            G = MAX(G, DABS(DI(I)-DI1(I)))
            G = MAX(G, DABS(DNI(I)-DNI1(I)))
         END DO
         IF (G.LE.EPS)  GO TO 2
      END DO
      IF (L.EQ.LMAX) WRITE (6,1305)
C     WRITE (6, 1302) ALB, XLSUM
   2  CONTINUE
C     WRITE (6, 1303)(DI(I)*WW(I), I=1,NG)
C     WRITE (6, 1304)(DNI(I)*WW(I), I=1,NG)
 1302 FORMAT('  ALB = ', D12.6, 
     +       '  XLSUM = ', D12.6 )
 1303 FORMAT('   DI = ', 5D12.6)
 1304 FORMAT('   DNI= ', 5D12.6)
 1305 FORMAT('WARNING:  I(MU) HAS NOT CONVERGED.')
 1306 FORMAT('   SUM= ',D12.6)
      SUM=0D0
      DO I=1,NG
         SUM=SUM+(DI(I)+DNI(I))*WW(I)
      ENDDO
      SUM=SUM*ALB*0.5D0
c     WRITE (6,1306) SUM
      RETURN
      END

C**************************************************************

      SUBROUTINE RENRM2                            
C__  REMORMALIZATION OF HIGHER-ORDER FOURIER COMONENTS OF
C__  THE PHASE FUNCTION
      PARAMETER (NG=100)
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 PH1(NG,NG),PH2(NG,NG),RNRM(NG)                         
      COMMON /PH/ PH1,PH2,RNRM                                               
      DO I=1,NG
         PH2(I,I)=PH2(I,I)*RNRM(I)                                          
      ENDDO   
      RETURN
      END

C**************************************************************

      SUBROUTINE FFR (ALB) 
      PARAMETER (NG=100)
      IMPLICIT REAL*8 (A-H,O-Z)           
      REAL*8 PH1(NG,NG),PH2(NG,NG),  
     *       R(NG,NG),F(NG,NG),XX(NG),XXX(NG,NG),B(NG)          
      COMMON /PH/ PH1,PH2
      COMMON /FR/ R,F,XX,XXX
      DO 100 I=1,NG                    
         XI=XX(I)    
         DO 20 J=1,NG                
            A=0D0                 
            DO 15 K=1,NG       
               A=A+R(I,K)*PH1(K,J)                
   15       CONTINUE                           
            B(J)=A                         
   20    CONTINUE                       
         DO 100 J=1,I               
            A1=0D0              
            A2=0D0          
            A3=0D0                  
            DO 30 K=1,NG           
               RKJ=R(K,J)         
               A1=A1+R(I,K)*PH2(K,J)                  
               A2=A2+PH2(I,K)*RKJ                  
               A3=A3+B(K)*RKJ                   
   30       CONTINUE                          
            F(I,J)=ALB*(XXX(I,J)*PH1(I,J)+XX(J)*A1+XI*A2+A3)  
  100 CONTINUE                              
      RETURN                              
      END                              

C**********************************************************************

      SUBROUTINE PHAS (M,MMAX1,AL1,X,WW) 
C__  COMPUTATION OF FOURIER COMPONENTS OF THE PHASE FUNCTION
C__  (SECTION 3.2)
      PARAMETER (NG=100, LMAX1=700)
      IMPLICIT REAL*8 (A-H,O-Z)                                    
      REAL*8 PH1(NG,NG),PH2(NG,NG),
     *       D1(LMAX1),AL1(LMAX1),    
     *       P(LMAX1,NG),PP(LMAX1,NG),    
     *       X(NG),WW(NG)           
      COMMON /PH/ PH1,PH2
      LMIN=M+1                                                 
      DO I=1,NG                                          
         XI=X(I)                                        
         XXI=-XI                                     
         WI=WW(I)                                
         CALL DD(XI,M,MMAX1,D1)
         DO L=LMIN,MMAX1                    
            DL=D1(L)*WI                      
            P(L,I)=DL    
         ENDDO                          
         CALL DD(XXI,M,MMAX1,D1)  
         DO L=LMIN,MMAX1           
            DL=D1(L)*WI            
            PP(L,I)=DL       
         ENDDO                  
      ENDDO                    
      DO J=1,NG
         DO I=1,NG
           PH1(I,J)=0D0
           PH2(I,J)=0D0
           DO L=LMIN,MMAX1
              ALL=AL1(L)
              PH1(I,J)=PH1(I,J)+ALL*PP(L,I)*P(L,J)
              PH2(I,J)=PH2(I,J)+ALL*P(L,I)*P(L,J)
           ENDDO
         ENDDO
      ENDDO
      RETURN   
      END   

C********************************************************************** 

C  CALCULATION OF GENERALIZED SPHERICAL FUNCTIONS THAT ARE NECESSARY
C  IN THE CALCULATION OF THE M-TH FOURIER COMPONENT OF THE PHASE MATRIX 
C  (SECTION 3.2)

      SUBROUTINE DD(X,M,MMAX1,D1)  
      PARAMETER (LMAX1=700)
      IMPLICIT REAL*8 (A-H,O-Z)          
      REAL*8 D1(LMAX1)
      LMAX=MMAX1-1                   
      DX=1D0-X*X                   
      DDX=DSQRT(DX)              
      IF (M.EQ.0) GO TO 10      
      IF (M.EQ.1) GO TO 60     
      IF (M.GE.2) GO TO 110   
   10 L=1
      D1(L)=1D0            
      IF (LMAX.LE.0) RETURN       
      L=L+1
      D1(L)=X                  
      IF (LMAX.LE.1) RETURN   
      L=L+1
      D1(L)=(3D0*X*X-1D0)*0.5D0          
      IF (LMAX.LE.2) RETURN           
      GO TO 160                     
   60 L=2
      D1(L)=DSQRT(0.5D0*DX)           
      IF (LMAX.LE.1) RETURN        
      L=L+1
      D1(L)=DSQRT(1.5D0*DX)*X   
      IF (LMAX.LE.2) RETURN    
      GO TO 160             
  110 A=DX*0.25D0*DSQRT(6D0)                 
      X1=1D0+X                          
      B=X1*X1*0.25D0                  
      X1=1D0-X                    
      C=X1*X1*0.25D0            
      IF (M.EQ.2) GO TO 140    
      DO L=3,M            
         A=A*DSQRT(DFLOAT(2*L-1)/DFLOAT(2*L))*DDX          
      ENDDO                                            

C************************************
      IF (A.LT.1D-302) A=0D0
C************************************

  140 M1=M+1                                         
      D1(M)=0D0                                  
      D1(M1)=A                                
      IF (LMAX.LE.2) RETURN                
  160 M1=MAX(M+1,3)                      
      DM=DFLOAT(M*M)                  
      DM2=DFLOAT(2*M)                
      DO L=M1,LMAX           
         DL=DFLOAT(L)        
         L1=L+1          
         L2=L-1         
         DL2=DFLOAT(L2)                       
         DL1=DFLOAT(2*L-1)              
         DL3=DL2*DL2                
         DDL=DL*DL                 
         A=DSQRT(DL3-DM)       
         B=1D0/DSQRT(DDL-DM)  
         C=DL1*X             
         D1(L1)=(C*D1(L)-A*D1(L2))*B              
      ENDDO                                   
      RETURN                               
      END                              

C**********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
C    N - NUMBER OF POINTS                                             *
C    Z - DIVISION POINTS                                              *
C    W - WEIGHTS                                                      *
C**********************************************************************
 
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END

C******************************************************************** 

      SUBROUTINE RMINF (M,X,WW,EP,ALB,MMAX1,AL1,IFLAG,NSEP)                  
C  ITERATIVE SOLUTION OF THE AMBARTSUMIAN'S EQUATION (SECTION 3.3)
      PARAMETER (NG=100, LMAX1=700)
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 RR(NG,NG),F3(NG,NG), 
     *       X(NG),XX(NG),XXX(NG,NG),WW(NG),A(NG),CC(NG,NG),
     *       PH1(NG,NG),PH2(NG,NG),DI(NG),DNI(NG),
     *       AL1(LMAX1)
      REAL*4 R(NG,NG),RZERO(NG,NG)                                    
      COMMON /RMATR/ R,RZERO                                                    
      COMMON /PH/ PH1,PH2
      COMMON /FR/ RR,F3,XX,XXX
      COMMON /IPARA/ DI,DNI
      M1=M+1                                                                    
      KPAR=0                                                            
      AAM1=0D0                                                          
      AAM2=0D0                                                          
      DO 1 I=M1,MMAX1                                                   
         AAM1=DMAX1(AAM1,AL1(I))                                        
         AAM2=DMAX1(AAM2,AL1(I)/DFLOAT(2*I-1))                          
    1 CONTINUE                                                          
      IF (AAM1*AAM2.GT.EP) KPAR=1                                        
      IF (KPAR.EQ.0) IFLAG=0
      DO I=1,NG                                                      
         XI=1D0/X(I)                                                    
         A(I)=XI                                                        
         XX(I)=XI*0.5D0                                                 
         DO J=1,I                                                    
            XJ=A(J)                                                     
            XIJ=XI*XJ*0.25D0                                            
            XXX(I,J)=XIJ                                                
            XXX(J,I)=XIJ                                                
            XIJ=1D0/(XI+XJ)                                             
            CC(I,J)=XIJ                                                
            CC(J,I)=XIJ                                                 
         ENDDO
      ENDDO

C****  INITIAL APPROXIMATION  *****        
                                              
      DO I=1,NG                         
         DO J=1,I                 
            RR(I,J)=ALB*XXX(I,J)*PH1(I,J)*CC(I,J)  
            RR(J,I)=RR(I,J)
         ENDDO
      ENDDO
      NITER=0             
      IF(KPAR.EQ.0) GO TO 100          

C****   ITERATIONS   *****                                            

   50 CALL FFR (ALB) 
      NITER=NITER+1                    
      IF (NITER.GT.150) GO TO 100  
      DELT=0D0              
      DO I=1,NG       
         WI=1D0/WW(I) 
         DO J=1,I    
            RRIJ=F3(I,J)*CC(I,J)     
            DELT=DMAX1(DELT,DABS(RRIJ-RR(I,J))*WI/WW(J))    
            RR(I,J)=RRIJ                     
            RR(J,I)=RRIJ
         ENDDO
      ENDDO
      IF (DELT.LE.EP) GO TO 100      
      IF (M.EQ.0.AND.ALB.GE.0.8D0) CALL RENRM4 (X,WW,ALB,EP) 
      GO TO 50                    
  100 CONTINUE
      IF (M.NE.0) GO TO 200
      DO J=1,NG             
         WJ=WW(J)             
         DO I=1,J     
            WIJ=RR(I,J)/(WJ*WW(I))          
            RZERO(I,J)=WIJ 
            RZERO(J,I)=WIJ 
         ENDDO
      ENDDO

C****  SUBTRACT THE FIRST-ORDER-SCATTERING CONTRIBUTION
                                              
  200 IF (NSEP.NE.1) GO TO 300
      DO I=1,NG                         
         DO J=1,I                 
            RRIJ=RR(I,J)-ALB*XXX(I,J)*PH1(I,J)*CC(I,J)  
            RR(I,J)=RRIJ  
            RR(J,I)=RRIJ 
         ENDDO
      ENDDO

C****  RETURN TO MAIN PROGRAM

  300 DO J=1,NG             
         WJ=WW(J)             
         DO I=1,J     
            WIJ=RR(I,J)/(WJ*WW(I))          
            R(I,J)=WIJ 
            R(J,I)=WIJ 
         ENDDO
      ENDDO
      WRITE (6,2000) M, NITER                           
      IF (NITER.GE.150) WRITE (6,2001)
      IF (NITER.GE.150) WRITE (6,2002)
 2000 FORMAT('M=',I4,' NITER=',I4)
 2001 FORMAT('WARNING:  NUMERICAL CONERGENCE HAS NOT BEEN REACHED.')
 2002 FORMAT('TRY A LARGER NG.')
      RETURN                   
      END                 

C****************************************************************

      SUBROUTINE RENRM4 (X,WW,ALB,EP)
C  SOBOLEV-VAN DE HULST RENORMALIZATION (SECTION 3.3)
      PARAMETER (NG=100)
      IMPLICIT REAL*8 (A-H,O-Z)         
      REAL*8 RR(NG,NG),RW(NG,NG),T(NG),
     &       W(NG),WW(NG),X(NG),DI(NG),DNI(NG),
     &       DIW(NG),DNIW(NG)
      COMMON /FR/ RR             
      COMMON /IPARA/ DI,DNI
      DELT=0.1*EP         
      FACTOR=0.5D0
      IF (ALB.LT.0.995D0) FACTOR=0.1D0
      IF (ALB.LT.0.95D0) FACTOR=0.05D0
      DO 10 J=1,NG       
         WWJ=WW(J)    
         W(J)=WWJ*WWJ         
         WWJ=1D0/WWJ         
         DIW(J) = DI(J)*WWJ
         DNIW(J) = DNI(J)*WWJ
         DO 10 I=1,J           
            RWIJ=RR(I,J)*WWJ/WW(I) 
            RW(I,J)=RWIJ 
            RW(J,I)=RWIJ 
   10 CONTINUE                 
      MRI=0                 
   20 MRI=MRI+1       
      IF (MRI.GT.200) WRITE (500)
  500 FORMAT ('RENORMALIZATION DIVERGES, REDUCE "FACTOR" IN
     & SUBROUTINE RENRM4')
      IF (MRI.GT.200) GO TO 60 
      U=0D0       
      DO 30 I=1,NG                
         TI=0D0             
         DO 25 J=1,NG    
            TI=TI+W(J)*X(J)*DIW(J)*RW(I,J)
   25    CONTINUE                 
         T(I)=DNIW(I)-2D0*TI          
         U=DMAX1(U,DABS(T(I)))
   30 CONTINUE
 1000 FORMAT('   NIT=',I3)
c     WRITE (6,*) MRI,U
      IF (U.LE.DELT) GO TO 60                
      DO 50 I=1,NG                        
         TI=T(I)                   
         DIWI=DIW(I)
         DO 50 J=1,I           
            RWIJ=RW(I,J)+FACTOR*(DIW(J)*TI+DIWI*T(J))
            RW(I,J)=RWIJ
            RW(J,I)=RWIJ
   50 CONTINUE                                
      GO TO 20                            
   60 DO 70 J=1,NG                      
         WWJ=WW(J)                    
         DO 70 I=1,NG              
            WWIJ=WW(I)*WWJ        
            RR(I,J)=RW(I,J)*WWIJ   
   70 CONTINUE               
      RETURN             
      END             

C**********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF THE MARKOV QUADRATURE       *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1], IF IND1 = 1 - ON      *
C    INTERVAL  (0,1].                                                 *
C    N - NUMBER OF POINTS                                             *
C    X - DIVISION POINTS                                              *
C    W - WEIGHTS                                                      *
C**********************************************************************

      SUBROUTINE MARK (N,X,W,IND1,IND2)                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      EXTERNAL F                                                        
      REAL*8 X(N),W(N)                                                  
      N=N-1                                                            
      BBB=0.01 D0                                                       
      IF (N.GE.30) BBB=0.001 D0                                         
      IF (N.GE.100) BBB=0.00015 D0                                      
      IF (N.GE.200) BBB=0.00003 D0                                      
      TOL=1D-15                                                         
      E=2D0                                                             
      DN=DFLOAT(N)                                                      
      DD=2D0*DN+1D0                                                     
      A=-1D0                                                            
      K=0                                                               
      DD=DD*DD                                                          
      YB=F(A,N)                                                         
   10 Y=YB                                                              
      B=A                                                               
   15 B=B+BBB                                                           
      IF (1D0-B.LE.1D-15) GO TO 100                                     
      YB=F(B,N)                                                         
      IF (Y/YB.GE.0D0) GO TO 15                                         
      K=K+1                                                             
      XX=ZEROIN(A,B,F,TOL,N)                                            
      IF (K.GE.2) BBB=(XX-X(K-1))*0.2D0                                 
      X(K)=XX                                                           
      D=FF(XX,N)                                                        
      C=D*2D0*(DN+1D0)*DN                                               
      WW=4D0*(1D0+XX)*DD/(C*C)                                          
      W(K)=WW                                                           
      E=E-WW                                                            
      A=B                                                               
      GO TO 10                                                          
  100 K=K+1                                                             
      W(K)=E                                                            
      X(K)=1D0                                                          
      IF (IND2.NE.1) GO TO 130                                          
      M=2*(K-1)                                                         
      A=0D0                                                             
      DM=DFLOAT(M)                                                      
      DO 120 I=1,K                                                      
  120     A=A+W(I)*DABS(X(I))**DM                                      
      B=0.5D0*DFLOAT(M+1)                                               
      A=A*B                                                             
  130 IF (IND1.NE.1) GO TO 140                                            
      DO 135 I=1,K                                                              
           W(I)=0.5D0*W(I)                                                      
  135      X(I)=0.5D0*(X(I)+1D0)                                                
  140 N=K                                                                       
      RETURN                                                                    
      END                                                                       

C***********************************************************************

      DOUBLE PRECISION FUNCTION F(X,N)                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      A=1D0                                                             
      B=X                                                               
      DO 10 I=1,N                                                       
           C=B                                                          
           DI=DFLOAT(I)                                                 
           B=((2D0*DI+1D0)*X*B-DI*A)/(DI+1D0)                           
   10      A=C                                                          
      DN=DFLOAT(N)                                                      
      F=(A-B)/(1D0-X)                                                   
      RETURN                                                            
      END                                                               

C***********************************************************************

      DOUBLE PRECISION FUNCTION FF(X,N)                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      M=N-1                                                             
      A=1D0                                                             
      B=X                                                               
      DO 10 I=1,M                                                       
           C=B                                                          
           DI=DFLOAT(I)                                                 
           B=((2D0*DI+1D0)*X*B-DI*A)/(DI+1D0)                           
   10      A=C                                                          
      DN=DFLOAT(N)                                                      
      FF=(A-B)/(1D0-X)                                                  
      RETURN                                                            
      END                                                               

C***********************************************************************

      DOUBLE PRECISION FUNCTION ZEROIN (AX,BX,F,TOL,N)                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      EPS=1D0                                                           
   10 EPS=0.5D0*EPS                                                     
      TOL1=1D0+EPS                                                      
      IF (TOL1.GT.1D0) GO TO 10                                         
   15 A=AX                                                              
      B=BX                                                              
      FA=F(A,N)                                                         
      FB=F(B,N)                                                         
   20 C=A                                                               
      FC=FA                                                             
      D=B-A                                                             
      E=D                                                               
   30 IF (DABS(FC).GE.DABS(FB)) GO TO 40                                
   35 A=B                                                               
      B=C                                                               
      C=A                                                               
      FA=FB                                                             
      FB=FC                                                             
      FC=FA                                                             
   40 TOL1=2D0*EPS*DABS(B)+0.5D0*TOL                                    
      XM=0.5D0*(C-B)                                                    
      IF (DABS(XM).LE.TOL1) GO TO 90                                    
   44 IF (FB.EQ.0D0) GO TO 90                                           
   45 IF (DABS(E).LT.TOL1) GO TO 70                                     
   46 IF (DABS(FA).LE.DABS(FB)) GO TO 70                                
   47 IF (A.NE.C) GO TO 50                                              
   48 S=FB/FA                                                           
      P=2D0*XM*S                                                        
      Q=1D0-S                                                           
      GO TO 60                                                          
   50 Q=FA/FC                                                           
      R=FB/FC                                                           
      S=FB/FA                                                           
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))                                
      Q=(Q-1D0)*(R-1D0)*(S-1D0)                                         
   60 IF (P.GT.0D0) Q=-Q                                                
      P=DABS(P)                                                         
      IF ((2D0*P).GE.(3D0*XM*Q-DABS(TOL1*Q))) GO TO 70                  
   64 IF (P.GE.DABS(0.5D0*E*Q)) GO TO 70                                
   65 E=D                                                               
      D=P/Q                                                             
      GO TO 80                                                          
   70 D=XM                                                              
      E=D                                                               
   80 A=B                                                               
      FA=FB                                                             
      IF (DABS(D).GT.TOL1) B=B+D                                        
      IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)                           
      FB=F(B,N)                                                         
      IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20                           
   85 GO TO 30                                                          
   90 ZEROIN=B                                                          
      RETURN                                                            
      END         
