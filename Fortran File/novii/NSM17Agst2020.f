      PROGRAM SNM        
      
c     ****************************************************
c     verified 17 August 2020 check "new parameter sets"
c     +Delta
c     ****************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      INTEGER I
      PI  = 3.14159265358979D0
      HC  = 197.327D0

 

C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c     FSU 
       KFO =1.306D0*HC

C------------------------------------------------------------------------------
      RHO = 2D0*KFO*KFO*KFO/(3D0*PI*PI)
 
     
      OPEN(unit=1,status='unknown',file='EOS_NSM_FSU2HZ2.dat')
      OPEN(unit=2,status='unknown',file='SOS_NSM_FSU2HZ2.dat')





c     DO 10 I=1,85
      DO 10 I=1,850
        RHBO=I * 0.01D0  
             
         RB=RHBO*RHO    

           
   

C     CALL ENERGY DENSITY AND BINDING ENERGY
C     

         CALL FED(RHBO,ED,EB,EL,PT)
        
c     CALL other properties
c         CALL FERMIM(RB,KFP,KFN,KFE,KFM,YE,YN,YP,YM)
c         CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP)
         CALL SPOSN(RHBO,DPDE,CL,SOS,ED,PRESSCR)
c         DEDEN=DE(RHBO)
c         DPRESS=DP(RHBO)
C     
C     CALL PRESSURE
C     
         PRESS= PT
c         ABC=func1(RHBO)

C     RESULTS
C     
C     ****************************************************
C     RB is Baryon density,  EB is binding energy, ED is energy density
C     PRESS is pressure
C     KNCM is incompresibility
c     RB in fm^-3, KF in fm^-1, EB in MeV, ED  and PRESS in MEV/fm^3

C     ****************************************************
       WRITE(1,*)RHBO,(RB/(HC*HC*HC)),PRESS,ED
       WRITE(2,*)RHBO,SOS,CL,PRESSCR,ED
c       WRITE(*,*)RHBO,CL,SOS,ED,PRESSCR,PRESS
       WRITE(*,*)RHBO
      
 10   CONTINUE

         STOP
         END
c-----------------------------------------------------------------
      INCLUDE "parset.f"
c------------------------------------------------------------------      
      FUNCTION gunc1(xa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      CALL FED(xa,ED,EB,EL,PT)
      gunc1=ED
      RETURN
      END
      FUNCTION gunc2(xa)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     
       CALL FED(xa,ED,EB,EL,PT)
      gunc2=PT
      RETURN
      END

      FUNCTION DE(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL gunc1
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      DE = (gunc1(x4)-8.D0*gunc1(x2)+8.D0*gunc1(x1)-gunc1(x3))/(12.D0*h)
      RETURN
      END

      FUNCTION DP(xa)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      EXTERNAL gunc2
      h  = 1.D-4
      x1 = xa+h
      x2 = xa-h 
      x3 = xa+2.D0*h
      x4 = xa-2.D0*h    
      DP = (gunc2(x4)-8.D0*gunc2(x2)+8.D0*gunc2(x1)-gunc2(x3))/(12.D0*h)
      RETURN
      END
      
      SUBROUTINE SPOSN(RHBO,DPDE,CL,SOS,ED,PRESSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL DE,DP,gunc2
      CL=1.D0/DSQRT(3.D0)
      DPRESS=DP(RHBO)
      DEDEN=DE(RHBO)
      DPDE=DPRESS/DEDEN
      
      IF (DPDE .GE. 0.0D0) THEN 
         SOS=DSQRT(DPDE)
      ELSE
         SOS=0.0D0
      ENDIF
c  Note: be carefull here token SOS=0 means unstable (imaginer speed of sound)!!!         
      IF (DSQRT((SOS-CL)*(SOS-CL)) .LT. 2.0D-3) THEN
         CALL FED(RHBO,ED,EB,EL,PT)
         EDX=ED
         PRX=gunc2(RHBO)
      END IF
      
      IF (SOS .GT. CL ) THEN
         PRESSCR=CL*(ED-EDX)+PRX
      ELSE   
         PRESSCR=gunc2(RHBO)
      END IF


      RETURN
      END
c------------------------------------------------------------------------------


C------------------------------------------------------------------------------
C      

      
C     ENERGY CALCULATIONS      
      SUBROUTINE FED(RHBO,ED,EB,EL,PT)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      HC  = 197.327D0


C------------------------------------------------------------------------------
C nuclear matter Saturation density !!!

c     FSU 
       KFO =1.306D0*HC

      
C-----------------------------------------------------------------------------
    
      RHO = 2.D0*KFO*KFO*KFO/(3.D0*PI*PI)   
      RB  = RHBO*RHO  
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)



C     ****************************************************
C      CHEK NUCLEAR MATTER
C     **************************************************** 
c      RP = 0.0*RB
c      RN = 1.0*RB
c      RE =0.0*RB
c      RM=0.0*RB

c      KFE=(3.D0*PI*PI*RE)**(1.D0/3.D0)
c      KFM=(3.D0*PI*PI*RM)**(1.D0/3.D0)
c      KFP=(3.D0*PI*PI*RP)**(1.D0/3.D0)
c      KFN=(3.D0*PI*PI*RN)**(1.D0/3.D0)
c       EE=0.0D0
c       EM=0.0D0
c     **************************************************** 

      CALL FERMIM(RB,KFP,KFN,KFE,KFM,YE,YN,YP,YM)

      CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP)

      RP = KFP*KFP*KFP/(3.D0*PI*PI)
      RN = KFN*KFN*KFN/(3.D0*PI*PI)   

      U = 0.5D0*MS*MS*SIG*SIG
      U = U + 0.5D0*MD*MD*DEL*DEL
      U = U - 0.5D0*MV*MV*V0*V0
      U = U - 0.5D0*MR*MR*B0*B0
      U = U + 0.333D0*B2*SIG*SIG*SIG
      U = U + 0.25D0*B3*SIG*SIG*SIG*SIG
      U = U - 0.25D0*C1*V0*V0*V0*V0
      U = U - D2*SIG*V0*V0  
      U = U - F2*SIG*B0*B0
      U = U - 0.5D0*D3*SIG*SIG*V0*V0
c---------------------------------------------------------------------------
c     should be modified      
      U = U - 0.5D0*G3*SIG*SIG*B0*B0
      U = U - 0.5D0*G4*V0*V0*B0*B0
      U = U - 0.25D0*G5*B0*B0*B0*B0
c-------------------------------------------------------------------------
      EP = KFP*DSQRT((KFP*KFP)+(MP*MP))
      EP = EP*(2.D0*KFP*KFP+MP*MP)
      TEMP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      EP = EP - MP*MP*MP*MP*TEMP
      EP = EP/(8.D0*PI*PI)

      EN = KFN * DSQRT((KFN*KFN)+(MN*MN))
      EN = EN * (2.D0*KFN*KFN+MN*MN)
      TEMN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      EN = EN - MN*MN*MN*MN*TEMN
      EN = EN/(8.D0*PI*PI)

      EE = KFE * DSQRT((KFE*KFE)+(ME*ME))
      EE = EE * (2.D0*KFE*KFE+ME*ME)
      TEME = DLOG((KFE+DSQRT((KFE*KFE)+(ME*ME)))/ME)
      EE = EE - ME*ME*ME*ME*TEME
      EE = EE /(8.D0*PI*PI)

      EM = KFM * DSQRT((KFM*KFM)+(MO*MO))
      EM = EM * (2.D0*KFM*KFM+MO*MO)
      TEMM = DLOG((KFM+DSQRT((KFM*KFM)+(MO*MO)))/MO)
      EM = EM - MO*MO*MO*MO*TEMM
      EM = EM /(8.D0*PI*PI)



      E = EP+EN+EE+EM+GV*V0*(RP+RN)+0.5D0*GR*B0*(RP-RN)

      ED = (E + U )/(HC*HC*HC)
      
      EB = (E + U )/RB

     
      EB = EB - MB

      EL= (EE+EM)/(RB*HC*HC*HC)
c------------------------------------------------------------------------
C Pressure "Analitic"     
       
      PP = KFP*DSQRT((KFP*KFP)+(MP*MP))
      PP = PP*(2.D0*KFP*KFP-3.0D0*MP*MP)
      TEMP2 = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      PP = PP + 3.0D0*MP*MP*MP*MP*TEMP2
      PP = PP/(8.D0*PI*PI)   

      PN = KFN * DSQRT((KFN*KFN)+(MN*MN))
      PN = PN * (2.D0*KFN*KFN-3.0D0*MN*MN)
      TEMN2 = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      PN = PN + 3.0D0*MN*MN*MN*MN*TEMN2
      PN = PN/(8.D0*PI*PI)

      PE = KFE * DSQRT((KFE*KFE)+(ME*ME))
      PE = PE * (2.D0*KFE*KFE-3.0D0*ME*ME)
      TEME2 = DLOG((KFE+DSQRT((KFE*KFE)+(ME*ME)))/ME)
      PE = PE + 3.0D0*ME*ME*ME*ME*TEME2
      PE = PE /(8.D0*PI*PI)

      PM = KFM * DSQRT((KFM*KFM)+(MO*MO))
      PM = PM * (2.D0*KFM*KFM-3.0D0*MO*MO)
      TEMM2 = DLOG((KFM+DSQRT((KFM*KFM)+(MO*MO)))/MO)
      PM = PM + 3.0D0*MO*MO*MO*MO*TEMM2
      PM = PM /(8.D0*PI*PI)

      P=(PP+PN+PE+PM)/3.0D0

      PT = (P - U )/(HC*HC*HC)

      RETURN
      END

 
c------------------------------------------------------------------------------ 

      SUBROUTINE  FERMIM(RB,KFPF,KFNF,KFEF,KFMF,YE,YN,YP,YM)
C
C   Fermi momenta calculations
C

      IMPLICIT DOUBLE PRECISION (A-H,K-Z)       

      PI  = 3.14159265358979D0
      MUE = 100.D0      
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)      
       
       YN0=0.8D0
       YP0=0.2D0
 
 50    RBN = YN0*RB
         RBP = YP0*RB
         KFP =(3.D0*PI*PI*RBP)**(1D0/3D0)
         KFN =(3.D0*PI*PI*RBN)**(1D0/3D0)
 
         CALL FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP)        
        
         PREYP = 0.0D0
 20      YP = PREYP-(1.D0/MUE)* FYP(PREYP,MN,MP,RB,ME,MO,SIG,DEL,V0,B0)
     
         IF(((YP-PREYP)*(YP-PREYP)).LT.1.D-24) GOTO 9
         PREYP = YP
         GOTO 20
 9       CONTINUE
             
C        
C     Fraction of each contituent
C
         A = (3.D0*PI*PI*(1-YP))**(2D0/3D0)
         B = (3.D0*PI*PI*YP)**(2D0/3D0)
         CN = (MN*MN)/RB**(2D0/3D0)
         CP = (MP*MP)/RB**(2D0/3D0)
         D = (ME*ME)/RB**(2D0/3D0)
         F = (MO*MO)/RB**(2D0/3D0)        
         G = 1.D0/(3.D0*PI*PI)
         MZ= (MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG+G5*B0*B0)
         R = 0.5D0*GR*GR*RB**(2.D0/3.D0)*(1.D0-2.D0*YP)/MZ
         FZ=(DSQRT(A+CN)-DSQRT(B+CP)+R)*(DSQRT(A+CN)-DSQRT(B+CP)+R)
         YN= (1.D0-YP)         
         YE= G*(FZ-D)**(3D0/2D0)
         IF ((FZ-F) .GT. 0D0) THEN
            YM = G*(FZ-F)**(3D0/2D0)
         ELSE 
            YM = 0D0
         ENDIF  
    

       IF(((YP-YP0)*(YP-YP0)).LT.1.D-6 ) GOTO 40
        YP0=YP
        YN0=YN
       GOTO 50
 40    CONTINUE
 
       KFEF= (3.D0*YE*RB*PI*PI)**(1D0/3D0)
       KFPF= (3.D0*YP*RB*PI*PI)**(1D0/3D0)
       KFNF= (3.D0*YN*RB*PI*PI)**(1D0/3D0)
       KFMF= (3.D0*YM*RB*PI*PI)**(1D0/3D0)
c      WRITE(*,*)YP,YN,YE,YM,B0
      RETURN
      END
c-------------------------------------------------------------------

      SUBROUTINE FRG2(KFP,KFN,SIG,DEL,V0,B0,MN,MP)
C
C   meson fields calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      MUE = 1.D+9      
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)  

       

       RP = KFP*KFP*KFP/(3.D0*PI*PI)
       RN = KFN*KFN*KFN/(3.D0*PI*PI)    

C     initial guess
       
c      PREDEL = (RP-RN)*GD/(MD*MD)
c      PRESIG = (RN+RP)*GS/(MS*MS)
c      PREV0  = (RN+RP)*GV/(MV*MV)
      
c     MR2EF=MR*MR+2.D0*F2*PRESIG+G4*PREV0*PREV0+G3*PRESIG*PRESIG
c      MR2EF=MR*MR
c      PREB0  = 0.5D0*(RP-RN)*GR/MR2EF

      PREDEL = 0.0D0
      PRESIG = 0.0D0
      PREV0  = 0.0D0
      PREB0  = 0.0D0
      
 20   SIG = PRESIG-(1.D0/MUE)*FS(PRESIG,PREV0,PREDEL,PREB0,MB,KFP,
     &     KFN,GS,GR,GD,MS,MR,B2,B3,D2,D3,F2,G3,G4)
      V0 = PREV0-(1.D0/MUE)*FS1(PREV0,PRESIG,PREB0,MV,MR,GR,GV,KFP,
     &     KFN,C1,D2,D3,F2,G3,G4)     
      B0 = PREB0-(1.D0/MUE)*FS3(PREB0,PREV0,PRESIG,MV,MR,GR,GV,KFP,KFN,
     &     F2,G3,G4,G5)
      DEL = PREDEL-(1.D0/MUE)*FS4(PREDEL,PRESIG,MD,GS,GD,KFP,KFN,MB)

      IF(((SIG-PRESIG)*(SIG-PRESIG)).LT.1.D-28 .AND.
     &   ((DEL-PREDEL)*(DEL-PREDEL)).LT.1.D-28 .AND.
     &   ((B0-PREB0)*(B0-PREB0)).LT.1.D-28     .AND. 
     &   ((V0-PREV0)*(V0-PREV0)).LT.1.D-28           ) GOTO 8

      PRESIG = SIG
      PREV0 = V0
      PREB0= B0
      PREDEL = DEL      
      GOTO 20
 8    CONTINUE 
               
      MN = MB+GS*SIG-GD*DEL
      MP = MB+GS*SIG+GD*DEL 

      RETURN
      END

c------------------------------------------------------------------------------

      FUNCTION FS(SIG,V0,DEL,B0,MB,KFP,KFN,GS,GR,GD,MS,MR,B2,B3,D2,
     &            D3,F2,G3,G4)
C
C  sigma meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI = 3.14159265358979D0
      RP = 1.D0*KFP*KFP*KFP/(3.D0*PI*PI)
      RN = 1.D0*KFN*KFN*KFN/(3.D0*PI*PI)
    
      MN= MB+GS*SIG-GD*DEL
      MP = MB +GS*SIG+GD*DEL
      FS = MS*MS*SIG
      FS = FS + B2*SIG*SIG
      FS = FS + B3*SIG*SIG*SIG
      FS = FS - D2*V0*V0
      FS = FS - D3*SIG*V0*V0
c-------------------------------------------------------------------------
c  should be modifed      
      FS = FS - F2*B0*B0
      FS = FS - G3*SIG*B0*B0
c------------------------------------------------------------------------
      VN = KFN*DSQRT(((KFN*KFN)+(MN*MN)))
      VIN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      VN = VN - MN*MN*VIN
      VN = VN*MN
      VN = VN/(2.D0*PI*PI)

      VP = KFP*DSQRT(((KFP*KFP)+(MP*MP)))
      VIP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      VP = VP - MP*MP*VIP
      VP = VP*MP
      VP = VP/(2.D0*PI*PI)

      V = (VP+VN)*GS
      FS= FS + V

      RETURN
      END
      
      
c------------------------------------------------------------------------------      
      
      FUNCTION FS1(V0,SIG,B0,MV,MR,GR,GV,KFP,KFN,C1,D2,D3,F2,G3,G4)
C
C omega meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      RP = 1.D0*KFP*KFP*KFP/(3.D0*PI*PI)
      RN = 1.D0*KFN*KFN*KFN/(3.D0*PI*PI)

      FS1 = MV*MV*V0
      FS1 = FS1 - GV*(RP+RN)
      FS1 = FS1 + 2.D0*D2*SIG*V0
      FS1 = FS1 + D3*SIG*SIG*V0
c--------------------------------------------------------------------
c should be modifed      
      FS1 = FS1 + C1*V0*V0*V0
      FS1 = FS1 + G4*V0*B0*B0
c-------------------------------------------------------------------      
      RETURN
      END


       FUNCTION FS3(B0,V0,SIG,MV,MR,GR,GV,KFP,KFN,F2,G3,
     %             G4,G5)
C
C rho meson calculation
C
c--------------------------------------------------------------------
c     should be modifed
       
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      RP = 1.D0*KFP*KFP*KFP/(3.D0*PI*PI)
      RN = 1.D0*KFN*KFN*KFN/(3.D0*PI*PI)
 
      FS3 = MR*MR*B0
      FS3 = FS3 - 0.5D0*GR*(RP-RN)
      FS3 = FS3 + 2.D0*F2*SIG*B0
      FS3 = FS3 + G3*SIG*SIG*B0
      FS3 = FS3 + G5*B0*B0*B0
      FS3 = FS3 + G4*V0*V0*B0
      RETURN
      END
c-----------------------------------------------------------------------------

      FUNCTION FS4(DEL,SIG,MD,GS,GD,KFP,KFN,MB)
C
C  Delta meson calculation
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      MN = MB+GS*SIG-GD*DEL
      MP = MB+GS*SIG+GD*DEL
 
      VN = KFN*DSQRT(((KFN*KFN)+(MN*MN)))
      VIN = DLOG((KFN+DSQRT((KFN*KFN)+(MN*MN)))/MN)
      VN = VN - MN*MN*VIN
      VN = VN*MN
      VN = VN /(2.D0*PI*PI)

      VP = KFP*DSQRT(((KFP*KFP)+(MP*MP)))
      VIP = DLOG((KFP+DSQRT((KFP*KFP)+(MP*MP)))/MP)
      VP = VP - MP*MP*VIP
      VP = VP*MP
      VP = VP/(2.D0*PI*PI)
 
      V = (VP-VN)*GD
      FS4 = MD*MD*DEL+V 
      RETURN
      END
c-----------------------------------------------------------------------------
 
      FUNCTION FYP(YP,MN,MP,RB,ME,MO,SIG,DEL,V0,B0)
C
C   Function of fraction
C
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PI  = 3.14159265358979D0
      CALL SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
 
      A = (3.D0*PI*PI*(1.D0-YP))**(2D0/3D0)
      B = (3.D0*PI*PI*YP)**(2D0/3D0)
      CN = (MN*MN)/RB**(2D0/3D0)
      CP = (MP*MP)/RB**(2D0/3D0)
      D = (ME*ME)/RB**(2D0/3D0)
      F = (MO*MO)/RB**(2D0/3D0)    
      G = 1.D0/(3.D0*PI*PI)
      MZ= (MR*MR+2.D0*F2*SIG+G4*V0*V0+G3*SIG*SIG+G5*B0*B0)
      R = 0.5D0*GR*GR*RB**(2.D0/3.D0)*(1.D0-2.D0*YP)/MZ   
      H = (DSQRT(A+CN)-DSQRT(B+CP)+R)*(DSQRT(A+CN)-DSQRT(B+CP)+R)
      IF ((H-F).GE.0.D0) THEN 
         FYP = -G**(2D0/3D0)*((H-D)**(3D0/2D0)
     &        +(H-F)**(3D0/2D0))**(2D0/3D0)+YP**(2D0/3D0)
      ELSE     
         FYP =-G**(2D0/3D0)*(H-D)+YP**(2D0/3D0)
      ENDIF
      RETURN
      END

c-----------------------------------------------------------------------------
