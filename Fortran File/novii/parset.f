c     List parameter used 1 July 2020
c-----------------------------------------------------------------------   
C
C Parameter BSP B
C  KFO =1.303
C
      SUBROUTINE QSET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0
      MB  = 939.20D0
    
      K2  = -1.0681D0
      K3  = 14.9857D0
      ET2 = -0.0872D0
      ET3 = 3.1265D0

      EP3 = 0.0D0
      EPR3= 0.0D0

      MS = 0.5384D0*MB
      MV  =0.8333D0*MB
      MR  =0.8200D0*MB
      MD  = 980D0

      ETR = 0.0D0   
      ETR3 = 0.0D0
      ETR4 =  53.7642D0
 

      GS = 0.8764D0*4.D0*PI
      GV = 1.1481D0*4.D0*PI
      GR = 1.0508D0*4.D0*PI
 
      GD = 0.0D0

      B2 = GS*MS*MS*K2/(2.D0*MB)
      B3 = GS*GS*MS*MS*K3/(6.D0*MB*MB)
      C1 = GV*GV*EP3/(6.D0)
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*MR*MR*ETR4/(2.D0*MB*MB)
      G5 = GR*GR*EPR3/(6.D0)     


c--------------------------------------------------------------
      RETURN
      END 
c---------------------------------------------------------------
C
C Parameter G3 A
C  KFO =1.30
C
      SUBROUTINE SET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0
      MB  = 939.0D0
    
      K2  = -2.606D0
      K3  =  1.694D0
      ET2 = -0.424D0
      ET3 =  0.114D0

      EP3 = 1.010D0
    
      MS = 0.559D0*MB
      MV  =0.832D0*MB
      MR  =0.820D0*MB
      MD  =1.043D0*MB 

      ETR = -0.645D0   
      ETR3 = 0.0D0
      ETR4 = 33.250D0
      EPR3=0.0D0

      GS = 0.782D0*4.D0*PI
      GV = 0.923D0*4.D0*PI
      GR = 0.962D0*4.D0*PI
 
      GD =0.160D0*4.D0*PI 
     
      B2 = GS*MS*MS*K2/(2.D0*MB)
      B3 = GS*GS*MS*MS*K3/(6.D0*MB*MB)
      C1 = GV*GV*EP3/(6.D0)
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*MR*MR*ETR4/(2.D0*MB*MB)
      G5 = GR*GR*EPR3/(6.D0)     


c--------------------------------------------------------------
      RETURN
      END 
c---------------------------------------------------------------

C   TM1e parameter set C
C  KFO =1.290
C
      SUBROUTINE CSET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0
      MB  = 938.0D0
  
   
      GS =  0.798070534D0*4.D0*PI
      GV =  1.00378128D0*4.D0*PI
      GR =1.11180874D0*4.D0*PI
      GD = 0.0D0


      K2  = -1.02149613D0
      K3  = 0.124465984D0
      ET2 = 0.0D0
      ET3 = 0.0D0

      EP3 = 2.68897577D0
      EPR3= 0.0D0

      MS  =  511.198D0
      MV  =  783.D0
      MR  =  770.D0
      MD  =  980.D0

      ETR = 0.0D0   
      ETR3 =0.0D0
C     perhatikan ada faktor 1/2 dari parameter karena penulisan lagrangian !!!
      
      ETR4 = 2.D0*0.0429D0
c-------------------------------------------------------------  
  
      B2 = GS*MS*MS*K2/(2.D0*MB)
      B3 = GS*GS*MS*MS*K3/(6.D0*MB*MB)
      C1 = GV*GV*EP3/(6.D0)
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*GR*GR*ETR4
      G5 = GR*GR*EPR3/6.D0     


c--------------------------------------------------------------
      RETURN
      END 
c---------------------------------------------------------------
C
C FSU Garnet  (kf=1.31 fm^-1)D
C

      SUBROUTINE XSET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD) 
  

   
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0
      MB  = 939.00D0
 
      MV  = 782.500D0
      MR  = 763.00D0
      MD  = 980D0
      MS  = 496.939D0


      GS = 10.5047228D0
      GV = 13.7001716D0

      GR = 13.8898305D0

      GD = 0.0D0

      K2  = -3.26018D0
      K3  = -0.003551D0
      ET2 = -2.D0*MB*GV*GV*HC/(MV*MV)*0.0D0
      ET3 = 2.D0*MB*MB*GV*GV/(MV*MV)*0.0D0
      EP3 = 0.02350D0
      EPR3= 0.0D0
      ETR = -2.D0*MB*GR*GR*HC/(MR*MR)*0.0D0
      ETR3 = 2.D0*MB*MB*GR*GR/(MR*MR)*0.0D0

C perhatikan ada faktor 1/2 dari parameter karena penulisan lagrangian !!!

      ETR4 = 2.D0*0.043377D0
            
c      write(*,*)GS,GV,GR,K2,K3,EP3,ETR4     

      

      B2 = GS*GS*GS*K2/2.D0
      B3 = GS*GS*GS*GS*K3/6.D0
      C1 = GV*GV*GV*GV*EP3/6.D0
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*GR*GR*ETR4
      G5= GR*GR*EPR3/(6.D0)     

      
      RETURN
      END 

c--------------------------------------------------------------------------
C
C FSU2H  (kf=1.306 fm^-1)E
C

      SUBROUTINE ESET(MB,MS,ME,MO,HC,ET2,ET3,EP3,ETR,MV,MR,GS,GV,GR,B2
     &     ,B3,C1,D2,D3,F2,G3,G4,G5,GD,MD) 
  

   
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
C     Initialisation
      PI  = 3.14159265358979D0
      HC  = 197.327D0
      ME  = 0.511D0
      MO  = 105.6D0
      MB  = 939.00D0
 
      MV  = 782.500D0
      MR  = 763.00D0
      MD  = 980D0
      MS  = 497.479D0


      GS = 9.9712888D0*DSQRT(102.7200D0/99.4266D0)
      GV = 13.032072D0*DSQRT(169.5315D0/169.8349D0)

      GR = 13.589985D0*DSQRT(247.3409D0/184.6877D0)

      GD = 0.0D0

      K2  = -4.0014D0
      K3  = -0.013298D0
      ET2 = -2.D0*MB*GV*GV*HC/(MV*MV)*0.0D0
      ET3 = 2.D0*MB*MB*GV*GV/(MV*MV)*0.0D0
      EP3 = 0.008D0
      EPR3= 0.0D0
      ETR = -2.D0*MB*GR*GR*HC/(MR*MR)*0.0D0
      ETR3 = 2.D0*MB*MB*GR*GR/(MR*MR)*0.0D0

C perhatikan ada faktor 1/2 dari parameter karena penulisan lagrangian !!!

      ETR4 = 2.D0*0.05D0
            
c      write(*,*)(GS/(4*PI)),(GV/(4*PI)),(GR/(4*PI)),K2,K3,EP3,ETR4,MS/MB     

        


      B2 = GS*GS*GS*K2/2.D0
      B3 = GS*GS*GS*GS*K3/6.D0
      C1 = GV*GV*GV*GV*EP3/6.D0
      D2 = GS*MV*MV*ET2/(2.D0*MB)
      D3 = GS*GS*MV*MV*ET3/(2.D0*MB*MB)  
      F2 = ETR*GS*MR*MR/(2.D0*MB)
      G3 = GS*GS*MR*MR*ETR3/(2.D0*MB*MB)
      G4 = GV*GV*GR*GR*ETR4
      G5= GR*GR*EPR3/(6.D0)     
      RETURN
      END 



