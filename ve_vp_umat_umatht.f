
!   ##################################################################################################################################
!   USER SUBROUTINE TO COMPUTE DRUCKER-PRAGER VISCOPLASTICITY WITH
!   CONSISTENT JACOBIAN OBTAINED BY LINEARIZING THE RETURN MAPPING
!   EQUATION (RESIDUAL VISCOPLASTIC POTENTIAL FUCNTION). 
!
!   WE USE DISLOCATION CREEP TO MODEL HIGH-TEMPERATURE DUCTILE DEFORMATION.
!   
!   WRITTEN BY EKEABINO MOMOH, HARSHA BHAT AND STEVE TAIT
!   ################################################################################################################################## 
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!     DOUBLE PRECISION FACTOR, STRESS
      CHARACTER*80 CMNAME
! 
      COMMON /CONDUCTIVITY_PARAMS/ TREF, ALPH_LE, RHOD_REF, TEMPERA, PRESS

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DIMENSION DEVSTRAN(NTENS),S(NTENS),STRIAL(NTENS),SOID(NTENS),
     1          DVPRJT(NTENS,NTENS),VSTRESS_DEV(NTENS), EPSILON_VP(NTENS), 
     2          DEPSILON_VP(NTENS), SPRJD(NTENS,NTENS),ELAS(NTENS), 
     3          HC(NTENS),DEPSILON_V(NTENS), EPSILON_V(NTENS), DDSDDE_V(NTENS,NTENS)
      PARAMETER (R0=0.0D0, P01=0.1D0, R1=1.0D0, R2=2.0D0, R3=3.0D0,  
     1           R4=4.0D0, R6=6.0D0, R7=7.0D0,R9=9.0D0, R12=12.0D0, 
     2           R40=40.0D0,P13=1.0D0/3.0D0, TOL=1E-6, INIT_GUESS=1E-21, 
     3           P05=0.5D0,P098=0.98D0,P23=2.0D0/3.0D0,P06=6.0D0/10.0D0,
     4           P07=7.0D0/10.0D0, A11=1085.7D0,A22=132.9D0,A33=-5.1D0,
     5           B11=1475.0D0,B22=80.D0,B33=-3.2D0, PWRFPT=1.5D0, TQC=1.0D0)
!   ##################################################################################################################################
!     INITIALIZE STATE VARIABLES VARIABLES
!   ################################################################################################################################## 
      EPSILON_V=R0
      EPSILON_VP=R0
      DEPSILON_V=R0
      DEPSILON_VP=R0
      PHI_Y=R0
      DGAMA_VISC=R0
      DGAMA_VISC_OLD=R0
      DGAMA_VP=R0
      c=R0
      FPT=R0
      T_SOLIDUS=R0
      T_LIQUIDUS=R0
      GAMMA_BAR=R0 
      RPL=R0
      RPL_VOL=R0
      RPL_SHEAR=R0
      DOT_VPSTRAN_1ST_INV=R0
      VPSTRAN_1ST_INV=R0
      VPSTRAN_2ND_INV=R0
      EPSILON_BAR_C=R0
      DTEMPERATURE_VOL=R0
      DTEMPERATURE_SHEAR=R0
      VISCOSITYY=R0
      PHI_YDP=R0
      DO I=1, NTENS
        EPSILON_VP(I)=STATEV(I)          		        !     SDV 1-4
      END DO
!     CALL HARDENING HISTORY FROM PREVIOUS TIME STEP 
      GAMMA_BAR=STATEV(NTENS+1)       			        !     SDV 5     ! DEFORMATION HISTORY 
      DGAMA_VISC_OLD=STATEV(NTENS+6)                    !     SDV 6     ! 
      VPSTRAN_1ST_INV=STATEV(NTENS+11)                  !     SDV 13    ! FIRST INVARIANT OF VISCOPLASTIC STRAIN
      VPSTRAN_2ND_INV=STATEV(NTENS+12)                  !     SDV 14    ! SECOND INVARIANT OF VISCOPLASTIC STRAIN 
!     INITIAL TEMPERATURE
      IF(KINC==1)THEN
        STATEV(NTENS+14)=TEMP
      END IF
      TEMP_INIT=STATEV(NTENS+14)
      
      EPSILON_V(1)=STATEV(NTENS+16)                    !     SDV 20    ! VISCOUS STRAINS(11)
      EPSILON_V(2)=STATEV(NTENS+17)                    !     SDV 21    ! VISCOUS STRAINS(22)
      EPSILON_V(3)=STATEV(NTENS+18)                    !     SDV 22    ! VISCOUS STRAINS(33)
      EPSILON_V(4)=STATEV(NTENS+19)                    !     SDV 23    ! VISCOUS STRAINS(12)  
!   ##################################################################################################################################
!     RETRIEVE RHEOLOGICAL PROPERTIES
!   ################################################################################################################################## 
        E=PROPS(1)                        		        !     YOUNG'S MODULUS 
        v=PROPS(2)                        		        !     POISSON'S RATIO
        c_0=PROPS(3)                      		        !     INITIAL COHESION
        H=PROPS(4)                        		        !     HARDENING MODULUS
        PWR=PROPS(5)                        	        !     STRESS EXPONENT 
        phi_i=PROPS(6)                      	        !     INITIAL FRICTION ANGLE
        phi_f=PROPS(7)      				            !     FINAL FRICTION ANGLE    
        RHOD_REF=PROPS(8)		        		        !     REFERENCE DENSITY
        C_P=PROPS(9)					                !     SPECIFIC HEAT CAPACITY
        E_A=PROPS(10)					                !     ACTIVATION ENERGY
        R=PROPS(11)					                    !     MOLECULAR GAS CONSTANT
        A=PROPS(12)                       		        !     CREEP PRE-EXPONENTIAL CONSTANT
        PLAST=PROPS(13)					                !     PLASTIC VISCOSITY
        psi=PROPS(14)                    		        !     DILATANCY ANGLE
        TREF=PROPS(15)					                !     REFERENCE TEMPERATURE
        ALPH_LE=PROPS(16)					            !     COEFFICIENT OF LINEAR EXPANSION
        EPSILON_BAR_C=PROPS(17)                         !     CRITICAL EFFECTIVE PLASTIC STRAIN
!   ##################################################################################################################################
!     COMPUTE QUANTITIES DEPENDENT ON THE RHEOLOGICAL PROPERTIES
!   ##################################################################################################################################      
!     SET FRICTION ANGLE AS A FUNCTION OF EFFECTIVE PLASTIC STRAIN RATE      
!      CONSTAN=R2*(SIND(phi_f)-SIND(phi_i))*SQRT(EPSILON_BAR_C)
!      phi=ASIND(SIND(phi_i)+((CONSTAN*VPSTRAN_EFF)/(VPSTRAN_EFF
!     1    +EPSILON_BAR_C)))
      EBULK=E/(R3*(R1-(R2*v)))          		! BULK MODULUS
      EBULK3=E/((R1-(R2*v)))
      G=(E/(R2*(R1+v)))                 		! SHEAR MODULUS
      R2G=R2*G
      ALPH_1=R3*TAND(phi_f)/(SQRT(R9+R12*(TAND(phi_f))**2))
      ALPH_2=R3/(SQRT(R9+R12*(TAND(phi_f))**2))
      ALPH_3=(R3*TAND(psi))/(SQRT(R9+R12*(TAND(psi))**2)) 
      DILATANCY_ANGLE=psi
!   ##################################################################################################################################
!     INITALIZE MATRICES AND ARRAYS
!   ##################################################################################################################################
!     C^e 	       ----> FOR THE ELASTIC STATE, DDSDDE IS THE ELASTICITY MATRIX
!     DVPRJT       ----> DEVIATORIC PROJECTION MATRIX
!     SOID         ----> SECOND ORDER IDENTITY MATRIX SAVED IN ARRAY FORM
!     SPRJD        ----> FOURTH ORDER SYMMETRIC IDENTITY MATRIX 
!     DEPSILON_VP  ----> VISCOPLASTIC STRAIN INCREMENT
!     TINC         ----> TEMPERATURE INCREMENT SAVED IN ARRAY FORM
      DO I=1,NTENS
        DO J=1,NTENS
        SPRJD(I,J)=R0
        DVPRJT(I,J)=R0
        END DO      
	    DEPSILON_VP(I)=R0
        DEPSILON_V(I)=R0
        SOID(I)=R0
        S(I)=R0
        ELAS(I)=R0
        STRIAL(I)=R0
        HC(I)=R0
      END DO   
!     INITIALIZE STRAIN INVARIANTS
      DOT_VPSTRAN_INV=R0
      DOT_VPSTRAN_1ST_INV=R0
      DOT_VPSTRAN_2ND_INV=R0
!     INITIALIZE TEMPERATURE CHANGE
      DT0=R0  
!   ##################################################################################################################################
!     ALLOCATE VALUES TO MATRICES AND ARRAYS
!     SPRJD        ----> FOURTH ORDER SYMMETRIC IDENTITY MATRIX 
!     SOID         ----> SECOND ORDER IDENTITY MATRIX SAVED IN ARRAY FORM
!   ##################################################################################################################################
      DO I=1,NDI
        SPRJD(I,I)=R1
        SOID(I)=R1
      END DO
      DO I=NDI+1,NTENS
        SPRJD(I,I)=P05
        SOID(I)=R0
      END DO  
!     DVPRJT       ----> DEVIATORIC PROJECTION MATRIX
      DO  M=1,NTENS
        DO  N=1,NTENS
        DVPRJT(M,N)=SPRJD(M,N)-(SOID(M)*SOID(N)/R3)
        END DO
      END DO
!   ##################################################################################################################################
!     DEFINE THE ELASTICITY MATRIC (C^e) I.E., THE CONSISTENT JACOBIAN (DDSDDE) FOR ELASTICITY PROBLEMS
!   ################################################################################################################################## 
      IF(NDI.EQ.3. AND. NSHR.EQ.1)THEN 
!     PLANE STRAIN/AXISYMMETRIC PROBLEMS
      	DDSDDE(1,1)=R1-v
      	DDSDDE(2,2)=R1-v
      	DDSDDE(3,3)=R1-v
      	DDSDDE(4,4)=P05*(R1-(R2*v))
      	DDSDDE(1,2)=v
      	DDSDDE(1,3)=v
      	DDSDDE(2,1)=v
      	DDSDDE(2,3)=v
      	DDSDDE(3,1)=v
      	DDSDDE(3,2)=v
      	DDSDDE=DDSDDE*E/((R1+v)*(R1-(R2*v)))
      ELSEIF(NDI.EQ.2 .AND. NSHR.EQ.1)THEN 
!     PLANE STRESS PROBLEMS  
      	DDSDDE(1,1)=R1
      	DDSDDE(2,2)=R1
      	DDSDDE(3,3)=P05*(R1-v)
      	DDSDDE(1,2)=v
      	DDSDDE(2,1)=v
      	DDSDDE=DDSDDE*E/(R1+(v*v))
      ELSE
!     3-D CONDITIONS
      	DDSDDE(1,1)=R1-v
      	DDSDDE(2,2)=R1-v
      	DDSDDE(3,3)=R1-v
      	DDSDDE(4,4)=P05*(R1-(R2*v))
      	DDSDDE(5,5)=P05*(R1-(R2*v))
      	DDSDDE(6,6)=P05*(R1-(R2*v))
      	DDSDDE(1,2)=v
      	DDSDDE(1,3)=v
      	DDSDDE(2,1)=v
      	DDSDDE(2,3)=v
      	DDSDDE(3,1)=v
      	DDSDDE(3,2)=v
      	DDSDDE=DDSDDE*E/((R1+v)*(R1-(R2*v)))
      END IF
!     ELASTIC STRESS STATE
!     DEPSILON        ----> DSTRAN
!     EPSILON^e       ----> ELAS 
      DO I = 1, NTENS
        DO J = 1, NTENS
        STRESS(I)=STRESS(I)+(DDSDDE(I,J)*(DSTRAN(J)))
        END DO
        ELAS(I)=STRAN(I)+DSTRAN(I)
      END DO      
      DOT_EPSTRAN_2ND_INV=R0
      DO M=1,NDI
         DOT_EPSTRAN_2ND_INV=DOT_EPSTRAN_2ND_INV+(DSTRAN(M)/DTIME)**R2
       END DO
      DO M=NDI+1,NTENS
        DOT_EPSTRAN_2ND_INV=DOT_EPSTRAN_2ND_INV+(R2*(DSTRAN(M)/DTIME)**R2)
      END DO
!   ##################################################################################################################################
!   ASSEMBLE THE TRIAL QUANTITIES (PRESSURE, DEVIATORIC STRESSES, J2 AND ELASTIC STRAIN MAGNITUDE)
!   ##################################################################################################################################
!     PTRIAL          ----> TRIAL ELASTIC PRESSURE
!     STRIAL          ----> TRIAL ELASTIC DEVIATORIC STRESS
!     DEVSTRAN        ----> TRIAL ELASTIC DEVIATORIC STRAIN
!     DEVSTRAN_NORM   ----> TRIAL ELASTIC DEVIATORIC STRAIN MAGNITUDE
!     DEVJ2           ----> J2
!     SQRTJ2          ----> SQUARE-ROOT OF J2
      PTRIAL=(STRESS(1)+STRESS(2)+STRESS(3))/R3
      DEVSTRAN_NORM=R0
      DEVSTRAN=R0
      DEVJ2=R0
      SQRTJ2=R0
      SQRTJ2T=R0
      DO I=1,NDI
        STRIAL(I)=STRESS(I)-(PTRIAL*SOID(I))
        DEVSTRAN(I)=STRIAL(I)/R2G
        DEVSTRAN_NORM=DEVSTRAN_NORM+DEVSTRAN(I)**2
        DEVJ2=DEVJ2+(STRIAL(I)*STRIAL(I))
      END DO
      DEVJ2=P05*DEVJ2
      DO I=NDI+1,NTENS
        STRIAL(I)=STRESS(I)-(PTRIAL*SOID(I))
        DEVSTRAN(I)=STRIAL(I)/G
        DEVSTRAN_NORM=DEVSTRAN_NORM+(R2*DEVSTRAN(I)**2)
        DEVJ2=DEVJ2+(STRIAL(I)*STRIAL(I))
      END DO
      DEVSTRAN_NORM=SQRT(DEVSTRAN_NORM)
      SQRTJ2=SQRT(DEVJ2)
      SQRTJ2T=SQRTJ2
      DOT_VPSTRAN_2ND_INV_ELAS=R0
      DO M=1,NDI
         DOT_VPSTRAN_2ND_INV_ELAS=DOT_VPSTRAN_2ND_INV_ELAS+(DSTRAN(M)/DTIME)**R2
      END DO
      DO M=NDI+1,NTENS
        DOT_VPSTRAN_2ND_INV_ELAS=DOT_VPSTRAN_2ND_INV_ELAS+(R2*(DSTRAN(M)/DTIME)**R2)
      END DO
      DOT_VPSTRAN_2ND_INV_ELAS=P05*DOT_VPSTRAN_2ND_INV_ELAS  
      P=PTRIAL
      DGAMA_VP=R0
      KK=1
!     c               ----> COHESION HARDENING
      c=c_0!+(H*GAMMA_BAR)  
!     PHI_Y           ----> DRUCKER-PRAGER YIELD CRITERION
!      PHI_Y=SQRTJ2+(ALPH_1*PTRIAL)-(ALPH_2*c)
!     CRP             ----> DISLOCATION CREEP
!      CRP=((A*EXP(-E_A/(R*TEMP))))
!      PWC=G+(ALPH_1*ALPH_3*EBULK)+(ALPH_2*ALPH_2*H)   
!   ##################################################################################################################################
!     CHECK IF THE MATERIAL YIELDS AND ACTIVATE PLASTICITY MODULE 
!     R_DGAMA         ----> RESIDUAL YIELD FUNCTION (R(DGAMA))
!     DR_DGAMA       ----> DERIVATIVE OF RESIDUAL YIELD FUNCTION 
!   ##################################################################################################################################
      IF(SQRTJ2.GT.1E1.AND.KSTEP.GT.R1)THEN
         CRP=A*EXP(-E_A/(R*TEMP))
!     INITIAL GUESSES FOR VISCOUS CREEP MULTIPLIERS
         DGAMA_OLD=R0
         DGAMA_VISC=R0
!     RESIDUAL VISCOUS CREEP EQUATION
      R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**PWR-DGAMA_VISC
      DR_DGAMA_VISC=-G*PWR*(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**(PWR-R1)-R1
      K=R1
      RESIDUALS=R1
!     ENTER NEWTON-RAPHSON LOOP TO COMPUTE THE VISCOUS CREEP MULTIPLIER
      	 DO WHILE(ABS(RESIDUALS).GT.TOL)
            R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**PWR-DGAMA_VISC
            DR_DGAMA_VISC=-G*PWR*(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**(PWR-R1)-R1
            DDGAMA=R_DGAMA_VISC/DR_DGAMA_VISC
            DGAMA_VISC=DGAMA_OLD-DDGAMA
            DJ_NEW=DGAMA_VISC
            IF(DGAMA_VISC.GE.SQRTJ2/G)THEN
              DJ_OLD=DGAMA_OLD
              DJ_NEW=DGAMA_VISC
              DO WHILE(RESIDUALS.GT.TOL .AND. DJ_NEW .GT. R0)
                DJ_NEW=(DJ_NEW-DGAMA_OLD)/R2
                R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DJ_NEW))**PWR-DJ_NEW
                RESIDUALS=ABS(R_DGAMA_VISC)
              END DO
            ELSEIF(DGAMA_VISC.LT.R0)THEN
              DJ_OLD=R0
              DJ_NEW=SQRTJ2/G
              DO WHILE (RESIDUALS.GT.TOL)
                DJ_NEW=(DJ_NEW-DJ_OLD)/R2
                R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DJ_NEW))**PWR-DJ_NEW
                RESIDUALS=ABS(R_DGAMA_VISC)
              END DO                
            END IF
            DGAMA_VISC=DJ_NEW
            RESIDUALS=ABS(R_DGAMA_VISC)
            DGAMA_OLD=DGAMA_VISC
            K=K+1 
         END DO
         IF(SQRTJ2.EQ.R0)THEN
            FACTOR=R1
         ELSE
            FACTOR=R1-(G*DGAMA_VISC/(SQRTJ2))
         END IF
         VSTRESS_DEV=R0
         VSTRESS_DEV_NORM=R0
         DEVJ2=R0
         DO I=1,NDI
           VSTRESS_DEV(I)=FACTOR*STRIAL(I)
           VSTRESS_DEV_NORM=VSTRESS_DEV_NORM+(VSTRESS_DEV(I))**R2
         END DO
         VSTRESS_DEV_NORM=P05*VSTRESS_DEV_NORM
         DO I=NDI+1,NTENS
           VSTRESS_DEV(I)=FACTOR*STRIAL(I)
           VSTRESS_DEV_NORM=VSTRESS_DEV_NORM+(VSTRESS_DEV(I))**R2
         END DO
         VSTRESS_DEV_NORM=SQRT(VSTRESS_DEV_NORM)
         SQRTJ2_V=VSTRESS_DEV_NORM
         DO I=1,NTENS
           STRESS(I)=VSTRESS_DEV(I)+PTRIAL*SOID(I)
           DEPSILON_V(I)=DGAMA_VISC*STRIAL(I)/(R2*SQRTJ2)
           EPSILON_V(I)=EPSILON_V(I)+DEPSILON_V(I)
         END DO
         DEVJ2=R0 
         SQRTJ2T=R0
         DO I=1,NDI
           DEVJ2=DEVJ2+(VSTRESS_DEV(I)*VSTRESS_DEV(I))
         END DO
         DEVJ2=P05*DEVJ2
         DO I=NDI+1,NTENS
           DEVJ2=DEVJ2+(VSTRESS_DEV(I)*VSTRESS_DEV(I))
         END DO
         SQRTJ2T=SQRT(DEVJ2) 
!     ASSEMBLE CONSISTENT JACOBIAN MATRIX FOR CREEP DEFORMATION
        b1=(SQRT(R2)*G*(DTIME*CRP)**(R1/PWR))/
     1     (((DGAMA_VISC)**((R1-PWR)/PWR))/PWR+G*(DTIME*CRP)**(R1/PWR))
         FVP1=R2G*(R1-DGAMA_VISC/(SQRT(R2)*DEVSTRAN_NORM))
         FVP2=SQRT(R2)*G*(DGAMA_VISC/DEVSTRAN_NORM-b1)
         DO M=1,NTENS
            DO N=1,NTENS
              DEV=FVP1*DVPRJT(M,N)+
     1        FVP2*DEVSTRAN(M)*DEVSTRAN(N)/(DEVSTRAN_NORM*DEVSTRAN_NORM)
              VOL=EBULK*SOID(M)*SOID(N)
              DDSDDE_V(M,N)=DEV
              DDSDDE(M,N)=DEV+VOL
            END DO
          END DO
!     CHECK FOR YIELD (VISCOPLASTICITY)
         PHI_YDP=SQRTJ2_V+(ALPH_1*PTRIAL)-(ALPH_2*c_0)
         c1=G+(ALPH_1*ALPH_3*EBULK)
         DEPSILON_VP=R0
         IF(PHI_YDP .GT. R0)THEN
!     INITIAL GUESSES FOR VISCOPLASTIC MULTIPLIERS
            DGAMA_OLD=INIT_GUESS
            DGAMA_VP=INIT_GUESS
!     RESIDUAL VISCOPLASTIC RETURN MAPPING EQUATION
            R_DGAMA_VP=(DTIME/PLAST)*((PHI_YDP-(c1*DGAMA_VP))/c_0)**PWR-DGAMA_VP               
            DR_DGAMA_VP=-c1*PWR*DTIME/(c_0*PLAST)*((PHI_YDP-(c1*DGAMA_VP))
     1                  /c_0)**(PWR-R1)-R1
            RESIDUALS=R1 
!     ENTER NEWTON-RAPHSON LOOP TO COMPUTE THE VISCOPLASTIC MULTIPLIER
!     BY SOLVING THE RESIDUAL VISCOPLASTIC FUNCTION
            DO WHILE(ABS(RESIDUALS).GT.TOL)
                R_DGAMA_VP=(DTIME/PLAST)*((PHI_YDP-(c1*DGAMA_VP))/c_0)**PWR-DGAMA_VP               
                DR_DGAMA_VP=-c1*PWR*DTIME/(c_0*PLAST)*((PHI_YDP-(c1*DGAMA_VP))
     1                  /c_0)**(PWR-R1)-R1
                DDGAMA=R_DGAMA_VP/DR_DGAMA_VP
                DGAMA_VP=DGAMA_OLD-DDGAMA
                DJ_NEW=DGAMA_VP
                IF(DJ_NEW.GE.PHI_YDP/c1)THEN
                    DJ_OLD=DGAMA_OLD
                    DJ_NEW=DGAMA_VP
                    DO WHILE(RESIDUALS.GT.TOL .AND. DJ_NEW .GT. R0)
                        DJ_NEW=(DJ_NEW-DGAMA_OLD)/R2
                        R_DGAMA_VP=(DTIME/PLAST)*((PHI_YDP-(c1*DJ_NEW))/c_0)**PWR-DJ_NEW
                        RESIDUALS=ABS(R_DGAMA_VP)
                    END DO
                ELSEIF(DGAMA_VP.LE.R0)THEN
                    DJ_OLD=R0
                    DJ_NEW=PHI_YDP/c1
                    DO WHILE(RESIDUALS.GT.TOL)
                        DJ_NEW=(DJ_NEW-DJ_OLD)/R2
                        R_DGAMA_VP=(DTIME/PLAST)*((PHI_YDP-(c1*DJ_NEW))/c_0)**PWR-DJ_NEW
                        RESIDUALS=ABS(R_DGAMA_VP)
                    END DO                
                END IF
                DGAMA_VP=DJ_NEW
                RESIDUALS=ABS(R_DGAMA_VP)
                DGAMA_OLD=DGAMA_VP
                K=K+1
            END DO
!     UPDATE STATE VARIABLES
            IF (SQRTJ2_V.EQ.R0)THEN
                FACTOR=R1
            ELSE
                FACTOR=R1-(G*DGAMA_VP/(SQRTJ2_V))
            END IF
!     UPDATE PRESSURE
            P=PTRIAL-(ALPH_3*EBULK*DGAMA_VP)
!     UPDATE DEVIATORIC STRESSES, STRESSES, INCREMENTS  
            DO I=1,NTENS
                S(I)=FACTOR*VSTRESS_DEV(I)
                STRESS(I)=S(I)+(P*SOID(I))        
                DEPSILON_VP(I)=DGAMA_VP*(VSTRESS_DEV(I)/(R2*SQRTJ2_V)
     1                         +(ALPH_3/R3*SOID(I)))
            END DO
            DEVJ2=R0 
            SQRTJ2T=R0
            DO I=1,NDI
                DEVJ2=DEVJ2+(S(I)*S(I))
            END DO
            DEVJ2=P05*DEVJ2
            DO I=NDI+1,NTENS
                DEVJ2=DEVJ2+(S(I)*S(I))
            END DO
            SQRTJ2T=SQRT(DEVJ2)
!     ASSEMBLE CONSISTENT JACOBIAN MATRIX FOR RETURN MAPPING TO THE SMOOTH PART OF THE DRUCKER-PRAGER CONE
            PHI_YDPR=SQRTJ2_V+(ALPH_1*PTRIAL)-(ALPH_2*c_0)-c1*DGAMA_VP
            b2=R1/c_0*(DTIME/PLAST)**(R1/PWR)/
     1        (((DGAMA_VP)**((R1-PWR)/PWR))
     2         /PWR+c1/c_0*(DTIME/PLAST)**(R1/PWR))
            FVP1=R1-G*DGAMA_VP/SQRTJ2_V
            FVP2=ALPH_1*b2*EBULK*G/(SQRTJ2_V)
            FVP3=G/SQRTJ2_V
            FVP4=SQRT(R2)*G*b2*(R1-b1/R2)
            FVP5=SQRT(R2)*DGAMA_VP/SQRTJ2_V*(b1/SQRT(R2)-R1)
            FVP6=R1-ALPH_1*ALPH_3*b2*EBULK
            FVP7=ALPH_3*b2*EBULK*(b1-SQRT(R2)*G)
            DO M=1,NTENS
                DO N=1,NTENS
                    DEV=FVP1*DDSDDE_V(M,N)-
     1                  FVP2*SOID(M)*VSTRESS_DEV(N)/(SQRTJ2_V)
     2                  -FVP3*(FVP4*DEVSTRAN(M)/DEVSTRAN_NORM
     3                  +(FVP5*DEVSTRAN(M)/DEVSTRAN_NORM))
     4                   *VSTRESS_DEV(N)/(SQRTJ2_V) 
                    VOL=EBULK*FVP6*SOID(M)*SOID(N)
     1                     +(FVP7*DEVSTRAN(M*SOID(M))/DEVSTRAN_NORM) 
                    DDSDDE(M,N)=DEV+VOL	           
            END DO
          END DO
!   UPDATE DEFORMATION AND HEAT GENERATION
      RPL=R0    
      RPL_SHEAR=R0  
      RPL_VOL=R0
      VPSTRAN_1ST_INV=R0 
      DOT_VPSTRAN_2ND_INV=R0
      VPSTRAN_2ND_INV=R0
      DOT_VPSTRAN_EFF=R0
      DO M=1,NTENS
        EPSILON_VP(M)=EPSILON_VP(M)+DEPSILON_VP(M)
        RPL=RPL+(STRESS(M)*(DEPSILON_VP(M))/DTIME)
        VPSTRAN_1ST_INV=VPSTRAN_1ST_INV+EPSILON_VP(M)*SOID(M)
        DOT_VPSTRAN_1ST_INV=DOT_VPSTRAN_1ST_INV+(DEPSILON_VP(M))
     1                      /DTIME*SOID(M)
      END DO    
      DO M=1,NDI
         DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+((DEPSILON_VP(M))
     1                      /DTIME)**R2
         VPSTRAN_2ND_INV=VPSTRAN_2ND_INV+(EPSILON_VP(M))**R2
      END DO
                
      DO M=NDI+1,NTENS
        DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+
     1                      (R2*((DEPSILON_VP(M))/DTIME)**R2)
        VPSTRAN_2ND_INV=VPSTRAN_2ND_INV+(R2*(EPSILON_VP(M))**R2)
      END DO
      DOT_VPSTRAN_2ND_INV=P05*DOT_VPSTRAN_2ND_INV
      VPSTRAN_2ND_INV=P05*VPSTRAN_2ND_INV
      DOT_VPSTRAN_2ND_INV=SQRT(DOT_VPSTRAN_2ND_INV)
      VPSTRAN_2ND_INV=SQRT(VPSTRAN_2ND_INV)        
      ELSE
!   UPDATE DEFORMATION AND HEAT GENERATION
      RPL=R0    
      RPL_SHEAR=R0  
      RPL_VOL=R0
      VPSTRAN_1ST_INV=R0 
      DOT_VPSTRAN_2ND_INV=R0
      VPSTRAN_2ND_INV=R0
      DOT_VPSTRAN_EFF=R0
      DO M=1,NTENS
        EPSILON_VP(M)=EPSILON_VP(M)+DEPSILON_V(M)
        RPL=RPL+(STRESS(M)*(DEPSILON_V(M))/DTIME)
        VPSTRAN_1ST_INV=VPSTRAN_1ST_INV+EPSILON_VP(M)*SOID(M)
        DOT_VPSTRAN_1ST_INV=DOT_VPSTRAN_1ST_INV+(DEPSILON_V(M))
     1                      /DTIME*SOID(M)
      END DO    
      DO M=1,NDI
         DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+((DEPSILON_V(M))
     1                      /DTIME)**R2
         VPSTRAN_2ND_INV=VPSTRAN_2ND_INV+(EPSILON_VP(M))**R2
      END DO
                
      DO M=NDI+1,NTENS
        DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+
     1                      (R2*((DEPSILON_V(M))/DTIME)**R2)
        VPSTRAN_2ND_INV=VPSTRAN_2ND_INV+(R2*(EPSILON_V(M))**R2)
      END DO
      DOT_VPSTRAN_2ND_INV=P05*DOT_VPSTRAN_2ND_INV
      VPSTRAN_2ND_INV=P05*VPSTRAN_2ND_INV
      DOT_VPSTRAN_2ND_INV=SQRT(DOT_VPSTRAN_2ND_INV)
      VPSTRAN_2ND_INV=SQRT(VPSTRAN_2ND_INV)   
      END IF
      END IF     
!   ##################################################################################################################################
!     END OF (VISCO)PLASTICITY ROUTINES  
!   ##################################################################################################################################      
      DOT_VPSTRAN_1ST_INV=R0
      DOT_VPSTRAN_2ND_INV=R0
      DO M=1,NDI
         DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+((DEPSILON_V(M))/
     1              DTIME*(DEPSILON_V(M))/DTIME)+
     2             ((DEPSILON_VP(M))/DTIME*(DEPSILON_VP(M))/DTIME)+
     3              (DSTRAN(M)/DTIME*DSTRAN(M)/DTIME)
        DOT_VPSTRAN_1ST_INV=DOT_VPSTRAN_1ST_INV+(DEPSILON_V(M))
     1               /DTIME+(DEPSILON_VP(M))/DTIME+(DSTRAN(M))/DTIME
      END DO                
      DO M=NDI+1,NTENS
        DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+
     1             R2*(((DEPSILON_V(M))/DTIME)*((DEPSILON_V(M))/DTIME)+
     2                 ((DEPSILON_VP(M))/DTIME*(DEPSILON_VP(M))/DTIME)+
     3                 DSTRAN(M)/DTIME*DSTRAN(M)/DTIME)
      END DO
      DOT_VPSTRAN_2ND_INV=P05*DOT_VPSTRAN_2ND_INV
      DOT_VPSTRAN_2ND_INV=SQRT(DOT_VPSTRAN_2ND_INV)
!   ##################################################################################################################################
      VPSTRAN_1ST_INV=R0 
      VPSTRAN_2ND_INV=R0
      DO M=1,NTENS
        VPSTRAN_1ST_INV=VPSTRAN_1ST_INV+EPSILON_VP(M)*SOID(M)
      END DO    
      DO M=1,NDI
         VPSTRAN_2ND_INV=VPSTRAN_2ND_INV+(EPSILON_VP(M))**R2
      END DO                
      DO M=NDI+1,NTENS
        VPSTRAN_2ND_INV=VPSTRAN_2ND_INV+(R2*(EPSILON_VP(M))**R2)
      END DO
      VPSTRAN_2ND_INV=SQRT(P05*VPSTRAN_2ND_INV)

      PRESS=-R1*(STRESS(1)+STRESS(2)+STRESS(3))/R3
      PRESS=PRESS/1E9 
      T_SOLIDUS=A11+(A22*PRESS)+(A33*PRESS**2)         ! FROM KATZ, SPIEGELMAN & LANGMUIR (2003), A NEW PARAMETERIZATION FOR HYDROUS MANTLE MELTING
                                                       ! GCUBED
      T_LIQUIDUS=B11+(B22*PRESS)+(B33*PRESS**2)        ! FROM KATZ, SPIEGELMAN & LANGMUIR (2003), A NEW PARAMETERIZATION FOR HYDROUS MANTLE MELTING
                                                       ! GCUBED                                                 
      T_SOLIDUS=T_SOLIDUS+273.15D0
      T_LIQUIDUS=T_LIQUIDUS+273.15D0       
!      COMPUTE FRACTION OF MELT AT TIME STEP       
       IF(T_SOLIDUS.LT.TEMP.AND.TEMP.LT.T_LIQUIDUS)THEN
          FPT=((TEMP-T_SOLIDUS)/(T_LIQUIDUS-T_SOLIDUS))**PWRFPT
       ELSE
          FPT=R0
       END IF  
!     UPDATE STATE VARIABLES
      MAXSTRESS=MAX(STRESS(1), STRESS(2), STRESS(3))
      MINSTRESS=MIN(STRESS(1), STRESS(2), STRESS(3))
      DIFF_STRESS=MAXSTRESS-MINSTRESS
      DO I=1,NTENS
        STATEV(I)=EPSILON_VP(I)                        !     SDV 1-4   ! VISCOPLASTIC STRAINS (11,22,33,12)
      END DO 
      STATEV(NTENS+1)=GAMMA_BAR                        !     SDV 5     ! HARDENING HISTORY
      STATEV(NTENS+2)=PHI_YDP       	               !     SDV 6     ! DRUCKER-PRAGER YIELD CRITERION
      STATEV(NTENS+3)=SQRTJ2                           !     SDV 7     ! ELASTIC J2
      STATEV(NTENS+4)=SQRTJ2T		                   !     SDV 8     ! VISCOPLASTIC J2
      STATEV(NTENS+5)=PRESS*1E9       		           !     SDV 9     ! PRESSURE
      STATEV(NTENS+6)=DGAMA_VP          		       !     SDV 10    ! VISCOPLASTIC SMOOTH CONE MULTIPLIER
      STATEV(NTENS+7)=DOT_VPSTRAN_1ST_INV              !     SDV 11    ! FIRST INVARIANT OF VISCOPLASTIC STRAIN RATE
      STATEV(NTENS+8)=LOG10(DOT_VPSTRAN_2ND_INV)       !     SDV 12    ! SECOND INVARIANT OF VISCOPLASTIC STRAIN RATE 
      STATEV(NTENS+9)=VPSTRAN_1ST_INV                  !     SDV 13    ! FIRST INVARIANT OF VISCOPLASTIC STRAIN
      STATEV(NTENS+10)=VPSTRAN_2ND_INV            	   !     SDV 14    ! SECOND INVARIANT OF VISCOPLASTIC STRAIN 
      STATEV(NTENS+11)=DIFF_STRESS                     !     SDV 15    ! DIFFERENTIAL STRESS      
      STATEV(NTENS+12)=RPL                             !     SDV 16    ! HEAT GENERATION PER UNIT TIME
      STATEV(NTENS+13)=DTEMP/DTIME                     !     SDV 17    ! CHECK FOR THERMAL STEADY STATE 
      STATEV(NTENS+15)=TEMP-TEMP_INIT                  !     SDV 19    ! THERMAL PERTURBATIONS
      STATEV(NTENS+16)=DEPSILON_V(1)                   !     SDV 20    ! VISCOUS STRAIN INCREMENT(11)
      STATEV(NTENS+17)=DEPSILON_V(2)                   !     SDV 21    ! VISCOUS STRAIN INCREMENT(22)
      STATEV(NTENS+18)=DEPSILON_V(3)                   !     SDV 22    ! VISCOUS STRAIN INCREMENT(33)
      STATEV(NTENS+19)=DEPSILON_V(4)                   !     SDV 23    ! VISCOUS STRAIN INCREMENT(12)
      STATEV(NTENS+20)=DEPSILON_VP(1)                  !     SDV 24    ! VISCOPLASTIC STRAIN INCREMENT(11)
      STATEV(NTENS+21)=DEPSILON_VP(2)                  !     SDV 25    ! VISCOPLASTIC STRAIN INCREMENT(22)
      STATEV(NTENS+22)=DEPSILON_VP(3)                  !     SDV 26    ! VISCOPLASTIC STRAIN INCREMENT(33)
      STATEV(NTENS+23)=DEPSILON_VP(4)                  !     SDV 27    ! VISCOPLASTIC STRAIN INCREMENT(12)
      STATEV(NTENS+24)=FPT                             !     SDV 28    ! FRACTION OF MELT      
      RETURN
      END      
!   ##################################################################################################################################
!   ##################################################################################################################################
!     END OF USER MATERIAL SUBROUTINE  
!   ##################################################################################################################################
!   ################################################################################################################################## 
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! 
      INCLUDE 'ABA_PARAM.INC'
      COMMON /CONDUCTIVITY_PARAMS/ TREF, ALPH_LE, RHOD_REF, TEMPERA, PRESS
      CHARACTER*80 CMNAME
       
       DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
       
       PARAMETER (EXMULT=0.0000004D0)
!       IF(CMNAME.EQ.'CRUST')THEN
!          COND=1.18D0+474.0D0/(TEMP+77.0D0)
!       ELSEIF(CMNAME.EQ.'LITHOSPHERE'.OR.CMNAME.EQ.'WEAKSHEARZONE')THEN
!          COND=0.73D0+1293.0D0/(TEMP+77.0D0)
!      ELSEIF(CMNAME.EQ.'Matrix'.OR.CMNAME.EQ.'Inclusion')THEN
!          COND=1.72D0+807.0D0/(TEMP+350.0D0)
!       ELSE
!          COND=0.73D0+1293.0D0/(TEMP+77.0D0)
!       END IF
       COND = PROPS(1)
       C_P = PROPS(2)
       RHOD = PROPS(3)
       DUDT = C_P
       DU = C_P*DTEMP
       U = U+DU
       DO I=1,NTGRD
         FLUX(I) = -COND*DTEMDX(I)
         DFDG(I,I) = -COND
       END DO
      RETURN 
      END
!   ##################################################################################################################################
!   ##################################################################################################################################
