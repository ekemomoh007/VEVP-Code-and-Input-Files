!   ##################################################################################################################################
!   USER SUBROUTINE TO COMPUTE DRUCKER-PRAGER VISCOPLASTICITY WITH
!   CONSISTENT JACOBIAN OBTAINED BY LINEARIZING THE RETURN MAPPING
!   EQUATION (RESIDUAL VISCOPLASTIC POTENTIAL FUNCTION)
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
      COMMON /CONDUCTIVITY_PARAMS/ TREF, ALPH_LE, RHOD_REF, TEMPERA, HL,
     1 PRESS, NCOMPS
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DIMENSION DEVSTRAN(NTENS),S(NTENS),STRIAL(NTENS),SOID(NTENS),
     1          DVPRJT(NTENS,NTENS),VSTRESS_DEV(NTENS),EPSILON_VP(NTENS),
     2          DEPSILON_VP(NTENS), SPRJD(NTENS,NTENS),ELAS(NTENS),
     3          HC(NTENS),DEPSILON_V(NTENS), EPSILON_V(NTENS), 
     4          DDSDDE_V(NTENS,NTENS)
      PARAMETER (R0=0.0D0, P01=0.1D0, R1=1.0D0, R2=2.0D0, R3=3.0D0,
     1           R4=4.0D0, R6=6.0D0, R7=7.0D0,R9=9.0D0, R10=10.0D0,
     2           R40=40.0D0,P13=1.0D0/3.0D0, TOL=1E-6, INIT_GUESS=1E-21,
     3           P05=0.5D0,P098=0.98D0,P23=2.0D0/3.0D0,P06=6.0D0/10.0D0,
     4           P07=7.0D0/10.0D0, A11=1085.7D0,A22=132.9D0,A33=-5.1D0,
     5           B11=1475.0D0,B22=80.D0,B33=-3.2D0, PWRFPT=1.5D0, 
     6           R12=12.0D0, TQC=1.0D0)
      DOUBLE PRECISION FPT
!   ##################################################################################################################################
!     INITIALIZE STATE VARIABLES VARIABLES
!   ##################################################################################################################################
       CREEP_COUNT=R0                                     
       PLAST_COUNT=R0
       NCOMPS=NTENS
       HGR=R0
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
       EPSILON_BAR=R0
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
         EPSILON_VP(I)=STATEV(I)          		      !     SDV 1-6
       END DO
!     CALL HARDENING HISTORY FROM PREVIOUS TIME STEP
        EPSILON_BAR=STATEV(NTENS+1)       			  !     SDV 7     ! DEFORMATION HISTORY
        VPSTRAN_1ST_INV=STATEV(NTENS+9)               !     SDV 17    ! FIRST INVARIANT OF VISCOPLASTIC STRAIN
        VPSTRAN_2ND_INV=STATEV(NTENS+10)              !     SDV 18    ! SECOND INVARIANT OF VISCOPLASTIC STRAIN
!     INITIAL TEMPERATURE
        IF(KINC==1)THEN
          STATEV(NTENS+14)=TEMP
        END IF
        TEMP_INIT=STATEV(NTENS+14)
!   ##################################################################################################################################
!     RETRIEVE RHEOLOGICAL PROPERTIES
!   ##################################################################################################################################
        E=PROPS(1)                        		      !     YOUNG'S MODULUS
        v=PROPS(2)                        		      !     POISSON'S RATIO
        c_0=PROPS(3)                      	 	      !     INITIAL COHESION (THIS IS CONSTANT IF COHESION HARDENING OR SOFTENING IS USED)
        H=PROPS(4)                        		      !     HARDENING MODULUS (IF HARDENING IS USED.)
        PWR=PROPS(5)                        	      !     STRESS EXPONENT
        phi_i=PROPS(6)                      	      !     INITIAL FRICTION ANGLE 
        phi_f=PROPS(7)      				          !     FINAL FRICTION ANGLE   
        RHOD_REF=PROPS(8)		                      !     REFERENCE DENSITY (IF BOUSSINESQ IS USED)
        C_P=PROPS(9)					              !     SPECIFIC HEAT CAPACITY
        E_A=PROPS(10)					              !     ACTIVATION ENERGY
        R=PROPS(11)					                  !     MOLECULAR GAS CONSTANT
        A=PROPS(12)                       		      !     CREEP PRE-EXPONENTIAL CONSTANT
        PLAST=PROPS(13)					              !     RELATIVE RATE OF VISCOPLASTIC STRAIN (INVERSE)
        psi=PROPS(14)                    		      !     DILATANCY ANGLE
        TREF=PROPS(15)					              !     REFERENCE TEMPERATURE
        ALPH_LE=PROPS(16)			                  !     COEFFICIENT OF LINEAR EXPANSION (USED HERE)
        EPSILON_BAR_C=PROPS(17)                       !     CRITICAL EFFECTIVE PLASTIC STRAIN (USED HERE)
        LH=PROPS(20)                                  !     LATENT HEAT CAPACITY
        RADIOGENIC=PROPS(21)                          !     RADIOGENIC HEAT PRODUCTION
!   ##################################################################################################################################
!     COMPUTE QUANTITIES DEPENDENT ON THE RHEOLOGICAL PROPERTIES
!   ##################################################################################################################################
!     SET FRICTION ANGLE AS A FUNCTION OF EFFECTIVE PLASTIC STRAIN RATE (FRICTION HARDENING)
      VPSTRAN_EFF=STATEV(NTENS+35)
      IF(VPSTRAN_EFF.GT.R0)THEN
         CONSTAN=R2*(SIND(phi_f)-SIND(phi_i))*SQRT(EPSILON_BAR_C)
          phi=ASIND(SIND(phi_i)+((CONSTAN*SQRT(VPSTRAN_EFF))/(VPSTRAN_EFF
     1        +EPSILON_BAR_C)))
      ELSE
        phi=phi_i
      END IF
      EBULK=E/(R3*(R1-(R2*v)))          		      ! BULK MODULUS
      EBULK3=E/((R1-(R2*v)))
      G=(E/(R2*(R1+v)))                   		      ! SHEAR MODULUS
      R2G=R2*G
!     THE DEFINITIONS FOR THE FRICTION-DEPENDENT PARAMETERS ARE FOUND IN GLERUM ET AL., 2017:
!     Implementing nonlinear viscoplasticity in ASPECT: benchmarking and applications to 3D subduction modeling, Solid Earth
!     EQUATION 6.122 OF DE SOUZA NETO'S COMPUTATIONAL PLASTICITY BOOK
      ALPH_1=R6*SIND(phi)/(SQRT(R3)*(R3-(SIND(phi))))
      ALPH_2=R6*COSD(phi)/(SQRT(R3)*(R3-(SIND(phi))))
      ALPH_3=R6*SIND(psi)/(SQRT(R3)*(R3-(SIND(psi))))
      DILATANCY_ANGLE=psi
      VISCOUS=R0
      PLASTIC=R0
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
      ETS=ELAS(1)+ELAS(2)+ELAS(3)+(R3*(ALPH_LE*(TEMP-TREF)))
      DOT_EPSTRAN_2ND_INV=R0
      DO M=1,NDI
         DOT_EPSTRAN_2ND_INV=DOT_EPSTRAN_2ND_INV+(DSTRAN(M)/DTIME)**R2
       END DO
      DO M=NDI+1,NTENS
        DOT_EPSTRAN_2ND_INV=DOT_EPSTRAN_2ND_INV+
     1      (R2*(DSTRAN(M)/DTIME)**R2)
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
      IF(KSTEP .GT. R1)THEN
      CRP=A*EXP(-E_A/(R*TEMP))
!     INITIAL GUESSES FOR VISCOUS CREEP MULTIPLIERS
!     RESIDUAL VISCOUS CREEP EQUATION
      DGAMA_VISC=R0
      DGAMA_OLD=DGAMA_VISC
      R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**PWR-DGAMA_VISC
      DR_DGAMA_VISC=-G*PWR*(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**(PWR-R1)-R1
      RESIDUALS=ABS(SQRTJ2)
      K=R0
       DO WHILE(ABS(RESIDUALS).GT.TOL)
        DR_DGAMA_VISC=-G*PWR*(DTIME*CRP)*(SQRTJ2-(G*DGAMA_VISC))**(PWR-R1)-R1
        DDGAMA=R_DGAMA_VISC/DR_DGAMA_VISC
        DGAMA_VISC=DGAMA_OLD-DDGAMA
        DJ_NEW=DGAMA_VISC
        IF(DGAMA_VISC.GT.SQRTJ2/G)THEN
            DJ_OLD=R0
            DJ_NEW=SQRTJ2/G!DTIME*CRP*SQRTJ2/(R1+(DTIME*G*CRP))
            R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DJ_NEW))**PWR-DJ_NEW
            RESIDUALS=ABS(R_DGAMA_VISC)
            DO WHILE(RESIDUALS.GT.TOL)
              DJ_NEW=(DJ_NEW-DGAMA_OLD)/R2
              R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DJ_NEW))**PWR-DJ_NEW
              RESIDUALS=ABS(R_DGAMA_VISC)
            END DO
            ELSEIF(DGAMA_VISC.LT.R0)THEN
              DJ_OLD=R0
              DJ_NEW=DTIME*CRP*SQRTJ2/(R1+(DTIME*G*CRP))
              R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DJ_NEW))**PWR-DJ_NEW
              RESIDUALS=ABS(R_DGAMA_VISC)
              DO WHILE (RESIDUALS.GT.TOL)
                DJ_NEW=(DJ_NEW-DJ_OLD)/R2
                R_DGAMA_VISC=(DTIME*CRP)*(SQRTJ2-(G*DJ_NEW))**PWR-DJ_NEW
                RESIDUALS=ABS(R_DGAMA_VISC)
              END DO
              DGAMA_VISC=DJ_NEW
         END IF
            R_DGAMA_VISC=(DTIME*CRP)*
     1                  (SQRTJ2-(G*DGAMA_VISC))**PWR-DGAMA_VISC
            RESIDUALS=ABS(R_DGAMA_VISC)
            DGAMA_OLD=DGAMA_VISC
            K=K+1
        END DO
        CREEP_COUNT=K
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
!   ##################################################################################################################################
!     CHECK IF THE MATERIAL YIELDS AND ACTIVATE PLASTICITY MODULE
!     R_DGAMA         ----> RESIDUAL YIELD FUNCTION (R(DGAMA))
!     DR_DGAMA        ----> DERIVATIVE OF RESIDUAL YIELD FUNCTION
!   ##################################################################################################################################
!     WE HAVE INCLUDED HARDENING FORMULATIONS HERE
!     c               ----> COHESION HARDENING
      c=c_0+(H*EPSILON_BAR)
!     PHI_Y           ----> DRUCKER-PRAGER YIELD CRITERION
!     PHI_YDP=SQRTJ2_V+(ALPH_1*PTRIAL)-(ALPH_2*c_0)
      PHI_YDP=SQRTJ2_V+(ALPH_1*PTRIAL)-(ALPH_2*c)
      c1=G+(ALPH_1*ALPH_3*EBULK)+(ALPH_2*ALPH_2*H)
      DEPSILON_VP=R0
      IF(PHI_YDP.GT.R0)THEN
!     INITIAL GUESSES FOR VISCOPLASTIC MULTIPLIERS
!     RESIDUAL VISCOPLASTIC RETURN MAPPING EQUATION
!        DGAMA_VP=(PHI_YDP/c_0)/((PLAST/DTIME)+(c1/c_0))
!        DGAMA_VP=(PHI_YDP*DTIME)/((PLAST*c_0)+(DTIME*c1))
         DGAMA_VP=R0
        R_DGAMA_VP=(DTIME/PLAST)*
     1              ((PHI_YDP-(c1*DGAMA_VP))/c_0)**PWR-DGAMA_VP
        DR_DGAMA_VP=-c1*PWR*DTIME/(c_0*PLAST)*((PHI_YDP-(c1*DGAMA_VP))
     1                  /c_0)**(PWR-R1)-R1
        RESIDUALS=ABS(PHI_YDP)
        DGAMA_OLD=DGAMA_VP
        K=R0
!     ENTER NEWTON-RAPHSON LOOP TO COMPUTE THE VISCOPLASTIC MULTIPLIER
!     BY SOLVING THE RESIDUAL VISCOPLASTIC FUNCTION
        DO WHILE(ABS(RESIDUALS).GT.TOL)
            DDGAMA=R_DGAMA_VP/DR_DGAMA_VP
            DGAMA_VP=DGAMA_OLD-DDGAMA
            IF(DGAMA_VP.LE.R0)THEN
               DJ_OLD=R0
               DJ_NEW=(PHI_YDP*DTIME)/((PLAST*c_0)+(DTIME*c1))!PHI_YDP/c1
               RESIDUALS=ABS((DTIME/PLAST)*
     1                      ((PHI_YDP-(c1*DJ_NEW))/c_0)**PWR-DJ_NEW)
               DO WHILE(RESIDUALS.GT.TOL)
                 DJ_NEW=(DJ_NEW-DJ_OLD)/R2
                 RESIDUALS=ABS((DTIME/PLAST)*
     1                       ((PHI_YDP-(c1*DJ_NEW))/c_0)**PWR-DJ_NEW)
               END DO
               DGAMA_VP=DJ_NEW
               R_DGAMA_VP=(DTIME/PLAST)*
     1                   ((PHI_YDP-(c1*DGAMA_VP))/c_0)**PWR-DGAMA_VP
            ELSE
               R_DGAMA_VP=(DTIME/PLAST)*
     1                     ((PHI_YDP-(c1*DGAMA_VP))/c_0)**PWR-DGAMA_VP
               RESIDUALS=ABS((DTIME/PLAST)*
     1                     ((PHI_YDP-(c1*DGAMA_VP))/c_0)**PWR-DGAMA_VP)
            END IF
            DR_DGAMA_VP=-c1*PWR*DTIME/(c_0*PLAST)*((PHI_YDP-
     1                  (c1*DGAMA_VP))/c_0)**(PWR-R1)-R1
            DGAMA_OLD=DGAMA_VP
            K=K+1
        END DO
        PLAST_COUNT=K
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
!    UPDATE ACCUMMULATED VISCOPLASTIC STRAINS
            EPSILON_BAR=EPSILON_BAR+(DGAMA_VP*ALPH_2) ! SEE EQUATION 6.205 ON PAGE 184 OF DE SOUZA NETO ET AL (2008): COMPUTATIONAL METHODS FOR PLASTICITY: THEORY AND APPLICATIONS
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
            PHI_YDPR=SQRTJ2_V+(ALPH_1*PTRIAL)-(ALPH_2*c)-c1*DGAMA_VP
            b2=R1/c*(DTIME/PLAST)**(R1/PWR)/
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
      END IF
!   ##################################################################################################################################
!    PLASTICITY ROUTINE ENDS HERE
!   ##################################################################################################################################      
!   ##################################################################################################################################
!    COMPUTE STRAIN RATE INVARIANTS
!   ##################################################################################################################################
      DOT_VPSTRAN_1ST_INV=R0
      DOT_VPSTRAN_2ND_INV=R0
      DOT_VPSTRAN_2ND_INV_IN=R0
      DO M=1,NDI
         DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+((DEPSILON_V(M))/
     1              DTIME*(DEPSILON_V(M))/DTIME)+
     2             ((DEPSILON_VP(M))/DTIME*(DEPSILON_VP(M))/DTIME)+
     3              (DSTRAN(M)/DTIME*DSTRAN(M)/DTIME)
        DOT_VPSTRAN_1ST_INV=DOT_VPSTRAN_1ST_INV+(DEPSILON_V(M))
     1               /DTIME+(DEPSILON_VP(M))/DTIME+(DSTRAN(M))/DTIME
        DOT_VPSTRAN_2ND_INV_IN=DOT_VPSTRAN_2ND_INV_IN+((DEPSILON_V(M))/
     1              DTIME*(DEPSILON_V(M))/DTIME)+
     2             ((DEPSILON_VP(M))/DTIME*(DEPSILON_VP(M))/DTIME)
      END DO
      DO M=NDI+1,NTENS
        DOT_VPSTRAN_2ND_INV=DOT_VPSTRAN_2ND_INV+
     1             R2*(((DEPSILON_V(M))/DTIME)*((DEPSILON_V(M))/DTIME)+
     2                 ((DEPSILON_VP(M))/DTIME*(DEPSILON_VP(M))/DTIME)+
     3                 DSTRAN(M)/DTIME*DSTRAN(M)/DTIME)
        DOT_VPSTRAN_2ND_INV_IN=DOT_VPSTRAN_2ND_INV_IN+
     1             R2*(((DEPSILON_V(M))/DTIME)*((DEPSILON_V(M))/DTIME)+
     2                 ((DEPSILON_VP(M))/DTIME*(DEPSILON_VP(M))/DTIME))
      END DO
      DOT_VPSTRAN_2ND_INV=P05*DOT_VPSTRAN_2ND_INV
      DOT_VPSTRAN_2ND_INV=SQRT(DOT_VPSTRAN_2ND_INV)
      DOT_VPSTRAN_2ND_INV_IN=P05*DOT_VPSTRAN_2ND_INV_IN
      DOT_VPSTRAN_2ND_INV_IN=SQRT(DOT_VPSTRAN_2ND_INV_IN)
!   UPDATE EFFECTIVE VISCOPLASTIC STRAIN
      VPSTRAN_EFF=R0
      DO M=1,NDI
       VPSTRAN_EFF=VPSTRAN_EFF+(R2/R3*(EPSILON_VP(M))*(EPSILON_VP(M)))
      END DO
      DO M=NDI+1,NTENS
       VPSTRAN_EFF=VPSTRAN_EFF+(R1/R3*(EPSILON_VP(M))*(EPSILON_VP(M)))
      END DO
      VPSTRAN_EFF=SQRT(VPSTRAN_EFF)
!   ##################################################################################################################################
!   UPDATE FINAL J2 AND ESTIMATE VISCOSITY
!   UPDATE DEVIATORIC STRESSES, STRESSES, INCREMENTS
!   ##################################################################################################################################
      DO I=1,NTENS
        S(I)=STRESS(I)-(P*SOID(I))
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
      VISCOSITY=R0
      IF(DOT_VPSTRAN_2ND_INV .GT. R0)THEN
         VISCOSITY=SQRTJ2T/DOT_VPSTRAN_2ND_INV
      END IF
!   ##################################################################################################################################
!    COMPUTE INEALSTIC STRAIN INVARIANTS
!   ##################################################################################################################################
      DO M=1,NTENS
        EPSILON_VP(M)=EPSILON_VP(M)+DEPSILON_V(M)+DEPSILON_VP(M)
      END DO
      VPSTRAN_1ST_INV=R0 ! 1ST INVARIANT
      VPSTRAN_2ND_INV=R0 ! 2ND INVARIANT
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
!   ##################################################################################################################################
!    COMPUTE FRACTION OF MELT AT TIME STEP AND ACCOUNT FOR LATENT HEAT EFFECTS
!   ##################################################################################################################################
        HL=R1
        IF(T_SOLIDUS.LE.TEMP.AND.TEMP.LE.T_LIQUIDUS)THEN
          FPT=(((TEMP-T_SOLIDUS)/(T_LIQUIDUS-T_SOLIDUS))**PWRFPT)
          FRAC=(((TEMP-T_SOLIDUS)/(T_LIQUIDUS-T_SOLIDUS)))
          HL=R1/(R1+(R3*LH/(R2*C_P*(T_LIQUIDUS-T_SOLIDUS))*SQRT(FRAC)))
        ELSEIF(TEMP.GT.T_LIQUIDUS)THEN
          FPT=R1
          FRAC=R1
          HL=R1/(R1+(R3*LH/(R2*C_P*(T_LIQUIDUS-T_SOLIDUS))*SQRT(FRAC)))
        ELSE
          FPT=R0
          FRAC=R0
          HL=R1/(R1+(R3*LH/(R2*C_P*(T_LIQUIDUS-T_SOLIDUS))*SQRT(FRAC)))
          HL=R1
       END IF
!   ##################################################################################################################################
!    COMPUTE HEAT GENERATION RATE (HGR) WRITTEN OUT AS RPL FOR VISCOPLASTIC DEFORMATION
!   ##################################################################################################################################
      HGR=R0
      DO M=1,NTENS
      HGR=HGR+(STRESS(M)*((DEPSILON_V(M)/DTIME)+(DEPSILON_VP(M)/DTIME)))
      END DO 
      IF(CMNAME.EQ.'LOWERCRUST'.OR.CMNAME.EQ.'LowerCrust'
     1 .OR.CMNAME.EQ.'CRUST'.OR.CMNAME.EQ.'Crust'.OR.CMNAME
     2 .EQ.'WEAKSHEARZONE'.OR.CMNAME.EQ.'WeakShearZone')THEN
        HGR=HGR+RADIOGENIC
      END IF
      RPL=(HGR)*HL ! ACCOUNTING FOR LATENT HEAT IN THE PRESENCE OF MELT
!   ##################################################################################################################################
!     UPDATE STATE VARIABLES
      MAXSTRESS=MAX(STRESS(1), STRESS(2), STRESS(3))
      MINSTRESS=MIN(STRESS(1), STRESS(2), STRESS(3))
      DIFF_STRESS=MAXSTRESS-MINSTRESS
      JI=R0
      PRESS_FOR_JI=R0
      PRESS_FOR_JI=(STRESS(1)+STRESS(2)+STRESS(3))/R3      
      END IF
!   ##################################################################################################################################
!     END OF (VISCO)PLASTICITY ROUTINES
!   ##################################################################################################################################
      DO I=1,NTENS
        STATEV(I)=EPSILON_VP(I)                                   !     SDV 1-6   ! VISCOUS + PLASTIC STRAINS (11,22,33,12,13,23)
      END DO
      STATEV(NTENS+1)=EPSILON_BAR                                 !     SDV 7     ! HARDENING HISTORY (IF COHESION HARDENING IS USED)
      STATEV(NTENS+2)=PHI_YDP       	                          !     SDV 8     ! DRUCKER-PRAGER YIELD CRITERION
      STATEV(NTENS+3)=SQRTJ2                                      !     SDV 9     ! INITIAL J2
      STATEV(NTENS+4)=SQRTJ2T		                              !     SDV 10    ! FINAL J2
      STATEV(NTENS+5)=PRESS*1E9       		                      !     SDV 11    ! PRESSURE
      STATEV(NTENS+6)=DGAMA_VP          	                      !     SDV 12    ! VISCOPLASTIC SMOOTH CONE MULTIPLIER
      STATEV(NTENS+7)=DOT_VPSTRAN_1ST_INV                         !     SDV 13    ! FIRST INVARIANT OF VISCOPLASTIC STRAIN RATE
      STATEV(NTENS+8)=LOG10(DOT_VPSTRAN_2ND_INV)                  !     SDV 14    ! SECOND INVARIANT OF VISCOPLASTIC STRAIN RATE
      STATEV(NTENS+9)=VPSTRAN_1ST_INV                             !     SDV 15    ! FIRST INVARIANT OF VISCOPLASTIC STRAIN
      STATEV(NTENS+10)=VPSTRAN_2ND_INV                            !     SDV 16    ! SECOND INVARIANT OF VISCOPLASTIC STRAIN
      STATEV(NTENS+11)=LOG10(VISCOSITY)                           !     SDV 17    ! EFFECTIVE VISCOSITY
      STATEV(NTENS+12)=HGR                                        !     SDV 18    ! HEAT GENERATION PER UNIT TIME
      STATEV(NTENS+13)=DTEMP/DTIME                                !     SDV 19    ! CHECK FOR THERMAL STEADY STATE
      STATEV(NTENS+15)=TEMP-TEMP_INIT                             !     SDV 21    ! THERMAL PERTURBATIONS
      STATEV(NTENS+16)=DEPSILON_V(1)/DTIME                        !     SDV 22    ! VISCOUS STRAIN RATE(11)
      STATEV(NTENS+17)=DEPSILON_V(2)/DTIME                        !     SDV 23    ! VISCOUS STRAIN RATE(22)
      STATEV(NTENS+18)=DEPSILON_V(3)/DTIME                        !     SDV 24    ! VISCOUS STRAIN RATE(33)
      STATEV(NTENS+19)=DEPSILON_V(4)/DTIME                        !     SDV 25    ! VISCOUS STRAIN RATE(12)
      STATEV(NTENS+20)=DEPSILON_V(5)/DTIME                        !     SDV 26    ! VISCOUS STRAIN RATE(13)
      STATEV(NTENS+21)=DEPSILON_V(6)/DTIME                        !     SDV 27    ! VISCOUS STRAIN RATE(23)
      STATEV(NTENS+22)=DEPSILON_VP(1)/DTIME                       !     SDV 28    ! VISCOPLASTIC STRAIN RATE(11)
      STATEV(NTENS+23)=DEPSILON_VP(2)/DTIME                       !     SDV 29    ! VISCOPLASTIC STRAIN RATE(22)
      STATEV(NTENS+24)=DEPSILON_VP(3)/DTIME                       !     SDV 30    ! VISCOPLASTIC STRAIN RATE(33)
      STATEV(NTENS+25)=DEPSILON_VP(4)/DTIME                       !     SDV 31    ! VISCOPLASTIC STRAIN RATE(12)
      STATEV(NTENS+26)=DEPSILON_VP(5)/DTIME                       !     SDV 32    ! VISCOPLASTIC STRAIN RATE(13)
      STATEV(NTENS+27)=DEPSILON_VP(6)/DTIME                       !     SDV 33    ! VISCOPLASTIC STRAIN RATE(23)
      STATEV(NTENS+28)=FPT                                        !     SDV 34    ! FRACTION OF MELT
      STATEV(NTENS+32)=DGAMA_VISC                                 !     SDV 38    ! DGAMA_VISC
      STATEV(NTENS+33)=DOT_VPSTRAN_1ST_INV_IN                     !     SDV 39    ! INELASTIC STRAIN RATE 1ST INVARIANT
      STATEV(NTENS+34)=LOG10(SQRTJ2T/DOT_VPSTRAN_2ND_INV_IN)      !     SDV 40    ! INELASTIC VISCOSITY
      STATEV(NTENS+35)=VPSTRAN_EFF                                !     SDV 41    ! EFFECTIVE VISCOPLASTIC STRAIN
      STATEV(NTENS+36)=CREEP_COUNT                                !     SDV 42    ! NEWTON-RAPHSON ITERATION COUNTS FOR CREEP
      STATEV(NTENS+37)=PLAST_COUNT                                !     SDV 43    ! NEWTON-RAPHSON ITERATION COUNTS FOR PLASTICITY
      STATEV(NTENS+38)=phi                                        !     SDV 44    ! FRICTION ANGLE
      STATEV(NTENS+39)=DOT_VPSTRAN_1ST_INV *DOT_VPSTRAN_2ND_INV   !     SDV 45    ! PRODUCT OF DEVIATORIC AND VOLUMETRIC STRAIN RATE TENSORS
      STATEV(NTENS+40)=VPSTRAN_1ST_INV *VPSTRAN_2ND_INV           !     SDV 46    ! PRODUCT OF DEVIATORIC AND VOLUMETRIC STRAIN TENSORS
      STATEV(NTENS+41)=ETS                                        !     SDV 47    ! ELASTIC AND THERMAL STRAINS
      RETURN
      END
!   ##################################################################################################################################
!   ##################################################################################################################################
!   ##################################################################################################################################
!     END OF USER MATERIAL SUBROUTINE
!   ##################################################################################################################################
!   ##################################################################################################################################
!   ##################################################################################################################################
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
      COMMON /CONDUCTIVITY_PARAMS/ TREF, ALPH_LE, RHOD_REF, TEMPERA, HL,
     1 PRESS, NCOMPS
      CHARACTER*80 CMNAME

       DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)

       PARAMETER (EXMULT=0.0000004D0,R0=0.0D0, R1=1.0D0)
       RPL=STATEV(NCOMPS+12)
       IF(KSTEP .GT. R1)THEN
        COND = PROPS(1)*HL ! ACCOUNTING FOR LATENT HEAT
       ELSE
        COND = PROPS(1)
       END IF
       C_P = PROPS(2)
       RHOD = PROPS(3)*(R1-(ALPH_LE*(TEMP-TREF)))
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
!   ##################################################################################################################################
!     END OF USER MATERIAL HEAT TRANSFER SUBROUTINE
!   ##################################################################################################################################
!   ##################################################################################################################################
!   ##################################################################################################################################
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,
     3 LACCFLA)
!
      INCLUDE 'ABA_PARAM.INC'
!
      COMMON /CONDUCTIVITY_PARAMS/ TREF, ALPH_LE, RHOD_REF, TEMPERA, HL,
     1 PRESS, NCOMPS
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(47)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(47),JARRAY(47),JMAC(*),JMATYP(*),
     1 COORD(*)
      PARAMETER (R0=0.0D0, R1=1.0D0)
      DOUBLE PRECISION FPT,CEQ,DENSITY_UPDATED,HFL_OLD,HFL_NEW,DHFL_DT
      CALL GETVRM('TEMP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
      TEMPERA=ARRAY(1)
      CEQ=TEMPERA-TREF
!      IF(KSTESP.GT.1)THEN
!      VOLPLAST=STATEV(NTENS+9)+STATEV(NTENS+41)
!      PRINT*,"ETS=",STATEV(NTENS+41)
!      ELSE
!      VOLPLAST=R0
!      ENDIF
      IF(CEQ.GT.R0)THEN
         DENSITY_UPDATED=RHOD_REF*(R1-(ALPH_LE*(TEMPERA-TREF)))
      ELSE
        DENSITY_UPDATED=RHOD_REF
      END IF
      CALL GETVRM('FV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
      FIELD(1)=DENSITY_UPDATED
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
      STATEV(35)=FIELD(1)
      IF(KINC==1)THEN
        CALL GETVRM('HFL',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
        HFL_OLD=ARRAY(1)
        HFL_NEW=ARRAY(1)
        DHFL_DT=R0
        CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
        STATEV(36)=HFL_OLD
        STATEV(37)=DHFL_DT
        CALL GETVRM('FV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
        FIELD(2)=HFL_NEW
      ELSE
        CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
        HFL_OLD=ARRAY(36)
        CALL GETVRM('HFL',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
        HFL_NEW=ARRAY(1)
        DHFL_DT=(HFL_NEW-HFL_OLD)/DTIME
        HFL_OLD=HFL_NEW
        STATEV(36)=HFL_OLD
        STATEV(37)=DHFL_DT
        CALL GETVRM('FV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
        FIELD(2)=HFL_NEW
        END IF
      RETURN
      END
!   ##################################################################################################################################
!   ##################################################################################################################################
!   ##################################################################################################################################
!     END OF USER DEFINED FIELD SUBROUTINE
!   ##################################################################################################################################
!   ##################################################################################################################################
!   ##################################################################################################################################
