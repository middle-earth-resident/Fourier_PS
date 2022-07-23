MODULE GEN_DATA
        REAL(8),PARAMETER :: L = 80
        INTEGER,PARAMETER :: N = 256,N2 = N/2, NN = N-1
        REAL(8),PARAMETER :: DT = 0.02 , TMAX = 20
        INTEGER,PARAMETER :: NMAX = TMAX/DT,IS = 10
        REAL(8),PARAMETER :: DX = L/N
       
        REAL(8),PARAMETER :: PI = 3.1415926358979D0,TWO_PI = 2* PI        
        COMPLEX(8),PARAMETER :: CI = (0.0D0,1.0D0)
       
        REAL(8),DIMENSION(0:N) :: X,K,K2
        COMPLEX(8),DIMENSION(0:N) :: U
END MODULE GEN_DATA


MODULE FFTW3
        USE,INTRINSIC :: ISO_C_BINDING
        IMPLICIT NONE
        INCLUDE 'fftw3.f03'
END MODULE FFTW3




PROGRAM FOURIER_PS
        USE GEN_DATA
        INTERFACE

                SUBROUTINE CALC_DU(LU,U,DU)
                        USE GEN_DATA,ONLY : N,NN
                        IMPLICIT NONE
                        COMPLEX(8),DIMENSION(0:N),INTENT(IN) :: LU,U 
                        COMPLEX(8),DIMENSION(0:N),INTENT(OUT) :: DU 
                END SUBROUTINE CALC_DU

                SUBROUTINE CALC_DER(U,LU)
                        USE GEN_DATA,ONLY : N
                        IMPLICIT NONE
                        COMPLEX(8),DIMENSION(0:N),INTENT(IN) :: U
                        COMPLEX(8),DIMENSION(0:N),INTENT(OUT) :: LU
                END SUBROUTINE CALC_DER

                SUBROUTINE INIT()
                        IMPLICIT NONE
                END SUBROUTINE INIT
                

                SUBROUTINE WRITE_DEN(FP,U2)
                        USE GEN_DATA,ONLY : X,N
                        IMPLICIT NONE
                        REAL(8),DIMENSION(0:N),INTENT(IN) :: U2
                        INTEGER,INTENT(IN) :: FP
                END SUBROUTINE WRITE_DEN
         END INTERFACE
        COMPLEX(8),DIMENSION(0:N) :: LU
        COMPLEX(8),DIMENSION(0:N) :: DU1,DU2,DU3,DU4,V

        INTEGER :: ITER,I,S,STEP
        CHARACTER(256) :: ITER_FILENAME

        ITER = NMAX/IS
        CALL INIT()
        STEP = 1
        
        DO S = 1,NMAX
                CALL CALC_DER(U,LU)
                CALL CALC_DU(LU,U,DU1)
                V = U + 0.5 * DU1 * DT
                
                CALL CALC_DER(V,LU)
                CALL CALC_DU(LU,V,DU2)
                V = U + 0.5 * DU2 * DT
        
                CALL CALC_DER(V,LU)
                CALL CALC_DU(LU,V,DU3)
                V = U +  DU3 * DT

                CALL CALC_DER(V,LU)
                CALL CALC_DU(LU,V,DU4)
                
                U = U + (DU1 + 2 * DU2 + 2 * DU3 + DU4) * DT/6
                
                
                IF(MOD(S,ITER)==0) THEN

                        OPEN(1,FILE="den.txt")
                                CALL WRITE_DEN(1,ABS(U))
                                WRITE(1,*)
                        STEP = STEP+1
                END IF
                

        END DO
        CLOSE(1)

END PROGRAM FOURIER_PS



SUBROUTINE CALC_DU(LU,U,DU)
        USE GEN_DATA,ONLY : N,NN,CI
        IMPLICIT NONE
        COMPLEX(8),DIMENSION(0:N),INTENT(IN) :: LU,U 
        COMPLEX(8),DIMENSION(0:N),INTENT(OUT) :: DU 

        INTEGER :: I

        DO I=0,N
                DU(I) = CI * (LU(I) + 2 * U(I) * U(I) * CONJG(U(I)))
        END DO
        RETURN 
END SUBROUTINE CALC_DU

SUBROUTINE CALC_DER(U,LU)
        USE GEN_DATA,ONLY : N,NN,K2,CI
        USE FFTW3,ONLY: C_PTR,FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE,FFTW_PLAN_DFT_1D,FFTW_DESTROY_PLAN,FFTW_EXECUTE_DFT
        IMPLICIT NONE
        COMPLEX(8),DIMENSION(0:N),INTENT(IN) :: U
        COMPLEX(8),DIMENSION(0:N),INTENT(OUT) :: LU

        COMPLEX(8),DIMENSION(0:N) :: UT,FU,K2U,IFU

        TYPE(C_PTR) :: PLANF,PLANB
        PLANF = FFTW_PLAN_DFT_1D(N,UT,FU,FFTW_FORWARD,FFTW_ESTIMATE)
        PLANB = FFTW_PLAN_DFT_1D(N,K2U,IFU,FFTW_BACKWARD,FFTW_ESTIMATE)
        
        UT = U 
        CALL FFTW_EXECUTE_DFT(PLANF,UT,FU)
        K2U = -1*K2*FU
        CALL FFTW_EXECUTE_DFT(PLANB,K2U,IFU)
        LU =  IFU/N
        
        CALL FFTW_DESTROY_PLAN(PLANF)
        CALL FFTW_DESTROY_PLAN(PLANB)
        RETURN
END SUBROUTINE CALC_DER


SUBROUTINE INIT()
        USE GEN_DATA,ONLY : L,X,K,K2,N,DX,N2,NN,TWO_PI,U,CI
        INTEGER :: I

        DO I = 0,N
                X(I) = -L/2.D0 + (I) * DX
        END DO

        DO I = 1,N2
                K(I-1) = (I-1) * TWO_PI/(N * DX)
        END DO
        DO I = 1,N2
                K(N2+I-1) = -(N2-I+1) * TWO_PI/(N * DX)
        END DO
        K2 = K * K 
        DO I = 0,N
                U(I) = 1.2 * 1/COSH(1.2 * (X(I) + 20)) * EXP(CI * X(I)) + 0.8 * 1/COSH(0.8 * X(I))
        END DO

END SUBROUTINE INIT

SUBROUTINE WRITE_DEN(FP,U2)
        USE GEN_DATA,ONLY : X,N
        IMPLICIT NONE
        REAL(8),DIMENSION(0:N),INTENT(IN) :: U2
        INTEGER,INTENT(IN) :: FP
        INTEGER :: I 
        DO I = 0,N
                WRITE(FP,20) U2(I)
        END DO

       20 FORMAT(E12.5) 
END SUBROUTINE WRITE_DEN
