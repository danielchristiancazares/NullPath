PROGRAM MINO_TIME_BATCH_TESTER
C*******************************************************************************
C     MINO TIME INTEGRATION BATCH TESTER
C     Tests GEOPHITIME subroutine for comprehensive parameter validation
C     
C     This is the CRITICAL missing component - Mino time parameterization
C     GEOPHITIME computes lambda (affine parameter) which is fundamental
C     to the semi-analytic method's accuracy and stability
C     
C     Input format: u0 uf mu0 muf a l l2 q2 tpm tpr su sm 
C     Output: All GEOPHITIME outputs including lambda (Mino time)
C*******************************************************************************
C     TODO: Per-case adaptive grid and ISO_C_BINDING sampler
C     - Expose per-λ sampling via geokerr_eval_lambda (shared lib libgeokerr_ref)
C     - Grid policy: lam0=0, S=4096, K=6, λmax = K*max(Tr,Tθ), Δλ = λmax/(S-1)
C     - Persist lam0_bits/dlam_bits for audit; unwrap φ; indices 0-based
C     - See reference/geokerr_eval_lambda.f90 for the C ABI stub (to be wired)
      IMPLICIT NONE
      
C     Parameters from GEOPHITIME subroutine signature
      DOUBLE PRECISION u0, uf, mu0, muf, a, l, l2, q2
      INTEGER tpm, tpr, su, sm
      
C     GEOPHITIME outputs (the missing piece!)
      DOUBLE PRECISION iu, h1, phimu, tmu, lambda
      DOUBLE PRECISION phiu, tu
      INTEGER ncase
      DOUBLE PRECISION u1, u2, u3, u4
      DOUBLE PRECISION rffu0, rffu1, rffmu1, rffmu2, rffmu3
      DOUBLE PRECISION rdc, rjc
      DOUBLE PRECISION tu01, tu02, tu03, tu04, tmu1, tmu3
      DOUBLE PRECISION phimu1, phimu3
      LOGICAL firstpt
      
C     Loop variables and file handling
      INTEGER i, ntests, ios
      CHARACTER*256 line
      
C     Constants
      DOUBLE PRECISION PI, ZERO, ONE, TWO
      PARAMETER (PI=3.141592653589793D0, ZERO=0.D0, ONE=1.D0, TWO=2.D0)
      
C     Error handling
      LOGICAL error_occurred
      
      WRITE(6,*) 'GEOKERR MINO TIME INTEGRATION BATCH TESTER'
      WRITE(6,*) '=========================================='
      WRITE(6,*) 'Testing GEOPHITIME subroutine - core Mino time integration'
      WRITE(6,*) ''
      
C     Read number of test cases
      READ(5,*,iostat=ios) ntests
      IF (ios .NE. 0) THEN
        WRITE(6,*) 'ERROR: Cannot read number of test cases'
        STOP 1
      ENDIF
      
      WRITE(6,*) 'Processing', ntests, ' test cases...'
      WRITE(6,*) ''
      
C     Output header
      WRITE(6,100) 
  100 FORMAT('#',1X,'case_id',5X,'u0',12X,'uf',12X,'mu0',11X,'muf',11X,
     &       'a',12X,'l',12X,'q2',11X,'lambda',7X,'phimu',8X,'tmu',9X,
     &       'phiu',9X,'tu',10X,'ncase',2X,'status')
      
C     Process each test case
      DO i = 1, ntests
        
        error_occurred = .FALSE.
        
C       Read input parameters
        READ(5,*,iostat=ios) u0, uf, mu0, muf, a, l, l2, q2, 
     &                       tpm, tpr, su, sm
        
        IF (ios .NE. 0) THEN
          WRITE(6,*) 'ERROR: Cannot read test case', i
          error_occurred = .TRUE.
        ELSE
        
C         Initialize outputs
          iu = ZERO
          h1 = ZERO  
          phimu = ZERO
          tmu = ZERO
          lambda = ZERO
          phiu = ZERO
          tu = ZERO
          ncase = 0
          u1 = ZERO; u2 = ZERO; u3 = ZERO; u4 = ZERO
          rffu0 = ZERO; rffu1 = ZERO  
          rffmu1 = ZERO; rffmu2 = ZERO; rffmu3 = ZERO
          rdc = ZERO; rjc = ZERO
          tu01 = ZERO; tu02 = ZERO; tu03 = ZERO; tu04 = ZERO
          tmu1 = ZERO; tmu3 = ZERO
          phimu1 = ZERO; phimu3 = ZERO
          firstpt = .TRUE.
          
C         ** CRITICAL CALL: This is what we failed to port **
C         GEOPHITIME is the core semi-analytic integration routine
C         lambda output is the Mino time parameter - fundamental to method
          CALL GEOPHITIME(u0, uf, mu0, muf, a, l, l2, q2, 
     &                    tpm, tpr, su, sm, iu, h1, 
     &                    phimu, tmu, ncase, u1, u2, u3, u4,
     &                    phiu, tu, lambda, rffu0, rffu1,
     &                    rffmu1, rffmu2, rffmu3, rdc, rjc,
     &                    tu01, tu02, tu03, tu04, tmu1, tmu3,
     &                    phimu1, phimu3, firstpt)
                    
C         Check for computational errors  
          IF (lambda .NE. lambda) THEN  ! NaN check
            error_occurred = .TRUE.
          ENDIF
          
        ENDIF
        
C       Output results (focus on lambda - the missing piece!)
        IF (error_occurred) THEN
          WRITE(6,200) i, u0, uf, mu0, muf, a, l, q2, 
     &                 ZERO, ZERO, ZERO, ZERO, ZERO, -1, 'ERROR'
        ELSE
          WRITE(6,200) i, u0, uf, mu0, muf, a, l, q2,
     &                 lambda, phimu, tmu, phiu, tu, ncase, 'OK'
        ENDIF
        
  200   FORMAT(I8,1X,8(ES12.5,1X),I5,1X,A5)
        
      ENDDO
      
      WRITE(6,*) ''
      WRITE(6,*) 'Batch processing complete!'
      WRITE(6,*) ''
      WRITE(6,*) 'KEY OUTPUT: lambda = Mino time parameter'
      WRITE(6,*) 'This is the fundamental quantity missing from CUDA port!'
      WRITE(6,*) 'All geodesic integration should use d/dlambda, not d/dt'
      
      END PROGRAM
      
C*******************************************************************************
C     The GEOPHITIME subroutine should already be included from geokerr_original.f
C     If not, it must be linked during compilation
C*******************************************************************************