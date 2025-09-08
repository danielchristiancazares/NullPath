      PROGRAM GEOPHITIME_BATCH_TESTER
C*******************************************************************************
C     GEOPHITIME MINO TIME INTEGRATION BATCH TESTER
C     Tests the extracted GEOPHITIME subroutine for comprehensive validation
C     
C     This validates the CRITICAL missing component - Mino time computation
C     GEOPHITIME computes lambda (affine parameter λ) from Equation (49)
C     λ = λ_u + t_μ  (radial + polar angle contributions)
C     
C     Input format from mino_time_inputs.dat:
C     u0 uf mu0 muf a l l2 q2 tpm tpr su sm iu h1 ncase u1 u2 u3 u4
C     
C     Output format:
C     lambda phimu tmu phiu tu (Mino time and coordinate components)
C*******************************************************************************
      IMPLICIT NONE
      
C     GEOPHITIME input parameters
      DOUBLE PRECISION u0, uf, mu0, muf, a, l, l2, q2
      INTEGER tpm, tpr, su, sm, ncase
      DOUBLE PRECISION iu, h1, u1, u2, u3, u4
      
C     GEOPHITIME outputs (the critical values!)
      DOUBLE PRECISION phimu, tmu, phiu, tu, lambda
      DOUBLE PRECISION rffu0, rffu1, rffmu1, rffmu2, rffmu3
      DOUBLE PRECISION rdc, rjc
      DOUBLE PRECISION tu01, tu02, tu03, tu04, tmu1, tmu3
      DOUBLE PRECISION phimu1, phimu3
      LOGICAL firstpt
      
C     Loop and file handling
      INTEGER i, ntests, ios
      CHARACTER*256 line
      
C     Constants
      DOUBLE PRECISION PI, ZERO, ONE, TWO
      PARAMETER (PI=3.141592653589793D0, ZERO=0.D0, ONE=1.D0, TWO=2.D0)
      
      WRITE(6,*) 'GEOKERR GEOPHITIME MINO TIME BATCH TESTER'
      WRITE(6,*) '========================================='
      WRITE(6,*) 'Validating Mino time computation λ = λ_u + t_μ'
      WRITE(6,*) ''
      
C     Read number of test cases
      READ(5,*,iostat=ios) ntests
      IF (ios .NE. 0) THEN
        WRITE(6,*) 'ERROR: Cannot read number of test cases'
        STOP 1
      ENDIF
      
      WRITE(6,*) 'Processing', ntests, 'test cases...'
      WRITE(6,*) ''
      WRITE(6,*) '# lambda phimu tmu phiu tu (Mino time outputs)'
      
C     Process each test case
      DO i = 1, ntests
        READ(5,*,iostat=ios) u0, uf, mu0, muf, a, l, l2, q2, 
     &                       tpm, tpr, su, sm, iu, h1, ncase,
     &                       u1, u2, u3, u4
        
        IF (ios .NE. 0) THEN
          WRITE(6,*) 'ERROR: Cannot read test case', i
          STOP 1
        ENDIF
        
C       Initialize outputs and call GEOPHITIME
        phimu = ZERO
        tmu = ZERO
        phiu = ZERO
        tu = ZERO
        lambda = ZERO
        firstpt = .TRUE.
        
C       Call the critical GEOPHITIME subroutine
        CALL GEOPHITIME(u0, uf, mu0, muf, a, l, l2, q2, 
     &                  tpm, tpr, su, sm, iu, h1, 
     &                  phimu, tmu, ncase, u1, u2, u3, u4,
     &                  phiu, tu, lambda, rffu0, rffu1, 
     &                  rffmu1, rffmu2, rffmu3, rdc, rjc,
     &                  tu01, tu02, tu03, tu04, tmu1, tmu3,
     &                  phimu1, phimu3, firstpt)
        
C       Output the critical Mino time results
        WRITE(6,100) lambda, phimu, tmu, phiu, tu
100     FORMAT(5(E23.15,1X))
        
      ENDDO
      
      WRITE(6,*) ''
      WRITE(6,*) 'Batch testing completed successfully!'
      WRITE(6,*) 'Key output: lambda = Mino time parameter'
      
      END
      
C*******************************************************************************
C     Include the complete GEOPHITIME subroutine and dependencies
C     This would normally be linked from the extracted file
C*******************************************************************************