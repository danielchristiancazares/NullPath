      PROGRAM MINO_TIME_VALIDATION
C*******************************************************************************
C     MINO TIME PARAMETER VALIDATION USING GEOKERR SUBROUTINES
C     
C     This program validates the missing Mino time computation by directly
C     calling the GEOKERR geodesic solver and extracting λ (lambda) values
C     
C     Strategy: Use existing GEOKERR infrastructure to compute complete
C     geodesics and extract the critical λ parameter for validation
C     
C     Input format: u0 uf mu0 muf a l q2 alpha beta su sm
C     Output format: lambda (Mino time parameter)
C*******************************************************************************
      IMPLICIT NONE
      
C     GEOKERR parameters
      DOUBLE PRECISION u0, uf, uout, mu0, muf, a, l, q2
      DOUBLE PRECISION alpha, beta
      INTEGER tpm, tpr, su, sm, nup, ncase, kext, npts
      DOUBLE PRECISION offset
      LOGICAL phit, usegeor, mufill
      
C     GEOKERR outputs - we need LAMBDAI array
      DOUBLE PRECISION, ALLOCATABLE :: ufi(:), mufi(:), dti(:)
      DOUBLE PRECISION, ALLOCATABLE :: dphi(:), lambdai(:)
      INTEGER, ALLOCATABLE :: tpmi(:), tpri(:)
      
C     Loop and file handling
      INTEGER i, ntests, ios
      DOUBLE PRECISION l2
      
      WRITE(6,*) 'GEOKERR MINO TIME PARAMETER VALIDATION'
      WRITE(6,*) '====================================='
      WRITE(6,*) 'Extracting λ (lambda) from complete geodesic solver'
      WRITE(6,*) ''
      
C     Read number of test cases
      read(5,*,iostat=ios) ntests
      IF (ios .NE. 0) THEN
        WRITE(6,*) 'ERROR: Cannot read number of test cases'
        STOP 1
      ENDIF
      
C     Allocate arrays for single point computation
      nup = 1
      npts = 0
      kext = 0
      ALLOCATE(ufi(nup), mufi(nup), dti(nup))
      ALLOCATE(dphi(nup), lambdai(nup))
      ALLOCATE(tpmi(nup), tpri(nup))
      
      WRITE(6,*) 'Processing', ntests, 'Mino time test cases...'
      WRITE(6,*) ''
      WRITE(6,*) '# lambda (Mino time parameter from GEOKERR)'
      
C     Process each test case
      DO i = 1, ntests
        read(5,*,iostat=ios) u0, uf, mu0, muf, a, l, q2, 
     &                       alpha, beta, su, sm
        
        IF (ios .NE. 0) THEN
          WRITE(6,*) 'ERROR: Cannot read test case', i
          STOP 1
        ENDIF
        
C       Set up GEOKERR parameters
        l2 = l * l
        uout = u0
        offset = 0.0D0
        phit = .TRUE.      ! We need time/phi computation for lambda
        usegeor = .FALSE.  ! Solve for muf given uf
        mufill = .FALSE.   ! No extra resolution needed
        tpm = 0           ! Will be computed
        tpr = 0           ! Will be computed
        ncase = 0         ! Will be computed
        
C       Call GEOKERR to compute complete geodesic with lambda
        CALL GEOKERR(u0, uf, uout, mu0, muf, a, l, q2, 
     &               alpha, beta, tpm, tpr, su, sm, nup, offset, 
     &               phit, usegeor, mufill, ncase, kext, npts,
     &               ufi, mufi, dti, dphi, tpmi, tpri, lambdai)
        
C       Output the critical Mino time parameter
        WRITE(6,100) lambdai(1)
100     FORMAT(E23.15)
        
      ENDDO
      
      DEALLOCATE(ufi, mufi, dti, dphi, lambdai, tpmi, tpri)
      
      WRITE(6,*) ''
      WRITE(6,*) 'Mino time validation completed!'
      WRITE(6,*) 'Lambda values extracted from full GEOKERR solver'
      
      END