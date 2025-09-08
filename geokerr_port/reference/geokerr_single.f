      SUBROUTINE GEOKERR_SINGLE_CASE(u0, uf, mu0, muf, a, l, q2,
     &                            su, sm, lambda_result, muf_result,
     &                            ncase_result, valid_result)
C*******************************************************************************
C     SINGLE CASE GEOKERR INTERFACE
C     
C     Clean interface for single geodesic computation using original GEOKERR
C     algorithm. Designed to be called from C/CUDA comparison tools.
C     
C     Inputs:
C       u0, uf     : Initial and final inverse radius
C       mu0, muf   : Initial and final cos(theta) 
C       a          : Kerr spin parameter
C       l          : Angular momentum (L/E)
C       q2         : Carter constant (Q/E^2)
C       su, sm     : Initial velocities
C     
C     Outputs:
C       lambda_result : Mino time parameter λ
C       muf_result    : Computed final cos(theta)
C       ncase_result  : Orbit classification (1-8)
C       valid_result  : Success flag (1=success, 0=failure)
C*******************************************************************************
      IMPLICIT NONE
      
C     Input parameters
      DOUBLE PRECISION u0, uf, mu0, muf, a, l, q2
      INTEGER su, sm
      
C     Output parameters  
      DOUBLE PRECISION lambda_result, muf_result
      INTEGER ncase_result, valid_result
      
C     Internal GEOKERR parameters
      DOUBLE PRECISION uout, alpha, beta, offset, l2
      INTEGER tpm, tpr, nup, ncase, kext, npts
      LOGICAL phit, usegeor, mufill
      
C     GEOKERR output arrays (single element)
      DOUBLE PRECISION ufi(1), mufi(1), dti(1), dphi(1), lambdai(1)
      INTEGER tpmi(1), tpri(1)
      
C     Detailed computation parameters
      DOUBLE PRECISION iu, h1, u1, u2, u3, u4
      DOUBLE PRECISION rffu0, rffu1, rffmu1, rffmu2, rffmu3
      DOUBLE PRECISION iu0, i1mu, i3mu
      DOUBLE PRECISION phimu, tmu, phiu, tu, lambda
      DOUBLE PRECISION rdc, rjc, tu01, tu02, tu03, tu04
      DOUBLE PRECISION tmu1, tmu3, phimu1, phimu3
      LOGICAL firstpt, pht
      
C     Initialize outputs
      lambda_result = 0.0D0
      muf_result = 0.0D0  
      ncase_result = 0
      valid_result = 0
      
C     Set up GEOKERR parameters
      l2 = l * l
      uout = u0
      offset = 0.0D0
      nup = 1
      npts = 0
      kext = 0
      phit = .TRUE.      ! Need phi/time computation for λ
      usegeor = .FALSE.  ! Solve for muf given uf
      mufill = .FALSE.   
      tpm = 0           ! Will be computed
      tpr = 0           ! Will be computed
      ncase = 0         ! Will be computed
      alpha = 0.0D0     ! Not needed for direct computation
      beta = 0.0D0      ! Not needed for direct computation
      
      ! Initialize all arrays to prevent undefined behavior
      ufi(1) = 0.0D0
      mufi(1) = 0.0D0
      dti(1) = 0.0D0
      dphi(1) = 0.0D0
      lambdai(1) = 0.0D0
      tpmi(1) = 0
      tpri(1) = 0
      
C     Call main GEOKERR subroutine
      CALL GEOKERR(u0, uf, uout, mu0, muf, a, l, q2, 
     &             alpha, beta, tpm, tpr, su, sm, nup, offset, 
     &             phit, usegeor, mufill, ncase, kext, npts,
     &             ufi, mufi, dti, dphi, tpmi, tpri, lambdai)
      
C     Check if computation succeeded
      IF (ncase .GT. 0 .AND. ncase .LE. 8) THEN
        valid_result = 1
        ncase_result = ncase
        muf_result = mufi(1)
        
C       Get detailed Mino time computation
        firstpt = .TRUE.
        pht = .TRUE.
        
C       Call GEOMU first to set up parameters
        CALL GEOMU(u0, uf, mu0, muf, a, l, l2, q2, iu, tpm, tpr, 
     &             su, sm, ncase, h1, u1, u2, u3, u4, rffu0, rffu1,
     &             rffmu1, rffmu2, rffmu3, iu0, i1mu, i3mu, pht, firstpt)
        
C       Call GEOPHITIME to get proper λ computation
        CALL GEOPHITIME(u0, uf, mu0, muf, a, l, l2, q2, 
     &                  tpm, tpr, su, sm, iu, h1, 
     &                  phimu, tmu, ncase, u1, u2, u3, u4,
     &                  phiu, tu, lambda, rffu0, rffu1, 
     &                  rffmu1, rffmu2, rffmu3, rdc, rjc,
     &                  tu01, tu02, tu03, tu04, tmu1, tmu3,
     &                  phimu1, phimu3, firstpt)
        
        lambda_result = lambda
        
      ELSE
        ! Computation failed
        valid_result = 0
        ncase_result = 0
        lambda_result = 0.0D0
        muf_result = mu0  ! Return initial value
      ENDIF
      
      RETURN
      END
      
C*******************************************************************************
C     WRAPPER PROGRAM FOR SINGLE CASE TESTING
C*******************************************************************************
      PROGRAM TEST_GEOKERR_SINGLE
      IMPLICIT NONE
      
      DOUBLE PRECISION u0, uf, mu0, muf, a, l, q2
      DOUBLE PRECISION lambda_result, muf_result
      INTEGER su, sm, ncase_result, valid_result
      
      ! Test case parameters
      u0 = 0.1D0
      uf = 0.01D0
      mu0 = 0.8D0
      muf = 0.9D0
      a = 0.5D0
      l = 2.0D0
      q2 = 4.0D0
      su = 1
      sm = 1
      
      WRITE(6,*) 'Testing single case GEOKERR interface...'
      WRITE(6,*) 'Input: u0=', u0, ' uf=', uf, ' mu0=', mu0
      WRITE(6,*) '       muf=', muf, ' a=', a, ' l=', l, ' q2=', q2
      
      CALL GEOKERR_SINGLE_CASE(u0, uf, mu0, muf, a, l, q2,
     &                         su, sm, lambda_result, muf_result,
     &                         ncase_result, valid_result)
      
      WRITE(6,*) ''
      WRITE(6,*) 'Results:'
      WRITE(6,*) '  Valid:', valid_result
      WRITE(6,*) '  NCASE:', ncase_result
      WRITE(6,*) '  Lambda:', lambda_result
      WRITE(6,*) '  Final mu:', muf_result
      
      END