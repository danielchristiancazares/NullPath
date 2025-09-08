      PROGRAM GEOKERR_BATCH_ELLIPTIC
C     Batch testing harness for Carlson elliptic integral functions
C     Reads test cases and outputs reference values for CUDA validation
      IMPLICIT NONE
      
      INTEGER MAXTEST
      PARAMETER (MAXTEST=10000)
      
      INTEGER NTEST, I
      DOUBLE PRECISION X, Y, Z, P, RESULT
      CHARACTER*2 TEST_TYPE
      
C     Declare Carlson functions
      DOUBLE PRECISION RF, RC, RD, RJ
      EXTERNAL RF, RC, RD, RJ
      
      WRITE(6,*) 'GEOKERR BATCH ELLIPTIC INTEGRAL TESTER'
      WRITE(6,*) 'Reading from elliptic_inputs.dat'
      
C     Open input file
      OPEN(UNIT=10, FILE='elliptic_inputs.dat', STATUS='OLD')
      READ(10,*) NTEST
      
      WRITE(6,*) 'Processing', NTEST, 'test cases'
      
C     Open output file  
      OPEN(UNIT=20, FILE='elliptic_outputs.dat', STATUS='UNKNOWN')
      WRITE(20,*) NTEST
      
      DO I = 1, NTEST
        READ(10,*) TEST_TYPE, X, Y, Z, P
        
        IF (TEST_TYPE .EQ. 'RF') THEN
          RESULT = RF(X, Y, Z)
        ELSE IF (TEST_TYPE .EQ. 'RC') THEN
          RESULT = RC(X, Y)
        ELSE IF (TEST_TYPE .EQ. 'RD') THEN
          RESULT = RD(X, Y, Z)
        ELSE IF (TEST_TYPE .EQ. 'RJ') THEN
          RESULT = RJ(X, Y, Z, P)
        ELSE
          WRITE(6,*) 'Unknown test type:', TEST_TYPE
          RESULT = 0.0D0
        ENDIF
        
        WRITE(20,'(A2,1X,4(E20.12,1X),E20.12)') TEST_TYPE, X, Y, Z, P, 
     &        RESULT
        
        IF (MOD(I, 1000) .EQ. 0) THEN
          WRITE(6,*) 'Processed', I, 'cases'
        ENDIF
      END DO
      
      CLOSE(10)
      CLOSE(20)
      
      WRITE(6,*) 'Results written to elliptic_outputs.dat'
      
      END
      
C     Include the Carlson elliptic integral functions from original geokerr
      double precision FUNCTION rf(x,y,z)
C     COMPUTES CARLSON'S ELLIPTIC INTEGRAL OF THE FIRST KIND
C     RF(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C                   -1/2     -1/2     -1/2
C         (1/2)(T+X)   (T+Y)   (T+Z)   DT,
C     WHERE X,Y,Z ARE NONNEGATIVE AND AT MOST ONE OF THEM IS ZERO.
      implicit double precision (a-h,o-z)
      parameter(errtol=.08,tiny=1.5e-38,big=3.e37,third=1./3.,
     *          c1=1./24.,c2=.1,c3=3./44.,c4=1./14.)
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.tiny.or.
     *   max(x,y,z).gt.big)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      zt=.25*(zt+alamb)
      ave=third*(xt+yt+zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.errtol)goto 1
      e2=delx*dely-delz*delz
      e3=delx*dely*delz
      rf=(1.+(c1*e2-c2-c3*e3)*e2+c4*e3)/sqrt(ave)
      return
      end

      double precision FUNCTION rc(x,y)
C     COMPUTES CARLSON'S DEGENERATE ELLIPTIC INTEGRAL
C     RC(X,Y) = INTEGRAL FROM ZERO TO INFINITY OF
C                         -1/2      -1
C               (1/2)(T+X)    (T+Y)  DT,
C     WHERE X IS NONNEGATIVE AND Y IS POSITIVE.
      implicit double precision (a-h,o-z)
      parameter(errtol=.04,tiny=1.69e-38,sqrtny=1.3e-19,big=3.e37,
     *          tnbg=tiny*big,comp1=2.236/sqrtny,comp2=tnbg*tnbg/25.,
     *          third=1./3.,c1=.3,c2=1./7.,c3=.375,c4=9./22.)
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.tiny.or.(x+
     *   abs(y)).gt.big.or.(y.lt.-comp1.and.x.gt.0..and.x.lt.comp2))
     *   pause 'invalid arguments in rc'
      if(y.gt.0.)then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
      alamb=2.*sqrt(xt)*sqrt(yt)+yt
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      ave=third*(xt+yt+yt)
      s=(yt-ave)/ave
      if(abs(s).gt.errtol)goto 1
      rc=w*(1.+s*s*(c1+s*(c2+s*(c3+s*c4))))/sqrt(ave)
      return
      end

      double precision FUNCTION rd(x,y,z)
C     COMPUTES CARLSON'S ELLIPTIC INTEGRAL OF THE SECOND KIND
C     RD(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C                           -1/2     -1/2     -3/2
C                 (3/2)(T+X)   (T+Y)   (T+Z)   DT,
C     WHERE X,Y ARE NONNEGATIVE, AT MOST ONE OF THEM IS ZERO, AND Z IS
C     POSITIVE.
      implicit double precision (a-h,o-z)
      parameter(errtol=.05,tiny=1.e-25,big=4.5e21,c1=3./14.,c2=1./6.,
     *          c3=9./22.,c4=3./26.,c5=.25*c3,c6=1.5*c4)
      if(min(x,y).lt.0..or.min(x+y,z).lt.tiny.or.max(x,y,z).gt.big)
     *   pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
1     sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      sum=sum+fac/(sqrtz*(zt+alamb))
      fac=.25*fac
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      zt=.25*(zt+alamb)
      ave=.2*(xt+yt+3.*zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.errtol)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd=3.*sum+fac*(1.+ed*(-c1+c5*ed-c6*delz*ee)
     *    +delz*(c2*ee+delz*(-c3*ec+delz*c4*ea)))/(ave*sqrt(ave))
      return
      end

      double precision FUNCTION rj(x,y,z,p)
C     COMPUTES CARLSON'S ELLIPTIC INTEGRAL OF THE THIRD KIND
C     RJ(X,Y,Z,P) = INTEGRAL FROM ZERO TO INFINITY OF
C                             -1/2     -1/2     -1/2     -1
C                   (3/2)(T+X)   (T+Y)   (T+Z)   (T+P)  DT,
C     WHERE X,Y,Z ARE NONNEGATIVE, AT MOST ONE OF THEM IS ZERO, AND P
C     IS POSITIVE.
      implicit double precision (a-h,o-z)
      parameter(errtol=.05,tiny=2.5e-13,big=9.e11,c1=3./14.,c2=1./3.,
     *          c3=3./22.,c4=3./26.,c5=.75*c3,c6=1.5*c4,c7=.5*c2,
     *          c8=c3+c3)
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,p).lt.tiny.or.
     *   max(x,y,z,p).gt.big)pause 'invalid arguments in rj'
      sum=0.
      fac=1.
      if(p.gt.0.)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1./(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      sqrtp=sqrt(pt)
      dnm=sqrtp*(sqrtp+sqrtx)*(sqrtp+sqrty)*(sqrtp+sqrtz)
      sum=sum+fac/dnm
      fac=.25*fac
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
      beta=pt*(pt+alamb)**2
      xt=.25*(xt+alamb)
      yt=.25*(yt+alamb)
      zt=.25*(zt+alamb)
      pt=.25*(pt+alamb)
      ave=.2*(xt+yt+zt+pt+pt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave
      delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.errtol)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp*delp
      ed=ea-3.*ec
      ee=eb+2.*delp*(ea-ec)
      rj=3.*sum+fac*(1.+ed*(-c1+c5*ed-c6*ee)+eb*(c7+delp*(-c8+delp*c4))
     *    +delp*ea*(c2-delp*c3)-c2*delp*ec)/(ave*sqrt(ave))
      if(p.le.0.)rj=a*(b*rj+3.*(rcx-rf(xt,yt,zt)))
      return
      end