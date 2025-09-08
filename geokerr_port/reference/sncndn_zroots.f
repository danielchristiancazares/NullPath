      SUBROUTINE SNCNDN(UU,EMMC,SN,CN,DN)
***********************************************************************
*     PURPOSE:  Compute Jacobi-elliptic functions SN,CN,DN.
*     ARGUMENTS:  Given the arguments U,EMMC=1-k^2 calculate sn(u,k), cn(u,k), dn(u,k).
*     REMARKS:  
*     AUTHOR:  Press et al (1992).
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
      IMPLICIT NONE
      INTEGER I,II,L
      DOUBLE PRECISION A,B,C,CA,CN,D,DN,EMC,EMMC,SN,SQRT,U,UU
      PARAMETER (CA=3.d-8)
      LOGICAL BO
      double precision EM(13),EN(13)
      EMC=EMMC
      U=UU
      IF(EMC.NE.0.d0)THEN
        BO=(EMC.LT.0.d0)
        IF(BO)THEN
          D=1.d0-EMC
          EMC=-EMC/D
          D=sqrt(D)
          U=D*U
        ENDIF
        A=1.d0
        DN=1.d0
        DO 11 I=1,13
          L=I
          EM(I)=A
          EMC=sqrt(EMC)
          EN(I)=EMC
          C=0.5d0*(A+EMC)
          IF(ABS(A-EMC).LE.CA*A)GO TO 1
          EMC=A*EMC
          A=C
11      CONTINUE
1       U=C*U
        SN=DSIN(U)
        CN=DCOS(U)
        IF(SN.EQ.0.)GO TO 2
        A=CN/SN
        C=A*C
        DO 12 II=L,1,-1
          B=EM(II)
          A=C*A
          C=DN*C
          DN=(EN(II)+A)/(B+A)
          A=C/B
12      CONTINUE
        A=1.d0/sqrt(C*C+1.d0)
        IF(SN.LT.0.)THEN
          SN=-A
        ELSE
          SN=A
        ENDIF
        CN=C*SN
2       IF(BO)THEN
          A=DN
          DN=CN
          CN=A
          SN=SN/D
        ENDIF
      ELSE
        CN=1.d0/DCOSH(U)
        DN=CN
        SN=DTANH(U)
      ENDIF
      RETURN
      END
      SUBROUTINE ZROOTS(A,M,ROOTS,POLISH)
***********************************************************************
*     PURPOSE:  Find all roots of a polynomial.
*     ARGUMENTS:  Given the degree M and the M+1 complex coefficients
*       A of the polynomial (with A(0) being the constant term), this
*       routine returns all M roots in the complex array ROOTS.  The
*       logical variable POLISH should be input as .TRUE. if polishing
*       (by Laguerre's method) is desired, .FALSE. if the roots will be
*       subsequently polished by other means.
*     ROUTINES CALLED:  LAGUER.
*     ALGORITHM: Laguerre's method.
*     ACCURACY:  The parameter EPS sets the desired accuracy.
*     REMARKS:  
*     AUTHOR:  Press et al (1992).
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
      INTEGER M
      LOGICAL POLISH
*       
      INTEGER I,J,JJ,ITS,MAXM,ncc
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1d-6,MAXM=7)
      COMPLEX*16 AD(MAXM),X,B,C,XSUM,XSTART
      COMPLEX*16 A(M+1),ROOTS(M)
*       
      IF(M.GT.MAXM-1) THEN
        WRITE(6,*) 'M too large in ZROOTS'
      ENDIF
*       Copy of coefficients for successive deflation.
      DO 10 J=1,M+1
        AD(J)=A(J)
   10 CONTINUE
*       Loop over each root to be found.
      XSUM=0.d0
       XSTART=DCMPLX(0.D0,0.D0)
* If ncc=1, the previous root is a complex conjugate of another root
       ncc=0
      DO 20 J=M,1,-1
* Start at zero to favour convergence to smallest remaining
* root, or if the previous root was complex, start at its complex
* conjugate (since the coefficients are real):
        X=XSTART
*         if(J.lt.M.and.dimag(ROOTS(J+1)).ne.0.d0.and.ncc.eq.0) then
         if(J.lt.M) then
            if(dimag(ROOTS(J+1)).ne.0.d0.and.ncc.eq.0) then
               XSTART=DCMPLX(dble(roots(J+1)),-dble(roots(J+1)))
* Since we have chosen the second root to start at the complex conjugate,
* we don't want to use its complex conjugate again as a starting root:
               ncc=1
            else
               XSTART=DCMPLX(0.D0,0.D0)
               ncc=0
            endif
         else
           XSTART=DCMPLX(0.D0,0.D0)
           ncc=0
         endif
*        if(J.NE.1.or.a(M+1).eq.0.d0) then
* Find the root.
          CALL LAGUER(AD,J,X,ITS)
*           XSUM=XSUM+X
*        else
*          X=-a(M)/a(M+1)-XSUM
*        endif
        IF(ABS(DIMAG(X)).LE.2.d0*EPS*EPS*ABS(DBLE(X))) 
     &       X=DCMPLX(DBLE(X),0.d0)
        ROOTS(J)=X
        B=AD(J+1)
*         Forward deflation.
        DO 15 JJ=J,1,-1
          C=AD(JJ)
          AD(JJ)=B
          B=X*B+C
   15   CONTINUE
   20 CONTINUE
      IF(POLISH) THEN
*         Polish the roots using the undeflated coefficients.
        DO 30 J=1,M
          CALL LAGUER(A,M,ROOTS(J),ITS)
   30   CONTINUE
      ENDIF
      DO 40 J=2,M 
*         Sort roots by their real parts by straight insertion.
        X=ROOTS(J)
        DO 35 I=J-1,1,-1
          IF(DBLE(ROOTS(I)).LE.DBLE(X)) GO TO 37
          ROOTS(I+1)=ROOTS(I)
   35   CONTINUE
        I=0
   37   ROOTS(I+1)=X
   40 CONTINUE
      RETURN
      END

************************************************************************
       double precision FUNCTION rf(x,y,z)
************************************************************************
*     PURPOSE: Compute Carlson fundamental integral RF
*              R_F=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2)
*     ARGUMENTS: Symmetric arguments x,y,z
*     ROUTINES CALLED:  None.
*     ALGORITHM: Due to B.C. Carlson.
*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
*     REMARKS:  
*     AUTHOR:  Press et al (1992).
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
       implicit none
       double precision x,y,z,ERRTOL,THIRD,a1,C2,C3,C4
       PARAMETER(ERRTOL=0.0025d0,THIRD=1.d0/3.d0,
     *              a1=1.d0/24.d0,C2=0.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
       double precision alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,
     *              sqrtz,sqrt,xt,yt,zt
       xt=x
       yt=y
       zt=z
1       continue
              sqrtx=sqrt(xt)
              sqrty=sqrt(yt)
              sqrtz=sqrt(zt)
              alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
              xt=0.25d0*(xt+alamb)
              yt=0.25d0*(yt+alamb)
              zt=0.25d0*(zt+alamb)
              ave=THIRD*(xt+yt+zt)
              if(ave.eq.0.d0) then
                delx=0.d0
                dely=0.d0
                delz=0.d0
              else
                delx=(ave-xt)/ave
                dely=(ave-yt)/ave
                delz=(ave-zt)/ave
              endif
       if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)go to 1
       e2=delx*dely-delz*delz
       e3=delx*dely*delz
       rf=(1.d0+(a1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
       return
       END
       
*********************************************************************************
       SUBROUTINE GAULEG(X1,X2,X,W,N,NS)
********************************************************************************* 
*  PURPOSE: Subroutine to calculate the abscissas and weights for the Gauss-Legendre-
*  Quadrature. The routine is based on the NUMERICAL RECIPES and uses an
*  algorithem of G.B. Rybicki.
*  Input: x1 ,x2: range of integration.
*          n: order of the orthogonal polynomials and the quadrature formula.
*  Output: x = x(n): array of the abscissas.
*          w = w(n): array of the weights. 
      INTEGER N,M,I,J
      REAL*8 X1,X2,X(NS),W(NS)
      REAL*8 PI,XM,XL,Z,P1,P2,P3,Z1,PP,EPS
      PARAMETER (PI = 3.14159265358979323846D0)
      PARAMETER (EPS=3.D-14)
 
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
         Z=COS(PI*(I-.25D0)/(N+.5D0))
 1       CONTINUE
         P1=1.D0
         P2=0.D0
         DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
 11      CONTINUE
         PP=N*(Z*P1-P2)/(Z*Z-1.D0)
         Z1=Z
         Z=Z1-P1/PP
         IF(ABS(Z-Z1).GT.EPS) GOTO 1
         X(I)=XM-XL*Z
         X(N+1-I)=XM+XL*Z
         W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
         W(N+1-I)=W(I)
 12   CONTINUE
      RETURN
      END
       
***********************************************************************
      SUBROUTINE LAGUER(A,M,X,ITS)
***********************************************************************
*     PURPOSE:  Find one root of a polynomial.
*     ARGUMENTS:
*     ROUTINES CALLED:
*     ALGORITHM:  
*     ACCURACY:
*     REMARKS:  I don't have the documentation for this routine!
*     AUTHOR:  Press et al (1992)
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
      IMPLICIT NONE
      INTEGER ITS,M
      COMPLEX*16 A(M+1),X
*       
      INTEGER ITER,J,MAXIT,MR,MT
      PARAMETER (MR=8,MT=10,MAXIT=MT*MR)
       DOUBLE PRECISION ABX,ABP,ABM,ERR,EPSS,FRAC(MR) 
       PARAMETER (EPSS=1.d-15)
      COMPLEX*16 DX,X1,B,D,F,G,H,SQ,GP,GM,G2
      DATA FRAC /0.5d0,0.25d0,0.75d0,0.13d0,0.38d0,0.62d0,0.88d0,1.d0/
*       Loop over iterations up to allowed maximum.
      DO 20 ITER=1,MAXIT
*         
        ITS=ITER
        B=A(M+1)
        ERR=ABS(B)
        D=DCMPLX(0.d0,0.d0)
        F=DCMPLX(0.d0,0.d0)
        ABX=ABS(X)
        DO 10 J=M,1,-1
*           Efficient computation of the polynomial and its first TWO
*           derivatives.
          F=X*F+D
          D=X*D+B
          B=X*B+A(J)
          ERR=ABS(B)+ABX*ERR
   10   CONTINUE
        ERR=EPSS*ERR
*         
        IF(ABS(B).LE.ERR) THEN
*           Special case: we are on the root.
          RETURN
        ELSE
*           The generic case; use Laguerre's formula.
          G=D/B
          G2=G*G
          H=G2-2.d0*F/B
          SQ=SQRT((M-1)*(M*H-G2))
          GP=G+SQ
          GM=G-SQ
          ABP=ABS(GP)
          ABM=ABS(GM)
          IF(ABP.LT.ABM) GP=GM
          IF (MAX(ABP,ABM).GT.0.d0) THEN
            DX=M/GP
          ELSE
            DX=CEXP(CMPLX(LOG(1.d0+ABX),DBLE(ITER)))
          ENDIF
        ENDIF
        X1=X-DX
*         Check if we've converged.
        IF(X.EQ.X1)RETURN
        IF (MOD(ITER,MT).NE.0) THEN
          X=X1
        ELSE
*           
          X=X-DX*FRAC(ITER/MT)
        ENDIF
   20 CONTINUE
      WRITE(6,*) 'Too many iterations'
      RETURN
      END
************************************************************************
      double precision FUNCTION rc(x,y)
************************************************************************
*     PURPOSE: Compute Carlson degenerate integral RC
*              R_C(x,y)=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1)
*     ARGUMENTS: x,y
*     ROUTINES CALLED:  None.
*     ALGORITHM: Due to B.C. Carlson.
*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
*     REMARKS:  
*     AUTHOR:  Press et al (1992)
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
      double precision x,y,ERRTOL,TINY,BIG,THIRD,a1,C2,
     *C3,C4
      PARAMETER (ERRTOL=0.0012d0,TINY=1.69d-38,BIG=3.d37,
     *THIRD=1.d0/3.d0,
     *a1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
      double precision alamb,ave,s,w,xt,yt
         if(y.gt.0.)then
           xt=x
           yt=y
           w=1.d0
         else
           xt=x-y
           yt=-y
           w=sqrt(x)/sqrt(xt)
         endif
1        continue
           alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
           xt=.25d0*(xt+alamb)
           yt=.25d0*(yt+alamb)
           ave=THIRD*(xt+yt+yt)
           s=(yt-ave)/ave
         if(abs(s).gt.ERRTOL)goto 1
         rc=w*(1.d0+s*s*(a1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END

************************************************************************
      double precision FUNCTION rd(x,y,z)
************************************************************************
*     PURPOSE: Compute Carlson degenerate integral RD
*              R_D(x,y,z)=3/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2)
*     ARGUMENTS: x,y,z
*     ROUTINES CALLED:  None.
*     ALGORITHM: Due to B.C. Carlson.
*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
*     REMARKS:  
*     AUTHOR:  Press et al (1992)
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************
      double precision x,y,z,ERRTOL,a1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=0.0015d0,a1=3.d0/14.d0,C2=1.d0/6.d0,
     *C3=9.d0/22.d0,C4=3.d0/26.d0,C5=0.25d0*C3,C6=1.5d0*C4)
      double precision alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     *sqrtz,sum,xt,yt,zt
      xt=x
      yt=y
      zt=z
      sum=0.d0
      fac=1.d0
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=.2d0*(xt+yt+3.d0*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.d0*eb
      ee=ed+ec+ec
      rd=3.d0*sum+fac*(1.d0+ed*(-a1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*
     *ec+delz*C4*ea)))/(ave*sqrt(ave))
      return
      END

************************************************************************
      double precision FUNCTION rj(x,y,z,p)
************************************************************************
*     PURPOSE: Compute Carlson fundamental integral RJ
*     RJ(x,y,z,p) = 3/2 \int_0^\infty dt
*                      (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1)
*     ARGUMENTS: x,y,z,p
*     ROUTINES CALLED:  RF, RC.
*     ALGORITHM: Due to B.C. Carlson.
*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
*     REMARKS:  
*     AUTHOR:  Press et al (1992)
*     DATE WRITTEN:  25 Mar 91.
*     REVISIONS:
***********************************************************************      
      double precision p,x,y,z,ERRTOL,a1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=0.0015d0,a1=3.d0/14.d0,C2=1.d0/3.d0,
     *C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3,C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)
      double precision a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,
     *fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf
        sum=0.d0
        fac=1.d0
        if(p.gt.0.d0)then
          xt=x
          yt=y
          zt=z
          pt=p
        else
          xt=min(x,y,z)
          zt=max(x,y,z)
          yt=x+y+z-xt-zt
          a=1.d0/(yt-p)
          b=a*(zt-yt)*(yt-xt)
          pt=yt+b
          rho=xt*zt/yt
          tau=p*pt/yt
          rcx=rc(rho,tau)
        endif
1       continue
          sqrtx=sqrt(xt)
          sqrty=sqrt(yt)
          sqrtz=sqrt(zt)
          alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
          beta=pt*(pt+alamb)**2
          sum=sum+fac*rc(alpha,beta)
          fac=.25d0*fac
          xt=.25d0*(xt+alamb)
          yt=.25d0*(yt+alamb)
          zt=.25d0*(zt+alamb)
          pt=.25d0*(pt+alamb)
          ave=.2d0*(xt+yt+zt+pt+pt)
          delx=(ave-xt)/ave
          dely=(ave-yt)/ave
          delz=(ave-zt)/ave
          delp=(ave-pt)/ave
        if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
        ea=delx*(dely+delz)+dely*delz
        eb=delx*dely*delz
        ec=delp**2
        ed=ea-3.d0*ec
        ee=eb+2.d0*delp*(ea-ec)
        rj=3.d0*sum+fac*(1.d0+ed*(-a1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+
     *  delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
        if (p.le.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
      return
      END
