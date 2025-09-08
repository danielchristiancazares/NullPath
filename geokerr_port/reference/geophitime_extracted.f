      subroutine geophitime(u0,uf,mu0,muf,a,l,l2,q2,tpm,tpr,su,sm,iu,h1,phimu,tmu,
     &      ncase,u1,u2,u3,u4,phiu,tu,lambda,rffu0,rffu1,rffmu1,rffmu2,rffmu3,
     &      rdc,rjc,tu01,tu02,tu03,tu04,tmu1,tmu3,phimu1,phimu3,firstpt)
*******************************************************************************
*     PURPOSE: Computes final polar angle muf at an inverse radius uf given initial inverse radius
*               and polar angle, and geodesic constants of the motion.
*
*     INPUTS:   U0 -- Starting u value. If U0=0, DTI and LAMBDAI values will blow up.
*               UF -- Final u value.
*               MU0 -- Starting mu=cos(theta) value.
*               MUF -- Final mu value.
*               A -- Black hole spin, on interval [0,1).
*               L -- Dimensionless z-component of angular momentum.
*               L2 -- L*L
*               Q2 -- Dimensionless Carter's constant.
*               TPM -- Number of mu turning points between MU0 and MUF.
*               TPR -- Number of u turning points reached between U0 and UF.
*               SU -- Initial du/dlambda, =1 (-1) for ingoing (outgoing) rays.
*               SM -- Initital dmu/dlambda.
*               IU -- Value of IU=IMU
*               H1 -- Value of h1 from Eq. (21) if NCASE=6
*               NCASE -- Case label corresponding to Table 1.
*               U1,2,3,4 -- (Real) roots of U(u) in increasing order.
*               RFFU0 -- Value of RF relevant for U0 piece of IU integral.
*               RFFU1 -- Value of RF relevant for UF piece of IU integral.
*               RFFMU1 -- Value of RF relevant for MU0 piece of IMU integral.
*               RFFMU2 -- Value of RF relevant for MUF piece of IMU integral.
*               RFFMU3 -- Value of RF relevant for turning point piece of IMU integral.
*               FIRSTPT -- Boolean variable. If .TRUE., RDC, RJC, TU01, TU02, TU03, TU04, TMU1, TMU3, PHIMU1 and
*                          PHIMU3 are computed. Otherwise, all are input.
*     OUTPUTS:  PHIMU -- Value of MU term from Eq. (15)
*               TMU -- Value of MU term from Eq. (14)
*               PHIU -- Value of U term from Eq. (15)
*               TU -- Value of U term from Eq. (14)
*               LAMBDA -- Affine parameter from Eq. (49)
*               RDC -- Value of RD for the complete elliptic integral of the 2nd kind.
*               RJC -- Value of RJ for the complete elliptic integral of the 3rd kind.
*               TU01 -- Value of the U0 part of the fourth term of Eq. (47). Computed if FIRSTPT=.TRUE., ignored
*                       if no physical turning points are present. Input otherwise.
*               TU02 -- Value of the U0 part of the third term of Eq. (47). Computed if FIRSTPT=.TRUE., ignored
*                       if no physical turning points are present. Input otherwise.
*               TU03 -- Value of the U0 part of the first term of Eq. (47). Computed if FIRSTPT=.TRUE., ignored if
*                       no physical turning points are present. Input otherwise.
*               TU04 -- Value of the U0 part of the second term of Eq. (47). Computed if FIRSTPT=.TRUE., ignored if
*                       no physical turning points are present. Input otherwise.
*     ROUTINES CALLED: ELLCUBICREAL, ELLCUBICCOMPLEX, ELLQUARTICREAL, ELLQUARTICCOMPLEX, ELLDOUBLECOMPLEX,
*                       TFNKERR, PHIFNKERR, CALCPHITMUSYM, CALCPHITMUASYM                   
*     ACCURACY:   Machine.
*     REMARKS: Based on Dexter & Agol (2009), and labeled equations refer to that paper unless noted
*              otherwise.   
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit none
      double precision a,u1,u2,u3,u4,l,l2,mneg,mpos,mu0,muf,one,uplus,yy,
     *      phimu,phiu0,phiu,pi,pi2,q2,ql2,s1,sm,su,tmu,tu0,tu,two,umi,ee,dd,
     *      u0,uf,tmu1,tmu2,tmu3,lambda,
     *      tu01,tu11,tu02,tu12,tu03,tu13,tu04,tu14,phiu01,phiu02,phiu11,phiu12,rffu0,rffu1,
     *      ellcubicreal,ellcubiccomplex,ellquarticreal,ellquarticcomplex,elldoublecomplex,
     *      f1,f2,g1,g2,h1,h2,muplus,a1,a2,a3,tfnkerr,phifnkerr,qs,ur,f,g,h,tu1,tu2,
     *      tu3,tu4,phiu1,phiu2,phimu1,phimu2,vfm,vfp,vsm,vsp,phimu3,iu,lambdau,
     *      rffmu1,rffmu2,rffmu3,rdc,rjc
      integer tpm,tpr,ncase,p(5)
c firstpt is a boolean variable indicating if we need to calculate integrals involving only the initial and turning points.
      logical firstpt
      PARAMETER ( one=1.d0, two=2.d0 )
      pi=acos(-one)
      pi2=two*pi
      p(1)=-1
      p(2)=-1
      p(3)=-1
      uplus=one/(one+sqrt(one-a*a))
c Only calculate 1/u_- since u_- blows up when a=0.
      umi=one-sqrt(one-a*a)
      ur=-one/(two*sqrt(one-a*a))
      qs=sign(1.d0,q2)
      dd=two*((a-l)**two+q2)
c Eq. (22)
      ee=-a*a*q2
      ql2=q2+l2
      a1=sm
      a2=sm*(-one)**tpm
c Eq. (35)
      a3=two*int((two*dble(tpm)+3.d0-sm)/4.d0)-one
      if(ncase.eq.0) then
        tu=0.d0
        phiu=0.d0
        phimu=0.d0
        tmu=0.d0
        lambdau=0.d0
      elseif(ncase.lt.3) then
c These are the cubic real roots cases with u1<0<u2<=u3
        p(5)=0
        if(abs(u3-u2).lt.1.D-12) then
c These are the equal roots cases
          tu=su*tfnkerr(u0,uf,u1,u2,l,a)/sqrt(dd)
          phiu=su*(phifnkerr(uf,u1,u2,l,a)-phifnkerr(u0,u1,u2,l,a))/sqrt(dd)   
          if(u0.ge.u3) then
            phiu=-phiu
            tu=-tu
          endif       
        elseif(u0.le.u2) then
c Table 1 Row 1
          if(firstpt.and.(u0.ne.u2)) then
            p(4)=-2
c First three u0 integrals in Eq. (47)
            phiu01=ellcubicreal(p,-u1,one,u2,-one,u3,-one,one,-one/uplus,rffu0,u0,u2)
            phiu02=ellcubicreal(p,-u1,one,u2,-one,u3,-one,one,-umi,rffu0,u0,u2)
            tu02=ellcubicreal(p,-u1,one,u2,-one,u3,-one,0.d0,one,rffu0,u0,u2)
            tu03=phiu01
            tu04=phiu02
            p(4)=-4
c Fourth u0 integral in Eq. (47)
            tu01=ellcubicreal(p,-u1,one,u2,-one,u3,-one,0.d0,one,rffu0,u0,u2)
          elseif(u0.eq.u2) then
            tu01=0.d0
            tu02=0.d0
            tu03=0.d0
            tu04=0.d0
            phiu01=0.d0
            phiu02=0.d0
          else
            phiu01=tu03
            phiu02=tu04
          endif
          if(abs(uf-u2).gt.1d-16) then
            p(4)=-2
c First three uf integrals in Eq. (47)
            phiu11=ellcubicreal(p,-u1,one,u2,-one,u3,-one,one,-one/uplus,rffu1,uf,u2)
            phiu12=ellcubicreal(p,-u1,one,u2,-one,u3,-one,one,-umi,rffu1,uf,u2)
            tu12=ellcubicreal(p,-u1,one,u2,-one,u3,-one,0.d0,one,rffu1,uf,u2)
            tu13=phiu11
            tu14=phiu12
            p(4)=-4
c Fourth uf integral in Eq. (47)
            tu11=ellcubicreal(p,-u1,one,u2,-one,u3,-one,0.d0,one,rffu1,uf,u2)
          else
            tu11=0.d0
            tu12=0.d0
            tu13=0.d0
            tu14=0.d0
            phiu11=0.d0
            phiu12=0.d0
          endif
c Eq. (47) for u0 piece
          tu0=(umi-one/uplus)*tu01+(umi**two-one/uplus**two)*tu02-
     *      (2.d0*a*(a-l)+a**two/uplus+one/uplus**3)*tu03+
     *      (2.d0*a*(a-l)+a**two*umi+umi**3)*tu04
c Eq. (47) for uf piece
          tu1=(umi-one/uplus)*tu11+(umi**two-one/uplus**two)*tu12-
     *      (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu13+
     *      (2.d0*a*(a-l)+a**2*umi+umi**3)*tu14
c Combine according to Eq. (54)
          tu=su*(tu0-(-one)**tpr*tu1)/sqrt(dd)*ur
c Eq. (48) for u0 piece
          phiu0=-(l/uplus+2.d0*(a-l))*phiu01+(l*umi+2.d0*(a-l))*phiu02
c Eq. (48) for uf piece
          phiu1=-(l/uplus+2.d0*(a-l))*phiu11+(l*umi+2.d0*(a-l))*phiu12
c Combine according to Eq. (54)
          phiu=su*(phiu0-(-one)**tpr*phiu1)/sqrt(dd)*ur
c Eq. (49)
          lambdau=su*(tu01-(-one)**tpr*tu11)/sqrt(dd)
        elseif(u0.ge.u3) then
c Table 1 Row 2
          if(firstpt.and.(u0.ne.u3)) then
            p(4)=-2
c First three u0 integrals in Eq. (47)
            phiu01=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,one,-one/uplus,rffu0,u3,u0)
            phiu02=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,one,-umi,rffu0,u3,u0)
            tu02=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,0.d0,one,rffu0,u3,u0)
            tu03=phiu01
            tu04=phiu02
c Fourth u0 integral in Eq. (47)
            p(4)=-4
            tu01=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,0.d0,one,rffu0,u3,u0)
          elseif(u0.eq.u3) then
            tu01=0.d0
            tu02=0.d0
            tu03=0.d0
            tu04=0.d0
            phiu01=0.d0
            phiu02=0.d0
          else
            phiu01=tu03
            phiu02=tu04
          endif
          if(uf.ne.u3) then
            p(4)=-2
c First three uf integrals in Eq. (47)
            phiu11=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,one,-one/uplus,rffu1,u3,uf)
            phiu12=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,one,-umi,rffu1,u3,uf)
            tu12=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,0.d0,one,rffu1,u3,uf)
            tu13=phiu11
            tu14=phiu12
            p(4)=-4
c Fourth uf integral in Eq. (47)
            tu11=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,0.d0,one,rffu1,u3,uf)
          else
            tu11=0.d0
            tu12=0.d0
            tu13=0.d0
            tu14=0.d0
            phiu11=0.d0
            phiu12=0.d0
          endif
c Eq. (47) for u0 part
          tu0=(umi-one/uplus)*tu01+(umi**2-one/uplus**2)*tu02-
     *      (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu03+
     *      (2.d0*a*(a-l)+a**2*umi+umi**3)*tu04
c Eq. (47) for uf part
          tu1=(umi-one/uplus)*tu11+(umi**2-one/uplus**2)*tu12-
     *      (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu13+
     *      (2.d0*a*(a-l)+a**2*umi+umi**3)*tu14
c Combine according to Eq. (54)
          tu=su*(tu0-(-one)**tpr*tu1)/sqrt(dd)*ur
c Eq. (48) for u0 part
          phiu0=-(l/uplus+2.d0*(a-l))*phiu01+(l*umi+2.d0*(a-l))*phiu02
c Eq. (48) for uf part
          phiu1=-(l/uplus+2.d0*(a-l))*phiu11+(l*umi+2.d0*(a-l))*phiu12
c Combine according to Eq. (54)
          phiu=su*(phiu0-(-one)**tpr*phiu1)/sqrt(dd)*ur
c Eq. (49)
          lambdau=su*(tu01-(-one)**tpr*tu11)/sqrt(dd)
        endif
      elseif(ncase.eq.3) then
c This is a cubic complex case with one real root.
c Table 1 Row 3
        ncase=3
        f=-one/dd/u1
        g=f/u1
        h=one
        if(u0.lt.uf) then
          p(4)=-4
c Fourth integral in Eq. (47)
          tu1=ellcubiccomplex(p,-u1,one,0.d0,one,f,g,h,rffu0,u0,uf)
          p(4)=-2
c First three integrals in Eq. (47)
          tu2=ellcubiccomplex(p,-u1,one,0.d0,one,f,g,h,rffu0,u0,uf)
          tu3=ellcubiccomplex(p,-u1,one,one,-one/uplus,f,g,h,rffu0,u0,uf)
          tu4=ellcubiccomplex(p,-u1,one,one,-umi,f,g,h,rffu0,u0,uf)
c Eq. (47)
          tu=su*ur/sqrt(dd)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Minus sign in phi components is from flipping the sign of the u/u_\pm-1 factor to keep arguments positive.
          phiu0=-tu3
          phiu1=-tu4
c Eq. (48)
          phiu=su*ur/sqrt(dd)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
          lambdau=su*tu1/sqrt(dd)
        elseif(u0.gt.uf) then
          p(4)=-4
c Fourth integral in Eq. (47)
          tu1=-ellcubiccomplex(p,-u1,one,0.d0,one,f,g,h,rffu0,uf,u0)
          p(4)=-2
c First three integrals in Eq. (47)
          tu2=-ellcubiccomplex(p,-u1,one,0.d0,one,f,g,h,rffu0,uf,u0)
          tu3=-ellcubiccomplex(p,-u1,one,one,-one/uplus,f,g,h,rffu0,uf,u0)
          tu4=-ellcubiccomplex(p,-u1,one,one,-umi,f,g,h,rffu0,uf,u0)
c Eq. (47)
          tu=su*ur/sqrt(dd)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Minus sign in phi components is from flipping the sign of the u/u_\pm-1 factor to keep arguments positive.
          phiu0=-tu3
          phiu1=-tu4
c Eq. (48)
          phiu=su*ur/sqrt(dd)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
c Eq. (49)
          lambdau=su*tu1/sqrt(dd)
        else
          tu=0.d0
          phiu=0.d0
          lambdau=0.d0
        endif
      endif
      if(q2.eq.0.d0) then
c Calculate mu components in the special case q2=0, where the t and phi integrals are elementary.
        s1=sign(1.d0,mu0)
        a1=s1*sm
        a2=s1*sm*(-one)**(tpm+1)
        if(abs(l).lt.abs(a)) then
          muplus=s1*sqrt(one-l2/a/a)
          phimu1=sign(1.d0,muplus*l/a)*(atan((muplus-tan(asin(mu0/muplus)/two))/sqrt(one-muplus**2))+
     *            atan((muplus+tan(asin(mu0/muplus)/two))/sqrt(one-muplus**2)))
          phimu2=sign(1.d0,muplus*l/a)*(atan((muplus-tan(asin(muf/muplus)/two))/sqrt(one-muplus**2))+
     *           atan((muplus+tan(asin(muf/muplus)/two))/sqrt(one-muplus**2)))
          phimu=a1*phimu1+a2*phimu2
          tmu=abs(muplus*a)*(a2*sqrt(1.d0-muf**2/muplus**2)+a1*sqrt(1.d0-mu0**2/muplus**2))
        else
          phimu=0.d0
          tmu=0.d0
        endif
      elseif(a.eq.0.d0) then
c Calculate mu components in the special case a=0. The t and phi integrals are again elementary.
        tmu=0.d0
        muplus=sqrt(q2/ql2)
        vfm= (muf-muplus**2)/(one-muf)/muplus
        vfp=-(muf+muplus**2)/(one+muf)/muplus
c This is code to suppress floating point errors in phimu calculation.
        if(abs(vfm).gt.one) vfm=sign(1.d0,vfm)*one
        if(abs(vfp).gt.one) vfp=sign(1.d0,vfp)*one
        if((mu0-muplus**2).eq.0.d0) then
          vsm=-one
        else
          vsm=(mu0-muplus**2)/(one-mu0)/muplus
        endif
        if((mu0+muplus**2).eq.0.d0) then
          vsp=-one
        else
          vsp=-(mu0+muplus**2)/(one+mu0)/muplus
        endif
        if(abs(vsp).gt.one) vsp=sign(1.d0,vsp)*one
        if(abs(vsm).gt.one) vsm=sign(1.d0,vsm)*one
        phimu1=pi-asin(vsm)+asin(vsp)
        phimu2=asin(vfm)-asin(vfp)+pi
        phimu3=two*pi
        phimu=-l*iu+sign(1.d0,l)*0.5d0*(a1*phimu1+a2*phimu2+a3*phimu3)
      endif
      if(ncase.eq.4) then
c This is the special case where q2=0 and l=a. U(u)=1 and the t and phi u components are elementary.
c Eq. (48)
        phiu=su*ur*(l*log((uf/uplus-one)/(u0/uplus-one))
     *       -l*umi*umi*log((uf*umi-one)/(u0*umi-one)))
c Eq. (47)
        tu=su*ur*((umi-one/uplus)*(one/u0-one/uf)+(umi**2-one/uplus**2)*log(uf/u0)+
     *     (a**2/uplus+one/uplus**3)*uplus*log((uf/uplus-one)/(u0/uplus-one))-
     *     (a**2*umi+umi**3)*umi*log((uf*umi-one)/(u0*umi-one)))
c Eq. (49)
        lambdau=su*(one/u0-one/uf)
      elseif(ncase.eq.5) then
c This is the quartic case with one pair of complex roots
        p(4)=-1
c Table 1 Row 5
        f=-qs*one/abs(ee)/u1/u4
        g=(u4+u1)/u1/u4*f
        h=one
        if(u0.lt.uf) then
          p(5)=-4
c Fourth integral in Eq. (47)
          tu1=ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,0.d0,one,f,g,h,rffu0,u0,uf)
          p(5)=-2
c First three integrals in Eq. (47)
          tu2=ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,0.d0,one,f,g,h,rffu0,u0,uf)
          tu3=ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,one,-one/uplus,f,g,h,rffu0,u0,uf)
          tu4=ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,one,-umi,f,g,h,rffu0,u0,uf)
c Eq. (47)
          tu=su*ur/sqrt(-qs*ee)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Eq. (49)
          lambdau=su*tu1/sqrt(abs(ee))
c Minus sign in phi components is from flipping the sign of the u/u_\pm-1 factor to keep arguments positive.
          phiu0=-tu3
          phiu1=-tu4
c Eq. (48)
          phiu=su*ur/sqrt(-qs*ee)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
        elseif(u0.gt.uf) then
          p(5)=-4
c Fourth integral in Eq. (47)
          tu1=-ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,0.d0,one,f,g,h,rffu0,uf,u0)
          p(5)=-2
c First three integrals in Eq. (47)
          tu2=-ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,0.d0,one,f,g,h,rffu0,uf,u0)
          tu3=-ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,one,-one/uplus,f,g,h,rffu0,uf,u0)
          tu4=-ellquarticcomplex(p,-u1,one,u4*qs,-one*qs,one,-umi,f,g,h,rffu0,uf,u0)
c Eq. (47)
          tu=su*ur/sqrt(-qs*ee)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Minus sign in phi components is from flipping the sign of the u/u_\pm-1 factor to keep arguments positive.
          phiu0=-tu3
          phiu1=-tu4
c Eq. (48)
          phiu=su*ur/sqrt(-qs*ee)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
c Eq. (49)
          lambdau=su*tu1/sqrt(abs(ee))
        else
          tu=0.d0
          phiu=0.d0
          lambdau=0.d0
        endif        
      elseif(ncase.eq.6) then
c This is the quartic complex case with no real roots.
        p(4)=-1
c Table 1 Row 6
        h2=one/h1
        g1=dd/ee/(h2-h1)
        g2=-g1
        f1=one/sqrt(ee)
        f2=f1
        if(u0.lt.uf) then
          p(5)=-4
c Fourth integral in Eq. (47)
          tu1=elldoublecomplex(p,f1,g1,h1,f2,g2,h2,0.d0,one,rffu0,u0,uf)
          p(5)=-2
c First three integrals in Eq. (47)
          tu2=elldoublecomplex(p,f1,g1,h1,f2,g2,h2,0.d0,one,rffu0,u0,uf)
          tu3=elldoublecomplex(p,f1,g1,h1,f2,g2,h2,one,-one/uplus,rffu0,u0,uf)
          tu4=elldoublecomplex(p,f1,g1,h1,f2,g2,h2,one,-umi,rffu0,u0,uf)
c Minus sign in phi components is from flipping the sign of the u/u_\pm-1 factor to keep arguments positive.
          phiu0=-tu3
          phiu1=-tu4
c Eq. (47)
          tu=su*ur/sqrt(-qs*ee)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Eq. (48)
          phiu=su*ur/sqrt(-qs*ee)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
c Eq. (49)
          lambdau=su*tu1/sqrt(-qs*ee)
        elseif(u0.gt.uf) then
          p(5)=-4
c Fourth integral in Eq. (47)
          tu1=-elldoublecomplex(p,f1,g1,h1,f2,g2,h2,0.d0,one,rffu0,uf,u0)
          p(5)=-2
c First three integrals in Eq. (47)
          tu2=-elldoublecomplex(p,f1,g1,h1,f2,g2,h2,0.d0,one,rffu0,uf,u0)
          tu3=-elldoublecomplex(p,f1,g1,h1,f2,g2,h2,one,-one/uplus,rffu0,uf,u0)
          tu4=-elldoublecomplex(p,f1,g1,h1,f2,g2,h2,one,-umi,rffu0,uf,u0)
c Minus sign in phi components is from flipping the sign of the u/u_\pm-1 factor to keep arguments positive.
          phiu0=-tu3
          phiu1=-tu4
c Eq. (47)
          tu=su*ur/sqrt(-qs*ee)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Eq. (48)
          phiu=su*ur/sqrt(-qs*ee)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
c Eq. (49)
          lambdau=su*tu1/sqrt(-qs*ee)
        else
          tu=0.d0
          phiu=0.d0
        endif
      elseif(ncase.gt.6) then
c These are the quartic cases with all real roots
        p(4)=-1
        if(abs(u3-u2).lt.1.D-12) then
c These are the equal roots quartic cases
          if(u0.lt.uf) then
            if(u0.lt.u2) then
c Table 1 Row 7
              p(5)=-2
c First three integrals in Eq. (47)
              phiu1=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/uplus,rffu0,u0,uf)
              phiu2=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-umi,rffu0,u0,uf)
              tu2=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu0,u0,uf)
              tu3=phiu1
              tu4=phiu2
              p(5)=-4
c Fourth integral in Eq. (47)
              tu1=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu0,u0,uf)     
            elseif(u0.gt.u3) then 
c Table 1 Row 8
c First three integrals in Eq. (47)
              phiu1=ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/uplus,rffu0,u0,uf)
              phiu2=ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-umi,rffu0,u0,uf)
              tu2=ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu0,u0,uf)
              tu3=phiu1
              tu4=phiu2
              p(5)=-4
c Fourth integral in Eq. (47)
              tu1=ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu0,u0,uf)
            else
              tu11=0.d0
              tu12=0.d0
              tu13=0.d0
              tu14=0.d0
              phiu11=0.d0
              phiu12=0.d0
            endif
c Eq. (47)
            tu=su*ur/sqrt(-qs*ee)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Eq. (48)
            phiu=su*ur/sqrt(-qs*ee)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
c Eq. (49)
            lambdau=tu1/sqrt(-qs*ee)
          elseif(u0.gt.uf) then
            if(u0.lt.u2) then
c Table 1 Row 7
              p(5)=-2
c First three integrals in Eq. (47)
              phiu1=-ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/uplus,rffu0,uf,u0)
              phiu2=-ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-umi,rffu0,uf,u0)
              tu2=-ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu0,uf,u0)
              tu3=-phiu1
              tu4=-phiu2
              p(5)=-4
c Fourth integral in Eq. (47)
              tu1=-ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu0,uf,u0)
            elseif(u0.gt.u3) then
c Table 1 Row 8
              p(5)=-2
c First three integrals in Eq. (47)
              phiu1=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/uplus,rffu0,uf,u0)
              phiu2=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-umi,rffu0,uf,u0)
              tu2=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu0,uf,u0)
              tu3=-phiu1
              tu4=-phiu2
              p(5)=-4
c Fourth integral in Eq. (47)
              tu1=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu0,uf,u0)
            else
              tu11=0.d0
              tu12=0.d0
              tu13=0.d0
              tu14=0.d0
              phiu11=0.d0
              phiu12=0.d0
            endif
c Eq. (47)
            tu=su*ur/sqrt(-qs*ee)*((umi-one/uplus)*tu1+(umi**2-one/uplus**2)*tu2-
     *       (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu3+
     *       (2.d0*a*(a-l)+a**2*umi+umi**3)*tu4)
c Eq. (48)
            phiu=su*ur/sqrt(-qs*ee)*((l/uplus+2.d0*(a-l))*phiu0-(l*umi+2.d0*(a-l))*phiu1)
c Eq. (49)
            lambdau=tu1/sqrt(-qs*ee)    
          else
            tu=0.d0
            phiu=0.d0
            lambdau=0.d0
          endif
        elseif(u0.le.u2) then
c This is the quartic case with distinct, real roots
c Table 1 Row 7
          if(firstpt.and.(u0.ne.u2)) then
c Only compute integrals involving u0 once per geodesic
            p(5)=-2
c First three u0 integrals in Eq. (47)
            phiu01=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/uplus,rffu0,u0,u2)
            phiu02=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-umi,rffu0,u0,u2)
            tu02=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu0,u0,u2)
            tu03=phiu01
            tu04=phiu02
            p(5)=-4
c Fourth u0 integral in Eq. (48)
            tu01=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu0,u0,u2)
          elseif(u0.eq.u2) then
            tu01=0.d0
            tu02=0.d0
            tu03=0.d0
            tu04=0.d0
            phiu01=0.d0
            phiu02=0.d0
          else
            phiu01=tu03
            phiu02=tu04
          endif
          if(abs(uf-u2).gt.1d-12) then
            p(5)=-2
c First three uf integrals in Eq. (47)
            phiu11=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-one/uplus,rffu1,uf,u2)
            phiu12=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,one,-umi,rffu1,uf,u2)
            tu12=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu1,uf,u2)
            tu13=phiu11
            tu14=phiu12
            p(5)=-4
c Fourth uf integral in Eq. (47)
            tu11=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,one,rffu1,uf,u2)
          else
            tu11=0.d0
            tu12=0.d0
            tu13=0.d0
            tu14=0.d0
            phiu11=0.d0
            phiu12=0.d0
          endif
c Eq. (47) for u0 piece
          tu0=(umi-one/uplus)*tu01+(umi**2-one/uplus**2)*tu02-
     *        (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu03+
     *        (2.d0*a*(a-l)+a**2*umi+umi**3)*tu04
c Eq. (47) for uf piece
          tu1=(umi-one/uplus)*tu11+(umi**2-one/uplus**2)*tu12-
     *        (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu13+
     *        (2.d0*a*(a-l)+a**2*umi+umi**3)*tu14
          tu=su*(tu0-(-one)**tpr*tu1)/sqrt(abs(ee))*ur
c Eq. (48)
          phiu0=-(l/uplus+2.d0*(a-l))*phiu01+(l*umi+2.d0*(a-l))*phiu02
          phiu1=-(l/uplus+2.d0*(a-l))*phiu11+(l*umi+2.d0*(a-l))*phiu12
c Combine according to Eq. (54)
          phiu=su*(phiu0-(-one)**tpr*phiu1)/sqrt(abs(ee))*ur
c Eq. (49)
          lambdau=su*(tu01-(-one)**tpr*tu11)/sqrt(abs(ee))
        elseif(u0.ge.u3) then
          if(firstpt.and.(u0.ne.u3)) then
            p(5)=-2
c First three u0 integrals in Eq. (47)
            phiu01=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/uplus,rffu0,u3,u0)
            phiu02=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-umi,rffu0,u3,u0)
            tu02=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu0,u3,u0)
            tu03=phiu01
            tu04=phiu02
            p(5)=-4
c Fourth u0 integral in Eq. (47)
            tu01=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu0,u3,u0)
          elseif(u0.eq.u3) then
            tu01=0.d0
            tu02=0.d0
            tu03=0.d0
            tu04=0.d0
            phiu01=0.d0
            phiu02=0.d0
          else
            phiu01=tu03
            phiu02=tu04
          endif
          if(uf.ne.u3) then
            p(5)=-2
c First three uf integrals in Eq. (47)
            phiu11=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-one/uplus,rffu1,u3,uf)
            phiu12=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,one,-umi,rffu1,u3,uf)
            tu12=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu1,u3,uf)
            tu13=phiu11
            tu14=phiu12
            p(5)=-4
c Fourth uf integral in Eq. (47)
            tu11=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,one,rffu1,u3,uf)
          else
            tu11=0.d0
            tu12=0.d0
            tu13=0.d0
            tu14=0.d0
            phiu11=0.d0
            phiu12=0.d0
          endif
c Eq. (47) for u0 piece
          tu0=(umi-one/uplus)*tu01+(umi**2-one/uplus**2)*tu02-
     *        (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu03+
     *        (2.d0*a*(a-l)+a**2*umi+umi**3)*tu04
c Eq. (47) for uf piece
          tu1=(umi-one/uplus)*tu11+(umi**2-one/uplus**2)*tu12-
     *        (2.d0*a*(a-l)+a**2/uplus+one/uplus**3)*tu13+
     *        (2.d0*a*(a-l)+a**2*umi+umi**3)*tu14
c Eq. (47)
          tu=su*(tu0-(-one)**tpr*tu1)/sqrt(abs(ee))*ur
c Eq. (48) for u0 piece
          phiu0=-(l/uplus+2.d0*(a-l))*phiu01+(l*umi+2.d0*(a-l))*phiu02
c Eq. (48) for uf piece
          phiu1=-(l/uplus+2.d0*(a-l))*phiu11+(l*umi+2.d0*(a-l))*phiu12
c Eq. (48)
          phiu=su*(phiu0-(-one)**tpr*phiu1)/sqrt(abs(ee))*ur
c Eq. (49) for u term
          lambdau=su*(tu01-(-one)**tpr*tu11)/sqrt(abs(ee))
        endif      
      endif
      if(ncase.gt.4) then
c Find roots of biquadratic M(mu).
        yy=-0.5d0*(a*a-ql2+sign(one,a*a-ql2)*sqrt((a*a-ql2)**2+4.d0*q2*a*a))
        if((a*a-ql2).lt.0.d0) then
          mneg=-yy/a/a
          mpos=q2/yy
        else
          mneg=q2/yy
          mpos=-yy/a/a
        endif
c Protect against rounding error in mpos:
        if(mpos.gt.1.d0) mpos=1.d0
        muplus=sqrt(mpos)
c NOTE: This formula uses a slightly different prescription for splitting up the mu integral.
c       All integrals are with respect to muplus, but the procedure of splitting into coefficients
c       and finding them by writing out specific cases is the same.
        a1=sm
        a2=sm*(-1.d0)**(tpm+1)
        a3=2.d0*int((2*tpm-sm+1)/4.d0)
        if(mneg.lt.0.d0) then
c Protect against rounding errors in roots:
          if(muplus.lt.mu0) then
            muplus=mu0
            mpos=muplus*muplus
          endif
c This is the symmetric roots case, where the orbit can cross the equatorial plane.
c Calculate phi, t mu component integrals:
          call calcphitmusym(a,mneg,mpos,mu0,muf,muplus,phimu1,phimu2,phimu3,tmu1,
     &                       tmu2,tmu3,rffmu1,rffmu2,rffmu3,rdc,rjc,firstpt)
c Eq. (45)
          tmu=a*a*mneg*iu+(a1*tmu1+a2*tmu2+a3*tmu3)
c Eq. (46)
          phimu=-l*iu+l*(a1*phimu1+a2*phimu2+a3*phimu3)
       else
c This is the asymmetric roots case
          if(sign(1.d0,mu0).eq.-1) muplus=-muplus
c Protect for rounding error when mu0 is a turning point:
          if(abs(muplus).lt.abs(mu0)) then
            muplus=mu0
            mpos=muplus*muplus
          endif
          if(abs(mneg).gt.mu0*mu0) then
            mneg=mu0*mu0
          endif
c Calculate phi, t mu component integrals:
          call calcphitmuasym(a,mneg,mpos,mu0,muf,muplus,phimu1,phimu2,phimu3,tmu1,
     &                        tmu2,tmu3,rffmu1,rffmu2,rffmu3,rdc,rjc,firstpt)
c Eq. (45)
          tmu=(a1*tmu1+a2*tmu2+a3*tmu3)
c Eq. (46)    
          phimu=-l*iu+l*(a1*phimu1+a2*phimu2+a3*phimu3)
        endif                
      endif
c Eq. (49)
      lambda=lambdau+tmu
!      write(6,*) 'geokerr geophitime lambda: ',lambdau,tmu,phimu,tu01,tu11,ee
!      write(6,*) 'geokerr geophitime u: ',uf,u2,u3,uf-u2
!      write(6,*) 'geokerr geophitime phi: ',phiu11,phiu,ur
      if(l.eq.0.d0) phiu=phiu-sign(1.d0,a)*pi*tpm
      return
      end
