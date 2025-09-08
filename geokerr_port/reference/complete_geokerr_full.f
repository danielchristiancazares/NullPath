      PROGRAM COMPLETE_GEOKERR_WRAPPER
C*******************************************************************************
C     COMPLETE GEOKERR VALIDATION WRAPPER
C     
C     This program provides a clean interface to the complete GEOKERR algorithm
C     for validation against our CUDA implementation. Takes input parameters
C     and outputs all computed values including lambda (Mino time).
C     
C     Input format: u0 uf mu0 muf a l q2 alpha beta su sm
C     Output format: All GEOKERR results for comprehensive validation
C*******************************************************************************
      IMPLICIT NONE
      
C     GEOKERR input parameters  
      DOUBLE PRECISION u0, uf, uout, mu0, muf, a, l, q2
      DOUBLE PRECISION alpha, beta
      INTEGER tpm, tpr, su, sm, nup, ncase, kext, npts
      DOUBLE PRECISION offset, l2
      LOGICAL phit, usegeor, mufill
      
C     GEOKERR outputs
      DOUBLE PRECISION ufi(1), mufi(1), dti(1), dphi(1), lambdai(1)
      INTEGER tpmi(1), tpri(1)
      
C     Additional geodesic computation parameters
      DOUBLE PRECISION iu, h1, u1, u2, u3, u4
      DOUBLE PRECISION rffu0, rffu1, rffmu1, rffmu2, rffmu3
      DOUBLE PRECISION iu0, i1mu, i3mu
      DOUBLE PRECISION phimu, tmu, phiu, tu, lambda
      DOUBLE PRECISION rdc, rjc, tu01, tu02, tu03, tu04
      DOUBLE PRECISION tmu1, tmu3, phimu1, phimu3
      LOGICAL firstpt, pht
      
C     File handling
      INTEGER i, ntests, ios
      
      WRITE(6,*) 'COMPLETE GEOKERR VALIDATION WRAPPER'
      WRITE(6,*) '=================================='
      WRITE(6,*) 'Reference implementation for CUDA validation'
      WRITE(6,*) ''
      
C     Read number of test cases
      read(5,*,iostat=ios) ntests
      IF (ios .NE. 0) THEN
        WRITE(6,*) 'ERROR: Cannot read number of test cases'
        STOP 1
      ENDIF
      
      WRITE(6,*) 'Processing', ntests, 'complete geodesic test cases...'
      WRITE(6,*) ''
      WRITE(6,*) '# Complete GEOKERR validation output:'
      WRITE(6,*) '# lambda phimu tmu phiu tu ncase uf_computed muf_computed'
      
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
        nup = 1
        npts = 0
        kext = 0
        phit = .TRUE.      ! Need phi/time computation
        usegeor = .FALSE.  ! Solve for muf given uf
        mufill = .FALSE.   
        tpm = 0           ! Will be computed
        tpr = 0           ! Will be computed
        ncase = 0         ! Will be computed
        
C       Call main GEOKERR subroutine
        CALL GEOKERR(u0, uf, uout, mu0, muf, a, l, q2, 
     &               alpha, beta, tpm, tpr, su, sm, nup, offset, 
     &               phit, usegeor, mufill, ncase, kext, npts,
     &               ufi, mufi, dti, dphi, tpmi, tpri, lambdai)
        
C       Extract individual geodesic components by calling GEOMU and GEOR
C       to get more detailed validation data
        firstpt = .TRUE.
        pht = .TRUE.
        
C       Call GEOMU to compute polar angle evolution
        CALL GEOMU(u0, uf, mu0, muf, a, l, l2, q2, iu, tpm, tpr, 
     &             su, sm, ncase, h1, u1, u2, u3, u4, rffu0, rffu1,
     &             rffmu1, rffmu2, rffmu3, iu0, i1mu, i3mu, pht, firstpt)
        
C       Call GEOPHITIME to get detailed Mino time breakdown
        CALL GEOPHITIME(u0, uf, mu0, muf, a, l, l2, q2, 
     &                  tpm, tpr, su, sm, iu, h1, 
     &                  phimu, tmu, ncase, u1, u2, u3, u4,
     &                  phiu, tu, lambda, rffu0, rffu1, 
     &                  rffmu1, rffmu2, rffmu3, rdc, rjc,
     &                  tu01, tu02, tu03, tu04, tmu1, tmu3,
     &                  phimu1, phimu3, firstpt)
        
C       Output comprehensive validation data
        WRITE(6,100) lambda, phimu, tmu, phiu, tu, ncase, 
     &               ufi(1), mufi(1)
100     FORMAT(6(E23.15,1X),I3,1X,2(E23.15,1X))
        
      ENDDO
      
      WRITE(6,*) ''
      WRITE(6,*) 'Complete GEOKERR validation completed!'
      WRITE(6,*) 'Use this data to validate CUDA implementation'
      
      END      subroutine geomu(u0,uf,mu0,muf,a,l,l2,q2,iu,tpm,tpr,su,sm,ncase,h1,
     &     u1,u2,u3,u4,rffu0,rffu1,rffmu1,rffmu2,rffmu3,iu0,i1mu,i3mu,pht,firstpt)
*******************************************************************************
*     PURPOSE: Computes final polar angle muf at an inverse radius uf given initial inverse radius
*               and polar angle, and geodesic constants of the motion.
*
*     INPUTS:   U0 -- Starting u value. If U0=0, DTI and LAMBDAI values will blow up.
*               UF -- Final u value.
*               MU0 -- Starting mu=cos(theta) value.
*               A -- Black hole spin, on interval [0,1).
*               L -- Dimensionless z-component of angular momentum.
*               L2 -- L*L
*               Q2 -- Dimensionless Carter's constant.
*               TPR -- Number of u turning points reached between U0 and UF.
*               SU -- Initial du/dlambda, =1 (-1) for ingoing (outgoing) rays.
*               SM -- Initital dmu/dlambda.
*               PHT -- Boolean variable. If .TRUE., RFFMU2 is computed for use in
*                      SUBROUTINE GEOPHITIME.
*               FIRSTPT -- Boolean variable. If .TRUE., roots of U(u) U1,2,3,4 and NCASE are computed as well as H1
*                           if NCASE=6
*     OUTPUTS:  MUF -- Final polar angle.
*               IU -- Value of IU integral between U0 and UF.
*               TPM -- Number of mu turning points reached between MU0 and MUF.
*               NCASE -- Case number corresponding to Table 1. Output if FIRSTPT=.TRUE., input otherwise.
*               H1 -- Value of h1 from Eq. (21) for given constants of motion if NCASE=6. Output if FIRSTPT=.TRUE.,
*                     input otherwise.
*               U1,U2,U3,U4 -- Increasing (real) roots. Output if FIRSTPT=.TRUE., input otherwise.
*               RFFU0 -- Value of RF relevant for U0 piece of IU integral. Computed if FIRSTPT=.TRUE., input otherwise.
*               RFFU1 -- Value of RF relevant for UF piece of IU integral. Computed if PHT=.TRUE.
*               RFFMU1 -- Value of RF relevant for MU0 piece of IMU integral. Computed if FIRSTPT=.TRUE., input
*                         otherwise.
*               RFFMU2 -- Value of RF relevant for MUF piece of IMU integral. Computed if PHT=.TRUE.
*               RFFMU3 -- Value of RF relevant for turning point piece of IMU integral. Computed if FIRSTPT=.TRUE.,
*                         input otherwise.
*               IU0 -- Value of IU integral between U0 and relevant turning point if one exists. Computed if
*                      u turning point present and FIRSTPT=.TRUE., input otherwise. Ignored if no turning point
*                      is present.
*               I1MU -- Value of IMU integral between MU0 and MUPLUS. Computed if FIRSTPT=.TRUE., input otherwise.
*               I3MU -- Value of IMU integral between MUMINUS and MUPLUS. Computed if FIRSTPT=.TRUE., input otherwise.
*     ROUTINES CALLED: SNCNDN, ZROOTS, CALCIMUSYM, CALCIMUSYMF, CALCIMUASYM, CALCIMUASYMF, ELLCUBICREAL, 
*                       ELLCUBICCOMPLEX, ELLQUARTICREAL, ELLQUARTICCOMPLEX, ELLDOUBLECOMPLEX
*                    
*     ACCURACY:   Machine.
*     REMARKS: Based on Dexter & Agol (2009), and labeled equations refer to that paper unless noted
*              otherwise.   
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit none
      integer i,nreal,p(5)
      double precision a,aa,a1,a2,a3,cn,dn,dis,rr,i1mu,i3mu,rffu0,rffu1,rffmu1,rffmu2,rffmu3,
     *      f,iu,qs,h1,h2,f1,f2,g1,g2,l,l2,m1,mneg,mpos,mu0,muf,muplus,uarg,
     *      one,pi,pi2,q2,ql2,s1,sn,sm,su,theta,third,two,bb,farg,sarg,i2mu,
     *      u0,uf,u1,u2,u3,u4,g,cc,dd,ee,iu0,iu1,ellcubicreal,ellcubiccomplex,ellquarticreal,
     *      ellquarticcomplex,elldoublecomplex,asech,muminus,iarg,qq,uplus,yy
      integer tpm,tpr,ncase
      logical pht,firstpt
      complex*16 c(5),root(4),coefs(7),hroots(6)
      PARAMETER ( ONE=1.D0, TWO=2.D0, THIRD = 0.3333333333333333D0 )
      pi=acos(-one)
      pi2=two*pi
      uplus=one/(one+sqrt(one-a*a))
      l2=l*l
      ql2=q2+l2
c  Eq. (13): Coefficients of quartic in U(u)
      cc=a*a-q2-l2
      dd=two*((a-l)**2+q2)
      ee=-a*a*q2
c  Determine if U is cubic or quartic and find roots.
      if((ee.eq.0.d0) .and. (dd.ne.0.d0)) then
        p(1)=-1
        p(2)=-1
        p(3)=-1
        p(4)=0
        p(5)=0
        qq=cc*cc/dd/dd/9.d0
        rr=(two*cc**3/dd**3+27.d0/dd)/54.d0
        dis=rr*rr-qq**3
        if(dis.lt.-1.d-16) then
c These are the cubic real roots cases with u1<0<u2<=u3
          theta=acos(rr/qq**1.5d0)
          u1=-two*sqrt(qq)*cos(theta/3.d0)-cc/dd/3.d0
          u2=-two*sqrt(qq)*cos((theta-two*pi)/3.d0)-cc/dd/3.d0
          u3=-two*sqrt(qq)*cos((theta+two*pi)/3.d0)-cc/dd/3.d0
          iu1=0.d0
80        continue
          if(u0.le.u2) then
            if(uf.gt.u2) uf=u2
            ncase=1
c Table 1 Row 1
            if(firstpt.and.(u0.ne.u2)) then
              iu0=ellcubicreal(p,-u1,one,u2,-one,u3,-one,0.d0,0.d0,rffu0,u0,u2)
            elseif(u0.eq.u2) then
              iu0=0.d0
            endif
            if(uf.ne.u2) iu1=ellcubicreal(p,-u1,one,u2,-one,u3,-one,0.d0,0.d0,rffu1,uf,u2)
            iu=su*(iu0-(-one)**tpr*iu1)/sqrt(dd)
          elseif(u0.ge.u3) then
            if(uf.lt.u3) uf=u3
            ncase=2
c Table 1 Row 2
            if(firstpt.and.(u0.ne.u3)) then
              iu0=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,0.d0,0.d0,rffu0,u3,u0)
            elseif(u0.eq.u3) then
              iu0=0.d0
            endif
            if(uf.ne.u3) iu1=-ellcubicreal(p,-u1,one,-u2,one,-u3,one,0.d0,0.d0,rffu1,u3,uf)
            iu=su*(iu0-(-one)**tpr*iu1)/sqrt(dd)
          else
c If uf is in the forbidden region u2 < uf < u3, issue warning and modify it.
            write(6,*) 'WARNING - Unphysical Cubic Real. Input modified.'
            if(su.eq.1) then
              u0=u3
            else
              u0=u2
            endif
c Try again with valid uf.
            goto 80
          endif
        elseif(abs(dis).lt.1.d-16) then
c This is a cubic case with equal roots. We could use Carlson routines here, but it's faster to use elementary functions.
          u1=-two*sqrt(qq)-cc/dd/3.d0
          u2=-two*sqrt(qq)*cos((two*pi)/3.d0)-cc/dd/3.d0
          u3=u2
          tpr=0
          if(u0.le.u2) then
            ncase=1
            if(uf.gt.u2) then 
              uf=u2
              iu=su*1.d300
            else
              farg=(sqrt(uf-u1)+sqrt(u2-u1))/abs(sqrt(uf-u1)-sqrt(u2-u1))
              sarg=(sqrt(u0-u1)+sqrt(u2-u1))/abs(sqrt(u0-u1)-sqrt(u2-u1))
              iu=su*(log(farg)-log(sarg))/sqrt((u2-u1)*dd)    
            endif           
          elseif(u0.ge.u2) then
            ncase=2
            if(uf.lt.u2) then
              uf=u2
              iu=su*1.d300
            else
              farg=(sqrt(uf-u1)+sqrt(u2-u1))/abs(sqrt(uf-u1)-sqrt(u2-u1))
              sarg=(sqrt(u0-u1)+sqrt(u2-u1))/abs(sqrt(u0-u1)-sqrt(u2-u1))
              iu=-su*(log(farg)-log(sarg))/sqrt((u2-u1)*dd)     
            endif
          endif
        else
c This is a cubic complex case with one real root.
          ncase=3
          aa=-sign(1.d0,rr)*(abs(rr)+sqrt(dis))**third
          if(aa.ne.0.d0) then
            bb=qq/aa
          else
            bb=0.d0
          endif
c Table 1 Row 3
          u1=(aa+bb)-cc/dd/3.d0
          f=-one/dd/u1
          g=f/u1
          if(uf.gt.u0) then
            iu=su*ellcubiccomplex(p,-u1,one,0.d0,0.d0,f,g,one,rffu0,u0,uf)/sqrt(dd)
          elseif(uf.lt.u0) then
            iu=-su*ellcubiccomplex(p,-u1,one,0.d0,0.d0,f,g,one,rffu0,uf,u0)/sqrt(dd)
          else
            iu=0.d0         
          endif 
        endif
        if(q2.eq.0.d0) then
c Find muf in the special case q2=0. In this case there can only be 0 or 1 mu turning points.
          s1=sign(1.d0,mu0)
          if(l2.lt.a*a.and.mu0.ne.0.d0) then
            muplus=s1*sqrt(one-l2/a/a)
c Eq. (43)
            muf=muplus/cosh(abs(a*muplus)*iu-s1*sm*asech(mu0/muplus))
          else
            muf=0.d0
          endif
        elseif(a.eq.0.d0) then
c Find muf in the special case a=0
          a1=sm
          muplus=sqrt(q2/ql2)
          if(mu0.gt.muplus) muplus=mu0
          i1mu=acos(mu0/muplus)/sqrt(ql2)
          i3mu=pi/sqrt(ql2)
          if(sm.eq.1) then
c Eq. (30)
            tpm=int((iu-i1mu)/i3mu)+int((1+sign(1.d0,(iu-i1mu)/i3mu))/2.d0)
          else
c Eq. (31)
            tpm=int((iu+i1mu)/i3mu)
          endif
          a2=sm*(-one)**tpm
c Eq. (35)
          a3=two*int((two*dble(tpm)+3.d0-sm)/4.d0)-one
c Eq. (44) using muminus=-muplus
          muf=-muplus*cos(sqrt(ql2)*(iu-a1*i1mu-a3*i3mu)/a2)
        endif
      elseif(ee.eq.0.d0 .and. dd.eq.0.d0) then
c This is the special case where q2=0 and l=a. In this case, U(u)=1.
        ncase=4
        tpr=0
        iu=su*(uf-u0)
        muf=0.d0
      else
        p(1)=-1
        p(2)=-1
        p(3)=-1
        p(4)=-1
        p(5)=0
c These are the quartic cases. First we find the roots.
        if(firstpt) then
          c(1)=dcmplx(one,0.d0)
          c(2)=dcmplx(0.d0,0.d0)
          c(3)=dcmplx(cc,0.d0)
          c(4)=dcmplx(dd,0.d0)
          c(5)=dcmplx(ee,0.d0)
          call zroots(c,4,root,.true.)
          nreal=0
          do i=1,4
            if(dimag(root(i)).eq.0.d0) nreal=nreal+1
          enddo
          if(nreal.eq.2) then
c This is the quartic complex case with 2 real roots.
            ncase=5
            u1=dble(root(1))
            if(dimag(root(2)).eq.0.d0) then
              u4=dble(root(2))
            else
              u4=dble(root(4))
            endif
          elseif(nreal.eq.0) then
            ncase=6
          else
            u1=dble(root(1))
            u2=dble(root(2))
            u3=dble(root(3))
            u4=dble(root(4))
90          continue
c There are a few cases where all roots are real but unphysical.
            if(u2.gt.uplus.and.u3.gt.uplus) then
              ncase=5
            elseif(u0.le.u2) then
              ncase=7
            elseif(u0.ge.u3) then
              ncase=8
            else
c If u2 < u0 < u3, issue warning and modify it.
              write(6,*) 'WARNING--Unphysical Quartic Real. Inputs modified.'
              if(su.eq.1) then
                u0=0.d0
              else
                u0=uplus
              endif
              goto 90
            endif             
          endif
        endif
        if(ncase.eq.5) then
c Cases with one pair of complex roots
c Table 1 Row 5
          qs=sign(1.d0,q2)
          f=-qs*one/abs(ee)/u1/u4
          g=(u4+u1)/u1/u4*f
          if(uf.gt.u0) then
            iu=su*ellquarticcomplex(p,-u1,one,qs*u4,-qs*one,0.d0,0.d0,f,g,one,rffu0,u0,uf)/sqrt(abs(ee))
          elseif(uf.lt.u0) then
            iu=-su*ellquarticcomplex(p,-u1,one,qs*u4,-qs*one,0.d0,0.d0,f,g,one,rffu0,uf,u0)/sqrt(abs(ee))
          else
            iu=0.d0
          endif
        elseif(ncase.eq.6) then
c This is the quartic complex case with no real roots.
c Table 1 Row 6
          ncase=6
c Solve for h1 from Eq. (21):
          coefs(1)=dcmplx(one,0.d0)
          coefs(2)=dcmplx(-cc/sqrt(ee),0.d0)
          coefs(3)=dcmplx(-one,0.d0)
          coefs(4)=dcmplx(sqrt(ee)*(two*cc/ee-(dd/ee)**2),0.d0)
          coefs(5)=dcmplx(-one,0.d0)
          coefs(6)=dcmplx(-cc/sqrt(ee),0.d0)
          coefs(7)=dcmplx(one,0.d0)
          call zroots(coefs,6,hroots,.true.)
          i=0
          h1=0.d0
10        continue
            i=i+1
            if(dimag(hroots(i)).eq.0.d0) h1=dble(hroots(i))
          if(h1.eq.0.d0) goto 10
c Given h1, we can solve for the other coefficients:
          h2=one/h1
          g1=dd/ee/(h2-h1)
          g2=-g1
          f1=one/sqrt(ee)
          f2=f1
          if(uf.gt.u0) then
            iu=su*elldoublecomplex(p,f1,g1,h1,f2,g2,h2,0.d0,0.d0,rffu0,u0,uf)/sqrt(abs(ee))
          elseif(uf.lt.u0) then
            iu=-su*elldoublecomplex(p,f1,g1,h1,f2,g2,h2,0.d0,0.d0,rffu0,uf,u0)/sqrt(abs(ee))
          else
            iu=0.d0
          endif
        else
c These are the quartic real roots cases
          if(abs(u3-u2).gt.1.d-12) then
            iu1=0.d0
            if(ncase.eq.7) then
c Table 1 Row 7
              if(uf.gt.u2) uf=u2
              if(firstpt.and.(u0.ne.u2)) then 
                iu0=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,0.d0,rffu0,u0,u2)
              elseif(u0.eq.u2) then
                iu0=0.d0
              endif
              if(uf.ne.u2) iu1=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,0.d0,rffu1,uf,u2)
              iu=su*(iu0-(-one)**tpr*iu1)/sqrt(abs(ee))
            elseif(ncase.eq.8) then
c Table 1 Row 8
              if(uf.lt.u3) uf=u3
              if(firstpt.and.(u0.ne.u3)) then 
                iu0=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,0.d0,rffu0,u3,u0)
              elseif(u0.eq.u3) then
                iu0=0.d0
              endif
              if(uf.ne.u3) iu1=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,0.d0,rffu1,u3,uf)
              iu=su*(iu0-(-one)**tpr*iu1)/sqrt(abs(ee))
            endif
          else
c These are the equal roots quartic cases.
            tpr=0
            if(ncase.eq.7) then
              if(uf.gt.u2) uf=u2
              iu0=ellquarticreal(p,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,0.d0,rffu0,u0,uf)
              iu=su*iu0/sqrt(abs(ee))
            elseif(ncase.eq.8) then
              if(uf.lt.u3) uf=u3
              iu0=-ellquarticreal(p,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,0.d0,rffu0,u0,uf)
              iu=su*iu0/sqrt(abs(ee))
            else
              write(6,*) 'ERROR - Unphysical Quartic Real'
              ncase=0
            endif
          endif            
        endif
c Find roots of M(mu) (Eq. (12)).
        yy=-0.5d0*(a*a-ql2+sign(one,a*a-ql2)*sqrt((a*a-ql2)**2+4.d0*q2*a*a))
        if((a*a-ql2).lt.0.d0) then
          mneg=-yy/a/a
          mpos=q2/yy
        else
          mneg=q2/yy
          mpos=-yy/a/a
        endif
        if(mpos.gt.1.d0) mpos=1.d0
        muplus=sqrt(mpos)
        if(mneg.lt.0.d0) then
c This is the symmetric roots case, where the orbit can cross the equatorial plane.
          if(muplus.lt.mu0) then
            muplus=mu0
            mpos=muplus*muplus
          endif
          muminus=-muplus
          if(firstpt) call calcimusym(a,mneg,mpos,mu0,muplus,i1mu,i3mu,rffmu1,rffmu3)
          a1=sm
          if(sm.eq.1) then
            tpm=int((iu-i1mu)/i3mu)+int((1+sign(1.d0,(iu-i1mu)/i3mu))/2.d0)
          else
            tpm=int((iu+i1mu)/i3mu)
          endif
          a2=sm*(-one)**tpm
c Eq. (35)
          a3=two*int((two*dble(tpm)+3.d0-sm)/4.d0)-one
c Eq. (42)
          iarg=abs(a)/a2*(iu-a1*i1mu-a3*i3mu)
          uarg=sqrt(mpos-mneg)*iarg
          m1=-mneg/(mpos-mneg)
          call sncndn(uarg,m1,sn,cn,dn)
          muf=muminus*cn
          if(pht) call calcimusymf(a,mneg,mpos,muf,muplus,i2mu,i3mu,rffmu2)
        else
          muplus=sign(1.d0,mu0)*muplus
          muminus=sign(1.d0,mu0)*sqrt(mneg)
          if(abs(muminus).gt.abs(mu0)) then
            muminus=mu0
            mneg=muminus*muminus
          endif
          if(abs(mu0).gt.abs(muplus)) then
            muplus=mu0
            mpos=muplus*muplus
          endif
c          write(6,*) 'abs: ',abs(mu0)-abs(muplus),mu0-muplus,mneg,mpos
          if(firstpt) call calcimuasym(a,mneg,mpos,mu0,muplus,i1mu,i3mu,rffmu1,rffmu3)
          a1=sm
          if(sm.eq.1) then
c Eq. (30)
            tpm=int((iu-i1mu)/i3mu)+int(abs((1+sign(1.d0,(iu-i1mu)/i3mu))/2.d0))
          else
c Eq. (33)
            tpm=int((iu+i1mu)/i3mu)
          endif
          a2=sm*(-one)**tpm
c Eq. (35)
          a3=two*int((two*dble(tpm)+3.d0-sm)/4.d0)-one
c Eq. (39)
          iarg=abs(a)/a2*(iu-a1*i1mu-a3*i3mu)
          uarg=abs(muplus)*iarg
          m1=mneg/mpos
          call sncndn(uarg,m1,sn,cn,dn)
          muf=muminus/dn
          if(pht) call calcimuasymf(a,mneg,mpos,muf,muplus,i2mu,i3mu,rffmu2)
        endif                
      endif
c      write(6,*) 'geomu debug: ',firstpt,muminus,muplus,mneg,mpos
c      write(6,*) 'geomu debug: ',iu0,iu1,muf
      return
      end

      double precision function ellcubicreal(p,a1,b1,a2,b2,a3,b3,a4,b4,rff,y,x)
*******************************************************************************
*     PURPOSE: Computes \int_y^x dt \Pi_{i=1}^4 (a_i+b_i t)^{p_i/2} for Table 1 Rows 1,2.
*     INPUTS:  Arguments for above integral.
*     OUTPUTS:  Value of integral, and RFF is the RF piece.
*     ROUTINES CALLED: RF,RJ,RC,RD
*     ACCURACY:   Machine.
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit NONE
      double precision one,half,two,three,ellcubic,d12,d13,d14,d24,d34,X1,X2,X3,X4,
     *                 Y1,Y2,Y3,Y4,U1c,U32,U22,W22,U12,Q22,P22,I1c,I3c,r12,r13,r24i,r34i,
     *                 I2c,K2c,a1,b1,a2,b2,a3,b3,a4,b4,y,x,rc,rd,rj,rf,rff
      integer p(5)
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0 )
      ellcubic=0.d0
c (2.1) Carlson (1989)
      d12=a1*b2-a2*b1
      d13=a1*b3-a3*b1
      d14=a1*b4-a4*b1
      d24=a2*b4-a4*b2
      d34=a3*b4-a4*b3
c (2.2) Carlson (1989)
      X1=sqrt(a1+b1*x)
      X2=sqrt(a2+b2*x)
      X3=sqrt(a3+b3*x)
      X4=sqrt(a4+b4*x)
      Y1=sqrt(a1+b1*y)
      Y2=sqrt(a2+b2*y)
      Y3=sqrt(a3+b3*y)
      Y4=sqrt(a4+b4*y)
c (2.3) Carlson (1989)
      U1c=(X1*Y2*Y3+Y1*X2*X3)/(x-y)
      U12=U1c**2
      U22=((X2*Y1*Y3+Y2*X1*X3)/(x-y))**2
      U32=((X3*Y1*Y2+Y3*X1*X2)/(x-y))**2
c (2.4) Carlson (1989)
      W22=U12
      W22=U12-b4*d12*d13/d14
c (2.5) Carlson (1989)
      Q22=(X4*Y4/X1/Y1)**2*W22 
      P22=Q22+b4*d24*d34/d14
c Now, compute the three integrals we need [-1,-1,-1],[-1,-1,-1,-2], and 
c  [-1,-1,-1,-4]:
      if(p(4).eq.0) then
c (2.21) Carlson (1989)
        rff=rf(U32,U22,U12)
        ellcubic=two*rff
      else
c (2.12) Carlson (1989)
        I1c=two*rff
        if(p(4).eq.-2) then
c (2.14) Carlson (1989)
          I3c=two*rc(P22,Q22)-two*b1*d12*d13/three/d14*rj(U32,U22,U12,W22)
c (2.49) Carlson (1989)
          ellcubic=(b4*I3c-b1*I1c)/d14
        else
c (2.1)  Carlson (1989)
          r12=a1/b1-a2/b2
          r13=a1/b1-a3/b3
          r24i=b2*b4/(a2*b4-a4*b2)
          r34i=b3*b4/(a3*b4-a4*b3)
c (2.13) Carlson (1989)
          I2c=two/three*d12*d13*rd(U32,U22,U12)+two*X1*Y1/U1c
c (2.59) & (2.6) Carlson (1989)
          K2c=b2*b3*I2c-two*b4*(X1*X2*X3/X4**2-Y1*Y2*Y3/Y4**2)
c (2.62) Carlson (1989)
          ellcubic=half*b4/d14/d24/d34*K2c
     *     +(b1/d14)**2*(one-half*r12*r13*r24i*r34i)*I1c
        endif
      endif
      ellcubicreal=ellcubic
      return
      end

      double precision function ellcubiccomplex(p,a1,b1,a4,b4,f,g,h,rff,y,x)
*******************************************************************************
*     PURPOSE: Computes \int_y^x dt \Pi_{i=1,4} (a_i+b_i t)^{p_i/2} (f+gt+ht^2)^{p_2/2} for
*              Table 1 Row 3.
*     INPUTS:  Arguments for above integral.
*     OUTPUTS:  Value of integral, and RFF is the RF piece.
*     ROUTINES CALLED: RF,RJ,RC,RD
*     ACCURACY:   Machine.
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      double precision a1,b1,a4,b4,f,g,h,y,x,X1,X4,Y1,Y4,d14,beta1,beta4,a11,c44,
     *                 a142,xi,eta,M2,Lp2,Lm2,I1c,U,U2,Wp2,W2,Q2,P2,rho,I3c,r24xr34,r12xr13,
     *                 N2c,K2c,ellcubic,one,two,four,three,half,six,rc,rd,rj,rf,rff
      integer p(5)
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0, FOUR=4.d0, SIX=6.d0 )
      ellcubic=0.d0
      X1=sqrt(a1+b1*x)
      X4=sqrt(a4+b4*x)
      Y1=sqrt(a1+b1*y)
      Y4=sqrt(a4+b4*y)
      d14=a1*b4-a4*b1
c (2.2) Carlson (1991)
      beta1=g*b1-two*h*a1
      beta4=g*b4-two*h*a4
c (2.3) Carlson (1991)
      a11=sqrt(two*f*b1*b1-two*g*a1*b1+two*h*a1*a1)
      c44=sqrt(two*f*b4*b4-two*g*a4*b4+two*h*a4*a4)
      a142=two*f*b1*b4-g*(a1*b4+a4*b1)+two*h*a1*a4
c (2.4) Carlson (1991)
      xi=sqrt(f+g*x+h*x*x)
      eta=sqrt(f+g*y+h*y*y)
c (3.1) Carlson (1991):
      M2=((X1+Y1)*sqrt((xi+eta)**two-h*(x-y)**two)/(x-y))**two
c (3.2) Carlson (1991):
      Lp2=M2-beta1+sqrt(two*h)*a11
      Lm2=M2-beta1-sqrt(two*h)*a11
      if(p(4).eq.0) then
        rff=rf(M2,Lm2,Lp2)
c (1.2)   Carlson (1991)
        ellcubic=four*rff
      else
c (3.8)  1991
        I1c=four*rff
c (3.3) 1991
        U=(X1*eta+Y1*xi)/(x-y)
        U2=U*U
        Wp2=M2-b1*(a142+a11*c44)/d14
        W2=U2-a11**two*b4/two/d14
c (3.4) 1991
        Q2=(X4*Y4/X1/Y1)**two*W2
        P2=Q2+c44**two*b4/two/d14
c (3.5) 1991
        rho=sqrt(two*h)*a11-beta1
c (3.9) 1991
        if(p(4).eq.-2) then
c (2.49) Carlson (1989)
          I3c=(two*a11/three/c44)*((-four*b1/d14)*(a142+a11*c44)
     *    *rj(M2,Lm2,Lp2,Wp2)-six*rff+three*rc(U2,W2))
     *    +two*rc(P2,Q2)
          ellcubic=(b4*I3c-b1*I1c)/d14
        else
c (2.19) Carlson (1991)
          r24Xr34=half*c44**two/h/b4**two
          r12Xr13=half*a11**two/h/b1**two
c (3.11) Carlson (1991)
          N2c=two/three*sqrt(two*h)/a11*(four*rho*rd(M2,Lm2,Lp2)
     *      -six*rff+three/U)+two/X1/Y1/U
c (2.5) & (3.12) Carlson (1991)
          K2c=half*a11**two*N2c-two*d14*(xi/X1/X4**two-eta/Y1/Y4**two)
c (2.62) Carlson (1989)
          ellcubic=half/d14/(h*b4*r24Xr34)*K2c+(b1/d14)**two*(one-half*r12Xr13/r24Xr34)*I1c
        endif
      endif
      ellcubiccomplex=ellcubic
      return
      end      

      double precision function asech(x)
*******************************************************************************
*     PURPOSE: Computes asech(x)
*     INPUTS:  x>0.
*     OUTPUTS:  asech(x)
*     ROUTINES CALLED: *
*     ACCURACY:   Machine.
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit NONE
      double precision x
      if(x.gt.0.d0) then
        asech=log((1.d0+sqrt(1.d0-x*x))/x)
      else
        asech=0.d0
      endif
      return
      end

      double precision function ellquarticreal(p,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,rff,y,x)
*******************************************************************************
*     PURPOSE: Computes \int_y^x dt \Pi_{i=1}^5 (a_i+b_i t)^{p_i/2} for Table 1 Rows 7,8.
*     INPUTS:  Arguments for above integral.
*     OUTPUTS:  Value of integral, and RFF is the RF piece.
*     ROUTINES CALLED: RF,RJ,RC,RD
*     ACCURACY:   Machine.    
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit NONE
      double precision a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,y,x,rff,
     *                 d12,d13,d14,d24,d34,d15,d25,d35,d45,X1,X2,X3,X4,Y1,Y2,Y3,Y4,
     *                 U122,U132,U142,I1,X52,Y52,W22,Q22,P22,I3,I2,r12,r13,
     *                 r25i,r35i,A111m1m2,one,half,two,three,
     *                 rc,rd,rj,rf
      integer p(5)
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0 )
c (2.1) Carlson (1988)
      d12=a1*b2-a2*b1
      d13=a1*b3-a3*b1
      d14=a1*b4-a4*b1
      d24=a2*b4-a4*b2
      d34=a3*b4-a4*b3
      d15=a1*b5-a5*b1
      d25=a2*b5-a5*b2
      d35=a3*b5-a5*b3
      d45=a4*b5-a5*b4
c (2.2) Carlson (1988)
      X1=sqrt(a1+b1*x)
      X2=sqrt(a2+b2*x)
      X3=sqrt(a3+b3*x)
      X4=sqrt(a4+b4*x) 
      Y1=sqrt(a1+b1*y)
      Y2=sqrt(a2+b2*y)
      Y3=sqrt(a3+b3*y)
      Y4=sqrt(a4+b4*y)
c (2.3) Carlson (1988)
      U122=((X1*X2*Y3*Y4+Y1*Y2*X3*X4)/(x-y))**two
      U132=((X1*X3*Y2*Y4+Y1*Y3*X2*X4)/(x-y))**two
      U142=((X1*X4*Y2*Y3+Y1*Y4*X2*X3)/(X-y))**two
c Now, compute the three integrals we need [-1,-1,-1,-1],[-1,-1,-1,-1,-2], 
c  [-1,-1,-1,-1,-4]:
      if(p(5).eq.0) then
        rff=rf(U122,U132,U142)
c (2.17) Carlson (1988)
        ellquarticreal=two*rff
        return
      else
c (2.13) Carlson (1988)
        I1=two*rff
        X52=a5+b5*x
        Y52=a5+b5*y
c (2.4) Carlson (1988)
        W22=U122-d13*d14*d25/d15
c (2.5) Carlson (1989)
        Q22=X52*Y52/(X1*Y1)**two*W22
        P22=Q22+d25*d35*d45/d15
c (2.15) Carlson (1988)
        if(p(5).eq.-2) then
c (2.35) Carlson (1988)
          I3=two*d12*d13*d14/three/d15*rj(U122,U132,U142,W22)+two*rc(P22,Q22)
          ellquarticreal=(b5*I3-b1*I1)/d15
          return
        else
          I2=two/three*d12*d13*rd(U122,U132,U142)+two*X1*Y1/X4/Y4/sqrt(U142)
c (2.1)  Carlson (1988)
          r12=a1/b1-a2/b2
          r13=a1/b1-a3/b3
          r25i=b2*b5/(a2*b5-a5*b2)
          r35i=b3*b5/(a3*b5-a5*b3)
c (2.48) Carlson (1988)
          A111m1m2=X1*X2*X3/X4/X52-Y1*Y2*Y3/Y4/Y52
          ellquarticreal=half*b5**two*d24*d34/d15/d25/d35/d45*I2
     *   + (b1/d15)**two*(one-half*r12*r13*r25i*r35i)*I1-b5**two/d15/d25/d35*A111m1m2
          return
        endif
      endif
      ellquarticreal=0.d0
      return
      end

      double precision function ellquarticcomplex(p,a1,b1,a4,b4,a5,b5,f,g,h,rff,y,x)
*******************************************************************************
*     PURPOSE: Computes \int_y^x dt \Pi_{i=1,4,5} (a_i+b_i t)^{p_i/2} (f+gt+ht^2)^{p_2/2} 
*              for Table 1 Row 5.
*     INPUTS:  Arguments for above integral.
*     OUTPUTS:  Value of integral, and RFF is the RF piece.
*     ROUTINES CALLED: RF,RJ,RC,RD
*     ACCURACY:   Machine.    
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit NONE
c This routine computes Carlson elliptic integrals for Table 1, Row 5 using
c Carlson (1991)
c JAD 2/20/2009
      double precision a1,b1,a4,b4,a5,b5,f,g,h,y,x,one,half,two,three,four,six,
     *                 X1,X4,Y1,Y4,d14,a11,c44,a142,xi,eta,M2,Lp2,Lm2,I1,
     *                 d15,X52,Y52,c552,c55,a152,U,U2,Wp2,W2,Q2,P2,I3,I2,A111m1m2,
     *                 d45,rc,rd,rj,rf,rff
      integer p(5)
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0, FOUR= 4.D0,
     &            SIX=6.d0 )
c (2.1) Carlson (1991)
      X1=sqrt(a1+b1*x)
      X4=sqrt(a4+b4*x)
      Y1=sqrt(a1+b1*y)
      Y4=sqrt(a4+b4*y)
c (2.3) Carlson (1991)
      a11=sqrt(two*f*b1*b1-two*g*a1*b1+two*h*a1*a1)
      c44=sqrt(two*f*b4*b4-two*g*a4*b4+two*h*a4*a4)
      a142=two*f*b1*b4-g*(a1*b4+a4*b1)+two*h*a1*a4
c (2.4) Carlson (1991)
      xi= sqrt(f+g*x+h*x*x)
      eta=sqrt(f+g*y+h*y*y)
c (2.6) Carlson (1991):
      M2=((X1*Y4+Y1*X4)*sqrt((xi+eta)**two-h*(x-y)**two)/(x-y))**two
c (2.7) Carlson (1991):
      Lp2=M2+a142+a11*c44
      Lm2=max(M2+a142-a11*c44,0.d0)
      if(p(5).eq.0.d0) then
c (1.2) Carlson (1991)
        rff=rf(M2,Lm2,Lp2)
        ellquarticcomplex=four*rff
        return
      else
c (2.14)  1991
        I1=four*rff
c (2.1) 1991
        d14=a1*b4-a4*b1
        d15=a1*b5-a5*b1
        d45=a4*b5-a5*b4
        X52=a5+b5*x
        Y52=a5+b5*y
c (2.3) 1991
        c552=two*f*b5*b5-two*g*a5*b5+two*h*a5*a5
        c55=sqrt(c552)
        a152=two*f*b1*b5-g*(a1*b5+a5*b1)+two*h*a5*a1
c (2.8) 1991
        U=(X1*X4*eta+Y1*Y4*xi)/(x-y)
        U2=U*U
        Wp2=M2+d14*(a152+a11*c55)/d15
        W2=U2-a11**two*d45/two/d15
c (2.9) 1991
        Q2=X52*Y52/(X1*Y1)**two*W2
        P2=Q2+c55**two*d45/two/d15
c (2.15) 1991
        if(p(5).eq.-2) then
c (2.35) Carlson (1988)
          I3=(two*a11/three/c55)*((four*d14/d15)*(a152+a11*c55)*rj(M2,Lm2,Lp2,Wp2)
     *        -six*rff+three*rc(U2,W2))+two*rc(P2,Q2)
          ellquarticcomplex=(b5*I3-b1*I1)/d15
          return
        else
c (2.15) Carlson (1991)
          I2=two*a11/three/c44*(four*(a142+a11*c44)*rd(M2,Lm2,Lp2)-six*rff)
     *         +two*a11/c44/U+two*X1*Y1/X4/Y4/U
c (2.5)  Carlson (1991)
          A111m1m2=X1*xi/X4/X52-Y1*eta/Y4/Y52
          ellquarticcomplex=b5**two/(two*d15*d45)*c44**two/c552*I2+
     *                    b1**two/d15**two*(one-half*b5**two*a11**two/b1**two/c552)*I1-two*b5**two/d15/c552*A111m1m2          
          return
          endif 
        endif
      ellquarticcomplex=0.d0
      return
      end

      double precision function elldoublecomplex(p,f1,g1,h1,f2,g2,h2,a5,b5,rff,y,x)
*******************************************************************************
*     PURPOSE: Computes \int_y^x dt (f_1+g_1t+h_1t^2)^{p_1/2} (f_2+g_2t+h_2t^2)^{p_2/2} (a_5+b_5t)^{p_5/2}
*              for Table 1 Row 6.
*     INPUTS:  Arguments for above integral.
*     OUTPUTS:  Value of integral, and RFF is the RF piece.
*     ROUTINES CALLED: RF,RJ,RC,RD
*     ACCURACY:   Machine.    
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit NONE
c This routine computes Carlson elliptic integrals for Table 1, Row 6
c using Carlson (1992).
c JAD 2/20/2009
      double precision one,half,two,three,four,six,f1,g1,h1,f2,g2,h2,a5,b5,y,x,xi1,xi2,eta1,eta2,
     *                 theta1,theta2,zeta1,zeta2,M,M2,delta122,delta112,delta222,delta,deltap,lp2,lm2,
     *                 deltam,rff,ellquartic,U,U2,alpha15,beta15,alpha25,beta25,lambda,omega2,psi,xi5,eta5,
     *                 gamma1,gamma2,Am111m1,A1111m4,XX,S,mu,T,V2,b2,a2,H,A1111m2,xi1p,B,G,Sigma,
     *                 S2,T2,eta1p,psi2,rc,rd,rj,rf
      integer p(5)
      one=1.d0
      half=0.5d0
      two=2.d0
      three=3.d0
      four=4.d0
      six=6.d0
c (2.1) Carlson (1992)
      xi1=sqrt(f1+g1*x+h1*x**two)
      xi2=sqrt(f2+g2*x+h2*x**two)
      eta1=sqrt(f1+g1*y+h1*y*y)
      eta2=sqrt(f2+g2*y+h2*y*y)
c (2.4) Carlson (1992)
      theta1=two*f1+g1*(x+y)+two*h1*x*y
      theta2=two*f2+g2*(x+y)+two*h2*x*y
c (2.5) Carlson (1992)
      zeta1=sqrt(2*xi1*eta1+theta1)
      zeta2=sqrt(2*xi2*eta2+theta2)
c (2.6) Carlson (1992)
      M=zeta1*zeta2/(x-y)
      M2=M*M
c (2.7) Carlson (1992)
      delta122=two*f1*h2+two*f2*h1-g1*g2
      delta112=four*f1*h1-g1*g1
      delta222=four*f2*h2-g2*g2
      Delta=sqrt(delta122*delta122-delta112*delta222)
c (2.8) Carlson (1992)
      Deltap=delta122+Delta
      Deltam=delta122-Delta
      Lp2=M2+Deltap
      Lm2=M2+Deltam
      if(p(5).eq.0) then
c (2.36) Carlson (1992)
        rff=rf(M2,Lm2,Lp2)
        ellquartic=four*rff
      else
c (2.6) Carlson (1992)
        U=(xi1*eta2+eta1*xi2)/(x-y)
        U2=U*U
c (2.11) Carlson (1992)
        alpha15=two*f1*b5-g1*a5 
        alpha25=two*f2*b5-g2*a5
        beta15=g1*b5-two*h1*a5 
        beta25=g2*b5-two*h2*a5
c (2.12) Carlson (1992)
        gamma1=half*(alpha15*b5-beta15*a5)
        gamma2=half*(alpha25*b5-beta25*a5)
c (2.13) Carlson (1992)
        Lambda=delta112*gamma2/gamma1
        Omega2=M2+Lambda
        psi=half*(alpha15*beta25-alpha25*beta15)
        psi2=psi*psi
c (2.15) Carlson (1992)
        xi5=a5+b5*x
        eta5=a5+b5*y
c (2.16) Carlson (1992)
        Am111m1=one/xi1*xi2-one/eta1*eta2
        A1111m4=xi1*xi2/xi5**two-eta1*eta2/eta5**two
c (2.17) Carlson (1992)
        XX = xi5*eta5*(theta1*half*Am111m1-xi5*eta5*A1111m4)/(x-y)**two
c (2.18) Carlson (1992)
        S=half*(M2+delta122)-U2
        S2=S*S
c (2.19) Carlson (1992)
        mu=gamma1*xi5*eta5/xi1/eta1
        T=mu*S+two*gamma1*gamma2
        T2=T*T
        V2=mu**two*(S2+Lambda*U2)
c (2.20) Carlson (1992)
        b2=Omega2**two*(S2/U2+Lambda)
        a2=b2+Lambda**two*psi2/gamma1/gamma2
c (2.22) Carlson (1992)
        H=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+half*rc(a2,b2))/gamma1**two-XX*rc(T2,V2)
        if (p(5).eq.-2) then
c (2.39) Carlson (1992)
          ellquartic=-two*(b5*H+beta15*rff/gamma1)
        else
          A1111m2=xi1*xi2/xi5-eta1*eta2/eta5
c (2.2) Carlson (1992)
          xi1p=half*(g1+two*h1*x)/xi1
          eta1p=half*(g1+two*h1*y)/eta1
c (2.3) Carlson (1992)
          B=xi1p*xi2-eta1p*eta2
c (2.9) Carlson (1992)
          G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+(delta122*theta1-delta112*theta2)/four/xi1/eta1/U  
c (2.10) Carlson (1992)  
          Sigma=G-Deltap*rff+B
c (2.41) Carlson (1992)
          ellquartic=b5*(beta15/gamma1+beta25/gamma2)*H+beta15**two*rff/gamma1**two+
     *                 b5**two*(Sigma-b5*A1111m2)/gamma1/gamma2
        endif
      endif
      elldoublecomplex=ellquartic
      return
      end
      subroutine geor(u0,uf,mu0,muf,a,l,l2,q2,imu,tpm,tpr,su,sm,ncase,h1,
     &     u1,u2,u3,u4,rffu0,rffu1,rffmu1,rffmu2,rffmu3,iu0,i1mu,i3mu,pht,firstpt)
*******************************************************************************
*     PURPOSE: Computes final radial coordinate corresponding to the final polar angle
*              MUF given initial coordinates U0, MU0 and constants of the motion.
*
*     INPUTS:   U0 -- Starting u value. If U0=0, DTI and LAMBDAI values will blow up.
*               MU0 -- Starting mu=cos(theta) value.
*               MUF -- Final mu value.
*               A -- Black hole spin, on interval [0,1).
*               L -- Dimensionless z-component of angular momentum.
*               L2 -- L*L
*               Q2 -- Dimensionless Carter's constant.
*               TPR -- Number of u turning points reached between U0 and UF.
*               SU -- Initial du/dlambda, =1 (-1) for ingoing (outgoing) rays.
*               SM -- Initital dmu/dlambda.
*               PHT -- Boolean variable. If .TRUE., RFFU1 is computed for use in
*                      SUBROUTINE GEOPHITIME.
*               FIRSTPT -- Boolean variable. If .TRUE., roots of U(u) U1,2,3,4 and NCASE are computed as well as H1
*                           if NCASE=6
*     OUTPUTS:  UF -- Final inverse radius.
*               IMU -- Value of IMU integral between MU0 and MUF.
*               TPM -- Number of mu turning points reached between MU0 and MUF.
*               NCASE -- Case number corresponding to Table 1. Output if FIRSTPT=.TRUE., input otherwise.
*               H1 -- Value of h1 from Eq. (21) for given constants of motion if NCASE=6. Output if FIRSTPT=.TRUE.,
*                     input otherwise.
*               U1,U2,U3,U4 -- Increasing (real) roots. Output if FIRSTPT=.TRUE., input otherwise.
*               RFFU0 -- Value of RF relevant for U0 piece of IU integral. Computed if FIRSTPT=.TRUE., input otherwise.
*               RFFU1 -- Value of RF relevant for UF piece of IU integral. Computed if PHT=.TRUE.
*               RFFMU1 -- Value of RF relevant for MU0 piece of IMU integral. Computed if FIRSTPT=.TRUE., input
*                         otherwise.
*               RFFMU2 -- Value of RF relevant for MUF piece of IMU integral. Computed if PHT=.TRUE.
*               RFFMU3 -- Value of RF relevant for turning point piece of IMU integral. Computed if FIRSTPT=.TRUE.,
*                         input otherwise.
*               IU0 -- Value of IU integral between U0 and relevant turning point if one exists. Computed if
*                      u turning point present and FIRSTPT=.TRUE., input otherwise. Ignored if no turning point
*                      is present.
*               I1MU -- Value of IMU integral between MU0 and MUPLUS. Computed if FIRSTPT=.TRUE., input otherwise.
*               I3MU -- Value of IMU integral between MUMINUS and MUPLUS. Computed if FIRSTPT=.TRUE., input otherwise.
*     ROUTINES CALLED: SNCNDN, ZROOTS, ASECH, CALCIMUSYM, CALCIMUSYMF, CALCIMUASYM, CALCIMUASYMF, ELLCUBICREAL, 
*                       ELLCUBICCOMPLEX, ELLQUARTICREAL, ELLQUARTICCOMPLEX, ELLDOUBLECOMPLEX                   
*     ACCURACY:   Machine.
*     REMARKS: Based on Dexter & Agol (2009), and labeled equations refer to that paper unless noted
*              otherwise.   
*     AUTHOR:     Dexter & Agol (2009)
*     DATE WRITTEN:  4 Mar 2009
*     REVISIONS: ***********************************************************************
      implicit none
c Given a starting position (u0,mu0) where u0=1/r0, mu0=cos(theta0),
c and the final polar angle, muf, this subroutine calculates the final
c radius, uf, for a geodesic with parameters (l,q2)
C JAD 2/20/2009
      integer i,nreal,parr(5)
      double precision a,aa,a1,a2,a3,c1,c2,c3,c4,c5,cn,dn,dis,rr,i1mu,i3mu,dummy,
     *      f,half,iu,qs,h1,h2,f1,f2,g1,g2,l,l2,m1,mneg,mpos,mu0,muf,muplus,
     *      one,pi,pi2,q2,ql2,s1,sn,sm,su,theta,third,two,bb,sarg,rffu0,rffu1,
     *      u0,uf,u1,u2,u3,u4,g,cc,dd,ee,iu0,iu1,ellcubicreal,ellcubiccomplex,ellquarticreal,
     *      ellquarticcomplex,elldoublecomplex,asech,qq,uplus,yy,
     *      i2mu,imu,jarg,cd2,dc2,sc,three,m,n,n2,p,r,ua,ub,mn2,pr2,temp,
     *      rffmu1,rffmu2,rffmu3,iut
      integer tpm,tpr,ncase
      logical pht,firstpt
      complex*16 c(5),root(4),coefs(7),hroots(6)
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5D0, THIRD = 0.3333333333333333D0, THREE=3.D0 )
      pi=acos(-one)
      pi2=two*pi
      uplus=one/(one+sqrt(one-a*a))
      ql2=q2+l2
c  Eq. (22): Coefficients of quartic in U(u)
      cc=a*a-q2-l2
      dd=two*((a-l)**2+q2)
      ee=-a*a*q2
      if(q2.eq.0.d0) then
c Calculate imu in the special case q2=0. In this case there can only be 0 or 1 mu turning points.
        if(l2.ge.a*a.or.mu0.eq.0.d0) then
          imu=0.d0
          uf=-1.d0
          ncase=0
          return
        else
          s1=sign(1.d0,mu0)
          a1=s1*sm
          a2=s1*sm*(-one)**(tpm+1)
          muplus=s1*sqrt(one-l2/a/a)
          i1mu=one/abs(a*muplus)*asech(mu0/muplus)
          i2mu=one/abs(a*muplus)*asech(muf/muplus)
          imu=a1*i1mu+a2*i2mu
        endif
      elseif(a.eq.0.d0) then
c Calculate imu in the special case a=0
        a1=sm
        muplus=sqrt(q2/ql2)
        if(mu0.gt.muplus) muplus=mu0
        i1mu=(half*pi-asin(mu0/muplus))/sqrt(ql2)
        i3mu=pi/sqrt(ql2)
        a2=sm*(-one)**tpm
        a3=two*int((two*dble(tpm)+3.d0-sm)/4.d0)-one
        i2mu=(half*pi+asin(muf/muplus))/sqrt(ql2)
        imu=a1*i1mu+a2*i2mu+a3*i3mu
      endif
c  Determine if U is cubic or quartic and find roots.
      if((ee.eq.0.d0) .and. (dd.ne.0.d0)) then
        parr(1)=-1
        parr(2)=-1
        parr(3)=-1
        parr(4)=0
        parr(5)=0
        qq=cc*cc/dd/dd/9.d0
        rr=(two*cc**3/dd**3+27.d0/dd)/54.d0
        dis=rr*rr-qq**3
        if(dis.lt.-1.d-16) then
c These are the cubic real roots cases with u1<0<u2<=u3
          theta=acos(rr/qq**1.5d0)
          u1=-two*sqrt(qq)*cos(theta/3.d0)-cc/dd/3.d0
          u2=-two*sqrt(qq)*cos((theta-two*pi)/3.d0)-cc/dd/3.d0
          u3=-two*sqrt(qq)*cos((theta+two*pi)/3.d0)-cc/dd/3.d0
80        continue
          if(u0.le.u2) then
            ncase=1
            if(firstpt.and.(u0.ne.u2)) then
              iu0=ellcubicreal(parr,-u1,one,u2,-one,u3,-one,0.d0,0.d0,rffu0,u0,u2)
            elseif(u0.eq.u2) then
              iu0=0.d0
            endif
c Table 2 Row 1
            m1=(u3-u2)/(u3-u1)
            jarg=sqrt((u3-u1)*dd)*half*(imu-su*iu0/sqrt(dd))
            call sncndn(jarg,m1,sn,cn,dn)
            cd2=cn**2/dn**2
            uf=u1+(u2-u1)*cd2
            tpr=(sign(1.d0,su*(imu-su*iu0/sqrt(dd)))+1.d0)/2.d0
            if(pht) iu1=ellcubicreal(parr,-u1,one,u2,-one,u3,-one,0.d0,0.d0,rffu1,uf,u2)/sqrt(dd)
          elseif(u0.ge.u3) then
            ncase=2
            if(firstpt.and.(u0.ne.u3)) then
              iu0=ellcubicreal(parr,-u1,one,-u2,one,-u3,one,0.d0,0.d0,rffu0,u3,u0)
            elseif(u0.ne.u3) then
              iu0=0.d0
            endif
c Table 2 Row 2
            m1=(u3-u2)/(u3-u1)
            jarg=sqrt((u3-u1)*dd)*half*(imu+su*iu0/sqrt(dd))
            call sncndn(jarg,m1,sn,cn,dn)
            dc2=dn**2/cn**2
            dummy=uf
            uf=u1+(u3-u1)*dc2
            tpr=(-sign(1.d0,su*(imu+iu0/sqrt(dd)))+1.d0)/2.d0
            if(pht) iu1=ellcubicreal(parr,-u1,one,-u2,one,-u3,one,0.d0,0.d0,rffu1,u3,uf)/sqrt(dd)
          else
            write(6,*) 'WARNING - Unphysical Cubic Real. Input modified.'
            if(su.eq.1) then
              u0=u3
            else
              u0=u2
            endif
            goto 80
          endif
        elseif(abs(dis).lt.1.d-16) then
          tpr=0
c This is a cubic case with equal roots.
          u1=-two*sqrt(qq)-cc/dd/3.d0
          u2=-two*sqrt(qq)*cos((two*pi)/3.d0)-cc/dd/3.d0
          u3=u2
          tpr=0
          if(u0.le.u2) then
            ncase=1
            if(uf.gt.u2) then 
              uf=u2
              iu=su*1.d300
            else
              sarg=sqrt((u0-u1)/(u2-u1))
              jarg=sqrt((u2-u1)*dd)*half*su*imu+half*log((one+sarg)/(one-sarg))
              uf=u1+(u2-u1)*tanh(jarg)**2
            endif           
          elseif(u0.ge.u2) then
            ncase=2
            if(uf.lt.u2) then
              uf=u2
              iu=su*1.d300
            else
              sarg=sqrt((u2-u1)/(u0-u1))
              jarg=-sqrt((u2-u1)*dd)*half*su*imu+half*log((one+sarg)/(one-sarg))
              uf=u1+(u2-u1)/tanh(jarg)**2    
            endif
          endif
        else
c This is a cubic complex case with one real root.
          ncase=3
c Table 2 Row 3
          tpr=0
          aa=-sign(1.d0,rr)*(abs(rr)+sqrt(dis))**third
          if(aa.ne.0.d0) then
            bb=qq/aa
          else
            bb=0.d0
          endif
          u1=(aa+bb)-cc/dd/3.d0
          f=-one/dd/u1
          g=f/u1
c Make sure there is a valid solution.
          if(su.gt.0.d0) then 
            iut=ellcubiccomplex(parr,-u1,one,0.d0,0.d0,f,g,one,dummy,u0,uplus)/sqrt(dd)
          else
            iut=ellcubiccomplex(parr,-u1,one,0.d0,0.d0,f,g,one,dummy,0.d0,u0)/sqrt(dd)
          endif
          if(imu.gt.iut) then
            uf=-1.d0
            ncase=0
            tpr=0
            return
          endif
          m=-g/2.d0
          if(firstpt) iu0=su*ellcubiccomplex(parr,-u1,one,0.d0,0.d0,f,g,one,dummy,u1,u0)/sqrt(dd)
          c3=-one
          if(a.ne.0.d0.or.l.ne.0.d0) c3=(a+l)/(a-l)
          c2=sqrt(u1*(three*u1+c3))
          c1=sqrt(c2*dd)
          m1=half+(6.d0*u1+c3)/(8.d0*c2)
          jarg=c1*(imu+iu0)
          call sncndn(jarg,m1,sn,cn,dn)
          uf=(c2+u1-(c2-u1)*cn)/(one+cn)
          if(pht) iu=ellcubiccomplex(parr,-u1,one,0.d0,0.d0,f,g,one,rffu0,u0,uf)/sqrt(dd)
        endif
      elseif(ee.eq.0.d0 .and. dd.eq.0.d0) then
c This is the special case where q2=0 and l=a and mu=0 at all times, so we can't invert.
        ncase=4
        uf=-1.D0
        tpr=0
      else
c Find roots of M(mu) in biquadratic case.
        yy=-0.5d0*(a*a-ql2+sign(one,a*a-ql2)*sqrt((a*a-ql2)**2+4.d0*q2*a*a))
        if((a*a-ql2).lt.0.d0) then
          mneg=-yy/a/a
          mpos=q2/yy
        else
          mneg=q2/yy
          mpos=-yy/a/a
        endif
        muplus=sqrt(mpos)
c NOTE: This formula uses a slightly different prescription for splitting up the mu integral.
c       All integrals are with respect to muplus.
        a1=sm
        a2=sm*(-1.d0)**(tpm+1)
        a3=2.d0*int((2*tpm-sm+1)/4.d0)
        if(mneg.lt.0.d0) then
c This is the symmetric roots case, where the orbit can cross the equatorial plane.
          if(mu0.gt.muplus) then
            muplus=mu0
            mpos=muplus*muplus
          endif
c Compute integrals involving u0 and turning points once per geodesic:
          if(firstpt) call calcimusym(a,mneg,mpos,mu0,muplus,i1mu,i3mu,rffmu1,rffmu3)
          call calcimusymf(a,mneg,mpos,muf,muplus,i2mu,i3mu,rffmu2)
          imu=a1*i1mu+a2*i2mu+a3*i3mu
        else
c This is the asymmetric roots case.
          if(abs(muf).lt.sqrt(mneg)) then
            uf=-1.d0
            ncase=0
            tpr=0
            return
          else
            if(sign(1.d0,mu0).eq.-1) muplus=-muplus
            if(abs(muplus).lt.abs(mu0)) then
              muplus=mu0
              mpos=muplus*muplus
            endif
            mneg=min(mu0*mu0,mneg)
            if(firstpt) call calcimuasym(a,mneg,mpos,mu0,muplus,i1mu,i3mu,rffmu1,rffmu3)
            call calcimuasymf(a,mneg,mpos,muf,muplus,i2mu,i3mu,rffmu2)
          endif
          imu=a1*i1mu+a2*i2mu+a3*i3mu
        endif                
c These are the quartic cases. First we find the roots.
        parr(1)=-1
        parr(2)=-1
        parr(3)=-1
        parr(4)=-1
        parr(5)=0
        if(firstpt) then
          c(1)=dcmplx(one,0.d0)
          c(2)=dcmplx(0.d0,0.d0)
          c(3)=dcmplx(cc,0.d0)
          c(4)=dcmplx(dd,0.d0)
          c(5)=dcmplx(ee,0.d0)
          call zroots(c,4,root,.true.)
          nreal=0
          do i=1,4
            if(dimag(root(i)).eq.0.d0) nreal=nreal+1
          enddo
          if(nreal.eq.2) then
c This is the quartic complex case with 2 real roots.
            ncase=5
            tpr=0
            u1=dble(root(1))
            if(dimag(root(2)).eq.0.d0) then
              u4=dble(root(2))
            else
              u4=dble(root(4))
            endif
          elseif(nreal.eq.0) then
            ncase=6
            tpr=0
          else
            u1=dble(root(1))
            u2=dble(root(2))
            u3=dble(root(3))
            u4=dble(root(4))
90          continue
            if(u2.gt.uplus.and.u3.gt.uplus) then
              ncase=5
            elseif(u0.le.u2) then
              ncase=7
            elseif(u0.ge.u3) then
              ncase=8
            else
              write(6,*) 'WARNING--Unphysical Quartic Real. 
     &                      Inputs modified.'
              if(su.eq.1) then
                u0=u3
              else
                u0=u2
              endif
              goto 90
            endif             
          endif
        endif
        if(ncase.eq.5) then
c Table 2 Row 5
          qs=sign(1.d0,q2)
          f=-qs*one/abs(ee)/u1/u4
          g=(u4+u1)/u1/u4*f
c Make sure there is a valid solution.
          if(su.gt.0.d0) then
            iut=ellquarticcomplex(parr,-u1,one,qs*u4,-qs*one,0.d0,0.d0,f,g,one,dummy,u0,uplus)/sqrt(abs(ee))
          else
            iut=ellquarticcomplex(parr,-u1,one,qs*u4,-qs*one,0.d0,0.d0,f,g,one,dummy,0.d0,u0)/sqrt(abs(ee))
          endif
          if(imu.gt.iut) then
            uf=-1.d0
            ncase=0
            tpr=0
            return
          endif
          if(qs.eq.one) then
            ua=u4
            ub=u1
          else
            ua=u1
            ub=u4
          endif
          m=-g*half
          n2=f-g**2/4.d0
          c4=sqrt((m-u4)**2+n2)
          c5=sqrt((m-u1)**2+n2)
          c1=sqrt(abs(ee*c4*c5))
          if(firstpt) iu0=su*ellquarticcomplex(parr,-u1,one,qs*u4,-qs*one,0.d0,0.d0,f,g,one,dummy,ub,u0)/sqrt(abs(ee))
          m1=qs*((c4+qs*c5)**2-(u4-u1)**2)/(4.d0*c4*c5)
          jarg=c1*(imu+iu0)
          call sncndn(jarg,m1,sn,cn,dn)
          uf=(u4*c5+qs*u1*c4-(qs*u4*c5-u1*c4)*cn)/((c4-qs*c5)*cn+(qs*c4+c5))
          if(pht) iu=ellquarticcomplex(parr,-u1,one,qs*u4,-qs*one,0.d0,0.d0,f,g,one,rffu0,u0,uf)/sqrt(abs(ee))
        elseif(ncase.eq.6) then
c This is the quartic complex case with no real roots. First we need to find the real arguments f,g,h.
          ncase=6
c Table 2 Row 6
          coefs(1)=dcmplx(one,0.d0)
          coefs(2)=dcmplx(-cc/sqrt(ee),0.d0)
          coefs(3)=dcmplx(-one,0.d0)
          coefs(4)=dcmplx(sqrt(ee)*(two*cc/ee-(dd/ee)**2),0.d0)
          coefs(5)=dcmplx(-one,0.d0)
          coefs(6)=dcmplx(-cc/sqrt(ee),0.d0)
          coefs(7)=dcmplx(one,0.d0)
          call zroots(coefs,6,hroots,.true.)
          i=0
          h1=0.d0
10        continue
            i=i+1
            if(dimag(hroots(i)).eq.0.d0) h1=dble(hroots(i))
          if(h1.eq.0.d0) goto 10
          h2=one/h1
          g1=dd/ee/(h2-h1)
          g2=-g1
          f1=one/sqrt(ee)
          f2=f1
c Make sure there is a valid solution.
          if(su.gt.0.d0) then
            iut=elldoublecomplex(parr,f1,g1,h1,f2,g2,h2,0.d0,0.d0,dummy,u0,uplus)/sqrt(abs(ee))
          else
            iut=elldoublecomplex(parr,f1,g1,h1,f2,g2,h2,0.d0,0.d0,dummy,0.d0,u0)/sqrt(abs(ee))
          endif
          if(imu.gt.iut) then
            uf=-1.d0
            ncase=0
            tpr=0
            return
          endif
C         Next we want to write the real and imaginary parts of the roots m+/-in, p+/-ir in terms
C         of real quantities.
          coefs(1)=ee**(-three)
          coefs(2)=-cc/ee**three
          coefs(3)=-ee**(-two)
          coefs(4)=-ee**(-two)*(dd**two/ee-two*cc)
          coefs(5)=-one/ee
          coefs(6)=-cc/ee
          coefs(7)=one
          call zroots(coefs,6,hroots,.true.)
          i=0
          mn2=0.d0
20        continue
            i=i+1
            if(dimag(hroots(i)).eq.0.d0) mn2=dble(hroots(i))
          if(mn2.eq.0.d0) goto 20
          p=dd/(two*ee**2*(mn2**2-one/ee))
          m=-half*dd/ee-p
          if(m.lt.p) then
            temp=p
            p=m
            m=temp
          endif
          pr2=one/(mn2*ee)
          n=sqrt(mn2-m*m)
          r=sqrt(pr2-p*p)
          c4=sqrt((m-p)**2+(n+r)**2)
          c5=sqrt((m-p)**2+(n-r)**2)
          c1=(c4+c5)*half*sqrt(abs(ee))
          c2=sqrt((4.d0*n*n-(c4-c5)**2)/((c4+c5)**2-4.d0*n*n))
          c3=m+c2*n
          if(firstpt) iu0=sign(1.d0,(u0-c3))*elldoublecomplex(parr,f1,g1,h1,f2,g2,h2,0.d0,0.d0,dummy,c3,u0)/sqrt(abs(ee))
          m1=((c4-c5)/(c4+c5))**2
          jarg=c1*(su*imu+iu0)
          call sncndn(jarg,m1,sn,cn,dn)
          sc=sn/cn
          uf=c3+(n*(one+c2**2)*sc)/(one-c2*sc)
          if(pht) iu=elldoublecomplex(parr,f1,g1,h1,f2,g2,h2,0.d0,0.d0,rffu0,u0,uf)/sqrt(-ee)
        else
c These are the quartic real roots cases
          if(ncase.eq.7) then
c Table 2 Row 7
            if(firstpt.and.(u0.ne.u2)) 
     &       iu0=ellquarticreal(parr,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,0.d0,rffu0,u0,u2)
            jarg=sqrt(abs(ee)*(u3-u1)*(u4-u2))*half*(imu-su*iu0/sqrt(-ee))
            m1=(u4-u1)*(u3-u2)/((u4-u2)*(u3-u1))
            call sncndn(jarg,m1,sn,cn,dn)
            dummy=uf
            uf=((u2-u1)*u3*sn**2-u2*(u3-u1))/((u2-u1)*sn**2-(u3-u1))
            tpr=(sign(1.d0,su*(imu-su*iu0/sqrt(abs(ee))))+1.d0)/2.d0
            if(pht) iu1=ellquarticreal(parr,-u1,one,u2,-one,u3,-one,u4,-one,0.d0,0.d0,rffu1,uf,u2)/sqrt(-ee)
          elseif(ncase.eq.8) then
c Table 2 Row 8
            if(firstpt.and.(u0.ne.u3)) 
     &      iu0=ellquarticreal(parr,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,0.d0,rffu0,u3,u0)
            jarg=sqrt(abs(ee)*(u3-u1)*(u4-u2))*half*(imu+su*iu0/sqrt(-ee))
            m1=(u4-u1)*(u3-u2)/((u4-u2)*(u3-u1))
            call sncndn(jarg,m1,sn,cn,dn)
            uf=((u4-u3)*u2*sn**2-u3*(u4-u2))/((u4-u3)*sn**2-(u4-u2)) 
            tpr=(-sign(1.d0,su*(imu+su*iu0/sqrt(-ee)))+1.d0)/2.d0
            if(pht) iu1=ellquarticreal(parr,-u1,one,-u2,one,-u3,one,u4,-one,0.d0,0.d0,rffu1,u3,uf)/sqrt(-ee)
          endif
        endif
      endif
      return
      end
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
