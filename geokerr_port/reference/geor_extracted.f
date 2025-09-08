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
