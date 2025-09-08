      subroutine geomu(u0,uf,mu0,muf,a,l,l2,q2,iu,tpm,tpr,su,sm,ncase,h1,
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
