

c-----------------------------------------------------------------------
      subroutine SLSINF(idist,itype,zlv,zrv,
     &f11,f12,f22,nrows,ifault,irow)
CDEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'slsinf_' :: SLSINF
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension zlv(1),zrv(1),f11(1),f12(1),f22(1)
      do 22 irow=1,nrows
      ifault=111
cdddd      call intpt('itype0',6,itype,1)
cdddd      call intpt('idist0',6,idist,1)
cdddd      call intpt('ifault0',7,ifault,1)
cdddd      call intpt('nrows0',6,nrows,1)
cdddd      call dblept('zlv0',4,zlv,nrows)
cdddd      call dblept('fout0',5,fout,nrows*3)
      call lsinf(idist,itype,zlv(irow),zrv(irow),f11(irow),
     &f12(irow),f22(irow),ifault)
cdddd      call intpt('idist1',6,idist,1)
cdddd      call intpt('ifault1',7,ifault,1)
cdddd      call intpt('nrows1',6,nrows,1)
cdddd      call dblept('zlv1',4,zlv,nrows)
cdddd      call dblept('fout1',5,fout,nrows*3)
cwww      write(6,433)ifault,idist,itype,zlv(irow),
cwww     &fout(irow,1),fout(irow,2),fout(irow,3)
cwww433   format(' if,idist,it,z,fout=',3i5,4g12.4)
      if(ifault.ge.4)return
22    continue
      return
      end

c-----------------------------------------------------------------------
      subroutine lsinf(idist,itype,zl,zr,f11,f12,f22,ifault)
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
c     computes:  fisher information matrix elements for time (type i)
c                or failure (type ii) censored units from the smallest
c                extreme value (sev), largest extreme value (lev),
c                normal, or logistic distribution
c
c    parameters:
c
c     idist - integer - input:  idist=1 if the distribution is sev
c                                    =2 if the distribution is lev
c                                    =3 if the distribution is normal
c                                    =4 if the distribution is logistic
c
c     itype - integer - input:  itype=1 for no censoring
c                                    =2 right censoring
c                                    =3 left censoring
c                                    =4 interval censoring
c
c        zl - real -    input:  standardized left censoring point. zl is
c                               defined by zl=q(pl), where pl is the
c                               (expected) proportion of left censored
c                               units and q(p) is the pth quantile of
c                               the standardized distribution specified
c                               by idist.  zl is not referenced when
c                               itype=2
c
c        zr - real -    input:  standardized right censoring point.
c                               zr=q(1-pr), where pr is the (expected)
c                               proportion of right censored units.
c                               zr is not referenced when itype=3
c
c       f11 - real -    output: c*(entry(1,1) of the fisher matrix)
c       f12 - real -    output: c*(entry(1,2) of the fisher matrix)
c       f22 - real -    output: c*(entry(2,2) of the fisher matrix)
c                               where c=sigma*sigma/n
c
c    ifault - integer - output: ifault=0 indicates successful completion
c                                     =1 indicates idist is other than
c                                        1, 2, 3, 4
c                                     =2 indicates itype is other than
c                                        1, 2, 3, or 4
c                                     =3 indicates itype=4 and zr.lt.zl
c                                     =4 indicates that the series
c                                        expansions in sevint failed to
c                                        converge
c                                     =5 indicates the euler's transfor-
c                                        mation in logint failed to
c                                        converge
c
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      double precision  big, eta, f11, f12, f22, one, small, thet0l,
     * thet0r, thet1l, thet1r, thet2l, thet2r, zero, zl, zlow, zr,
     * zupp
c
      data  big/1.0d10/, one/1.0d0/, small/-1.0d10/, zero/0.0d0/
c
c       check for illegal idist
c
cwww      write(6,555)idist
cwww555   format(' idist at top=',i5)
cwwwcdddd      call intpt('idist2',6,idist,1)
cwwwcdddd      call intpt('ifault2',7,ifault,1)
      if (idist.lt.1.or.idist.gt.4) goto 60
c
c       check for illegal itype
c
      if (itype.lt.1.or.itype.gt.4) goto 70
c
c       set censoring bounds for itype = 1, 2, 3, 4
c       corresponding to complete, right, left, and
c       interval cases, respectively
c
      goto (10,20,30,40), itype
c
c       uncensored case
c
   10 zupp=big
      zlow=small
      goto 50
c
c       right censored at zr
c
   20 zupp=zr
      zlow=small
      goto 50
c
c       left censored at zl
c
   30 zupp=big
      zlow=zl
      goto 50
c
c       interval censored - first check that the left
c       censoring bound is not greater than the right bound
c
   40 if (zl.gt.zr) goto 80
      zupp=zr
      zlow=zl
c
c       fisher matrix elements
c
   50 call lsint(idist,zupp,thet0r,thet1r,thet2r,eta,ifault)
      if (ifault.eq.4.or.ifault.eq.5) return
      call lsint(idist,zlow,thet0l,thet1l,thet2l,eta,ifault)
      if (ifault.eq.4.or.ifault.eq.5) return
      f11=thet0r-thet0l+eta
      f12=thet1r-thet1l+zlow*eta
      f22=thet2r-thet2l+zlow*zlow*eta
      return
c
c       set ifault
c
   60 ifault=1
cwwwcdddd      call intpt('idist3',6,idist,1)
cwwwcdddd      call intpt('ifault3',7,ifault,1)
      return
   70 ifault=2
      return
   80 ifault=3
      return
      end

c=======================================================================
c
      subroutine lsint(idist,z,theta0,theta1,theta2,eta,ifault)
c
c-----------------------------------------------------------------------
c
c       routes the computation of theta0, theta1, theta2, and eta to
c       sevint, levint, norint, or logint according  with the value of
c       idist
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      double precision eta, theta0, theta1, theta2, z
      go to (10,20,30,40),idist
   10 call sevint(z,theta0,theta1,theta2,eta,ifault)
      return
   20 call levint(z,theta0,theta1,theta2,eta,ifault)
      return
   30 call norint(z,theta0,theta1,theta2,eta,ifault)
      return
   40 call logint(z,theta0,theta1,theta2,eta,ifault)
      return
      end

c=======================================================================
c
      subroutine sevint(z,theta0,theta1,theta2,eta,ifault)
c
c-----------------------------------------------------------------------
c
c     computes: quantities needed to obtain the elements of the fisher
c               information matrix from a smallest extreme value
c               distribution and censored data
c
c       theta0, theta1, and theta2 are the integrals over (-infinity,z)
c       of the functions: g0(x)=h(x)*h(x)*g(x), g1(x)=(1+x*h(x))*g(x),
c       g2(x)=g(x)*(1+x*h(x))**2, where h(x)=1.  g(x) and bg(x) are
c       pdf and the cdf for a smallest extreme value distribution
c
c       theta0 has a closed form formula. g1(x) and g2(x) are integrated
c       by using power series expansions when z.le.1, and by power
c       series expansions through z=1 plus a gaussian quadrature from
c       1 to z when z.gt.1
c
c       eta=g(x)*g(x)/((1-bg(x))*bg(x))
c
c       ifault= 4 if the series expansions to compute s1 and s2 do not
c               converge within tol of the true values after including
c               jmax terms in the expansion.  tol and jmax are set to 25
c               and 10**(-11) in the data statements
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      double precision  asymp1, asymp2, bgz, cdfsev, denom, eta, ez,
     * factor, gz, gzs, g1, g1xl, g1xr, g2xl, g2xr, half, one, pdfsev,
     * sum1, sum2, sz, s1, s2, theta0, theta1, theta2, thet11, thet21,
     * tol, two, x, x1, x2, x3, x4, x5, z, zero, zmax, zmin
c
      double precision  p(8), w(8)
c
c       p and w are constants used in the gaussian quadratures
c
      data p(1)/0.4947004674958250d0 /, p(2)/0.4722875115366163d0 /,
     *     p(3)/0.4328156011939159d0 /, p(4)/0.3777022041775015d0 /,
     *     p(5)/0.3089381222013219d0 /, p(6)/0.2290083888286137d0 /,
     *     p(7)/0.1408017753896295d0 /, p(8)/0.4750625491881872d-1/
      data w(1)/0.1357622970587705d-1/, w(2)/0.3112676196932395d-1/,
     *     w(3)/0.4757925584124639d-1/, w(4)/0.6231448562776694d-1/,
     *     w(5)/0.7479799440828837d-1/, w(6)/0.8457825969750127d-1/,
     *     w(7)/0.9130170752246179d-1/, w(8)/0.9472530522753425d-1/
c
c       asymp1, asymp2, thet11, and thet21 are the integrals
c       for g1 and g2 when z=+infinity and z=1, respectively
c
      data asymp1/0.4227843350984671d0/, asymp2/0.1823680660852879d1/,
     *     half/0.5d0/, jmax/25/, one/1.0d0/, tol/1.0d-11/,
     *     thet11/0.2720757938345342d0/, thet21/0.1475933122158450d1/,
     *     two/2.0d0/, zero/0.0d0/, zmax/3.6d0/, zmin/-34.0d0/
c
c       smallest extreme value pdf and cdf
c
      pdfsev(x)=dexp(x-dexp(x))
      cdfsev(x)=one-dexp(-dexp(x))
c
c       function needed in the gaussian quadrature
c
      g1(x)=(one+x)*dexp(x-dexp(x))
c
c       initial values
c
      ifault=0
      eta=zero
c
c       check for extreme z value
c
      if (z.ge.zmax) goto 50
      if (z.le.zmin) goto 60
c
c       compute theta0=integral of g0(x)=dexp(x-dexp(x)) and eta
c       over (-infinity,z)
c
      ez=dexp(z)
      sz=dexp(-ez)
      gz=pdfsev(z)
      gzs=gz*gz
      bgz=cdfsev(z)
      theta0=bgz
      denom=bgz*sz
      if (denom.gt.zero) eta=gzs/denom
c
c       computation of theta1 and theta2
c
c       select integration by power series expansions or
c       by power series expansions plus a gaussian quadrature
c
      if (z.lt.one) then
c
c       z is less than 1 - integration by power series expansions
c
            factor=-one
            s1=zero
            s2=zero
            do 10 j=1,jmax
c
c       terms to evaluate theta1=integral of g1(x)=(1+x)g(x)
c            over (-infinity,z)
c
            x5=dfloat(j)
            factor=-factor*ez/x5
            x1=z-one/x5
            x3=factor*x1
            s1=s1+x3
c
c       terms to evaluate theta2=integral of g2(x)=(1+x)g1(x)
c                over (-infinity,z)
c
            x2=x1*x1+one/(x5*x5)
            x4=factor*x2
            s2=s2+x4
c
c       tests for convergence of the series -  the summations stop when
c       the absolute value of the last added term is smaller than tol
c       for both series - a fault is declared if convergence is not
c       reached in a maximum of jmax terms
c
            x5=dmax1(dabs(x3),dabs(x4))
            if (x5.lt.tol) goto 20
   10    continue
            ifault=4
            return
c
c       add terms to obtain integrals
c
   20       theta1=theta0+s1
            theta2=two*theta1-theta0+s2
            return
      else
c
c       z is between 1 and 3.6.
c       thet11 and thet12 contain the integrals of g1(x) and g2(x)
c       over (-infinity,1) and they are defined in the data statements
c
c       gaussian quadrature to integrate g1(x) and g2(x) over (1,z)
c       sums are storaged in sum1 and sum2, respectively
c
   30    sum1=zero
         sum2=zero
         x1=half*(z+one)
         x2=z-one
         do 40 iquad=1,8
            x3=x1-p(iquad)*x2
            x4=x1+p(iquad)*x2
            g1xl=g1(x3)
            g1xr=g1(x4)
            g2xl=(one+x3)*g1xl
            g2xr=(one+x4)*g1xr
            sum1=sum1+w(iquad)*(g1xl+g1xr)
            sum2=sum2+w(iquad)*(g2xl+g2xr)
   40    continue
         sum1=sum1*x2
         sum2=sum2*x2
c
c       add terms to obtain integrals
c
         theta1=sum1+thet11
         theta2=sum2+thet21
      endif
      return
c
c       asymptotic values.  z is greater than or equal to zmax
c
   50 theta0=one
      theta1=asymp1
      theta2=asymp2
      return
c
c       asymptotic values.  z is smaller than or equal to zmin
c
   60 theta0=zero
      theta1=zero
      theta2=zero
      return
      end
c=======================================================================
c
      subroutine levint(z,theta0,theta1,theta2,eta,ifault)
c
c-----------------------------------------------------------------------
c
c     computes: quantities needed to obtain the elements of the fisher
c               information matrix from a largest extreme value
c               distribution and censored data
c
c       theta0, theta1, and theta2 are up to an additive constant the
c       integrals over (-infinity,z) of the functions:
c       g0(x)=h(x)*h(x)*g(x), g1(x)=(1+x*h(x))*g(x),
c       and g2(x)=g(x)*(1+x*h(x))**2, where h(x)=1.  g(x) and bg(x) are
c       pdf and the cdf for a largest extreme value distribution
c
c       theta0, theta1, theta2, and eta all are computed using
c       sevint and the relationship between a largest and a smallest
c       extreme value distribution
c
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      double precision  asymp1, asymp2,
     * eta, etasev, one, theta0, theta1, theta2, t0sev,
     * t1sev, t2sev, z, zsev
c
c      asymp1 and asymp2 are the integrals for g1 and g2
c      when z=+infinity
c
      data asymp1/-0.4227843350984671d0/, asymp2/0.1823680660852879d1/,
     *     one/1.0d0/
c
c       call sevint to compute the fisher matrix elements
c
        zsev=-z
        call sevint(zsev,t0sev,t1sev,t2sev,etasev,ifault)
c
c       compute thetas and eta for the lev
c
      theta0=one-t0sev+etasev
      theta1=asymp1+t1sev+z*etasev
      theta2=asymp2-t2sev+z*z*etasev
      eta=etasev
      return
      end
c=======================================================================
c
      subroutine norint(z,theta0,theta1,theta2,eta,ifault)
c
c-----------------------------------------------------------------------
c
c     computes: quantities needed to obtain the elements of the fisher
c               information matrix from a normal distribution and
c               censored data
c
c        theta0, theta1, and theta2 are the integrals over (-infinity,z)
c        of the following functions: g0(x)=h(x)*h(x)*g(x), g1(x)=(one+
c        x*h(x))*g(x), and g2(x)=g(x)*(one+x*h(x))**2, where h(x)=psi(x)
c        +g(x)/(one-bg(x)) and psi(x)=-x. g(x) and bg(x) are the pdf and
c        cdf for a standard normal
c
c        theta0, theta1, theta2, and eta depend all on bg(z) and they
c        are obtained using the complementary error function
c
c        eta=g(x)*g(x)/((1-bg(x))*bg(x))
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      double precision  bgz, cdfnor, cval, denom, eta, gz, gzs, half,
     * one, pdfnor, root, sz, theta0, theta1, theta2, two, x, z, zero,
     * zmax, zmin, alnorm,dxerc
c
c       cval=1/sqrt(2*pi) and root=1/sqrt(2), where pi=3.14159....
c
      data cval/0.3989422804014326d0/, half/0.5d0/, one/1.0d0/,
     *     root/0.7071067811865475d0/, two/2.0d0/,  zero/0.0d0/,
     *     zmax/8.0d0/, zmin/-8.0d0/
c
c       normal pdf and cdf
c
      pdfnor(x)=cval*dexp(-half*x*x)
      cdfnor(x)=half*dxerc(-x*root)
c     cdfnor(x)=alnorm(x,.false.)
c
c       initial values
c
      ifault=0
      eta=zero
c
c       check for extreme z value
c
      if (z.ge.zmax) goto 10
      if (z.le.zmin) goto 20
c
c       compute theta0=integral of g0(x) over (-infinity,z)
c
      sz=cdfnor(-z)
      gz=pdfnor(z)
      gzs=gz*gz
      bgz=cdfnor(z)
      theta0=bgz-z*gz+gzs/sz
      denom=bgz*sz
      if (denom.gt.zero) eta=gzs/denom
c
c       compute theta1=integral of g1(x) over (-infinity,z)
c
      theta1=z*theta0-z*bgz-gz
c
c       compute theta2=integral of g2(x) over (-infinity,z)
c
      theta2=z*theta1+two*bgz
      return
c
c       asymptotic values.  z is greater than or equal to zmax
c
   10 theta0=one
      theta1=zero
      theta2=two
      return
c
c       asymptotic values.  z is smaller than or equal to zmin
c
   20 theta0=zero
      theta1=zero
      theta2=zero
      return
      end
c=======================================================================
c
      subroutine logint(z,theta0,theta1,theta2,eta,ifault)
c
c-----------------------------------------------------------------------
c
c     computes: quantities needed to obtain the elements of the fisher
c               information matrix from a logistic distribution and
c               censored data
c
c       theta0, theta1, and theta2 are the integrals over (-infinity,z)
c       of the following functions: g0(x)=h(x)*h(x)*g(x), g1(x)=(one+
c       x*h(x))*g(x), and g2(x)=g(x)*(one+x*h(x))**2, where
c       h(x)=1-bg(x).  g(x) and bg(x) are the pdf and cdf for a standard
c       logistic
c
c       theta0 and theta1 have closed form formulas. theta2 is computed
c       using euler's transformation on a powers series expansion of the
c       integral over (-infinity,z) of log(1+exp(x))
c
c       eta=g(x)*g(x)/((1-bg(x))*bg(x))
c
c       ifault= 5 if the euler's transformation used to accelerate the
c               convergence of s3 does not converge within tol of the
c               true value after including jmax terms in the expansion.
c               tol and jmax are set to 40 and  10**(-11) in the data
c               statements
c-----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      double precision  anuz, asymp0, asymp2, bgz, cdflog, const0, dum,
     * emabsz, eta, ez, gz, g3z, half, one, pdflog, sz, s3, s3old,
     * theta0, theta1, theta2, three, tmp, tol, two, x, x1, x2, x3, x4,
     * z, zero, zmax, zmin
c
      double precision  wksp(40)
c
c       asymp0, asymp2 are the integrals for g0 and g2 when
c       z=+infinity. const0=pi*pi/6, where pi=3.14159....
c
      data asymp0/0.3333333333333333d0/, asymp2/0.14299560445654842d1/,
     *     const0/0.1644934066848226d1/, half/0.5d0/, jmax/40/,
     *     one/1.0d0/, three/3.0d0/, tol/1.0d-11/, two/2.0d0/,
     *     zero/0.0d0/, zmax/34.0d0/, zmin/-34.0d0/
c
c       logistic pdf and cdf
c
      pdflog(x)=one/(dexp(x)*(one+dexp(-x))**2)
      cdflog(x)=one/(one+dexp(-x))
c
c       initial values
c
      ifault=0
      eta=zero
c
c       check for extreme z value
c
      if (z.ge.zmax) goto 40
      if (z.le.zmin) goto 50
c
c       compute theta0=integral of g0(x) over (-infinity,z)
c
      ez=dexp(z)
      sz=one/(one+ez)
      theta0=(one-sz**3)/three
c
c       compute theta1=integral of g1(x) over (-infinity,z) and eta
c
      bgz=cdflog(z)
      gz=pdflog(z)
      g3z=dlog(one+ez)
      theta1=z*theta0+(gz-g3z)/three
      eta=gz
c
c       compute theta2=integral of g2(x) over (-infinity,z)
c       s3 is computed using a power series expansion
c       with accelerated convergence from euler transformation
c
      emabsz=dexp(-dabs(z))
      x1=emabsz
      s3=zero
c
c       euler partial sum based on the first term
c
      n=1
      wksp(1)=x1
      s3=half*x1
      s3old=s3
      do 20 jterm=2,jmax
         x2=dfloat(jterm)**2
         x3=dfloat(jterm-1)**2
         x1=-x1*emabsz*x3/x2
c
c       euler partial sum based on two or more terms
c
         tmp=wksp(1)
         wksp(1)=x1
         do 10 j=1,n
            if (j.lt.n)dum=wksp(j+1)
            wksp(j+1)=half*(wksp(j)+tmp)
            if (j.lt.n)tmp=dum
   10    continue
c
c       euler improved partial sums
c
        if (dabs(wksp(n+1)).gt.dabs(wksp(n))) then
           s3=s3+wksp(n+1)
        else
           s3=s3+half*wksp(n+1)
           n=n+1
        endif
c
c       tests for convergence of the series -  the summation stops when
c       the absolute difference between two consecutive partial sums is
c       less than tol - a fault is declared if convergence is not
c       reached in a maximum of jmax terms
c
         x4=dabs(s3old-s3)
         s3old=s3
         if (x4.lt.tol) goto 30
   20 continue
      ifault=5
      return
c
c       add terms to obtain the integral
c
   30 anuz=s3
      if (z.gt.zero) anuz=const0+z*z/two-s3
      theta2=z*theta1+(bgz-z*g3z+z*gz+two*anuz)/three
      return
c
c       asymptotic values.  z is greater than or equal to zmax
c
   40 theta0=asymp0
      theta1=zero
      theta2=asymp2
      return
c
c       asymptotic values.  z is smaller than or equal to zmin
c
   50 theta0=zero
      theta1=zero
      theta2=zero
      return
c=======last card in lsinf=======
      end

c-----------------------------------------------------------------------
      double precision function dxerc(x)
c-----------------------------------------------------------------------
c  #
c
c  function to compute the complimentary error function
c
c   author - w. j. cody
c
c   date - january 8, 1985
c
c  #
      integer jint
      double precision x, result
c  #
      jint = 1
      call calerf(x,result,jint)
      dxerc = result
      return
c---------- last card of dxerc ----------
      end

c-----------------------------------------------------------------------
      subroutine calerf(arg,result,jint)
c-----------------------------------------------------------------------
c  #
c
c this packet computes the error and complimentary error functions
c   for real arguments  arg.  it contains two function type
c   subprograms,  erf  and  erfc  (or  derf  and  dxerc),  and one
c   subroutine type subprogram,  calerf.  the calling statements
c   for the primary entries are
c
c                   y=erf(x)     (or   y=derf(x) )
c   and
c                   y=erfc(x)    (or   y=dxerc(x) ).
c
c   the routine  calerf  is intended for internal packet use only,
c   all computations within the packet being concentrated in this
c   routine.  the function subprograms invoke  calerf  with the
c   statement
c          call calerf(arg,result,jint)
c   where the parameter usage is as follows
c
c      function                     parameters for calerf
c       call              arg                  result          jint
c     erf(arg)      any real argument         erf(arg)          0
c     erfc(arg)     abs(arg) .lt. xmax        erfc(arg)         1
c
c   the main computation evaluates near minimax approximations
c   from 'rational chebyshev approximations for the error function'
c   by w. j. cody, math. comp., 1969, pp. 631-638.  this
c   transportable program uses rational functions that theoretically
c   approximate  erf(x)  and  erfc(x)  to at least 18 significant
c   decimal digits.  the accuracy achieved depends on the arithmetic
c   system, the compiler, the intrinsic functions, and proper
c   selection of the machine-dependent constants.
c
c*******************************************************************
c*******************************************************************
c
c explanation of machine-dependent constants
c
c   xsmall = argument below which erf(x) may be represented
c            by   2*x/sqrt(pi)  and above which  x*x  will
c            not underflow.  a conservative value is the
c            largest x such that   1.0 + x = 1.0   to machine
c            precision.
c   xmax   = largest argument acceptable to  erfc;  solution to
c            equation:  w(x) * (1-0.5/x**2) = xmin,  where
c            w(x) = exp(-x*x)/(x*sqrt(pi)),  and xmin is the
c            smallest positive machine number (see table below).
c
c     approximate values for some important machines are:
c
c                          xsmall     xmax     xmin
c
c    ibm 195     (d.p.)   1.39d-17   13.306   5.40d-79
c    cdc 7600    (s.p.)   7.11e-15   25.922   3.13e-294
c    cray-1      (s.p.)   7.11e-15   75.326   4.58e-2467
c    univac 1108 (d.p.)   1.73d-18   26.582   2.78d-309
c    vax 11/780  (s.p.)   5.96e-8     9.269   2.94e-39
c    vax 11/780  (d.p.)   1.39d-17    9.269   2.94d-39
c    ibm pc      (s.p.)   5.96e-8     9.194   1.18e-38
c    ibm pc      (d.p.)   1.11d-16   26.543   2.23d-308
c
c*******************************************************************
c*******************************************************************
c
c error returns
c
c  the program returns  erfc = 0  for  arg .gt. xmax.
c
c
c other subprograms required (single precision version)
c
c     abs, exp
c
c other subprograms required (double precision version)
c
c     dabs, dexp
c
c
c  author: w. j. cody
c          mathematics and computer science division
c          argonne national laboratory
c          argonne, il 60439
c
c  latest modification: january 8, 1985
c
c  #
      integer i,jint
      double precision a,arg,b,c,d,four,half,p,one,q,result,sqrpi,
     1               two,thresh,x,xmax,xden,xnum,xsmall,y,ysq,zero
      dimension a(5),b(4),c(9),d(8),p(6),q(5)
c  #
c  mathematical constants
c  #
      data four,one,half,two,zero/4.0d0,1.0d0,0.5d0,2.0d0,0.0d0/
      data sqrpi/5.6418958354775628695d-1/,thresh/0.46875d0/
c  #
c  machine-dependent parameters
c  #
      data xsmall/4.2d-16/, xmax/9.269d0/
c  #
c  coefficients for approximation to derf in first interval
c  #
      data a/3.16112374387056560d00,1.13864154151050156d02,
     1       3.77485237685302021d02,3.20937758913846947d03,
     2       1.85777706184603153d-1/
      data b/2.36012909523441209d01,2.44024637934444173d02,
     1       1.28261652607737228d03,2.84423683343917062d03/
c  #
c  coefficients for approximation to dxerc in second interval
c  #
      data c/5.64188496988670089d-1,8.88314979438837594d0,
     1       6.61191906371416295d01,2.98635138197400131d02,
     2       8.81952221241769090d02,1.71204761263407058d03,
     3       2.05107837782607147d03,1.23033935479799725d03,
     4       2.15311535474403846d-8/
      data d/1.57449261107098347d01,1.17693950891312499d02,
     1       5.37181101862009858d02,1.62138957456669019d03,
     2       3.29079923573345963d03,4.36261909014324716d03,
     3       3.43936767414372164d03,1.23033935480374942d03/
c  #
c  coefficients for approximation to dxerc in third interval
c  #
      data p/3.05326634961232344d-1,3.60344899949804439d-1,
     1       1.25781726111229246d-1,1.60837851487422766d-2,
     2       6.58749161529837803d-4,1.63153871373020978d-2/
      data q/2.56852019228982242d00,1.87295284992346047d00,
     1       5.27905102951428412d-1,6.05183413124413191d-2,
     2       2.33520497626869185d-3/
c  #
      x = arg
      y = dabs(x)
      if (y .gt. four) go to 200
      if (y .gt. thresh) go to 100
c  #
c  evaluate erf for abs(x) .le. 0.46875
c  #
      ysq = zero
      if (y .gt. xsmall) ysq = y * y
      xnum = a(5)*ysq
      xden = ysq
      do 20 i = 1, 3
         xnum = (xnum + a(i)) * ysq
         xden = (xden + b(i)) * ysq
   20 continue
      result = x * (xnum + a(4)) / (xden + b(4))
      if (jint .ne. 0) result = one - result
      go to 800
c  #
c  evaluate erfc for 0.46875 .lt. abs(x) .le. 4.0
c  #
  100 ysq = y * y
      xnum = c(9)*y
      xden = y
      do 120 i = 1, 7
         xnum = (xnum + c(i)) * y
         xden = (xden + d(i)) * y
  120 continue
      result = dexp(-ysq) * (xnum + c(8)) / (xden + d(8))
      go to 300
c  #
c  evaluate erfc for abs(x) .gt. 4.0
c  #
  200 result = zero
      if (y .ge. xmax) go to 300
  220 ysq = one / (y * y)
      xnum = p(6)*ysq
      xden = ysq
      do 240 i = 1, 4
         xnum = (xnum + p(i)) * ysq
         xden = (xden + q(i)) * ysq
  240 continue
      result = ysq *(xnum + p(5)) / (xden + q(5))
      result = (dexp(-y*y) / y) * (sqrpi - result)
c  #
c  fix up for neg. arg., erf, etc.
c  #
  300 if (jint .eq. 0) go to 350
      if (x .lt. zero) result = two - result
      go to 800
  350 result = (half - result) + half
      if (x .lt. zero) result = -result
c  #
  800 return
c---------- last card of calerf ----------
      end
