      subroutine diagm(ellow,elup,invs,all,eonly) 
c********************************************************************* 
c     solves the secular equation for the general hermitian 
c     eigenvalue problem.  for systems with inversion (i.e. that 
c     reduce to a real symmetric matrix), real versions of the 
c     routines are used. for all=.true., all eigenvalues and 
c     and eigenvectors are determined. for all=.false., only 
c     those between ellow and elup are determined. if eonly=.true., 
c     then no eigenvectors are determined (only has an effect if 
c     all=.false.). 
c                  m. weinert 
c********************************************************************* 
      implicit real*8(a-h,o-z) 
      parameter (nvd=1500,neigd=1500) 
      logical all,eonly,invs 
c--->    task commons 
      common /eigv/ e(neigd),zr(nvd,neigd),zi(nvd,neigd),ne,nblw 
      common /secl/ s(nvd,nvd),h(nvd,nvd),nmat 
c--->    work arrays 
      common d(nvd),e1(nvd),e2(nvd),tau(2,nvd),t1(nvd),t2(nvd),t3(nvd), 
     +       t4(nvd),t5(nvd),index(nvd) 
      common /io/ iin,iout
      nm=nvd 
c--->    set energy range: if elup<ellow set all=.true. 
      eu=elup 
      el=ellow 
      if(eu.lt.el) all=.true. 
c--->    reduce the general problem to the standard problem 
      if(invs) then 
      call reducr(nm,nmat,h,s) 
      else 
      call reducc(nm,nmat,h,s) 
      endif 
c--->    reduce the standard problem to real tridiagonal form 
      if(invs) then 
      call tredr(nm,nmat,h,d,e1,e2) 
      else 
      call htrid3(nm,nmat,h,d,e1,e2,tau) 
      endif 
c--->    find eigenvalues and eigenvectors of tridiagonal system 
      if(all) then 
c--->    obtain all eigenvalues and eigenvectors 
         if(nmat.gt.neigd) stop 'neigd' 
      do 2 i=1,nmat 
      do 1 j=1,nmat 
    1 zr(j,i)=0.0 
    2 zr(i,i)=1.0 
      call tql2(nm,nmat,d,e1,zr,ierr) 
      ne=nmat 
      if(ierr.ne.0) then 
      write(16,1000) ierr 
 1000 format(//,' **** tql2:',i3,'th eigenvalue didnot converge') 
      stop 'diagm1' 
      endif 
      do 3 i=1,nmat 
    3 e(i)=d(i) 
      nblw=0 
      else 
c--->    determine eigenvalues (and eigenvectors) in the range (el,eu) 
c--->    reduce the energy window to guarantee that dimensions do 
c--->    not overflow 
      nblw=neblw1(d,e1,e2,el,nmat) 
    4 ne=neblw1(d,e1,e2,eu,nmat)-nblw 
      if(ne.le.neigd) go to 5 
      write( 6,1001) eu,ne,neigd 
      write(16,1001) eu,ne,neigd 
 1001 format(' $$$ too many eigenvalues: eu, ne1,neigd=',f10.5,2i6, 
     +       ' eu will be reduced by 0.025 a.u.') 
      eu=eu-0.025 
      if(eu.lt.el) stop 'diagm2' 
      go to 4 
    5 if(ne.le.0) return 
c--->    obtain eigenvalues 
      q=-1.0 
      call bisect(nmat,q,d,e1,e2,el,eu,neigd,ne,e,index,ierr,t1,t2) 
      if(ierr.ne.0) stop 'diagm3' 
      if( (ne.le.0) .or. eonly) return 
c--->    obtain eigenvectors of tridiagonal matrix 
      call tinvit(nm,nmat,d,e1,e2,ne,e,index,zr,ierr,t1,t2,t3,t4,t5) 
      if(ierr.ne.0) then 
      write(16,1002) ierr 
 1002 format(1x//,' **** tinvit:',i4,'th eigenvector didnot converge') 
      stop 'diagm4' 
      endif 
      endif 
c--->    back transform the eigenvectors 
      if(invs) then 
      call trbakr(nm,nmat,h,ne,zr) 
      call rebakr(nm,nmat,ne,s,zr) 
      do 10 j=1,ne 
      do 10 k=1,nmat 
   10 zi(k,j)=0.0 
      else 
      call htrib3(nm,nmat,h,tau,ne,zr,zi) 
      call rebakc(nm,nmat,ne,s,zr,zi) 
      endif 
      return 
      end 
      subroutine bisect(n,eps1,d,e,e2,lb,ub,mm,m,w,ind,ierr,rv4,rv5) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
      integer i,j,k,l,m,n,p,q,r,s,ii,mm,m1,m2,tag,ierr,isturm 
      dimension d(n),e(n),e2(n),w(mm),rv4(n),rv5(n) 
c      real u,v,lb,t1,t2,ub,xu,x0,x1,eps1,machep 
      real*8 lb,machep 
      integer ind(mm) 
c     this subroutine is a translation of the bisection technique 
c     in the algol procedure tristurm by peters and wilkinson. 
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971). 
c     this subroutine finds those eigenvalues of a tridiagonal 
c     symmetric matrix which lie in a specified interval, 
c     using bisection. 
c     on input- 
c        n is the order of the matrix, 
c        eps1 is an absolute error tolerance for the computed 
c          eigenvalues.  if the input eps1 is non-positive, 
c          it is reset for each submatrix to a default value, 
c          namely, minus the product of the relative machine 
c          precision and the 1-norm of the submatrix, 
c        d contains the diagonal elements of the input matrix, 
c        e contains the subdiagonal elements of the input matrix 
c          in its last n-1 positions.  e(1) is arbitrary, 
c        e2 contains the squares of the corresponding elements of e. 
c          e2(1) is arbitrary, 
c        lb and ub define the interval to be searched for eigenvalues. 
c          if lb is not less than ub, no eigenvalues will be found, 
c        mm should be set to an upper bound for the number of 
c          eigenvalues in the interval.  warning- if more than 
c          mm eigenvalues are determined to lie in the interval, 
c          an error return is made with no eigenvalues found. 
c     on output- 
c        eps1 is unaltered unless it has been reset to its 
c          (last) default value, 
c        d and e are unaltered, 
c        elements of e2, corresponding to elements of e regarded 
c          as negligible, have been replaced by zero causing the 
c          matrix to split into a direct sum of submatrices. 
c          e2(1) is also set to zero, 
c        m is the number of eigenvalues determined to lie in (lb,ub), 
c        w contains the m eigenvalues in ascending order, 
c        ind contains in its first m positions the submatrix indices 
c          associated with the corresponding eigenvalues in w -- 
c          1 for eigenvalues belonging to the first submatrix from 
c          the top, 2 for those belonging to the second submatrix, etc., 
c        ierr is set to 
c          zero       for normal return, 
c          3*n+1      if m exceeds mm, 
c        rv4 and rv5 are temporary storage arrays. 
c     the algol procedure sturmcnt contained in tristurm 
c     appears in bisect in-line. 
c     note that subroutine tql1 or imtql1 is generally faster than 
c     bisect, if more than n/4 eigenvalues are to be found. 
c     questions and comments should be directed to b. s. garbow, 
c     applied mathematics division, argonne national laboratory 
c     ------------------------------------------------------------------ 
c     ********** machep is a machine dependent parameter specifying 
c                the relative precision of floating point arithmetic. 
c                ********** 
      machep = 2.**(-47) 
      ierr = 0 
      tag = 0 
      t1 = lb 
      t2 = ub 
c     ********** look for small sub-diagonal entries ********** 
      do 40 i = 1, n 
         if (i .eq. 1) go to 20 
         if (abs(e(i)) .gt. machep * (abs(d(i)) + abs(d(i-1)))) 
     x      go to 40 
   20    e2(i) = 0.0 
   40 continue 
c     ********** determine the number of eigenvalues 
c                in the interval ********** 
      p = 1 
      q = n 
      x1 = ub 
      isturm = 1 
      go to 320 
   60 m = s 
      x1 = lb 
      isturm = 2 
      go to 320 
   80 m = m - s 
      if (m .gt. mm) go to 980 
      q = 0 
      r = 0 
c     ********** establish and process next submatrix, refining 
c                interval by the gerschgorin bounds ********** 
  100 if (r .eq. m) go to 1001 
      tag = tag + 1 
      p = q + 1 
      xu = d(p) 
      x0 = d(p) 
      u = 0.0 
      do 120 q = p, n 
         x1 = u 
         u = 0.0 
         v = 0.0 
         if (q .eq. n) go to 110 
         u = abs(e(q+1)) 
         v = e2(q+1) 
  110    xu = dmin1(d(q)-(x1+u),xu) 
         x0 = dmax1(d(q)+(x1+u),x0) 
         if (v .eq. 0.0) go to 140 
  120 continue 
  140 x1 = dmax1(abs(xu),abs(x0)) * machep 
      if (eps1 .le. 0.0) eps1 = -x1 
      if (p .ne. q) go to 180 
c     ********** check for isolated root within interval ********** 
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940 
      m1 = p 
      m2 = p 
      rv5(p) = d(p) 
      go to 900 
  180 x1 = x1 * float(q-p+1) 
      lb = dmax1(t1,xu-x1) 
      ub = dmin1(t2,x0+x1) 
      x1 = lb 
      isturm = 3 
      go to 320 
  200 m1 = s + 1 
      x1 = ub 
      isturm = 4 
      go to 320 
  220 m2 = s 
      if (m1 .gt. m2) go to 940 
c     ********** find roots by bisection ********** 
      x0 = ub 
      isturm = 5 
      do 240 i = m1, m2 
         rv5(i) = ub 
         rv4(i) = lb 
  240 continue 
c     ********** loop for k-th eigenvalue 
c                for k=m2 step -1 until m1 do -- 
c                (-do- not used to legalize computee-go-to) ********** 
      k = m2 
  250    xu = lb 
c     ********** for i=k step -1 until m1 do -- ********** 
         do 260 ii = m1, k 
            i = m1 + k - ii 
            if (xu .ge. rv4(i)) go to 260 
            xu = rv4(i) 
            go to 280 
  260    continue 
  280    if (x0 .gt. rv5(k)) x0 = rv5(k) 
c     ********** next bisection step ********** 
  300    x1 = (xu + x0) * 0.5 
         if ((x0 - xu) .le. (2.0 * machep * 
     x      (abs(xu) + abs(x0)) + abs(eps1))) go to 420 
c     ********** in-line procedure for sturm sequence ********** 
  320    s = p - 1 
         u = 1.0 
         do 340 i = p, q 
            if (u .ne. 0.0) go to 325 
            v = abs(e(i)) / machep 
            go to 330 
  325       v = e2(i) / u 
  330       u = d(i) - x1 - v 
         if(u .lt. 0.0) s=s+1 
  340    continue 
         go to (60,80,200,220,360), isturm 
c     ********** refine intervals ********** 
  360    if (s .ge. k) go to 400 
         xu = x1 
         if (s .ge. m1) go to 380 
         rv4(m1) = x1 
         go to 300 
  380    rv4(s+1) = x1 
         if (rv5(s) .gt. x1) rv5(s) = x1 
         go to 300 
  400    x0 = x1 
         go to 300 
c     ********** k-th eigenvalue found ********** 
  420    rv5(k) = x1 
      k = k - 1 
      if (k .ge. m1) go to 250 
c     ********** order eigenvalues tagged with their 
c                submatrix associations ********** 
  900 s = r 
      r = r + m2 - m1 + 1 
      j = 1 
      k = m1 
      do 920 l = 1, r 
         if (j .gt. s) go to 910 
         if (k .gt. m2) go to 940 
         if (rv5(k) .ge. w(l)) go to 915 
         do 905 ii = j, s 
            i = l + s - ii 
            w(i+1) = w(i) 
            ind(i+1) = ind(i) 
  905    continue 
  910    w(l) = rv5(k) 
         ind(l) = tag 
         k = k + 1 
         go to 920 
  915    j = j + 1 
  920 continue 
  940 if (q .lt. n) go to 100 
      go to 1001 
c     ********** set error -- underestimate of number of 
c                eigenvalues in interval ********** 
  980 ierr = 3 * n + 1 
 1001 lb = t1 
      ub = t2 
      return 
c     ********** last card of bisect ********** 
      end 
      subroutine htrib3(nm,n,a,tau,m,zr,zi) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
c***refer to  eisdoc 
c     this subroutine is a translation of a complex analogue of 
c     the algol procedure trbak3, num. math. 11, 181-195(1968) 
c     by martin, reinsch, and wilkinson. 
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). 
c     this subroutine forms the eigenvectors of a complex hermitian 
c     matrix by back transforming those of the corresponding 
c     real symmetric tridiagonal matrix determined by  htrid3. 
c     on input 
c        nm must be set to the row dimension of two-dimensional 
c          array parameters as declared in the calling program 
c          dimension statement. 
c        n is the order of the matrix. 
c        a contains information about the unitary transformations 
c          used in the reduction by  htrid3. 
c        tau contains further information about the transformations. 
c        m is the number of eigenvectors to be back transformed. 
c        zr contains the eigenvectors to be back transformed 
c          in its first m columns. 
c     on output 
c        zr and zi contain the real and imaginary parts, 
c          respectively, of the transformed eigenvectors 
c          in their first m columns. 
c     note that the last component of each returned vector 
c     is real and that vector euclidean norms are preserved. 
c     questions and comments should be directed to b. s. garbow, 
c     applied mathematics division, argonne national laboratory 
c     ------------------------------------------------------------------ 
c***routines called  (none) 
c***end prologue  htrib3 
      integer i,j,k,l,m,n,nm 
      dimension a(nm,n),tau(2,n),zr(nm,m),zi(nm,m) 
cc      real h,s,si 
c***first executable statement  htrib3 
      if (m .eq. 0) go to 200 
c     .......... transform the eigenvectors of the real symmetric 
c                tridiagonal matrix to those of the hermitian 
c                tridiagonal matrix. .......... 
      do 50 k = 1, n 
         do 50 j = 1, m 
            zi(k,j) = -zr(k,j) * tau(2,k) 
            zr(k,j) = zr(k,j) * tau(1,k) 
   50 continue 
      if (n .eq. 1) go to 200 
c     .......... recover and apply the householder matrices .......... 
      do 140 i = 2, n 
         l = i - 1 
         h = a(i,i) 
         if (h .eq. 0.0e0) go to 140 
         do 130 j = 1, m 
            s = 0.0e0 
            si = 0.0e0 
            do 110 k = 1, l 
               s = s + a(i,k) * zr(k,j) - a(k,i) * zi(k,j) 
               si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j) 
  110       continue 
c     .......... double divisions avoid possible underflow .......... 
            s = (s / h) / h 
            si = (si / h) / h 
            do 120 k = 1, l 
               zr(k,j) = zr(k,j) - s * a(i,k) - si * a(k,i) 
               zi(k,j) = zi(k,j) - si * a(i,k) + s * a(k,i) 
  120       continue 
  130    continue 
  140 continue 
  200 return 
      end 
      subroutine htrid3(nm,n,a,d,e,e2,tau) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
c***refer to  eisdoc 
c     this subroutine is a translation of a complex analogue of 
c     the algol procedure tred3, num. math. 11, 181-195(1968) 
c     by martin, reinsch, and wilkinson. 
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). 
c     this subroutine reduces a complex hermitian matrix, stored as 
c     a single square array, to a real symmetric tridiagonal matrix 
c     using unitary similarity transformations. 
c     on input 
c        nm must be set to the row dimension of two-dimensional 
c          array parameters as declared in the calling program 
c          dimension statement. 
c        n is the order of the matrix. 
c        a contains the lower triangle of the complex hermitian input 
c          matrix.  the real parts of the matrix elements are stored 
c          in the full lower triangle of a, and the imaginary parts 
c          are stored in the transposed positions of the strict upper 
c          triangle of a.  no storage is required for the zero 
c          imaginary parts of the diagonal elements. 
c     on output 
c        a contains information about the unitary transformations 
c          used in the reduction. 
c        d contains the diagonal elements of the the tridiagonal matrix. 
c        e contains the subdiagonal elements of the tridiagonal 
c          matrix in its last n-1 positions.  e(1) is set to zero. 
c        e2 contains the squares of the corresponding elements of e. 
c          e2 may coincide with e if the squares are not needed. 
c        tau contains further information about the transformations. 
c     questions and comments should be directed to b. s. garbow, 
c     applied mathematics division, argonne national laboratory 
c     ------------------------------------------------------------------ 
c***end prologue  htrid3 
      integer i,j,k,l,n,ii,nm,jm1,jp1 
      dimension a(nm,n),d(n),e(n),e2(n),tau(2,n) 
c      real f,g,h,fi,gi,hh,si,scale 
c***first executable statement  htrid3 
      tau(1,n) = 1.0e0 
      tau(2,n) = 0.0e0 
c     .......... for i=n step -1 until 1 do -- .......... 
      do 300 i=n,1,-1 
         h = 0.0e0 
         scale = 0.0e0 
         if (i-1.lt. 1) go to 130 
c     .......... scale row (algol tol then not needed) .......... 
         do 120 k = 1, i-1 
  120    scale = scale + abs(a(i,k)) + abs(a(k,i)) 
         if (scale .ne. 0.0e0) go to 140 
         tau(1,i-1) = 1.0e0 
         tau(2,i-1) = 0.0e0 
  130    e(i) = 0.0e0 
         e2(i) = 0.0e0 
         go to 290 
  140    do 150 k = 1, i-1 
            a(i,k) = a(i,k) / scale 
            a(k,i) = a(k,i) / scale 
            h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i) 
  150    continue 
         e2(i) = scale * scale * h 
         g = sqrt(h) 
         e(i) = scale * g 
         f = sqrt( a(i,i-1)**2 + a(i-1,i)**2 ) 
c     .......... form next diagonal element of matrix t .......... 
         if (f .eq. 0.0e0) go to 160 
         tau(1,i-1) = (a(i-1,i) * tau(2,i) - a(i,i-1) * tau(1,i)) / f 
         si = (a(i,i-1) * tau(2,i) + a(i-1,i) * tau(1,i)) / f 
         h = h + f * g 
         g = 1.0e0 + g / f 
         a(i,i-1) = g * a(i,i-1) 
         a(i-1,i) = g * a(i-1,i) 
         if (i-1 .eq. 1) go to 270 
         go to 170 
  160    tau(1,i-1) = -tau(1,i) 
         si = tau(2,i) 
         a(i,i-1) = g 
  170    f = 0.0e0 
         do 240 j = 1, i-1 
            g = 0.0e0 
            gi = 0.0e0 
c     .......... form element of a*u .......... 
            do 180 k = 1, j-1 
               g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i) 
               gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k) 
  180       continue 
  190       g = g + a(j,j) * a(i,j) 
            gi = gi - a(j,j) * a(j,i) 
            do 200 k = j+1,i-1 
               g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i) 
               gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k) 
  200       continue 
c     .......... form element of p .......... 
  220       e(j) = g / h 
            tau(2,j) = gi / h 
            f = f + e(j) * a(i,j) - tau(2,j) * a(j,i) 
  240    continue 
         hh = f / (h + h) 
c     .......... form reduced a .......... 
         do 260 j = 1, i-1 
            f = a(i,j) 
            g = e(j) - hh * f 
            e(j) = g 
            fi = -a(j,i) 
            gi = tau(2,j) - hh * fi 
            tau(2,j) = -gi 
            a(j,j) = a(j,j) - 2.0e0 * (f * g + fi * gi) 
            if (j .eq. 1) go to 260 
            do 250 k = 1, j-1 
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k) 
     x                         + fi * tau(2,k) + gi * a(k,i) 
               a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i) 
     x                         - fi * e(k) - gi * a(i,k) 
  250       continue 
  260    continue 
  270    do 280 k = 1, i-1 
            a(i,k) = scale * a(i,k) 
            a(k,i) = scale * a(k,i) 
  280    continue 
         tau(2,i-1) = -si 
  290    d(i) = a(i,i) 
         a(i,i) = scale * sqrt(h) 
  300 continue 
      return 
      end 
      function neblw1(d,e,e2,ub,n) 
      implicit real*8 (a-h,o-z) 
c*********************************************************************** 
c     counts the number of eigenvalues below ub in the tridiagonal 
c     matrix of order n. 
c     d contains the diagonal elements, e the off-diagonal elements, 
c     and e2=e**2  (e(1), e2(1) are arbitary). 
c     this routine is based on the algol routine sturmcnt given in 
c     wilkenson and reinsch, linear algebra, p. 423 (1971). 
c           m. weinert  july 1982 
c*********************************************************************** 
c 
      real*8 machep 
      dimension d(n),e(n),e2(n) 
c--->    (relative machine precision)**-1 
      machep=2.**(47) 
      m=0 
      x=d(1)-ub 
      if(x.lt. 0.0) m=m+1 
      do 3 i=2,n 
      if(x.ne.0.0) go to 1 
      u=abs(e(i))*machep 
      go to 2 
    1 u=e2(i)/x 
    2 x=d(i)-ub-u 
      if(x.lt. 0.0) m=m+1 
    3 continue 
      neblw1=m 
      return 
      end 
      subroutine rebakc(nm,n,m,b,zr,zi) 
      implicit real*8 (a-h,o-z) 
c********************************************************************** 
c     complex version of the algol procedure rebaka, linear 
c     algebra, vol. ii, 1971 by wilkinson and reinsch. 
c     forms the eigenvectors of the generalized hermitian 
c     eigensystem by back-transforming those of the derived 
c     standard matrix determined by reducc. 
c     input: 
c      nm     row dimension of the 2-d arrays 
c      n      order of the matrix system 
c      m      number of eigenvectors to back-transform 
c      b      contains the cholesky decomposition obtained 
c             in reducc in compact storage mode. 
c      zr,zi  contain the real and imaginary parts of the 
c             eigenvectors to be back-transformed in the 
c             first m columns. 
c     output: 
c      zr,zi  contain the back-transformed eigenvectors 
c                 m. weinert   july 1983 
c********************************************************************** 
      dimension b(nm,n),zr(nm,n),zi(nm,n) 
      do 2 j=1,m 
      do 2 ii=1,n 
      i=n+1-ii 
      xr=zr(i,j) 
      xi=zi(i,j) 
      do 1 k=i+1,n 
      xr=xr - b(k,i)*zr(k,j) - b(i,k)*zi(k,j) 
      xi=xi - b(k,i)*zi(k,j) + b(i,k)*zr(k,j) 
    1 continue 
      zr(i,j)=xr/b(i,i) 
      zi(i,j)=xi/b(i,i) 
    2 continue 
      return 
      end 
      subroutine rebakr(nm,n,m,b,zr) 
      implicit real*8 (a-h,o-z) 
c********************************************************************** 
c     real version of the algol procedure rebaka, linear 
c     algebra, vol. ii, 1971 by wilkinson and reinsch. 
c     forms the eigenvectors of the generalized real symmetric 
c     eigensystem by back-transforming those of the derived 
c     standard matrix determined by reducr. 
c     input: 
c      nm     row dimension of the 2-d arrays 
c      n      order of the matrix system 
c      m      number of eigenvectors to back-transform 
c      b      contains the cholesky decomposition obtained 
c             in reducc (lower triangle only: b(i,j), i.gej) 
c      zr     contains the eigenvectors to be back-transformed in the 
c             first m columns. 
c     output: 
c      zr     contains the back-transformed eigenvectors 
c                 m. weinert   jan. 1987 
c********************************************************************** 
      dimension b(nm,n),zr(nm,n) 
      do 2 j=1,m 
      do 2 ii=1,n 
      i=n+1-ii 
      xr=zr(i,j) 
      do 1 k=i+1,n 
      xr=xr - b(k,i)*zr(k,j) 
    1 continue 
      zr(i,j)=xr/b(i,i) 
    2 continue 
      return 
      end 
      subroutine reducc(nm,n,a,b) 
      implicit real*8 (a-h,o-z) 
c********************************************************************** 
c     complex version of the algol procedure reduc1, linear 
c     algebra, vol. ii, wilkinson and reinsch, 1971. 
c     reduction of the general hermitian eigenvalue problem 
c     a*x=lamda*b*x to the equivalent problem p*z=lamda*z 
c     using cholesky decomposition. 
c     the procedure will fail if b, perhaps due to rounding 
c     errors, is not positive definite. 
c     input: 
c      nm     row dimension of 2-d arrays a and b as declared in 
c             call routine 
c      n      order of the eigensystem 
c      a,b    the lower triangle of the hermitian matrices stored 
c             in compact mode as:  (i.ge.j) 
c                a(i,j) = real ( (h(i,j) ) 
c                a(j,i) = imag ( (h(i,j) ) 
c     output: 
c      a      contains the lower triangle of the reduced problem 
c             stored in compact mode 
c      b      contains the lower triangular cholesky decomposition 
c             stored in compact mode 
c                     m. weinert     july 1983 
c     note that a is in the form required in htrid3 and b is in the 
c     form required in rebakc. 
c********************************************************************** 
      dimension a(nm,n),b(nm,n) 
ccc   almost absolute zero was defined by Mike 
ccc      data zero/0.0e0/ 
ccc   change it to the same as old aplmat.f 
      data zero/1.0d-55/ 
c--->    form l in lower triangle of b 
      do 3 j=1,n 
      xr=b(j,j) 
      do 1 k=1,j-1 
    1 xr=xr - b(j,k)*b(j,k) - b(k,j)*b(k,j) 
      if(xr.le.zero) go to 30 
      y=sqrt(xr) 
      b(j,j)=y 
      do 3 i=j+1,n 
      xr=b(i,j) 
      xi=b(j,i) 
      do 2 k=1,j-1 
      xr=xr - b(i,k)*b(j,k) - b(k,i)*b(k,j) 
      xi=xi - b(k,i)*b(j,k) + b(i,k)*b(k,j) 
    2 continue 
      b(j,i)=xi/y 
      b(i,j)=xr/y 
    3 continue 
c--->    form hermitian conjugate of inv(l)*a 
      do 13 j=1,n 
      y=b(j,j) 
      xr=a(j,j) 
      do 11 k=1,j-1 
   11 xr=xr - b(j,k)*a(j,k) - b(k,j)*a(k,j) 
      a(j,j)=xr/y 
      do 13 i=j+1,n 
      xr=a(i,j) 
      xi=a(j,i) 
      do 12 k=1,j-1 
      xr=xr - b(j,k)*a(i,k) - b(k,j)*a(k,i) 
      xi=xi + b(k,j)*a(i,k) - b(j,k)*a(k,i) 
   12 continue 
      a(j,i)=xi/y 
      a(i,j)=xr/y 
   13 continue 
c--->    premultiply by inv(l) 
      do 26 i=1,n 
      y=b(i,i) 
      do 24 j=1,i-1 
      xr=a(i,j) 
      xi=a(j,i) 
      xr=xr - b(i,j)*a(j,j) 
      xi=xi - b(j,i)*a(j,j) 
      do 21 k=j+1,i-1 
      xr=xr - b(i,k)*a(k,j) + b(k,i)*a(j,k) 
   21 xi=xi - b(k,i)*a(k,j) - b(i,k)*a(j,k) 
      do 23 k=1,j-1 
      xr=xr - b(i,k)*a(j,k) - b(k,i)*a(k,j) 
   23 xi=xi - b(k,i)*a(j,k) + b(i,k)*a(k,j) 
      a(j,i)=xi/y 
      a(i,j)=xr/y 
   24 continue 
      xr=a(j,j) 
      do 25 k=1,j-1 
   25 xr=xr - b(j,k)*a(j,k) - b(k,j)*a(k,j) 
      a(j,j)=xr/y 
   26 continue 
      return 
c--->    error section 
   30 write(6,1000) j,n,xr 
      stop 
 1000 format(//,' $$$reducc: singularity in row',i3,' of the',i3, 
     +          ' dimensional matrix. diagonal element=',e15.6) 
      end 
      subroutine reducr(nm,n,a,b) 
      implicit real*8 (a-h,o-z) 
c********************************************************************** 
c     real version of the algol procedure reduc1, linear 
c     algebra, vol. ii, wilkinson and reinsch, 1971. 
c     reduction of the general real symmetric eigenvalue problem 
c     a*x=lamda*b*x to the equivalent problem p*z=lamda*z 
c     using cholesky decomposition. 
c     the procedure will fail if b, perhaps due to rounding 
c     errors, is not positive definite. 
c     input: 
c      nm     row dimension of 2-d arrays a and b as declared in 
c             call routine 
c      n      order of the eigensystem 
c      a,b    the lower triangle of the symmetric matrices stored 
c             in standard mode: a(i,j), b(i,j), i.ge.j 
c     output: 
c      a      contains the lower triangle of the reduced problem 
c      b      contains the lower triangular cholesky decomposition 
c                     m. weinert     jan. 1987 
c     note that a is in the form required in tredr and b is in the 
c     form required in rebakr. 
c********************************************************************** 
      dimension a(nm,n),b(nm,n) 
ccc   almost absolute zero was defined by Mike 
ccc      data zero/0.0e0/ 
ccc   change it to the same as old aplmat.f 
      data zero/1.0d-55/ 
c--->    form l in lower triangle of b 
      do 3 j=1,n 
      xr=b(j,j) 
      do 1 k=1,j-1 
    1 xr=xr - b(j,k)*b(j,k) 
      if(xr.le.zero) go to 30 
      y=sqrt(xr) 
      b(j,j)=y 
      do 3 i=j+1,n 
      xr=b(i,j) 
      do 2 k=1,j-1 
    2 xr=xr - b(i,k)*b(j,k) 
      b(i,j)=xr/y 
    3 continue 
c--->    form hermitian conjugate of inv(l)*a 
      do 13 j=1,n 
      y=b(j,j) 
      xr=a(j,j) 
      do 11 k=1,j-1 
   11 xr=xr - b(j,k)*a(j,k) 
      a(j,j)=xr/y 
      do 13 i=j+1,n 
      xr=a(i,j) 
      do 12 k=1,j-1 
   12 xr=xr - b(j,k)*a(i,k) 
      a(i,j)=xr/y 
   13 continue 
c--->    premultiply by inv(l) 
      do 23 i=1,n 
      y=b(i,i) 
      do 23 j=1,i 
      xr=a(i,j) 
      do 21 k=1,j-1 
   21 xr=xr - b(i,k)*a(j,k) 
      do 22 k=j,i-1 
   22 xr=xr - b(i,k)*a(k,j) 
      a(i,j)=xr/y 
   23 continue 
      return 
c--->    error section 
   30 write(6,1000) j,n,xr 
      stop 
 1000 format(//,' $$$reducr: singularity in row',i3,' of the',i3, 
     +          ' dimensional matrix. diagonal element=',e15.6) 
      end 
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z, 
     x                  ierr,rv1,rv2,rv3,rv4,rv6) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group 
      dimension d(n),e(n),e2(n),w(m),z(nm,m), 
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n) 
cc      real u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,machep 
      real*8 norm,machep 
c     real sqrt,abs,float 
      integer ind(m) 
c     this subroutine is a translation of the inverse iteration tech- 
c     nique in the algol procedure tristurm by peters and wilkinson. 
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971). 
c     this subroutine finds those eigenvectors of a tridiagonal 
c     symmetric matrix corresponding to specified eigenvalues, 
c     using inverse iteration. 
c     on input- 
c        nm must be set to the row dimension of two-dimensional 
c          array parameters as declared in the calling program 
c          dimension statement, 
c        n is the order of the matrix, 
c        d contains the diagonal elements of the input matrix, 
c        e contains the subdiagonal elements of the input matrix 
c          in its last n-1 positions.  e(1) is arbitrary, 
c        e2 contains the squares of the corresponding elements of e, 
c          with zeros corresponding to negligible elements of e. 
c          e(i) is considered negligible if it is not larger than 
c          the product of the relative machine precision and the sum 
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain 
c          0.0 if the eigenvalues are in ascending order, or 2.0 
c          if the eigenvalues are in descending order.  if  bisect, 
c          tridib, or  imtqlv  has been used to find the eigenvalues, 
c          their output e2 array is exactly what is expected here, 
c        m is the number of specified eigenvalues, 
c        w contains the m eigenvalues in ascending or descending order, 
c        ind contains in its first m positions the submatrix indices 
c          associated with the corresponding eigenvalues in w -- 
c          1 for eigenvalues belonging to the first submatrix from 
c          the top, 2 for those belonging to the second submatrix, etc. 
c     on output- 
c        all input arrays are unaltered, 
c        z contains the associated set of orthonormal eigenvectors. 
c          any vector which fails to converge is set to zero, 
c        ierr is set to 
c          zero       for normal return, 
c          -r         if the eigenvector corresponding to the r-th 
c                     eigenvalue fails to converge in 5 iterations, 
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays. 
c     questions and comments should be directed to b. s. garbow, 
c     applied mathematics division, argonne national laboratory 
c     ------------------------------------------------------------------ 
c     ********** machep is a machine dependent parameter specifying 
c                the relative precision of floating point arithmetic. 
c                ********** 
      machep = 2.**(-47) 
      ierr = 0 
      if (m .eq. 0) go to 1001 
      tag = 0 
      order = 1.0 - e2(1) 
      q = 0 
c     ********** establish and process next submatrix ********** 
  100 p = q + 1 
      do 120 q = p, n 
         if (q .eq. n) go to 140 
         if (e2(q+1) .eq. 0.0) go to 140 
  120 continue 
c     ********** find vectors by inverse iteration ********** 
  140 tag = tag + 1 
      s = 0 
      do 920 r = 1, m 
         if (ind(r) .ne. tag) go to 920 
         its = 1 
         x1 = w(r) 
         if (s .ne. 0) go to 510 
c     ********** check for isolated root ********** 
         xu = 1.0 
         if (p .ne. q) go to 490 
         rv6(p) = 1.0 
         go to 870 
  490    norm = abs(d(p)) 
         ip = p + 1 
      i=q-ip+1 
      if(i.lt.1) go to 500 
      norm=norm+sasum(i,d(ip),1)+sasum(i,e(ip),1) 
  500 continue 
c     ********** eps2 is the criterion for grouping, 
c                eps3 replaces zero pivots and equal 
c                roots are modified by eps3, 
c                eps4 is taken very small to avoid overflow ********** 
         eps2 = 1.0e-3 * norm 
         eps3 = machep * norm 
         uk = float(q-p+1) 
         eps4 = uk * eps3 
         uk = eps4 / sqrt(uk) 
         s = p 
  505    group = 0 
         go to 520 
c     ********** look for close or coincident roots ********** 
  510    if (abs(x1-x0) .ge. eps2) go to 505 
         group = group + 1 
         if (order * (x1 - x0) .le. 0.0) x1 = x0 + order * eps3 
c     ********** elimination with interchanges and 
c                initialization of vector ********** 
  520    v = 0.0 
         do 580 i = p, q 
            rv6(i) = uk 
            if (i .eq. p) go to 560 
            if (abs(e(i)) .lt. abs(u)) go to 540 
c     ********** warning -- a divide check may occur here if 
c                e2 array has not been specified correctly ********** 
            xu = u / e(i) 
            rv4(i) = xu 
            rv1(i-1) = e(i) 
            rv2(i-1) = d(i) - x1 
            rv3(i-1) = 0.0 
            if (i .ne. q) rv3(i-1) = e(i+1) 
            u = v - xu * rv2(i-1) 
            v = -xu * rv3(i-1) 
            go to 580 
  540       xu = e(i) / u 
            rv4(i) = xu 
            rv1(i-1) = u 
            rv2(i-1) = v 
            rv3(i-1) = 0.0 
  560       u = d(i) - x1 - xu * v 
            if (i .ne. q) v = e(i+1) 
  580    continue 
         if (u .eq. 0.0) u = eps3 
         rv1(q) = u 
         rv2(q) = 0.0 
         rv3(q) = 0.0 
c     ********** back substitution 
c                for i=q step -1 until p do -- ********** 
  600    do 620 ii = p, q 
            i = p + q - ii 
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i) 
            v = u 
            u = rv6(i) 
  620    continue 
c     ********** orthogonalize with respect to previous 
c                members of group ********** 
         if (group .eq. 0) go to 700 
         j = r 
         do 680 jj = 1, group 
  630       j = j - 1 
            if (ind(j) .ne. tag) go to 630 
            xu = 0.0 
            do 640 i = p, q 
  640       xu = xu + rv6(i) * z(i,j) 
            do 660 i = p, q 
  660       rv6(i) = rv6(i) - xu * z(i,j) 
  680    continue 
  700    norm = 0.0 
         do 720 i = p, q 
  720    norm = norm + abs(rv6(i)) 
         if (norm .ge. 1.0) go to 840 
c     ********** forward substitution ********** 
         if (its .eq. 5) go to 830 
         if (norm .ne. 0.0) go to 740 
         rv6(s) = eps4 
         s = s + 1 
         if (s .gt. q) s = p 
         go to 780 
  740    xu = eps4 / norm 
         do 760 i = p, q 
  760    rv6(i) = rv6(i) * xu 
c     ********** elimination operations on next vector 
c                iterate ********** 
  780    do 820 i = ip, q 
            u = rv6(i) 
c     ********** if rv1(i-1) .eq. e(i), a row interchange 
c                was performed earlier in the 
c                triangularization process ********** 
            if (rv1(i-1) .ne. e(i)) go to 800 
            u = rv6(i-1) 
            rv6(i-1) = rv6(i) 
  800       rv6(i) = u - rv4(i) * rv6(i-1) 
  820    continue 
         its = its + 1 
         go to 600 
c     ********** set error -- non-converged eigenvector ********** 
  830    ierr = -r 
         xu = 0.0 
         go to 870 
c     ********** normalize so that sum of squares is 
c                1 and expand to full order ********** 
  840    u = 0.0 
         do 860 i = p, q 
  860    u = u + rv6(i)**2 
         xu = 1.0 / sqrt(u) 
  870    do 880 i = 1, n 
  880    z(i,r) = 0.0 
         do 900 i = p, q 
  900    z(i,r) = rv6(i) * xu 
         x0 = x1 
  920 continue 
      if (q .lt. n) go to 100 
 1001 return 
c     ********** last card of tinvit ********** 
      end 
      subroutine tql2(nm,n,d,e,z,ierr) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
c***refer to  eisdoc 
c     this subroutine is a translation of the algol procedure tql2, 
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and 
c     wilkinson. 
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971). 
c     this subroutine finds the eigenvalues and eigenvectors 
c     of a symmetric tridiagonal matrix by the ql method. 
c     the eigenvectors of a full symmetric matrix can also 
c     be found if  tred2  has been used to reduce this 
c     full matrix to tridiagonal form. 
c     on input 
c        nm must be set to the row dimension of two-dimensional 
c          array parameters as declared in the calling program 
c          dimension statement. 
c        n is the order of the matrix. 
c        d contains the diagonal elements of the input matrix. 
c        e contains the subdiagonal elements of the input matrix 
c          in its last n-1 positions.  e(1) is arbitrary. 
c        z contains the transformation matrix produced in the 
c          reduction by  tred2, if performed.  if the eigenvectors 
c          of the tridiagonal matrix are desired, z must contain 
c          the identity matrix. 
c      on output 
c        d contains the eigenvalues in ascending order.  if an 
c          error exit is made, the eigenvalues are correct but 
c          unordered for indices 1,2,...,ierr-1. 
c        e has been destroyed. 
c        z contains orthonormal eigenvectors of the symmetric 
c          tridiagonal (or full) matrix.  if an error exit is made, 
c          z contains the eigenvectors associated with the stored 
c          eigenvalues. 
c        ierr is set to 
c          zero       for normal return, 
c          j          if the j-th eigenvalue has not been 
c                     determined after 30 iterations. 
c     questions and comments should be directed to b. s. garbow, 
c     applied mathematics division, argonne national laboratory 
c     ------------------------------------------------------------------ 
c***end prologue  tql2 
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr 
      dimension d(n),e(n),z(nm,n) 
ccc      real b,c,c2,c3,dl1,el1,f,g,h,p,r,s,s2 
c***first executable statement  tql2 
      ierr = 0 
      if (n .eq. 1) go to 1001 
      do 100 i = 2, n 
  100 e(i-1) = e(i) 
      f = 0.0e0 
      b = 0.0e0 
      e(n) = 0.0e0 
      do 240 l = 1, n 
         j = 0 
         h = abs(d(l)) + abs(e(l)) 
         if (b .lt. h) b = h 
c     .......... look for small sub-diagonal element .......... 
         do 110 m = l, n 
            if (b + abs(e(m)) .eq. b) go to 120 
c     .......... e(n) is always zero, so there is no exit 
c                through the bottom of the loop .......... 
  110    continue 
  120    if (m .eq. l) go to 220 
  130    if (j .eq. 30) go to 1000 
         j = j + 1 
c     .......... form shift .......... 
         l1 = l + 1 
         l2 = l1 + 1 
         g = d(l) 
         p = (d(l1) - g) / (2.0e0 * e(l)) 
         r = sqrt( 1.0 + p**2 ) 
         d(l) = e(l) / (p + sign(r,p)) 
         d(l1) = e(l) * (p + sign(r,p)) 
         dl1 = d(l1) 
         h = g - d(l) 
         if (l2 .gt. n) go to 145 
         do 140 i = l2, n 
  140    d(i) = d(i) - h 
  145    f = f + h 
c     .......... ql transformation .......... 
         p = d(m) 
         c = 1.0e0 
         c2 = c 
         el1 = e(l1) 
         s = 0.0e0 
         mml = m - l 
c     .......... for i=m-1 step -1 until l do -- .......... 
         do 200 ii = 1, mml 
            c3 = c2 
            c2 = c 
            s2 = s 
            i = m - ii 
            g = c * e(i) 
            h = c * p 
            if (abs(p) .lt. abs(e(i))) go to 150 
            c = e(i) / p 
            r = sqrt(c*c+1.0e0) 
            e(i+1) = s * p * r 
            s = c / r 
            c = 1.0e0 / r 
            go to 160 
  150       c = p / e(i) 
            r = sqrt(c*c+1.0e0) 
            e(i+1) = s * e(i) * r 
            s = 1.0e0 / r 
            c = c * s 
  160       p = c * d(i) - s * g 
            d(i+1) = h + s * (c * g + s * d(i)) 
c     .......... form vector .......... 
            do 180 k = 1, n 
               h = z(k,i+1) 
               z(k,i+1) = s * z(k,i) + c * h 
               z(k,i) = c * z(k,i) - s * h 
  180       continue 
  200    continue 
         p = -s * s2 * c3 * el1 * e(l) / dl1 
         e(l) = s * p 
         d(l) = c * p 
         if (b + abs(e(l)) .gt. b) go to 130 
  220    d(l) = d(l) + f 
  240 continue 
c     .......... order eigenvalues and eigenvectors .......... 
      do 300 ii = 2, n 
         i = ii - 1 
         k = i 
         p = d(i) 
         do 260 j = ii, n 
            if (d(j) .ge. p) go to 260 
            k = j 
            p = d(j) 
  260    continue 
         if (k .eq. i) go to 300 
         d(k) = d(i) 
         d(i) = p 
         do 280 j = 1, n 
            p = z(j,i) 
            z(j,i) = z(j,k) 
            z(j,k) = p 
  280    continue 
  300 continue 
      go to 1001 
c     .......... set error -- no convergence to an 
c                eigenvalue after 30 iterations .......... 
 1000 ierr = l 
 1001 return 
      end 
      subroutine trbakr(nm,n,a,m,z) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
c     translation of the algol procedure trbak1, linear algebra, 
c     vol. ii, 1971 by wilkinson and reinsch. 
c     forms the eigenvectors of a real symmetric matrix by back 
c     transforming those of the corresponding symmetric tridiagonal 
c     matrix determined by tredr. 
c     input: 
c       nm    row dimension of the 2-d arrays 
c       n     order of the matrix system 
c       a     contains information about the orthogonal transformations 
c             obtained in tredr 
c       m     number of eigenvectors to be back transformed 
c       z     contains the eigenvectors to be back transformed in 
c             the first m columns 
c     output: 
c       z     contains the transformed eigenvectors in the first 
c             m columns 
c                  m. weinert  jan. 1987 
c********************************************************************* 
      dimension a(nm,n),z(nm,m) 
      if(m.eq.0 .or. n.eq.1) return 
      do 4 i=2,n 
      h=a(i,i) 
      if(h.ne.0.0) then 
      do 3 j=1,m 
      s=0.0 
      do 1 k=1,i-1 
    1 s=s+a(i,k)*z(k,j) 
      s=(s/h)/h 
      do 2 k=1,i-1 
    2 z(k,j)=z(k,j)-s*a(i,k) 
    3 continue 
      endif 
    4 continue 
      return 
      end 
      subroutine tredr(nm,n,a,d,e,e2) 
      implicit real*8 (a-h,o-z) 
c********************************************************************* 
c     real version of the algol procedure tred1, linear algebra, 
c     vol. ii, 1971 by wilkinson and reinsch. 
c     reduces the real symmetric matrix to tridiagonal form using 
c     householder's reduction. 
c     input: 
c       nm    row dimension of 2-d arrays 
c       n     order of the matix system 
c       a     lower triangle of the real symmetric matrix, 
c             only the elements a(i,j), i.ge.j, need be given 
c     output: 
c       a     contains information about the unitary transformations 
c       d     diagonal elements of the tridiagonal matrix 
c       e     subdiagonal elements of the tridiagonal matrix in its 
c             last n-1 positions. e(1)=0.0 
c       e2    squares of the corresponding elements of e. e2 may 
c             coincide with e if the squares are not needed. 
c                       m. weinert  jan. 1987 
c********************************************************************* 
      dimension a(nm,n),d(n),e(n),e2(n) 
      do 10 i=n,2,-1 
      do 1 k=1,i-1 
    1 d(k)=a(i,k) 
c--->    scale rows 
      scale=0.0 
      do 2 k=1,i-1 
    2 scale=max(scale,abs(d(k))) 
c--->    if scale is too small to guarantee orthogonality, 
c--->    transformation is skipped 
      if(scale.eq.0.0) then 
      e(i)=0.0 
      e2(i)=0.0 
      d(i)=a(i,i) 
      a(i,i)=0.0 
      else 
      do 3 k=1,i-1 
    3 d(k)=d(k)/scale 
      h=0.0 
      do 4 k=1,i-1 
    4 h=h+d(k)*d(k) 
      e2(i)=h*scale*scale 
      f=d(i-1) 
      g=-sign(sqrt(h),f) 
      e(i)=scale*g 
      h=h-f*g 
      d(i-1)=f-g 
      a(i,i-1)=scale*d(i-1) 
      f=0.0 
c--->    form element of a*u 
      do 7 j=1,i-1 
      g=0.0 
      do 5 k=1,j 
    5 g=g+d(k)*a(j,k) 
      do 6 k=j+1,i-1 
    6 g=g+d(k)*a(k,j) 
c--->    form element of p 
      e(j)=g/h 
      f=f+e(j)*d(j) 
    7 continue 
c--->    form k 
      hh=f/(h+h) 
c--->    form reduced a 
      do 9 j=1,i-1 
      f=d(j) 
      g=e(j)-hh*f 
      e(j)=g 
      do 8 k=1,j 
    8 a(j,k)=a(j,k)-f*e(k)-g*d(k) 
    9 continue 
      d(i)=a(i,i) 
      a(i,i)=scale*sqrt(h) 
      endif 
   10 continue 
      e2(1)=0.0 
      e(1)=0.0 
      d(1)=a(1,1) 
      a(1,1)=0.0 
      return 
      end 
      function sasum (n,sx,incx) 
      implicit real*8 (a-h,o-z) 
c                                 specifications for arguments 
      integer            n,incx 
      dimension          sx(1) 
c                                  specifications for local variables 
      integer            i,m,mp1,ns 
c                                  first executable statement 
      sasum = 0.00 
      if (n.le.0) return 
      if (incx.eq.1) go to 10 
c                                  code for increments not equal to 1. 
      ns = n*incx 
      do 5 i=1,ns,incx 
         sasum = sasum+abs(sx(i)) 
    5 continue 
      return 
c                                  code for increments equal to 1. 
c                                    clean-up loop so remaining vector 
c                                    length is a multiple of 6. 
   10 m = n-(n/6)*6 
      if (m.eq.0) go to 20 
      do 15 i=1,m 
         sasum = sasum+abs(sx(i)) 
   15 continue 
      if (n.lt.6) return 
   20 mp1 = m+1 
      do 25 i=mp1,n,6 
      sasum = sasum+abs(sx(i))+abs(sx(i+1)) 
     &    +abs(sx(i+2))+abs(sx(i+3)) 
     &    +abs(sx(i+4))+abs(sx(i+5)) 
   25 continue 
      return 
      end 
