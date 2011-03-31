c  Kernel Estimation of PDF 
c  2004-07-19
c
       subroutine kepdf(x, h, xeval, nx, ndim, neval, f)
       implicit double precision (a-h,o-z)
       dimension x(nx,ndim), h(nx,ndim), xeval(neval, ndim), f(neval)
       dimension x0(ndim)
       data zero/0.0d0/, one/1.0d0/, two/2.0d0/
       pi = 4.0d0 * datan(one)
       const = dsqrt((two*pi)**ndim)*nx
       do 100 ieval=1,neval
         s = zero
         do 105 j=1,ndim
           x0(j) = xeval(ieval,j)
  105    continue
         do 110 i=1,nx
           p=one
           do 120 j=1,ndim
             u = (x(i,j)-x0(j))/h(i,j)
             p = p * dexp(-u*u/two)/h(i,j)
  120      continue
           s = s + p
  110    continue
         f(ieval)= s/const
  100  continue
       return
       end


