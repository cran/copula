c
c     common declarations for stable programs 
c
c     Return codes for stable programs 
c
c          0 = no error
c          1 = parameter error on some input value (e.g. alpha .le. TOL(9),
c		 alpha>2, |beta|>1, asked for p-th quantile where p <= 0 
c                or p >= 1,etc.)     
c          2 = parameter outside of tabulated values in QKSPDF
c          3 = too many data points for internal array
c          4 = error computing the likelihood, e.g. pdf(x(i))=0 at some
c              data point x(i) (can occur during search when alpha<1,beta=+1/-1.
c		 5 = possible approximation error in qkspdf.  The current
c		 approximation to log likelihood is inaccurate for highly skewed
c		 cases.
c		 6 = possible error in sfitmleci because some parameter estimate is at/near
c			the boundary
c		 7 = alpha and/or beta near a special case and parameters were rounded 
c			See subroutine scheck for details.
c
      integer noerr,errpar,errtab,errlik,errsiz,errapr,errbnd, warnrnd
      parameter (  noerr = 0, errpar = 1, errtab = 2, errsiz = 3,
     1  errlik = 4, errapr = 5, errbnd=6, warnrnd=7)


c     Global constants
      real*8 pi,piby2,logpi2
      parameter (pi=3.141592653589793d0,piby2=1.570796326794896d0, 
     1  logpi2=0.4515827052894548d0)

	integer ntol
	parameter (ntol=11)

c     Global variables
      real*8 alpha0, beta0, zeta, zcos, theta0, theta1, a1, a2, xtol,
     1   coef, c1ab, za1ln, tol(ntol), qalpha,qbeta,qgamma,qdelta0,qp
      integer qparam     
      logical onesid, gauss, cauchy, levy, anear1, debug
      common /scom/ alpha0,beta0,zeta,zcos,theta0,theta1,a1,a2,xtol,
     1   coef, c1ab, za1ln, tol, 
     1   qalpha,qbeta,qgamma,qdelta0,qp,qparam,
     1   onesid,gauss,cauchy,levy,anear1,debug
