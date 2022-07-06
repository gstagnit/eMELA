************************************************************************
*
*     xDistributions.f:
*
*     Routine that returns the evolved x-space distributions.
*
************************************************************************
      subroutine xDistributions(x,Q,xfph)
*
      implicit none
**
*     Input Variables
*
      double precision x
      double precision Q
**
*     Internal Variables
*
      double precision xfev(19)
**
*     Output Variables
*      
      double precision xfph(-9:9)
*
      call xDistributionsEv(x,Q,xfev)
*
*     Rotate to the physical basis
*      
      call evln2lha(xfev,xfph)
*
      return
      end
************************************************************************ 
      subroutine xDistributionsEv(x,Q,xfev)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/ns.h"
      include "../commons/tecparam.h"
**
*     Input Variables
*
      double precision x
      double precision Q
**
*     Internal Variables
*
      integer nfi,nff
      integer ip
      double precision integrand
      external integrand
      integer j
      double precision r
      double precision theta,sigma,t
      double complex s,tmp(19)
      double complex xfN(19)
      double precision q2
**
*     Output Variables
*
      double precision xfev(19)
*
*     The initial scale is always the electron mass
      nfi = 1
*
      q2 = Q**2d0      
      do nff = 11, 1, -1
         if (q2.ge.q2th(nff)) exit
      enddo
*     
*     trapezioidal integration with Talbot path
      do ip=1,19
         tmp(ip) = (0d0,0d0)
      enddo      
      t = - dlog(x)
      r = 2d0 * rinvmel / 5d0 / t
      do j=1,minvmel-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(minvmel) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         s     = r * theta * dcmplx(1d0/tan(theta),1d0) + 1d0
         call NDistributionsNF(s,nfi,nff,q2,xfN)
         do ip=1,19
            tmp(ip) = tmp(ip)
     1              + exp( t * s ) * dcmplx(1d0,sigma) * xfN(ip)
         enddo
      enddo
      do ip=1,19
         xfev(ip) = x * r * dble(tmp(ip)) / minvmel
      enddo
*     
      return
      end
