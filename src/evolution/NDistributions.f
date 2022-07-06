************************************************************************
*
*     NDistributions.f:
*
*     Routine that returns the evolved N-space distributions.
*
************************************************************************
      subroutine NDistributions(N,Q,fnph)
*
      implicit none
**
*     Input Variables
*
      double complex N
      double precision Q
**
*     Internal Variables
*
      double complex fnev(19)
**
*     Output Variables
*      
      double complex fnph(-9:9)
*
      call NDistributionsEv(N,Q,fnev)
*
*     Rotate to the physical basis
*      
      call evln2lhac(fnev,fnph)
*
      return
      end
************************************************************************
      subroutine NDistributionsEv(N,Q,fn)
      implicit none
*      
      include "../commons/massthrs.h"      
**
*     Input Variables
*
      double complex N
      double precision Q
**      
*     Internal Variables
*
      integer nfi,nff
      double precision q2
*      
*     Output Variables
*
      double complex fn(19)
*      
*     The initial scale is always the electron mass
      nfi = 1
*
      q2 = Q**2d0      
      do nff = 11, 1, -1
         if (q2.ge.q2th(nff)) exit
      enddo
*
      call NDistributionsNF(N,nfi,nff,q2,fn)
*
      return
      end
************************************************************************      
      subroutine NDistributionsNF(N,nfi,nff,q2,fn)
*
      implicit none
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
      include "../commons/renscheme.h"
      include "../commons/tecparam.h"
      include "../commons/ns.h"      
**
*     Input Variables
*
      integer nfi,nff
      double complex N
      double precision q2
**
*     Internal Variables
*
      integer ifl
      double precision mu2,mu20
      double precision aQED
      double precision aq2,amu20,amu2
**
*     Output Variables
*
      double complex fn(19)
*
      double precision rangeA,rangeQ
      common/nstepMELA/rangeA,rangeQ
*      
*     Call N-space distributions
*     
      call electronPDFsn(N, fn) ! Electron PDFs
*
*     Run over sub-intervals
*
      aq2 = aQED(q2)
*
*     Adjust the number of iterations in the path ordering
*     (if there are thresholds, we need more points where
*     the range of variation is larger).
*     No less than NMINSTEP points for step.
      rangeA = aq2-ath(nfi)
      rangeQ = dlog(q2)-dlog(q2th(nfi))
*
*      write(*,*)"*************************************************"
*      write(*,*)"q2th",q2th
*      write(*,*) "nff",nff
*      write(*,*) "rangeQ",rangeQ
*      write(*,*) "rangeA",rangeA
*     Evolve from q2th(nfi) to q2th(nff)
      do ifl = nfi, nff-1         
         amu20 = ath(ifl)
         amu2  = ath(ifl+1)                  
         mu20  = q2th(ifl)
         mu2   = q2th(ifl+1)
         call nevolve(N,mu20,amu20,mu2,amu2,ifl,fn)
*        write(*,*) "Evolving from mu0=",sqrt(mu20),"to mu=",sqrt(mu2)
      enddo
*      
*     Add evolution between q2th(nff) and q2
      amu20 = ath(nff)
      amu2  = aq2
      mu20  = q2th(nff)
      mu2   = q2
      call nevolve(N,mu20,amu20,mu2,amu2,nff,fn)
*      write(*,*) "Evolving from mu0=",sqrt(mu20),"to mu=",sqrt(mu2)
*      write(*,*)"*************************************************"      
*
      return
      end
************************************************************************
      subroutine nevolve(N,mu20,amu20,mu2,amu2,nf,fn)
      implicit none
*      
      include "../commons/renscheme.h"
      include "../commons/tecparam.h"      
**
*     Input Variables
*
      double complex N
      double precision mu20,amu20
      double precision mu2,amu2
      integer nf
**
*     Internal Variables
*
      double precision perc
      integer npstep      
      double complex evf(19,19)      
**
*     Input/Output Variables
*      
      double complex fn(19)      
*
      double precision rangeA,rangeQ
      common/nstepMELA/rangeA,rangeQ
*      
      if(renscheme.eq."MSBAR")then
         perc  = dble(nphotot) * (amu2-amu20)/rangeA
      else
         perc  = dble(nphotot) * (dlog(mu2)-dlog(mu20))/rangeQ
      endif
      if (perc.lt.nphomin) then
         npstep = nphomin
      else
         npstep = int(perc)
      endif
*      write(*,*) "perc at nf",nf,"is",perc
*      write(*,*) "npstep at nf",nf,"is",npstep      
*         
      if(renscheme.eq."MSBAR")then
         call path_ordering(N,amu20,amu2,nf,evf,npstep)
      elseif(renscheme.eq."FIXED")then
         call alpha_fixed(N,mu20,mu2,nf,evf)
c         write(*,*) "q0=",sqrt(mu20),"q=",sqrt(mu2),evf
      elseif((renscheme.eq."ALPMZ").or.(renscheme.eq."ALGMU")
     .        .or.(renscheme.eq."AFAKE"))then
         if(afixsol.eq."PATHOR")then
            call alphaXX_pathordering(N,mu20,mu2,nf,evf,npstep)
         elseif(afixsol.eq."MAGNUS")then
            call alphaXX_magnus(N,mu20,mu2,nf,evf)
         endif
      endif
*      
      call mvmult(evf,19,19,fn,19)
*
      return
      end
