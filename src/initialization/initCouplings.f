************************************************************************
*
*     initCouplings.f:
*
*     Initialize the values of alpha at the initial scale Q20 and at 
*     the lepton thresholds.
*
************************************************************************
      subroutine initCouplings
*
      implicit none
*
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
**
*     Internal Variables
*
      integer i
      double precision aQED
*
      do i=1,11
         ath(i) = aQED(q2th(i))
      enddo
*
      return
      end
