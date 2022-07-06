************************************************************************
*
*     initialconditions.f:
*
*     These routines return PDFs at the intitial scale in N space
*     assuming that the initial scale is the mass of the electron.
*
************************************************************************
      subroutine electronPDFsn(N, npdf)
*
      implicit none
*
      include "../commons/alpha.h"
      include "../commons/ipt.h"
      include "../commons/facscheme.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      integer ipdf
      double complex deN, dgmN
      double complex fdeN, fdgmN
      double complex sg, gm
      double precision amu0
**
*     Output Variables
*
      double complex npdf(19)
*
*     Set the value of alpha(mu0)
*
      amu0 = ath(1)
*      
      deN  = (0d0, 0d0)
      dgmN = (0d0, 0d0)
      if (facscheme.eq."MSBAR") then
         deN  = - fdeN(N)
         dgmN = - fdgmN(N)
      endif
*
      sg = 1d0
      gm = 0d0
      if (ipt.gt.0) then
         sg = sg + amu0 * deN
         gm = gm + amu0 * dgmN
      endif
*
*     Set electron and photon PDFs in the evolution basis
*
*      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     gm Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
      npdf(1) = gm
      do ipdf = 2, 7
         npdf(ipdf) = sg
      enddo
*
*     Set the rest to zero
*      
      do ipdf = 8, 19
         npdf(ipdf) = (0d0,0d0)
      enddo
*
      return
      end
*
************************************************************************
      function fdeN(N)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      double complex Nm, N1
      double complex s1, s2
      double complex psi, dpsi
**
*     Output Variables
*
      double complex fdeN
*
      Nm = N - 1d0
      N1 = N + 1d0
      s1 = emc + psi(N1)
      s2 = zeta2 - dpsi(N1, 1)
*
      fdeN = - 2d0 * ( - ( 1d0 / N - 1d0 / N1 - 2d0 * s1 + 3d0 / 2d0 )
     1     - 2d0 * ( s1**2 + s2 - s1 / N + s1 / N1 + 1d0 / N1**2
     2     - 7d0 / 4d0 ) )
*
      return
      end
*
************************************************************************
      function fdgmN(N)
*
      implicit none
**
*     Input Variables
*
      double complex N
**
*     Internal Variables
*
      double complex Nm, N1
**
*     Output Variables
*
      double complex fdgmN
*
      Nm = N - 1d0
      N1 = N + 1d0
*
      fdgmN = - 2d0 * ( - ( 2d0 / Nm - 2d0 / N + 1d0 / N1 )
     1     - 2d0 * ( - 2d0 / Nm**2 + 2d0 / N**2 - 1d0 / N1**2 ) )
*
      return
      end
