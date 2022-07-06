*************************************************************
*
*     File: andim_lo.f
*     
*     Returns the LO anomalous dimensions in the 
*     a=alpha/4pi notation
*
*     Input variables: N  : Mellin variable
*                      NF : Number of flavours
*
*************************************************************
      SUBROUTINE ANDIM_LO(N,NF,P0NS,P0SG)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/charges.h"
      include "../commons/activeflavours.h"
      include "../commons/nfsum.h"      
*
*     ---------------------------------------------------------------------
*
*     Internal variables
*
      DOUBLE PRECISION NFS2
      DOUBLE PRECISION EL2T,EU2T,ED2T
      DOUBLE COMPLEX NS,N1,N2,NM
      DOUBLE COMPLEX PSI,S1
      DOUBLE COMPLEX PFFA,PFGA,PGFA,PGGB
*
*     ---------------------------------------------------------------------
*     
*     Input variables
*     
      DOUBLE COMPLEX N
      INTEGER NF
*     
*     Output variables  
*     
      DOUBLE COMPLEX P0NS(3)
      DOUBLE COMPLEX P0SG(4,4)
*     
*     ---------------------------------------------------------------------
*     
      NS = N * N
      N1 = N + 1D0
      N2 = N + 2D0
      NM = N - 1D0
*     
      S1 = EMC + PSI(N1)
*     
      PFFA = 3D0 - 4D0 * S1 + 2D0 / ( N * N1 )
      PFGA = 2D0 * ( NS + N + 2D0 ) / ( N * N1 * N2 )
      PGFA = 2D0 * ( NS + N + 2D0 ) / ( N * N1 * NM )
      PGGB = - 4D0 / 3D0
*
*     Sum of charges
*
      NFS2 = NFSUM2(NF)
*
      EL2T = EL2
      EU2T = EU2
      ED2T = ED2
      IF (NL(NF).EQ.0) EL2T = 0D0
      IF (NU(NF).EQ.0) EU2T = 0D0
      IF (ND(NF).EQ.0) ED2T = 0D0
*
*     Non-singlet
*
      P0NS(1) = EL2T * PFFA
      P0NS(2) = EU2T * PFFA
      P0NS(3) = ED2T * PFFA
*
*     Singlet
*
      P0SG(1,1) = NFS2 * PGGB
      P0SG(1,2) = EL2T * PGFA
      P0SG(1,3) = EU2T * PGFA
      P0SG(1,4) = ED2T * PGFA

      P0SG(2,1) = 2D0 * NL(NF) * EL2T * PFGA
      P0SG(2,2) = P0NS(1)
      P0SG(2,3) = (0D0, 0D0)
      P0SG(2,4) = (0D0, 0D0)

      P0SG(3,1) = 2D0 * NC * NU(NF) * EU2T * PFGA
      P0SG(3,2) = (0D0, 0D0)
      P0SG(3,3) = P0NS(2)
      P0SG(3,4) = (0D0, 0D0)

      P0SG(4,1) = 2D0 * NC * ND(NF) * ED2T * PFGA
      P0SG(4,2) = (0D0, 0D0)
      P0SG(4,3) = (0D0, 0D0)
      P0SG(4,4) = P0NS(3)
*     
      RETURN
      END
*     
*----------------------------------------------------------------------

