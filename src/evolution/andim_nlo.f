********************************************************************************
*                                                                              *
*     andim_nlo.f                                                              *
*                                                                              *
********************************************************************************
*                                                                              *
*     Returns the NLO space-like anomalous dimensions in the a=alpha/4pi       *
*     notation                                                                 *
*                                                                              *
*     Input variables:                                                         *
*     N  : Mellin variable                                                     *
*     NF : Number of flavours                                                  *
*                                                                              *
********************************************************************************
      SUBROUTINE ANDIM_NLO(N,NF,P1NS,P1SG)
*
      IMPLICIT NONE
*
      include "../commons/consts.h"
      include "../commons/charges.h"
      include "../commons/activeflavours.h"
      include "../commons/nfsum.h"
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
      DOUBLE PRECISION NFS2
      DOUBLE PRECISION NFS4
      DOUBLE PRECISION EL2T,EU2T,ED2T
      DOUBLE PRECISION EL4T,EU4T,ED4T
      DOUBLE COMPLEX NS,NT,NFO,NFI,NSI,NSE,NE,NN
      DOUBLE COMPLEX N1,N2,NM,NMS,N1S,N1T,N2S,N2T
      DOUBLE COMPLEX N3,N4,N5,N6
      DOUBLE COMPLEX S11,S12,S13,S14,S15,S16
      DOUBLE COMPLEX SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP
      DOUBLE COMPLEX SSTR2P,SSTR3P
      DOUBLE COMPLEX PPSA,PGFA,PFGB,PGFC,PGGC
      DOUBLE COMPLEX PNPA,PNMA,PNSC
      DOUBLE COMPLEX DPSI,PSI,S1,S2
*
* ---------------------------------------------------------------------
*
*     Input variables
*
      DOUBLE COMPLEX N
      INTEGER NF
*      
*     Output variables  
*
      DOUBLE COMPLEX P1NS(2,3)
      DOUBLE COMPLEX P1SG(4,4)
*
* ---------------------------------------------------------------------
*
      S1 = EMC + PSI(N+1d0)
      S2 = ZETA2 - DPSI(N+1d0,1)
*     
      NS = N * N
      NT = NS * N
      NFO = NT * N
      NFI = NFO * N
      NSI = NFI * N
      NSE = NSI * N
      NE = NSE * N
      NN = NE * N
*     
      NM = N - 1d0
      N1 = N + 1d0
      N2 = N + 2d0
      NMS = NM * NM
      N1S = N1 * N1
      N1T = N1S * N1
      N2S = N2 * N2
      N2T = N2S * N2
*
* ---------------------------------------------------------------------
*
* ..Analytic continuations of the occuring sums as given in GRV (1990) 
*   (with an improved parametrization of the moments of  Sp(x)/(1+x).)
*
      N3 = N + 3D0
      N4 = N + 4D0
      N5 = N + 5D0
      N6 = N + 6D0
      S11 = S1  + 1D0/N1
      S12 = S11 + 1D0/N2
      S13 = S12 + 1D0/N3
      S14 = S13 + 1D0/N4
      S15 = S14 + 1D0/N5
      S16 = S15 + 1D0/N6
      SPMOM = 1.0000D0 * (ZETA2 - S1 / N ) / N  -
     1     0.9992D0 * (ZETA2 - S11/ N1) / N1 +
     2     0.9851D0 * (ZETA2 - S12/ N2) / N2 -
     3     0.9005D0 * (ZETA2 - S13/ N3) / N3 +
     4     0.6621D0 * (ZETA2 - S14/ N4) / N4 -
     5     0.3174D0 * (ZETA2 - S15/ N5) / N5 +
     6     0.0699D0 * (ZETA2 - S16/ N6) / N6  
*     
      SLC = - 5D0/8D0 * ZETA3
      SLV = - ZETA2/2D0* (PSI(N1/2D0) - PSI(N/2D0)) 
     1     + S1/NS + SPMOM
      SSCHLM = SLC - SLV
      SSTR2M = ZETA2 - DPSI (N1/2D0,1)
      SSTR3M = 0.5D0 * DPSI (N1/2D0,2) + ZETA3

      SSCHLP = SLC + SLV
      SSTR2P = ZETA2 - DPSI(N2/2D0,1)
      SSTR3P = 0.5D0 * DPSI(N2/2D0,2) + ZETA3      
*     
*     ---------------------------------------------------------------------
*
*     ..The contributions to P1NS as given in Gonzalez-Arroyo et al. (1979) 
*     (Note that the anomalous dimensions in the literature often differ 
*     from these moments of the splitting functions by factors -1 or -2,
*     in addition to possible different normalizations of the coupling)
*
      PNMA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1     16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2M ) +
     2     64D0* SSCHLM + 24D0* S2 - 3D0 - 8D0* SSTR3M -
     3     8D0* (3D0* NT + NS -1D0) / (NT * N1T) +
     4     16D0* (2D0* NS + 2D0* N +1D0) / (NT * N1T) ) * (-0.5D0)
      PNPA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1     16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2P ) +
     2     64D0* SSCHLP + 24D0* S2 - 3D0 - 8D0* SSTR3P -
     3     8D0* (3D0* NT + NS -1D0) / (NT * N1T) -
     4     16D0* (2D0* NS + 2D0* N +1D0)/(NT * N1T) ) * (-0.5D0)
      
      PNSC = ( -160D0/9D0* S1 + 32D0/3.* S2 + 4D0/3D0 +
     1     16D0*(11D0*NS+5D0*N-3D0)/(9D0* NS * N1S))*(-0.5D0)
*     
*     ---------------------------------------------------------------------
*     
*     ..The contributions to P1SG as given in Floratos et al. (1981) 
*     ..Pure singlet (PS) and FG
*     
      PPSA = (5d0* NFI + 32d0* NFO + 49d0* NT+38d0* NS + 28d0* N + 8d0) 
     1     / (NM * NT * N1T * N2S) * 2d0
*     
      PFGB = (2d0* S1 * S1 - 2d0* S2 + 5d0) * (NS + N + 2d0)
     1     / (N * N1 * N2) - 4d0* S1 / NS
     2     + (11d0* NFO + 26d0* NT + 15d0* NS + 8d0* N + 4d0)
     3     / (NT * N1T * N2)
*     
*     ---------------------------------------------------------------------
*     
*     ..GF and GG
*     
      PGFA = (- S1 * S1 + 5d0* S1 - S2) * (NS + N + 2d0) 
     1     / (NM * N * N1)  -  2d0* S1 / N1S
     2     - (12d0* NSI + 30d0* NFI + 43d0* NFO + 28d0* NT - NS
     3     - 12d0* N - 4d0) / (2d0* NM * NT * N1T) 
      PGFC = (S1 - 8d0/3d0) * (NS + N + 2d0) / (NM * N * N1) + 1d0/ N1S
      PGFC = 4d0/3d0* PGFC
*     
      PGGC = (2d0* NSI + 4d0* NFI + NFO - 10d0* NT - 5d0* NS - 4d0* N
     1     - 4d0) * (-2d0) / (NM * NT * N1T * N2)  -  1d0
*
*     Sum of charges
*
      NFS2 = NFSUM2(NF)
      NFS4 = NFSUM4(NF)

      EL2T = EL2
      EU2T = EU2
      ED2T = ED2
      IF (NL(NF).EQ.0) EL2T = 0D0
      IF (NU(NF).EQ.0) EU2T = 0D0
      IF (ND(NF).EQ.0) ED2T = 0D0

      EL4T = EL4
      EU4T = EU4
      ED4T = ED4
      IF (NL(NF).EQ.0) EL4T = 0D0
      IF (NU(NF).EQ.0) EU4T = 0D0
      IF (ND(NF).EQ.0) ED4T = 0D0
*
*     Non-singlet
*
      P1NS(1,1) = (0D0,0D0)
      P1NS(1,2) = (0D0,0D0)
      P1NS(1,3) = (0D0,0D0)
      IF (NL(NF).GT.0) P1NS(1,1) = EL4 * PNPA + EL2T * NFS2 * PNSC
      IF (NU(NF).GT.0) P1NS(1,2) = EU4 * PNPA + EU2T * NFS2 * PNSC
      IF (ND(NF).GT.0) P1NS(1,3) = ED4 * PNPA + ED2T * NFS2 * PNSC

      P1NS(2,1) = (0D0,0D0)
      P1NS(2,2) = (0D0,0D0)
      P1NS(2,3) = (0D0,0D0)
      IF (NL(NF).GT.0) P1NS(2,1) = EL4 * PNMA + EL2T * NFS2 * PNSC
      IF (NU(NF).GT.0) P1NS(2,2) = EU4 * PNMA + EU2T * NFS2 * PNSC
      IF (ND(NF).GT.0) P1NS(2,3) = ED4 * PNMA + ED2T * NFS2 * PNSC
*
*     Singlet
*
      P1SG(1,1) = 4D0 * NFS4 * PGGC
      P1SG(1,2) = 4D0 * ( EL4 * PGFA + EL2T * NFS2 * PGFC )
      P1SG(1,3) = 4D0 * ( EU4 * PGFA + EU2T * NFS2 * PGFC )
      P1SG(1,4) = 4D0 * ( ED4 * PGFA + ED2T * NFS2 * PGFC )

      P1SG(2,1) = 4D0 * NL(NF) * EL4 * PFGB
      P1SG(2,2) = 4D0 * NL(NF) * EL2T * EL2T * PPSA + P1NS(1,1)
      P1SG(2,3) = 4D0 * NL(NF) * EL2T * EU2T * PPSA
      P1SG(2,4) = 4D0 * NL(NF) * EL2T * ED2T * PPSA

      P1SG(3,1) = 4D0 * NC * NU(NF) * EU4 * PFGB
      P1SG(3,2) = 4D0 * NC * NU(NF) * EU2T * EL2T * PPSA
      P1SG(3,3) = 4D0 * NC * NU(NF) * EU2T * EU2T * PPSA + P1NS(1,2)
      P1SG(3,4) = 4D0 * NC * NU(NF) * EU2T * ED2T * PPSA

      P1SG(4,1) = 4D0 * NC * ND(NF) * ED4 * PFGB
      P1SG(4,2) = 4D0 * NC * ND(NF) * ED2T * EL2T * PPSA
      P1SG(4,3) = 4D0 * NC * ND(NF) * ED2T * EU2T * PPSA
      P1SG(4,4) = 4D0 * NC * ND(NF) * ED2T * ED2T * PPSA + P1NS(1,3)
*     
      RETURN
      END
*     
*     =================================================================av==
