***********************************************************************
*
*     alpha_fixed.f:
*
*     This subroutine returns the QED evolution factors, in N space, from
*     Q2I to Q2F with NF active flavours, for Non Singlet and Singlet-Gluon
*     when alpha is fixed with value AREF in commons/alpha.h.
*
***********************************************************************
      SUBROUTINE ALPHA_FIXED(ZN,Q2I,Q2F,NF,EVF)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/activeflavours.h"
      include "../commons/facscheme.h"
      include "../commons/tecparam.h"
      include "../commons/alpha.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION Q2I,Q2F
      DOUBLE COMPLEX ZN
**
*     Internal Variables
*
      INTEGER I,J,K,NA
      INTEGER MP(0:3)
      DOUBLE PRECISION LN2
      DOUBLE PRECISION AFPI
      DOUBLE COMPLEX G0(4,4),G1(4,4)
      DOUBLE COMPLEX G0NS(3),G1NS(2,3)
      DOUBLE COMPLEX JLL,JGL,FDEN,FDGMN
      DOUBLE COMPLEX JM(4,4),JS(4,4),JSINV(4,4)
      DOUBLE COMPLEX SPSG(4,4),SPNS(2,3)
      DOUBLE COMPLEX SGTMP1(4,4)
      DOUBLE COMPLEX SGTMP2(4,4)
      DOUBLE COMPLEX SG(4,4)
      DOUBLE COMPLEX EFNS(2,3),EFSG(4,4)
      DATA MP / 1, 2, 8, 14 /
**
*     Output Variables
*
      DOUBLE COMPLEX EVF(19,19)
*
*     LO splitting functions
*
      CALL ANDIM_LO(ZN,NF,G0NS,G0)
*
*     NLO splitting functions
*
      IF(IPT.GE.1)THEN
         CALL ANDIM_NLO(ZN,NF,G1NS,G1)
      ENDIF
*
*     Log of the scales
      LN2 = DLOG(Q2F/Q2I)
*
*     Our splitting functions require alpha/(4*pi)
      AFPI = AREF/4D0/PI
*      
*     Singlet
*
*     Solution at LO
*     
      IF(IPT.EQ.0)THEN
         DO I=1,4
            DO J=1,4
               SPSG(I,J) = AFPI * G0(I,J) * LN2
            ENDDO
         ENDDO
*
*     Now we need to exponentiate SPSG that is a 4x4 matrix. In this
*     case, analytical formulas are too involved, therefore we adopt an
*     expanded approach truncating the series to the first NMATEXP+1 terms.
*
         CALL MATRIXEXP(NMATEXP,4,SPSG,EFSG)
*     
*     Non Singlet
*     
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = AFPI * G0NS(J) * LN2
            ENDDO
         ENDDO
*     
*     Solution at NLO
*
      ELSE
         DO I=1,4
            DO J=1,4
               EFSG(I,J) = (0D0,0D0)
               JM(I,J)   = (0D0,0D0)
            ENDDO
            EFSG(I,I) = (1D0,0D0)
         ENDDO
*
*     Delta-scheme corrections
*
         JLL = (0D0, 0D0)
         JGL = (0D0, 0D0)
         IF (FACSCHEME.EQ."DELTA") THEN
            JLL = FDEN(ZN)
            JGL = FDGMN(ZN)
            JM(1,2) = JGL
            JM(2,2) = JLL
         ENDIF
*     
         DO I=1,4
            DO J=1,4
               SGTMP1(I,J) = G0(I,J) + AFPI * G1(I,J)
               JS(I,J) = AFPI * JM(I,J)
            ENDDO
            JS(I,I) = 1D0 + JS(I,I)
         ENDDO
*
         DO I=1,4
            DO J=1,4
               JSINV(I,J) = (0D0, 0D0)
            ENDDO
            JSINV(I,I) = (1D0, 0D0)
         ENDDO
         JSINV(1,2) = JSINV(1,2) - AFPI * JGL / ( 1D0 + AFPI * JLL )
         JSINV(2,2) = JSINV(2,2) / ( 1D0 + AFPI * JLL )
*
         CALL MMULT(JS,4,4,SGTMP1,4,4,SGTMP2)
         CALL MMULT(SGTMP2,4,4,JSINV,4,4,SG)
*
         DO I=1,4
            DO J=1,4
               SPSG(I,J) = AFPI * SG(I,J) * LN2
            ENDDO
         ENDDO
*
*     Now we need to exponentiate SPSG that is a 4x4 matrix. In this
*     case, analytical formulas are too involved, therefore we adopt an
*     expanded approach truncating the series to the first NMATEXP+1 terms.
*
         CALL MATRIXEXP(NMATEXP,4,SPSG,EFSG)
*
*     
*     Non Singlet (equal in MSbar and Delta scheme for alpha fixed)
*     
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = AFPI * ( G0NS(J) + AFPI * G1NS(I,J) ) * LN2
            ENDDO
         ENDDO
*         
      ENDIF
*
*     Exponentiate the non-singlet (here because same at LO and NLO)
*      
      DO I=1,2
         DO J=1,3
            EFNS(I,J) = ZEXP( SPNS(I,J) )
         ENDDO
      ENDDO
*     
*     Contruct evolution matrix according to nf
*
*     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     g Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
*     Initialise to zero
*     
      DO I=1,19
         DO J=1,19
            EVF(I,J) = (0D0, 0D0)
         ENDDO
      ENDDO
*
*     Singlet matrix
*
      DO I=1,4
         DO J=1,4
            EVF(MP(I-1),MP(J-1)) = EFSG(I,J)
         ENDDO
      ENDDO
*
*     Total valence
*
      DO I=1,3
         EVF(MP(I)+3,MP(I)+3) = EFNS(2,I)
      ENDDO
*
      DO K=1,3
         IF (K.EQ.1) NA = NL(NF)
         IF (K.EQ.2) NA = NU(NF)
         IF (K.EQ.3) NA = ND(NF)
         DO I=1,2
            IF (NA.GT.I) THEN
               EVF(MP(K)+I,  MP(K)+I)   = EFNS(1,K)
               EVF(MP(K)+I+3,MP(K)+I+3) = EFNS(2,K)
            ELSE
               DO J=1,4
                  EVF(MP(K)+I,MP(J-1)) = EFSG(K+1,J)
               ENDDO
               EVF(MP(K)+I+3,MP(K)+3)  = EFNS(2,K)
            ENDIF
         ENDDO
      ENDDO
*
      RETURN
      END
 
