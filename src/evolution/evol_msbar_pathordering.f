***********************************************************************
*
*     path_ordering.f:
*
*     This subroutine returns the QED evolution factors, in N space, from
*     Q2I to Q2F with NF active flavours, for Non Singlet and Singlet-Gluon
*     using the path ordering method.
*
***********************************************************************
      SUBROUTINE PATH_ORDERING(ZN,AII,AFF,NF,EVF,NINTS)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/beta.h"
      include "../commons/activeflavours.h"
      include "../commons/facscheme.h"
      include "../commons/ns.h"
      include "../commons/tecparam.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION AII,AFF
      DOUBLE COMPLEX ZN
      INTEGER NINTS
**
*     Internal Variables
*
      INTEGER I,J,K,NA
      INTEGER MP(0:3)
      DOUBLE PRECISION BT0,BT1
      DOUBLE PRECISION AF,AI,DA,AK
      DOUBLE COMPLEX G0(4,4),G1(4,4)
      DOUBLE COMPLEX G0NS(3),G1NS(2,3)
      DOUBLE COMPLEX JLL,JGL,FDEN,FDGMN
      DOUBLE COMPLEX JM(4,4),JSI(4,4),JSF(4,4),JSIINV(4,4),JSFINV(4,4)
      DOUBLE COMPLEX SPSG(4,4),SPNS(2,3)
      DOUBLE COMPLEX SGITMP1(4,4),SGFTMP1(4,4)
      DOUBLE COMPLEX SGITMP2(4,4),SGFTMP2(4,4)
      DOUBLE COMPLEX SGI(4,4),SGF(4,4)
      DOUBLE COMPLEX DEFSG(4,4)
      DOUBLE COMPLEX EFNS(2,3),EFSG(4,4),EFSGTMP(4,4)
      DATA MP / 1, 2, 8, 14 /
**
*     Output Variables
*
      DOUBLE COMPLEX EVF(19,19)
*
      BT0 = BETA0(NF)
      BT1 = 0D0
*
*     LO splitting functions
*
      CALL ANDIM_LO(ZN,NF,G0NS,G0)
*
*     NLO splitting functions
*
      IF(IPT.GE.1)THEN
         CALL ANDIM_NLO(ZN,NF,G1NS,G1)
         BT1 = BETA1(NF)
c$$$*     The B1 coefficient is different from zero,
c$$$*     only when alpha is running at two-loops
c$$$         IF(IPTALPHA.GT.0)THEN
c$$$            BT1 = BETA1(NF)
c$$$         ENDIF
      ENDIF
*
*     Singlet
*
*     Solution at LO
*
      DA = ( AFF - AII ) / DBLE(NINTS)
*
      IF(IPT.EQ.0)THEN
*
         DO I=1,4
            DO J=1,4
               EFSG(I,J) = (0D0,0D0)
            ENDDO
            EFSG(I,I) = (1D0,0D0)
         ENDDO
*     
         AI = AII
         DO K=1,NINTS
            AF = AI + DA
            DO I=1,4
               DO J=1,4
                  SGI(I,J) = - G0(I,J) / AI / BT0
                  SGF(I,J) = - G0(I,J) / AF / BT0
                  SPSG(I,J) = DA * ( SGI(I,J) + SGF(I,J) ) / 2D0
               ENDDO
            ENDDO
*     
            CALL MATRIXEXP(NMATEXP,4,SPSG,DEFSG)
            CALL MMULT(DEFSG,4,4,EFSG,4,4,EFSGTMP)
            CALL EQUATEM(EFSGTMP,4,4,EFSG)
*     
            AI = AI + DA
         ENDDO                  
*     
*     Non Singlet
*     
         AI = AII
         AF = AFF
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = - DA * ( G0NS(J) / BT0 / AI         
     1              + G0NS(J) / BT0 / AF ) / 2D0
            ENDDO
         ENDDO
*     
         AK = AI
         DO K=1,NINTS-1
            AK = AK + DA
*     
            DO I=1,2
               DO J=1,3
                  SPNS(I,J) = SPNS(I,J)
     1                 - DA * G0NS(J) / BT0 / AK
               ENDDO
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
         AI = AII
         DO K=1,NINTS
            AF = AI + DA
            DO I=1,4
               DO J=1,4
                  SGITMP1(I,J) = - ( G0(I,J) + AI * ( G1(I,J)
     1                 - G0(I,J) * BT1 / BT0 ) ) / BT0 / AI
                  SGFTMP1(I,J) = - ( G0(I,J) + AF * ( G1(I,J)
     1                 - G0(I,J) * BT1 / BT0 ) ) / BT0 / AF
                  JSI(I,J) = AI * JM(I,J)
                  JSF(I,J) = AF * JM(I,J)
               ENDDO
               JSI(I,I) = 1D0 + JSI(I,I)
               JSF(I,I) = 1D0 + JSF(I,I)
            ENDDO
*
            DO I=1,4
               DO J=1,4
                  JSIINV(I,J) = (0D0, 0D0)
                  JSFINV(I,J) = (0D0, 0D0)
               ENDDO
               JSIINV(I,I) = (1D0, 0D0)
               JSFINV(I,I) = (1D0, 0D0)
            ENDDO
            JSIINV(1,2) = JSIINV(1,2) - AI * JGL / ( 1D0 + AI * JLL )
            JSIINV(2,2) = JSIINV(2,2) / ( 1D0 + AI * JLL )
            JSFINV(1,2) = JSFINV(1,2) - AF * JGL / ( 1D0 + AF * JLL )
            JSFINV(2,2) = JSFINV(2,2) / ( 1D0 + AF * JLL )
*
            CALL MMULT(JSI,4,4,SGITMP1,4,4,SGITMP2)
            CALL MMULT(SGITMP2,4,4,JSIINV,4,4,SGI)
            CALL MMULT(JSF,4,4,SGFTMP1,4,4,SGFTMP2)
            CALL MMULT(SGFTMP2,4,4,JSFINV,4,4,SGF)
*
            DO I=1,4
               DO J=1,4
                  SGI(I,J) = SGI(I,J) + JM(I,J) / ( 1D0 + AI * JLL )
                  SGF(I,J) = SGF(I,J) + JM(I,J) / ( 1D0 + AF * JLL )
                  SPSG(I,J) = DA * ( SGI(I,J) + SGF(I,J) ) / 2D0
               ENDDO
            ENDDO
*
*     Now we need to exponentiate SPSG that is a 4x4 matrix. In this
*     case, analytical formulas are too involved, therefore we adopt an
*     expanded approach truncating the series to the first NMATEXP+1 terms.
*
            CALL MATRIXEXP(NMATEXP,4,SPSG,DEFSG)
            CALL MMULT(DEFSG,4,4,EFSG,4,4,EFSGTMP)
            CALL EQUATEM(EFSGTMP,4,4,EFSG)
*
            AI = AI + DA
         ENDDO
*     
*     Non Singlet
*     
         AI = AII
         AF = AFF
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = DA * (
     1              JLL / ( 1D0 + AI * JLL)
     2              - ( G0NS(J) + AI * ( G1NS(I,J)
     3              - G0NS(J) * BT1/BT0 ) ) / BT0 / AI +
     4              JLL / ( 1D0 + AF * JLL)
     5              - ( G0NS(J) + AF * ( G1NS(I,J)
     6              - G0NS(J) * BT1/BT0 ) ) / BT0 / AF
     7              ) / 2D0
            ENDDO
         ENDDO
*     
         AK = AI
         DO K=1,NINTS-1
            AK = AK + DA
            DO I=1,2
               DO J=1,3
                  SPNS(I,J) = SPNS(I,J)
     1                 + DA * ( JLL / ( 1D0 + AK * JLL) - ( G0NS(J)
     2                 + AK * ( G1NS(I,J) - G0NS(J) * BT1/BT0 ) )
     3                 / BT0 / AK )
               ENDDO
            ENDDO
         ENDDO
      ENDIF
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
 
