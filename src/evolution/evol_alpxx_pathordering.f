***********************************************************************
*
***********************************************************************
      SUBROUTINE ALPHAXX_PATHORDERING(ZN,Q2I,Q2F,NF,EVF,NINTS)
*
      IMPLICIT NONE
*
      include "../commons/ipt.h"
      include "../commons/beta.h"
      include "../commons/activeflavours.h"
      include "../commons/facscheme.h"
      include "../commons/ns.h"
      include "../commons/tecparam.h"
      include "../commons/consts.h"
      include "../commons/alpha.h"
      include "../commons/massthrs.h"
**
*     Input Variables
*
      INTEGER NF
      DOUBLE PRECISION Q2I,Q2F
      DOUBLE COMPLEX ZN
      INTEGER NINTS     
**
*     Internal Variables
*
      INTEGER I,J,K,NA
      INTEGER MP(0:3)
      DOUBLE PRECISION AONE2PI
      DOUBLE PRECISION TTI,TTF,DTT,TTK,TMID
      DOUBLE PRECISION LN2
      DOUBLE PRECISION B0, DTILDE, MF2
      DOUBLE PRECISION DK
      DOUBLE PRECISION BT0
      DOUBLE PRECISION MZ2
      DOUBLE COMPLEX VERYTMP
      DOUBLE COMPLEX G0(4,4),G1(4,4)
      DOUBLE COMPLEX G0NS(3),G1NS(2,3)
      DOUBLE COMPLEX JLL,JGL,FDEN,FDGMN
      DOUBLE COMPLEX JM(4,4),JS(4,4),JSINV(4,4)
      DOUBLE COMPLEX SPSG(4,4),SPNS(2,3)
      DOUBLE COMPLEX SPNSAN(2,3)!,SPNSAN2(2,3)
      DOUBLE COMPLEX SGTMP1(4,4),SGTMP2(4,4),SG(4,4)
c      DOUBLE COMPLEX SGITMP1(4,4),SGITMP2(4,4),SGI(4,4)
c      DOUBLE COMPLEX SGFTMP1(4,4),SGFTMP2(4,4),SGF(4,4)
      DOUBLE COMPLEX DEFSG(4,4)
      DOUBLE COMPLEX EFNS(2,3),EFSG(4,4),EFSGTMP(4,4)
      DATA MP / 1, 2, 8, 14 /
***
*     Output Variables
*
      DOUBLE COMPLEX EVF(19,19)
*
      BT0 = BETA0(NF)
      B0 = - BT0 / 4d0 / PI
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
      AONE2PI = AREF / 2D0 / PI
      LN2 = DLOG(Q2F/Q2I)
*
*     Solution at LO
*
      IF(IPT.EQ.0)THEN
*     
*     Singlet
*     
         DO I=1,4
            DO J=1,4
               SPSG(I,J) = AONE2PI * G0(I,J) / 2D0 * LN2
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
               SPNS(I,J) = AONE2PI * G0NS(J) / 2D0 * LN2
            ENDDO
         ENDDO
*
         DO I=1,2
            DO J=1,3
               EFNS(I,J) = ZEXP( SPNS(I,J) )
            ENDDO
         ENDDO
*         
*     Solution at NLO
*
      ELSE
*
*     Initialize evolution operator to unity
*
         DO I=1,4
            DO J=1,4
               EFSG(I,J) = (0D0,0D0)
            ENDDO
            EFSG(I,I) = (1D0,0D0)
         ENDDO
*     
*     Delta-scheme corrections
*     
         DO I=1,4
            DO J=1,4
               JM(I,J) = (0D0,0D0)
            ENDDO
         ENDDO
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
*     Path ordering step
*
         MZ2 = Q2TH(10)
         IF (NF.GE.10) THEN
            MF2 = MZ2
         ELSEIF (NF.LT.10) THEN
            MF2 = Q2TH(NF+1)
         ENDIF
*      
         TTF = AONE2PI * DLOG(Q2F/MF2)
         TTI = AONE2PI * DLOG(Q2I/MF2)
*         TTF = AONE2PI * DLOG(Q2F/MZ2)
*         TTI = AONE2PI * DLOG(Q2I/MZ2)         
         DTT = ( TTF - TTI ) / DBLE(NINTS)
*         WRITE(6,*)"TTI",TTI,"TTF",TTF,"DTT",DTT
*
         CALL GETDK(NF,DK)
         DTILDE = DK
*        DTILDE = DK - 2D0*PI*B0 * DLOG(MF2/MZ2)
*        WRITE(6,*)"DTILDE",DTILDE
*
         DO K=1,NINTS
            TTF = TTI + DTT
            TMID = ( TTI + TTF ) / 2D0
            DO I=1,4
               DO J=1,4
                  SGTMP1(I,J) = G0(I,J) / 2D0 * (
     1                 1D0 + AONE2PI * DTILDE
     2                 + 2D0*PI*B0 * TMID )
     3                 + AONE2PI * ( G1(I,J) / 4D0 )
c                  SGITMP1(I,J) = G0(I,J) / 2D0 * (
c     1                 1D0 + AONE2PI * DTILDE
c     2                 + 2D0*PI*B0 * TTI )
c     3                 + AONE2PI * ( G1(I,J) / 4D0 )
c                  SGFTMP1(I,J) = G0(I,J) / 2D0 * (
c     1                 1D0 + AONE2PI * DTILDE
c     2                 + 2D0*PI*B0 * TTF )
c     3                 + AONE2PI * ( G1(I,J) / 4D0 )                                    
                  JS(I,J) = AONE2PI/2D0 * JM(I,J)
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
            JSINV(1,2) = JSINV(1,2)
     1           - AONE2PI/2D0 * JGL / ( 1D0 + AONE2PI/2D0 * JLL )
            JSINV(2,2) = JSINV(2,2) / ( 1D0 + AONE2PI/2D0 * JLL )
*
            CALL MMULT(JS,4,4,SGTMP1,4,4,SGTMP2)
            CALL MMULT(SGTMP2,4,4,JSINV,4,4,SG)
c            CALL MMULT(JS,4,4,SGITMP1,4,4,SGITMP2)
c            CALL MMULT(SGITMP2,4,4,JSINV,4,4,SGI)
c            CALL MMULT(JS,4,4,SGFTMP1,4,4,SGFTMP2)
c            CALL MMULT(SGFTMP2,4,4,JSINV,4,4,SGF)            
*
            DO I=1,4
               DO J=1,4
                  SPSG(I,J) = DTT * SG(I,J)
c                  SPSG(I,J) = DTT * ( SGI(I,J) + SGF(I,J) ) / 2D0
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
            TTI = TTI + DTT
         ENDDO
*     
*     Non Singlet (equal for Delta and MSbar factorization scheme)
*
*     reset TTF and TTI
         TTF = AONE2PI * DLOG(Q2F/MF2)
         TTI = AONE2PI * DLOG(Q2I/MF2)         
c$$$         TTI = AONE2PI * DLOG(Q2I/MZ2)         
c$$$         TTF = AONE2PI * DLOG(Q2F/MZ2)
*         
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = DTT * ( G0NS(J) / 2D0 * (
     1              1D0 + AONE2PI * DTILDE )
     2              + AONE2PI * ( G1NS(I,J) / 4D0 )
     3              + PI*B0 * (TTI+TTF) * G0NS(J) / 2D0 )
            ENDDO
         ENDDO
*     
         TTK = TTI
         DO K=1,NINTS-1
            TTK = TTK + DTT
            DO I=1,2
               DO J=1,3
                  VERYTMP = DTT * (
     1                 G0NS(J) / 2D0 * (
     2                 1D0 + AONE2PI * DTILDE
     3                 + 2D0*PI*B0 * TTK )
     4                 + AONE2PI * ( G1NS(I,J) / 4D0 ) )
*                  IF ((I.EQ.1).AND.(J.EQ.1)) THEN
*                     WRITE(6,*)"TTK",TTK,VERYTMP
*                  ENDIF
                  SPNS(I,J) = SPNS(I,J) + VERYTMP
               ENDDO
            ENDDO
         ENDDO
*
*     Write the analytic solution for the non-singlet ev.op. 
         DO I=1,2
            DO J=1,3
               SPNSAN(I,J) = ( ( AONE2PI * G0NS(J)/2D0
     1              + AONE2PI**2 *
     2              ( G1NS(I,J)/4D0 + DK*G0NS(J)/2D0 ) ) * DLOG(Q2F/Q2I)
     3              + AONE2PI**2 * PI*B0*G0NS(J)/2D0 *
     4              (DLOG(Q2F/MF2)**2 - DLOG(Q2I/MF2)**2) )            
            ENDDO
         ENDDO
**
*     Same as above, just for check
*      DO I=1,2
*         DO J=1,3
*            SPNSAN2(I,J) = ( AONE2PI * DLOG(Q2F/Q2I) *
*     1           ( G0NS(J)/2D0 + AONE2PI * (
*     2           G1NS(I,J)/4D0 + DTILDE * G0NS(J)/2D0 ) )
*     3           + AONE2PI**2 * PI*B0*G0NS(J)/2D0 *
*     4           (DLOG(Q2F/MZ2)**2 - DLOG(Q2I/MZ2)**2) )
*         ENDDO
*      ENDDO
*      
*      WRITE(6,*)"SPNS",SPNS
*     WRITE(6,*)"SPNSAN",SPNSAN
*      WRITE(6,*)(SPNS-SPNSAN)/SPNSAN
*      WRITE(6,*)(SPNSAN2-SPNSAN)/SPNSAN
*
*     Try to use the analytic solution
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = SPNSAN(I,J)
            ENDDO
         ENDDO
*     
         DO I=1,2
            DO J=1,3
               EFNS(I,J) = ZEXP( SPNS(I,J) )
            ENDDO
         ENDDO
*         
      ENDIF      
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
************************************************************************
*     calculate d^{(k)}, needed in the alphamz scheme
      subroutine getdk(k,dk)
      implicit none
      include "../commons/consts.h"
      include "../commons/massthrs.h"
      include "../commons/beta.h"
      include "../commons/nfsum.h"
      include "../commons/activeflavours.h"
      include "../commons/renscheme.h"      
*     input
      integer k
*     output
      double precision dk
*     internal
      integer i
      double precision b0
      double precision mz2,mw2
      double precision deltaGMU
      double precision wbos
*
      mz2 = q2th(10)
      mw2 = q2th(9)
*
*     FORMAL REPLACEMENT TO HAVE GMU FROM ALPHA(MZ):
*        5/(9*Pi)*C2M + WBOSON -> DELTA(MSb->GMU) 
*     or
*        10/9*C2M -> 2*Pi*DELTA(MSb->GMU)
*
      if(RENSCHEME.eq."ALPMZ")then
         wbos = 1d0/pi*(7d0/4d0*dlog(mw2/mz2)-1d0/6d0)*waem
         dk = 2d0*pi * (5d0/9d0/pi * nfsum2(10) + wbos)
      endif
c      
      if(RENSCHEME.eq."ALGMU")then
cccccccccccccccccccccccccccccccccccccccccccccccccc         
c     WARNING: ASSOCIATED TO VALUE OF ALPHA = 1/132.18289853516
         deltaGMU = 4.09293586102687773d0
cccccccccccccccccccccccccccccccccccccccccccccccccc         
      endif
      if(RENSCHEME.eq."AFAKE")then
cccccccccccccccccccccccccccccccccccccccccccccccccc         
c     WARNING: ASSOCIATED TO VALUE OF ALPHA = 1/130.50453291347856
         deltaGMU = 2.5d0
cccccccccccccccccccccccccccccccccccccccccccccccccc         
      endif      
      if((RENSCHEME.eq."ALGMU").or.(RENSCHEME.eq."AFAKE"))then
         dk = 2d0*pi * deltaGMU
      endif
*      
      if (k.ge.9) then
         return
      endif
*
*     ADDING PERTURBATIVE EFFECTS
      do i=k+1,9
         b0 = -beta0(i)/(4d0*pi)
         dk = dk + 2d0*pi*b0*dlog(q2th(i)/q2th(i+1))
*         write(6,*) "i",i,"dlog of",q2th(i),"over",q2th(i+1),dk,"b0",b0
      enddo
*     
      return
      end
