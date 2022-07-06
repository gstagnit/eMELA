***********************************************************************
*
***********************************************************************
      SUBROUTINE ALPHAXX_MAGNUS(ZN,Q2I,Q2F,NF,EVF)
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
**
*     Internal Variables
*
      INTEGER I,J,K,NA
      INTEGER MP(0:3)
      DOUBLE PRECISION AONE2PI
      DOUBLE PRECISION LN2
      DOUBLE PRECISION B0, MF2
      DOUBLE PRECISION DK, DCAL
      DOUBLE PRECISION BT0
      DOUBLE PRECISION MZ2
      DOUBLE COMPLEX G0(4,4),G1(4,4)
      DOUBLE COMPLEX G0NS(3),G1NS(2,3)
      DOUBLE COMPLEX JLL,JGL,FDEN,FDGMN
      DOUBLE COMPLEX JM(4,4),JS(4,4),JSINV(4,4)
      DOUBLE COMPLEX SPSG(4,4),SPNS(2,3)
      DOUBLE COMPLEX EFNS(2,3),EFSG(4,4)
      DOUBLE COMPLEX SGATMP1(4,4),SGATMP2(4,4),SGA(4,4)      
      DOUBLE COMPLEX SGBTMP1(4,4),SGBTMP2(4,4),SGB(4,4)
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
*     Prepare for Magnus
*     
         LN2 = DLOG(Q2F/Q2I)
*
         MZ2 = Q2TH(10)         
         IF (NF.GE.10) THEN
            MF2 = MZ2
         ELSEIF (NF.LT.10) THEN
            MF2 = Q2TH(NF+1)
         ENDIF
*         
         CALL GETDK(NF,DK)         
         DCAL = DK + 2D0*PI*B0 * DLOG(Q2I/MF2)
*
         DO I=1,4
            DO J=1,4
               SGATMP1(I,J) = AONE2PI**2 * G0(I,J)/2D0 * 2D0*PI*B0
               SGBTMP1(I,J) = AONE2PI * ( G0(I,J)/2D0 + AONE2PI * (
     1              G1(I,J)/4D0 + DCAL * G0(I,J)/2D0 ) )
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
     1        - AONE2PI/2D0 * JGL / ( 1D0 + AONE2PI/2D0 * JLL )
         JSINV(2,2) = JSINV(2,2) / ( 1D0 + AONE2PI/2D0 * JLL )
*
         CALL MMULT(JS,4,4,SGATMP1,4,4,SGATMP2)
         CALL MMULT(SGATMP2,4,4,JSINV,4,4,SGA)
         CALL MMULT(JS,4,4,SGBTMP1,4,4,SGBTMP2)
         CALL MMULT(SGBTMP2,4,4,JSINV,4,4,SGB)         
*
c         WRITE(6,*) "------------"
c         CALL MAGNUS(LN2,SGA,SGB,4,1,EFSG)
c         WRITE(6,*) EFSG
C          CALL MAGNUS(LN2,SGA,SGB,4,2,EFSG)
c         WRITE(6,*) EFSG
c         CALL MAGNUS(LN2,SGA,SGB,4,3,EFSG)
c         WRITE(6,*) EFSG
          CALL MAGNUS(LN2,SGA,SGB,4,4,EFSG)
c         WRITE(6,*) EFSG
*     
*     Non Singlet (equal for Delta and MSbar factorization scheme)
*     
*     Write the analytic solution for the non-singlet ev.op. 
         DO I=1,2
            DO J=1,3
               SPNS(I,J) = ( ( AONE2PI * G0NS(J)/2D0
     1              + AONE2PI**2 *
     2              ( G1NS(I,J)/4D0 + DK*G0NS(J)/2D0 ) ) * DLOG(Q2F/Q2I)
     3              + AONE2PI**2 * PI*B0*G0NS(J)/2D0 *
     4              (DLOG(Q2F/MF2)**2 - DLOG(Q2I/MF2)**2) )            
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
************************************************************************
*
*     Returns the Magnus solution of the matrix equation
*           dY/dt = M(t) * Y(t)
*     and with the special case M(t) = A*t + B
*     The solution can be written as
*            MSOL(t) = exp(sum_{k=1}^{NMAG} Omega_k(t))
*     such that Y(t) = MSOL(t)*Y(0)
*     see e.g. https://en.wikipedia.org/wiki/Magnus_expansion
*
      SUBROUTINE MAGNUS(T,AMAT,BMAT,N,NMAG,MSOL)
*
      IMPLICIT NONE
*      
      include "../commons/tecparam.h"            
**
*     Input Variables
*
      INTEGER N,NMAG
      DOUBLE PRECISION T
      DOUBLE COMPLEX AMAT(N,N),BMAT(N,N) 
**
*     Internal Variables
*
      INTEGER I,J
      DOUBLE COMPLEX OMEGASUM(N,N)
      DOUBLE COMPLEX AB(N,N),BA(N,N),AAB(N,N),ABA(N,N),BAA(N,N)
      DOUBLE COMPLEX ABB(N,N),BAB(N,N),BBA(N,N)
      DOUBLE COMPLEX AAAB(N,N),AABA(N,N),ABAA(N,N),BAAA(N,N)
      DOUBLE COMPLEX AABB(N,N),ABAB(N,N),BABA(N,N),BBAA(N,N)
      DOUBLE COMPLEX ABBB(N,N),BABB(N,N),BBAB(N,N),BBBA(N,N)
**
*     Output Variables
*
      DOUBLE COMPLEX MSOL(N,N)
*
      IF (NMAG.LT.1) THEN
         WRITE(6,*) "NMAG SHOULD BE GREATER THAN 1"
         CALL EXIT(-10)
      ENDIF
*
c     counting based on B = order alpha, A = order alpha^2
      DO I=1,N
         DO J=1,N
c           order alpha
            IF (NMAG.GE.1) THEN
               OMEGASUM(I,J) = 1D0/2D0*T**2 * AMAT(I,J) + T * BMAT(I,J)
            ENDIF
c           order alpha^3
            IF (NMAG.GE.2) THEN
               CALL MMULT(AMAT,N,N,BMAT,N,N,AB)
               CALL MMULT(BMAT,N,N,AMAT,N,N,BA)
               OMEGASUM(I,J) = OMEGASUM(I,J) +
     1              1D0/12D0*T**3 * (AB(I,J) - BA(I,J))
            ENDIF
c           order alpha^5
            IF (NMAG.GE.3) THEN
               CALL MMULT(AMAT,N,N,AB,N,N,AAB)
               CALL MMULT(AB,N,N,AMAT,N,N,ABA)
               CALL MMULT(BA,N,N,AMAT,N,N,BAA)
               OMEGASUM(I,J) = OMEGASUM(I,J) + 1D0/240D0*T**5 *
     1              (AAB(I,J) - 2D0*ABA(I,J) + BAA(I,J))
            ENDIF
c           order alpha^5            
            IF (NMAG.GE.4) THEN
               CALL MMULT(AB,N,N,BMAT,N,N,ABB)
               CALL MMULT(BA,N,N,BMAT,N,N,BAB)
               CALL MMULT(BMAT,N,N,BA,N,N,BBA)
               CALL MMULT(AMAT,N,N,AAB,N,N,AAAB)
               CALL MMULT(AMAT,N,N,ABA,N,N,AABA)
               CALL MMULT(ABA,N,N,AMAT,N,N,ABAA)
               CALL MMULT(BAA,N,N,AMAT,N,N,BAAA)
               CALL MMULT(AMAT,N,N,ABB,N,N,AABB)
               CALL MMULT(ABA,N,N,BMAT,N,N,ABAB)
               CALL MMULT(BMAT,N,N,ABA,N,N,BABA)
               CALL MMULT(BMAT,N,N,BAA,N,N,BBAA)
               CALL MMULT(ABB,N,N,BMAT,N,N,ABBB)
               CALL MMULT(BMAT,N,N,ABB,N,N,BABB)
               CALL MMULT(BMAT,N,N,BAB,N,N,BBAB)               
               CALL MMULT(BMAT,N,N,BBA,N,N,BBBA)
               OMEGASUM(I,J) = OMEGASUM(I,J) + 1D0/5040D0 * (
     1              T**7 * ( -AAAB(I,J) + 3D0*AABA(I,J)
     2              - 3D0*ABAA(I,J) + BAAA(I,J) )
     3              + T**6 * ( 7D0*AABB(I,J) - 14D0*ABAB(I,J)
     4              + 14D0*BABA(I,J) - 7D0*BBAA(I,J) )
     5              + T**5 * ( -7D0*ABBB(I,J) + 21D0*BABB(I,J)
     6              - 21D0*BBAB(I,J) + 7D0*BBBA(I,J) ) )
            ENDIF
         ENDDO
      ENDDO
*
      CALL MATRIXEXP(NMATEXP,N,OMEGASUM,MSOL)
*
      RETURN
      END        
