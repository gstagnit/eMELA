************************************************************************
*
*      mmult.f:
*
*      Compute exponential of square matrix.
*
************************************************************************
      SUBROUTINE MATRIXEXP(NT,N,M,EXPM)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER NT,N
      DOUBLE COMPLEX M(N,N)
**
*     Internal Variables
*
      INTEGER I,J,K
      DOUBLE COMPLEX TKM1(N,N),TK(N,N)
**
*     Output Variables
*
      DOUBLE COMPLEX EXPM(N,N)
*
*     Initialize the output matrix to unit matrix
*
      DO I=1,N
         DO J=1,N
            TKM1(I,J) = (0D0,0D0)
         ENDDO
         TKM1(I,I) = (1D0,0D0)
      ENDDO
*
      CALL EQUATEM(TKM1,N,N,EXPM)
      DO K=1,NT
         CALL MMULTF(TKM1,N,N,M,N,N,TK,1D0/DBLE(K))
         CALL SUMM(EXPM,N,N,TK)
         CALL EQUATEM(TK,N,N,TKM1)
      ENDDO
*
      RETURN
      END
