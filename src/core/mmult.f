************************************************************************
*
*      mmult.f:
*
*      Matrix multiplication utilities
*
************************************************************************
      SUBROUTINE MMULT(A,ROWSA,COLSA,B,ROWSB,COLSB,C)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER ROWSA,COLSA,ROWSB,COLSB
      DOUBLE COMPLEX A(ROWSA,COLSA),B(ROWSB,COLSB)
**
*     Internal Variables
*
      INTEGER I,J,K
**
*     Output Variables
*
      DOUBLE COMPLEX C(ROWSA,COLSB)
*
*     Check that the matrices to multiply have the correct dimensions.
*
      IF(COLSA.NE.ROWSB)THEN
          WRITE(6,*) "In mmult:"
          WRITE(6,*) "Input matrices do not match"
          CALL EXIT(-10)
      ENDIF
*
      DO I=1,ROWSA
         DO J=1,COLSB
            C(I,J) = (0D0,0D0)
            DO K=1,COLSA
	       C(I,J) = C(I,J) + A(I,K) * B(K,J)
	    ENDDO
	 ENDDO
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE MVMULT(M,ROWSM,COLSM,V,ROWSV)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER ROWSM,COLSM,ROWSV
      DOUBLE COMPLEX M(ROWSM,COLSM),V(ROWSV)
**
*     Internal Variables
*
      INTEGER I,J
      DOUBLE COMPLEX W(ROWSM)
*
*     Check that dimenstions match.
*
      IF(COLSM.NE.ROWSV)THEN
          WRITE(6,*) "In mvmult:"
          WRITE(6,*) "Input matrices do not match"
          CALL EXIT(-10)
      ENDIF
*
      DO I=1,ROWSM
         W(I) = (0D0,0D0)
         DO J=1,COLSM
            W(I) = W(I) + M(I,J) * V(J)
	 ENDDO
      ENDDO
      DO I=1,COLSM
         V(I) = W(I)
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE MMULTF(A,ROWSA,COLSA,B,ROWSB,COLSB,C,F)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER ROWSA,COLSA,ROWSB,COLSB
      DOUBLE COMPLEX A(ROWSA,COLSA),B(ROWSB,COLSB)
      DOUBLE PRECISION F
**
*     Internal Variables
*
      INTEGER I,J,K
**
*     Output Variables
*
      DOUBLE COMPLEX C(ROWSA,COLSB)
*
*     Check that the matrices to multiply have the  correct dimensions.
*
      IF(COLSA.NE.ROWSB)THEN
          WRITE(6,*) "In mmultf:"
          WRITE(6,*) "Input matrices do not match"
          CALL EXIT(-10)
      ENDIF
*
      DO I=1,ROWSA
         DO J=1,COLSB
            C(I,J) = (0D0,0D0)
            DO K=1,COLSA
	       C(I,J) = C(I,J) + A(I,K) * B(K,J)
	    ENDDO
            C(I,J) = C(I,J) * F
	 ENDDO
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE EQUATEM(A,ROWS,COLS,B)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER ROWS,COLS
      DOUBLE COMPLEX A(ROWS,COLS)
**
*     Internal Variables
*
      INTEGER I,J
**
*     Output Variables
*
      DOUBLE COMPLEX B(ROWS,COLS)
*
*     Initialize the output matrix to zero.
*
      DO I=1,ROWS
         DO J=1,COLS
            B(I,J) = A(I,J)
         ENDDO
      ENDDO
*
      RETURN
      END
*
************************************************************************
      SUBROUTINE SUMM(A,ROWS,COLS,B)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER ROWS,COLS
      DOUBLE COMPLEX A(ROWS,COLS),B(ROWS,COLS)
**
*     Internal Variables
*
      INTEGER I,J
*
*     Initialize the output matrix to zero.
*
      DO I=1,ROWS
         DO J=1,COLS
            A(I,J) = A(I,J) + B(I,J)
         ENDDO
      ENDDO
*
      RETURN
      END
