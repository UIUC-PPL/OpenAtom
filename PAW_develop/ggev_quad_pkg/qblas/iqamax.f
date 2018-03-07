      INTEGER FUNCTION IQAMAX(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL*16 DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL*16 DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
      IQAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IQAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = ABS(DX(1))
         DO I = 2,N
            IF (ABS(DX(I)).GT.DMAX) THEN
               IQAMAX = I
               DMAX = ABS(DX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = ABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (ABS(DX(IX)).GT.DMAX) THEN
               IQAMAX = I
               DMAX = ABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
