      REAL*16 FUNCTION QLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      REAL*16   X, Y
*     ..
*
*  Purpose
*  =======
*
*  QLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) REAL*16
*  Y       (input) REAL*16
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*16   ZERO
      PARAMETER          ( ZERO = 0.0Q0 )
      REAL*16   ONE
      PARAMETER          ( ONE = 1.0Q0 )
*     ..
*     .. Local Scalars ..
      REAL*16   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         QLAPY2 = W
      ELSE
         QLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of QLAPY2
*
      END
