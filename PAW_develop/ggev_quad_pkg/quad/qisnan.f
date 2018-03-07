      LOGICAL FUNCTION QISNAN( DIN )
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2010
*
*     .. Scalar Arguments ..
      REAL*16   DIN
*     ..
*
*  Purpose
*  =======
*
*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN     (input) REAL*16
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL QLAISNAN
      EXTERNAL QLAISNAN
*  ..
*  .. Executable Statements ..
      QISNAN = QLAISNAN(DIN,DIN)
      RETURN
      END
