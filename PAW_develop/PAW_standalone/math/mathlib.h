//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  The function prototypes for mathlib.C
//
//==========================================================================
// PINY constants, inline functions and typedefs

#ifndef _MATHLIB_
#define _MATHLIB_

//==========================================================================
// constants for approx to erfc
#ifndef M_PI_QI
#define M_PI_QI 3.14159265358979323846264338327950288419716939937510
#endif
#define PERFC  ( 0.3614)
#define CERFC1 ( 0.20414220964220030)
#define CERFC2 ( 0.19975359569614810)
#define CERFC3 ( 0.22131765964055760)
#define CERFC4 ( 0.03360430734640255)
#define CERFC5 ( 0.47325925787217550) 
#define CERFC6 (-0.50907852006973500)
#define CERFC7 ( 0.67726314919476460)
#define CERFC8 (-0.36991297909221700)
#define CERFC9 ( 0.06965131976970335)
#define DCERFC1 (1.0*CERFC1)
#define DCERFC2 (2.0*CERFC2)
#define DCERFC3 (3.0*CERFC3)
#define DCERFC4 (4.0*CERFC4)
#define DCERFC5 (5.0*CERFC5)
#define DCERFC6 (6.0*CERFC6)
#define DCERFC7 (7.0*CERFC7)
#define DCERFC8 (8.0*CERFC8)
#define DCERFC9 (9.0*CERFC9)
#define PRE_ERFC (2.0/sqrt(M_PI_QI))
//==========================================================================
 
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Function prototypes
void gethinv(double *, double *, double *, int );
double dist(double, double, double);	
inline double erfc_a_r_over_r(double, double, double, double);
inline double gerfc(double, double, double *);
//==========================================================================

#endif
