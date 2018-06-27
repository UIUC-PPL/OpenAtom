//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//             Test the Gaussian integral of the i_3(x^2)
//==========================================================================
// Standard include files
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <time.h>
using namespace std;

//==========================================================================
// General constants and functions
#define M_PI_G 3.14159265358979323846264338327950288419716939937510
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

//==========================================================================
// Local funcitons
int main (int , char *[]);
inline double gauss_i_3_integral(double , double );
inline double gauss_i_3_integrand(double , double );

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : test integral my comparing numerical derivate to integrand
//==========================================================================
int main (int argc, char *argv[]){
//==========================================================================  
// Constants and inputs
  double delta = 1.0e-5;
  double a     = 1.24751;
  double x     = atof(argv[1]);
//==========================================================================  
// Check the integral at the current value of x
  double xp          = x+delta;
  double xm          = x-delta;
  double deriv_exact = gauss_i_3_integrand(a, x);
  double yp          = gauss_i_3_integral(a, xp);
  double ym          = gauss_i_3_integral(a, xm);
  double deriv_apprx = 0.5*(yp-ym)/delta;
  printf("results %g : %.10g %.10g\n",x,deriv_exact,deriv_apprx);
//==========================================================================  
  return 1;
 }//end routine main
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//   Gaussian integral of i_3(x^2)
//==========================================================================
inline double gauss_i_3_integral(double a, double x){
//==========================================================================
   double temp;
   double pre    = sqrt(M_PI_G)/14.0;
   double a2     = a*a;
   double a4     = a2*a2;
   double a6     = a2*a4;
   double ap2    = a2 + 1.0;
   double am2    = a2 - 1.0;
   double ap     = sqrt(ap2);
   double am     = sqrt(am2);
   double x2     = x*x;
   double x4     = x2*x2;
   double x6     = x2*x4;
   double x5     = x*x4;
   double ch_x2  = cosh(x2);
   double shc_x2 = sinh(x2)/x2;
   double gH     = exp(-a2*x2);
   if(x>0){
     temp = pre*( (8.0*a6 + 4.0*a4 - 4.0*a2 - 1.0)*am*erfc(am*x)
                 -(8.0*a6 - 4.0*a4 - 4.0*a2 + 1.0)*ap*erfc(ap*x))
            -gH*(-(15.0  - 6.0*a2*x2 + (4.0*a4 + 4.0)*x4 +(-8.0*a6 + 4.0*a2)*x6)*shc_x2
                 +(15.0  - 6.0*a2*x2 + (4.0*a4 - 1.0)*x4                       )*ch_x2
                )/(7.0*x5);
   }else{
     temp = pre*( (+8.0*a6 + 4.0*a4 - 4.0*a2 - 1.0)*am*erfc(am*x)
                 -(+8.0*a6 - 4.0*a4 - 4.0*a2 + 1.0)*ap*erfc(ap*x));
   }//endif
//==========================================================================
   return temp;
}// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//              e^(-a^2x^2) i_3(x^2) = Gaussian \times i_3(x^2)
//==========================================================================
inline double gauss_i_3_integrand(double a, double x){
//==========================================================================
   double temp;
//==========================================================================
   if(x>0){
     double a2     = a*a;
     double x2     = x*x;
     double x4     = x2*x2;
     double x6     = x4*x2;
     double ch_x2  = cosh(x2);
     double shc_x2 = sinh(x2)/x2;
     temp = exp(-a2*x2)*((x4 + 15.0)*ch_x2 - (6.0*x4 + 15.0)*shc_x2)/x6;
   }else{
     temp = 0.0;
   }//endif
//==========================================================================
   return temp;
}// end routine
//==========================================================================
