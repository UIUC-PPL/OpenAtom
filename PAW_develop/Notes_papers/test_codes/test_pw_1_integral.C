//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//             Test the Gaussian integral of the i_1(x^2)
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
inline double gauss_i_1_integral(double , double );
inline double gauss_i_1_integrand(double , double );

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
  double deriv_exact = gauss_i_1_integrand(a, x);
  double yp          = gauss_i_1_integral(a, xp);
  double ym          = gauss_i_1_integral(a, xm);
  double deriv_apprx = 0.5*(yp-ym)/delta;
  printf("results %g : %.10g %.10g\n",x,deriv_exact,deriv_apprx);
//==========================================================================  
  return 1;	
 }//end routine main
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//   Gaussian integral of i_1(x^2)
//==========================================================================
inline double gauss_i_1_integral(double a, double x){
//==========================================================================
   double temp;
   double pre = sqrt(M_PI_G)/6.0;
   double a2  = a*a;
   double ap2 = a*a + 1.0;
   double am2 = a*a - 1.0;
   double ap = sqrt(ap2);
   double am = sqrt(am2);
   double x2  = x*x;
   if(x>0){
      temp = pre*( (2.0*a2 + 1.0)*am*erfc(am*x) - (2.0*a2 - 1.0)*ap*erfc(ap*x) )
            -exp(-a2*x2)/(3.0*x)*((2.0*a2*x2 - 1.0)*(sinh(x2)/x2) + cosh(x2));
   }else{
      temp = pre*( (2.0*a2 + 1.0)*am*erfc(am*x) - (2.0*a2 - 1.0)*ap*erfc(ap*x) );
   }//endif
//==========================================================================
   return temp;
}// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//              e^(-a^2x^2) i_l(x^2) = Gaussian \times i_l(x^2)
//==========================================================================
inline double gauss_i_1_integrand(double a, double x){
//==========================================================================
   double temp;
//==========================================================================
   double a2 = a*a;
   double x2 = x*x;
   double x4 = x2*x2;
   if(x>0){
     temp = exp(-a2*x2)* (x2*cosh(x2) - sinh(x2))/x4;
   }else{
     temp = 0.0;
   }//endif
//==========================================================================
   return temp;
}// end routine
//==========================================================================
