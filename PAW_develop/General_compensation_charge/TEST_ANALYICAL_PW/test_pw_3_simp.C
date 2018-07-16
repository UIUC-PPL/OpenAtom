//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files
//==========================================================================
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <time.h>

//==========================================================================
// General constants and Functions
using namespace std;
#define M_PI_G 3.14159265358979323846264338327950288419716939937510
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

//==========================================================================
// Local funcitons
int main (int , char *[]);
inline double gauss_i_3_integral(double , double );
inline double pw_erf_3(double , double, double);

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller to manage compensation charge energy
//==========================================================================
int main (int argc, char *argv[]){
//==========================================================================  
// Constants and inputs
  double beta  = 1.375715;
  double rlt   = atof(argv[1]);
  double rgt   = atof(argv[2]);
  
  double gmean = sqrt(2.0*rgt*rlt);
  double a     = sqrt(rgt*rgt+rlt*rlt)/gmean;
  double xl    = 0.0;
  double xu    = beta*gmean;
  double pre   = 14.0/(sqrt(M_PI_G)*gmean);

//==========================================================================  
// Check the integral at the current value of x
  double fromInt = (gauss_i_3_integral(a, xu) - gauss_i_3_integral(a, xl))*pre;
  double limit   = -gauss_i_3_integral(a, xl)*pre;
  double simp    =  pw_erf_3(rlt,rgt,beta);
  printf("pw_erf_3(%g %g): %.10g %.10g\n",rlt,rgt,fromInt,simp);
  printf("beta_lim(%g %g): %.10g %.10g\n",rlt,rgt,limit,rlt*rlt*rlt/(rgt*rgt*rgt*rgt));

//==========================================================================
  return 1;
}//end routine
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
//                prefactor*Integral evaluated at limits and simplified
//==========================================================================
inline double pw_erf_3(double rlt, double rgt, double beta){
//==========================================================================
   double sqrtPi  = sqrt(M_PI_G);
   double pre     = 2.0*beta/sqrtPi;
   double beta2   = beta*beta;
   double beta4   = beta2*beta2;
   double beta6   = beta2*beta4;
   double rgt2    = rgt*rgt;
   double rlt2    = rlt*rlt;
   double rgt4    = rgt2*rgt2;
   double rlt4    = rlt2*rlt2;
   double rgt6    = rgt2*rgt4;
   double rlt6    = rlt2*rlt4;
   double argP    = beta*(rgt+rlt);
   double argM    = beta*(rgt-rlt);
   double argH    = 2.0*beta2*rgt*rlt;
   double argE    = beta2*(rgt2+rlt2);
   double erfP    = erf(argP);
   double erfM    = erf(argM);
   double derfP   = (erfP+erfM)*0.5;
   double derfM   = (erfP-erfM)*0.5;
   double gaussH  = pre*exp(-argE);
   double sinchH  = sinh(argH)/argH;
   double coshH   = cosh(argH);
   double poly    = (1.0 + beta2*rlt2*(1.0+beta2*rgt2)*2.0/3.0);
   double derfM_srlt =  gaussH*rlt*poly;
//---------------------------------------------------------------------------
// check each step of the simplification 
   double temp =  (derfP)*(rlt2*rlt/rgt4)  +  (derfM-derfM_srlt)*(rgt2*rgt/rlt4)
                - (gaussH)*( (-15.0 + 6.0*argE - 4.0*beta4*(rgt4 + rlt4 + 6.0*rgt2*rlt2)
                                               + 8.0*beta6*(rgt6 + rlt4*rgt2 + rlt2*rgt4 + rlt6) )*sinchH 
                            +(+15.0 - 6.0*argE + 4.0*beta4*(rgt4 + rlt4 + 1.0*rlt2*rgt2))*coshH
                            -8.0*poly*rgt6*beta6
                           )/(argH*argH*argH);
//==========================================================================
   return temp;
 }// end routine
//==========================================================================
