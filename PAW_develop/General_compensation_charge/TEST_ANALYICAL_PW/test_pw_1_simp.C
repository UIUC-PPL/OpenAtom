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
inline double gauss_i_1_integral(double , double );
inline double pw_erf_1(double , double, double);

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
  double pre   = 6.0/(sqrt(M_PI_G)*gmean);

//==========================================================================  
// Check the integral at the current value of x
  double fromInt = (gauss_i_1_integral(a, xu) - gauss_i_1_integral(a, xl))*pre;
  double simp    = pw_erf_1(rlt,rgt,beta);
  printf("pw_erf_1(%g %g): %.10g %.10g\n",rlt,rgt,fromInt,simp);

//==========================================================================
  return 1;	
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                Gaussian integral of i_1(x^2)
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
//                prefactor*Integral evaluated at limits and simplified
//==========================================================================
inline double pw_erf_1(double rlt, double rgt, double beta){
//==========================================================================
 double sqrtPi = sqrt(M_PI_G);
 double pre    = 2.0*beta/sqrtPi;
 double beta2  = beta*beta;
 double rgt2   = rgt*rgt;
 double rlt2   = rlt*rlt;
 double argP   = beta*(rgt+rlt);
 double argM   = beta*(rgt-rlt);
 double argH   = 2.0*beta2*rgt*rlt;
 double argE   = beta2*(rgt2+rlt2);
 double erfP   = erf(argP);
 double erfM   = erf(argM);
 double derfP  = (erfP+erfM)*0.5;
 double derfM  = (erfP-erfM)*0.5;
 double gaussH = pre*exp(-argE);
 double sinchH = sinh(argH)/argH;
 double coshH  = cosh(argH);
 double temp   = (derfP)*(rlt/rgt2) 
               + (derfM-rlt*gaussH)*(rgt/rlt2)
               - gaussH*( (2.0*argE-1.0)*sinchH + coshH - 2.0*beta2*rgt2)/argH;
 double part1  = (derfP)*(rlt/rgt2) 
               + (derfM-rlt*gaussH)*(rgt/rlt2);
 double part2  = - gaussH*( (2.0*argE-1.0)*sinchH + coshH - 2.0*beta2*rgt2)/argH;
 printf("parts : %.10g %.10g %.10g %.10g\n",part1,part2,part1+part2,temp);
//==========================================================================
   return temp;
 }// end routine
//==========================================================================
