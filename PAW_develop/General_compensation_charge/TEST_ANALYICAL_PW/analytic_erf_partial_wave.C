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
#define RT_PI_G (sqrt(M_PI_G))
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

//==========================================================================
// Local funcitons
int main (int , char *[]);
inline double gauss_i_0_integrand(double , double );
inline double gauss_i_1_integrand(double , double );
inline double gauss_i_2_integrand(double , double );
inline double gauss_i_3_integrand(double , double );
inline double gauss_i_0_integral(double , double );
inline double gauss_i_1_integral(double , double );
inline double gauss_i_2_integral(double , double );
inline double gauss_i_3_integral(double , double );
inline double pw_erf_0(double , double, double);
inline double pw_erf_1(double , double, double);
inline double pw_erf_2(double , double, double);
inline double pw_erf_3(double , double, double);
inline void pw_erf_all(int , double *, double , double , double );
void readtoendofline(FILE *);

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller to manage compensation charge energy
//==========================================================================
int main (int argc, char *argv[]){
//==========================================================================  
// Constants 
  int lmax     = 3;
  double beta  = 1.375715;
  double delta = 1.0e-5;
  double pw_erf_tst[10], pw_erf[10];
//==========================================================================  
// Inputs and angular moment invariant quantities
  printf("\n");
  if(argc<4){
    printf("=======================================\n");
    printf(" Command line argument: code.exe r< r> nx\n");
    printf("=======================================\n");
    exit(1);
  }//endif
  double rlt   = atof(argv[1]);
  double rgt   = atof(argv[2]);
  if(rlt>rgt){
    printf("Correcting r< greater than r> error\n\n");
    double t=rlt; rlt = rgt; rgt = t;
  }//endif
  double gmean = sqrt(2.0*rgt*rlt);
  double a     = sqrt(rgt*rgt+rlt*rlt)/gmean;
  double xl    = 0.0;
  double xu    = beta*gmean;
  double x     = xu;
  double xp    = x+delta;
  double xm    = x-delta;
//==========================================================================  
// Super stable global function
  pw_erf_all(lmax,pw_erf,beta,rlt,rgt);

//==========================================================================  
// Check the integrals of the modified spherical bessel functions
  printf("===============================================================\n");
  printf("       Checking the indefinite integrals for partial waves\n");
  for(int il=0;il<=lmax;il++){
    double deriv_exact,yp,ym;
    switch(il){
       case 0: deriv_exact = gauss_i_0_integrand(a, x);
               yp          = gauss_i_0_integral(a, xp);
               ym          = gauss_i_0_integral(a, xm);
               break;
       case 1: deriv_exact = gauss_i_1_integrand(a, x);
               yp          = gauss_i_1_integral(a, xp);
               ym          = gauss_i_1_integral(a, xm);
               break;
       case 2: deriv_exact = gauss_i_2_integrand(a, x);
               yp          = gauss_i_2_integral(a, xp);
               ym          = gauss_i_2_integral(a, xm);
               break;
       case 3: deriv_exact = gauss_i_3_integrand(a, x);
               yp          = gauss_i_3_integral(a, xp);
               ym          = gauss_i_3_integral(a, xm);
               break;
    }//end switch
    double deriv_apprx = 0.5*(yp-ym)/delta;
    printf("---------------------------------------------------------------\n");
    printf("  Integral_%d (%g)  : %.10g %.10g\n",il,x,deriv_exact,deriv_apprx);
  }//endif
  printf("===============================================================\n");
  printf("\n");

//==========================================================================  
// Check the partial waves

  printf("===============================================================\n");
  printf("       Checking the 1st simplification of the partial waves\n");
  for(int il=0;il<=lmax;il++){
    double dl  = (double)il;
    double dl1 = dl+1.0;
    double pre   = 2.0*(2.0*dl+1.0)/(sqrt(M_PI_G)*gmean);
    double lower, upper, simp;
    switch(il){
       case 0: lower   =  gauss_i_0_integral(a, xl);
               upper   =  gauss_i_0_integral(a, xu);
               simp    =  pw_erf_0(rlt,rgt,beta);
             break;
       case 1: lower   =  gauss_i_1_integral(a, xl);
               upper   =  gauss_i_1_integral(a, xu);
               simp    =  pw_erf_1(rlt,rgt,beta);
             break;
       case 2: lower   =  gauss_i_2_integral(a, xl);
               upper   =  gauss_i_2_integral(a, xu);
               simp    =  pw_erf_2(rlt,rgt,beta);
             break;
       case 3: lower   =  gauss_i_3_integral(a, xl);
               upper   =  gauss_i_3_integral(a, xu);
               simp    =  pw_erf_3(rlt,rgt,beta);
             break;
    }//end switch
    double fromInt = (upper-lower)*pre;
    pw_erf_tst[il] = simp;
    printf("---------------------------------------------------------------\n");
    printf("  pw_erf_%d(%g %g; %g)  : %.10g %.10g %.10g\n",il,rlt,rgt,beta,fromInt,simp, fabs(fromInt-simp));
    printf("  beta_lim_%d(%g %g; %g): %.10g %.10g\n",il,rlt,rgt,beta,-lower*pre,(pow(rlt,dl)/pow(rgt,dl1)));
  }//endfor : il
  printf("===============================================================\n");
  printf("\n");

//==================================================================================
// stuff to evaluate mod spherical bessels 1st kind and gaussian integrals over them
  int lmax1 = lmax+1;
  double *pw_erf_num = new double [lmax1];
  double *i_n        = new double [lmax1];
  double *i_n_0      = new double [lmax1];
  double *i_n_1      = new double [lmax1];
  double *i_n_2      = new double [lmax1];
  double *i_n_3      = new double [lmax1];
  double *i_n_1r     = new double [lmax1];
  double *i_n_2r     = new double [lmax1];
  double *i_n_3r     = new double [lmax1];
  double *pref       = new double [lmax1];

  for(int il=0;il<=lmax;il++){
    double al = (double) il;
    pref[il]  = 2.0*(2.0*al+1.0)/sqrt(M_PI_G);
    i_n_1[il] = 0.5/(2.0*al+3.0);
    i_n_2[il] = 0.125/((2.0*al+3.0)*(2.0*al+5.0));
  }//endfor
  i_n_0[0] = 1.0;
  for(int il=1;il<=lmax;il++){
    double al = (double) il;
    i_n_0[il] = i_n_0[il-1]/(2*al+1.0);
  }//endfor

  for(int il=0;il<=lmax;il++){
    double nu   = 0.5 + ((double) il); // need bessel order not spherical bessel order so add 1/2
    double fnu2 = 4.0*nu*nu;
    i_n_1r[il]  = -(fnu2-1.0)/8.0; 
    i_n_2r[il]  = -i_n_1r[il]*(fnu2-9.0)/(8.0*2.0);
    i_n_3r[il]  = -i_n_2r[il]*(fnu2-25.0)/(8.0*3.0);
  }//endfor

//==========================================================================  
// Trapezoidal rule integration

  int nx = atoi(argv[3]);
  double nxd = (double) nx;
  double dx  = beta/nxd;

  pw_erf_num[0] = 0.5*dx*pref[0]; 
  for(int il=1;il<=lmax;il++){pw_erf_num[il] = 0.0; }
  for(int ix=1;ix<=nx;ix++){
    double x   = dx*((double)ix);
    double x2  = x*x;
    double a2  = (rlt*rlt+rgt*rgt);
    double eee = exp(-a2*x2);
    double arg = 2.0*x2*rlt*rgt;
    double arg2 = arg*arg;
    double arg4 = arg2*arg2;
    double sinch_a = sinh(arg)/arg; 
    double cosh_a  = cosh(arg);
    //  use general recurion relation at large enough argument
    if(arg<30.0){
      if(arg>0.01) {i_n[0] = sinch_a;}
      if(arg2>0.01){i_n[1] = (cosh_a - sinch_a)/arg;}
      double argl_p1 = arg2;
      for(int il=1;il<=lmax-1;il++){
        double al = (double) il;
        if(argl_p1*arg > 0.01){i_n[il+1] = i_n[il-1] - (2.0*al+1.0)*i_n[il]/arg;}
        argl_p1 *= arg;
     }//endfor
     //  use expansion at small enough argument
     double argl = 1.0;
     for(int il=0;il<=lmax;il++){
        if(argl*arg<=0.01){i_n[il] = (1.0 + arg2*i_n_1[il] + arg4*i_n_2[il])*argl*i_n_0[il];}
        argl *= arg;
     }//endfor
    }else{
      //  Use asymptotic expansion at large argument to compute renormalized bessel: i_l*2*z*exp(-arg)
      eee = 0.5*exp(-a2*x2+arg)/arg;      // absorb renormalization into a renormalized eee
      double argi = 1.0/arg; double argi2 = argi*argi; double argi3 = argi2*argi;
      for(int il=0;il<=lmax;il++){i_n[il] = (1.0 + argi*i_n_1r[il] + argi2*i_n_2r[il] + argi3*i_n_3r[il]);}
    }//endif
    double w = (ix==nx ? 0.5 : 1.0);
    for(int il=0;il<=lmax;il++){pw_erf_num[il] += pref[il]*eee*dx*i_n[il]*w;}
  }//endfor

  printf("===============================================================\n");
  printf("       Checking the numerical trap integration method\n");
  printf("---------------------------------------------------------------\n");
  for(int il=0;il<=lmax;il++){
    double rat = rlt/rgt;
    double al  = (double)il;
    printf("  pw_erf_%d(%g %g; %g)  : %.10g %.10g %.10g %.10g\n",
           il,rlt,rgt,beta,pw_erf_num[il],pw_erf[il],fabs(pw_erf_num[il]-pw_erf[il]),pow(rat,al)/rgt);  
  }//endfor
  printf("===============================================================\n");
  printf("\n");

//==========================================================================  
// Gauss-Legendre

  int ng;
  FILE *fp = fopen("gl_w_x_200.dat","r");
  if(fp==NULL){
    printf("=============================================================\n");
    printf(" File gl_w_x_200.dat not found\n");
    printf("=============================================================\n");
    exit(1);
  }//endif
  int iii = fscanf(fp,"%d",&ng); 
  if(iii!=1){
    printf("=============================================================\n");
    printf(" 1st line of file gl_w_x_200.dat invalid\n");
    printf("=============================================================\n");
    exit(1);
  }//endif
  readtoendofline(fp);
  double *gx = new double [ng];
  double *wx = new double [ng];
  for(int ig=0;ig<ng;ig++){
    int jjj = fscanf(fp,"%lg %lg",&gx[ig],&wx[ig]);
    if(jjj!=2){
      printf("=============================================================\n");
      printf(" %d th line of file gl_w_x_200.dat invalid\n",ig+2);
      printf("=============================================================\n");
      exit(1);
    }//endif
    readtoendofline(fp);
    gx[ig] = (gx[ig]+1.0)*0.5;
    wx[ig] *= 0.5;
  }//endfor
  fclose(fp);

  for(int il=0;il<=lmax;il++){pw_erf_num[il] = 0.0; }
  for(int ig=0;ig<ng;ig++){
    double sig_inv  = 2.0*(rgt-rlt);
    double beta_now = 1.0/MAX(1.0/beta,sig_inv/10.0);
    double x       = gx[ig]*beta_now;
    double x2      = x*x;
    double a2      = (rlt*rlt+rgt*rgt);
    double eee     = exp(-a2*x2);
    double arg     = 2.0*x2*rlt*rgt;
    double arg2    = arg*arg;
    double arg4    = arg2*arg2;
    double arg6    = arg2*arg4;
    double sinch_a = (arg > 0.0 ? sinh(arg)/arg : 1.0);
    double cosh_a  = cosh(arg);
    if(arg<30.0){
     //  use general recurion relation at intermediate argument
      if(arg>0.01) {i_n[0] = sinch_a;}
      if(arg2>0.01){i_n[1] = (cosh_a - sinch_a)/arg;}
      double argl_p1 = arg2;
      for(int il=1;il<=lmax-1;il++){
        double al = (double) il;
        if(argl_p1*arg > 0.01){i_n[il+1] = i_n[il-1] - (2.0*al+1.0)*i_n[il]/arg;}
        argl_p1 *= arg;
      }//endfor
     //  Use expansion at small argument
      double argl = 1.0;
      for(int il=0;il<=lmax;il++){
         if(argl*arg<=0.01){
           i_n[il] = (1.0 + arg2*i_n_1[il] + arg4*i_n_2[il] + arg6*i_n_3[il])*argl*i_n_0[il];
         }//endif
         argl *= arg;
      }//endfor
    }else{
      //  Use asymptotic expansion at large argument to compute renormalized bessel: i_l*2*z*exp(-arg)
      eee = 0.5*exp(-a2*x2+arg)/arg;      // absorb renormalization into a renormalized eee
      double argi = 1.0/arg; double argi2 = argi*argi; double argi3 = argi2*argi;
      for(int il=0;il<=lmax;il++){i_n[il] = (1.0 + argi*i_n_1r[il] + argi2*i_n_2r[il] + argi3*i_n_3r[il]);}
    }//endif
    for(int il=0;il<=lmax;il++){pw_erf_num[il] += pref[il]*eee*i_n[il]*wx[ig]*beta_now;}
  }//endfor: GL points

  printf("===============================================================\n");
  printf("       Checking the numerical GL integration method\n");
  printf("---------------------------------------------------------------\n");
  for(int il=0;il<=lmax;il++){
    double rat = rlt/rgt;
    double al  = (double)il;
    printf("  pw_erf_%d(%g %g; %g)  : %.10g %.10g %.10g %.10g\n",il,rlt,rgt,beta,pw_erf_num[il],pw_erf[il],
                       fabs(pw_erf_num[il]-pw_erf[il]),pow(rat,al)/rgt);  
  }//endfor
  printf("===============================================================\n");
  printf("\n");
  printf("\n");

//==========================================================================  
// Check the global function

  printf("===============================================================\n");
  printf("       Checking the final form of analytical partial waves\n");
  printf("---------------------------------------------------------------\n");
  pw_erf_all(lmax,pw_erf,beta,rlt,rgt);
  for(int il=0;il<=lmax;il++){
    printf("  pw_erf_%d(%g %g; %g)  : %.10g %.10g %.10g\n",il,rlt,rgt,beta,
                  pw_erf_tst[il],pw_erf[il],fabs(pw_erf_tst[il]-pw_erf[il]));  
  }//endfor
  printf("===============================================================\n");
  printf("\n");

//==========================================================================
  return 1;
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//              e^(-a^2x^2) i_2(x^2) = Gaussian \times i_2(x^2)
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
   double temp =  (derfP)*(rlt2*rlt/rgt4)  
               +  (derfM-derfM_srlt)*(rgt2*rgt/rlt4)
                - (gaussH)*( (-15.0 + 6.0*argE - 4.0*beta4*(rgt4 + rlt4 + 6.0*rgt2*rlt2)
                                               + 8.0*beta6*(rgt6 + rlt4*rgt2 + rlt2*rgt4 + rlt6) )*sinchH 
                            +(+15.0 - 6.0*argE + 4.0*beta4*(rgt4 + rlt4 + 1.0*rlt2*rgt2))*coshH
                            -8.0*poly*rgt6*beta6
                           )/(argH*argH*argH);
//==========================================================================
   return temp;
 }// end routine
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//              e^(-a^2x^2) i_2(x^2) = Gaussian \times i_2(x^2)
//==========================================================================
inline double gauss_i_2_integrand(double a, double x){
//==========================================================================
   double temp;
//==========================================================================
   if(x>0){
     double a2     = a*a;
     double x2     = x*x;
     double x4     = x2*x2;
     double ch_x2  = cosh(x2);
     double shc_x2 = sinh(x2)/x2;
     temp = exp(-a2*x2)* ((x4+3.0)*shc_x2 - 3.0*ch_x2)/x4;
   }else{
     temp = 0.0;
   }//endif
//==========================================================================
   return temp;
}// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//   Gaussian integral of i_2(x^2)
//==========================================================================
inline double gauss_i_2_integral(double a, double x){
//==========================================================================
   double temp;
   double pre    = sqrt(M_PI_G)/10.0;
   double a2     = a*a;
   double a4     = a2*a2;
   double ap2    = a2 + 1.0;
   double am2    = a2 - 1.0;
   double ap     = sqrt(ap2);
   double am     = sqrt(am2);
   double x2     = x*x;
   double x3     = x2*x;
   double x4     = x2*x2;
   double ch_x2  = cosh(x2);
   double shc_x2 = sinh(x2)/x2;
   double gH     = exp(-a2*x2);
   if(x>0){
      temp = pre*( (4.0*a4 + 2.0*a2 - 1.0)*am*erfc(am*x) 
                  -(4.0*a4 - 2.0*a2 - 1.0)*ap*erfc(ap*x) )
           + (gH/(5.0*x3))*( ((1.0-4.0*a4)*x4 + 2.0*a2*x2 - 3.0)*shc_x2 
                                             -(2.0*a2*x2 - 3.0)*ch_x2 + 4.0*a4*x4 )
           - 4.0*x*a4*gH/5.0;
   }else{
      temp =  pre*( (4.0*a4 + 2.0*a2 - 1.0)*am*erfc(am*x) 
                   -(4.0*a4 - 2.0*a2 - 1.0)*ap*erfc(ap*x) );
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
inline double pw_erf_2(double rlt, double rgt, double beta){
//==========================================================================
   double sqrtPi  = sqrt(M_PI_G);
   double pre     = 2.0*beta/sqrtPi;
   double beta2   = beta*beta;
   double beta4   = beta2*beta2;
   double rgt2    = rgt*rgt;
   double rlt2    = rlt*rlt;
   double rgt4    = rgt2*rgt2;
   double rlt4    = rlt2*rlt2;
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
   double combo   = rgt4 + rlt4 + rgt2*rlt2;
   double bombo   = rgt4 + rlt4;
   double mombo   = beta2*rlt2*(1.0+rgt2*beta2)*2.0/3.0;
//---------------------------------------------------------------------------
   double temp3  = (rlt2/(rgt2*rgt))*(derfP-gaussH*rgt)  
                 + (rgt2/(rlt2*rlt))*(derfM-gaussH*rlt*(1.0+mombo))
                 - gaussH*((4.0*beta4*combo -2.0*argE + 3.0)*sinchH 
                                           +(2.0*argE - 3.0)*coshH - 4.0*beta4*(bombo+mombo*rgt4))/(argH*argH);
//==========================================================================
   double temp = temp3;
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
   if(x>0){
     double a2 = a*a;
     double x2 = x*x;
     double x4 = x2*x2;
     temp = exp(-a2*x2)* (x2*cosh(x2) - sinh(x2))/x4;
   }else{
     temp = 0.0;
   }//endif
//==========================================================================
   return temp;
}// end routine
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
 double gaussH  = pre*exp(-argE);
 double sinchH  = sinh(argH)/argH;
 double coshH   = cosh(argH);

 double temp   = (derfP-rgt*gaussH)*(rlt/rgt2) 
               + (derfM-rlt*gaussH)*(rgt/rlt2)
               - gaussH*( (2.0*argE-1.0)*sinchH + coshH - 2.0*argE)/argH;
//==========================================================================
   return temp;
 }// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//              e^(-a^2x^2) i_l(x^2) = Gaussian \times i_l(x^2)
//==========================================================================
inline double gauss_i_0_integrand(double a, double x){
//==========================================================================
   double temp;
//==========================================================================
   if(x>0){
     double a2 = a*a;
     double x2 = x*x;
     temp = exp(-a2*x2)*(sinh(x2)/x2);
   }else{
     temp = 0.0;
   }//endif
//==========================================================================
   return temp;
}// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//   Gaussian integral of i_1(x^2)
//==========================================================================
inline double gauss_i_0_integral(double a, double x){
//==========================================================================
   double temp;
   double pre = sqrt(M_PI_G)/2.0;
   double a2  = a*a;
   double ap2 = a*a + 1.0;
   double am2 = a*a - 1.0;
   double ap = sqrt(ap2);
   double am = sqrt(am2);
   double x2  = x*x;
   if(x>0){
      temp = pre*( am*erfc(am*x) - ap*erfc(ap*x) )
            -exp(-a2*x2)*x*(sinh(x2)/x2);
   }else{
      temp = pre*( am*erfc(am*x) - ap*erfc(ap*x) );
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
inline double pw_erf_0(double rlt, double rgt, double beta){
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
 double sinchH  = sinh(argH)/argH;

 double temp   = derfP/rgt
               + (derfM-rlt*gaussH)/rlt
               - gaussH*(sinchH-1.0);
//==========================================================================
   return temp;
 }// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                prefactor*Integral evaluated at limits and simplified
//==========================================================================
inline void pw_erf_all(int lmax, double *pw_erf, double beta, double rlt, double rgt){
//==========================================================================
   double rgt2    = rgt*rgt;           double rlt2     = rlt*rlt;
   double rgt4    = rgt2*rgt2;         
   double brgt    = beta*rgt;          double brlt     = beta*rlt;
   double brgt2   = brgt*brgt;         double brlt2    = brlt*brlt;
   double brgt4   = brgt2*brgt2;       double brlt4    = brlt2*brlt2;
   double brgt6   = brgt2*brgt4;       double brlt6    = brlt2*brlt4;
   double brgt8   = brgt4*brgt4;       double brlt8    = brlt4*brlt4;

   double argP    = brgt+brlt;         double argM     = brgt-brlt;
   double argH    = 2.0*brgt*brlt;     double argE     = brgt2+brlt2;

   double erfP    = erf(argP);         double erfM     = erf(argM);
   double derfP   = (erfP+erfM)*0.5;   double derfM    = (erfP-erfM)*0.5;

   double sinchH  = sinh(argH)/argH;   double coshH   = cosh(argH);
   double norm    = 2.0*beta/RT_PI_G;  double gaussH  = norm*exp(-argE);

   double poly    = 1.0 + 2.0*brlt2*(1.0+brgt2)/3.0;
   double derfM_rlt_S = gaussH;        double derfM_rlt_F = gaussH*poly;

   double combo   = brgt4 + brlt4 + brgt2*brlt2;
//==========================================================================
// Treat derfM at small beta*r< with 6th order Taylor expansion (good to O(brlt^8))
   double dderfM_rlt_S;
   double dderfM_rlt_F;
   double brlt_cut = 0.085; 
   if(brlt>brlt_cut){
     dderfM_rlt_S = (derfM-rlt*derfM_rlt_S)/rlt;
     dderfM_rlt_F = (derfM-rlt*derfM_rlt_F)/rlt;
   }else{
     // remove gaussian diff for precision improvement 
     double poly_m1      = 2.0*brlt2*(1.0+brgt2)/3.0;
     double derfM_rlt_S_G= 0.0;        double derfM_rlt_F_G = gaussH*poly_m1;
     double eeegt        = norm*exp(-brgt2);
     double diffG        = eeegt*(1.0-exp(-brlt2));
     double derfM_rlt_sm = eeegt*(+ ( 2.0*brgt2 -   1.0)*(brlt2/3.0) 
                                  + ( 4.0*brgt4 -  12.0*brgt2 + 3.0)*(brlt4/30.0)
                                  + ( 8.0*brgt6 -  60.0*brgt4 +  90.0*brgt2 - 15.0)*(brlt6/630.0)  
                                  + (16.0*brgt8 - 224.0*brgt6 + 840.0*brgt4 - 840.0*brgt2 + 105.0)*(brlt8/22680.0)
                                  );
     dderfM_rlt_S = derfM_rlt_sm - derfM_rlt_S_G + diffG;
     dderfM_rlt_F = derfM_rlt_sm - derfM_rlt_F_G + diffG;
   }//endif
//==========================================================================
// The first 4 partial waves from analytic work
  //-----------------------------------------------------------------------
  //Intermediate argH :
   if(argH<20.0 && argH>0.2){ 
     pw_erf[0] = (derfP)/rgt + (dderfM_rlt_S) - (gaussH)*(sinchH-1.0);

     if(lmax>=1){
       pw_erf[1] = (derfP)*(rlt/rgt2) + (dderfM_rlt_S)*(rgt/rlt)
                 - (gaussH)*( (2.0*argE-1.0)*sinchH + coshH - 2.0*brgt2)/argH;
     }//endif

     if(lmax>=2){
       pw_erf[2] = (derfP)*(rlt2/(rgt2*rgt)) + (dderfM_rlt_F)*(rgt2/(rlt2))
                 - (gaussH)*((4.0*combo -2.0*argE + 3.0)*sinchH 
                                     +(2.0*argE - 3.0)*coshH - 4.0*brgt4*poly)/(argH*argH);
     }//endif

     if(lmax>=3){
       pw_erf[3] = (derfP)*(rlt2*rlt/rgt4)  + (dderfM_rlt_F)*(rgt2*rgt/(rlt2*rlt))
                 - (gaussH)*( (-15.0 + 6.0*argE - 4.0*(brgt4 + brlt4 + 6.0*brgt2*brlt2)
                                     + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6) )*sinchH 
                             +(+15.0 - 6.0*argE + 4.0*(brgt4 + brlt4 + 1.0*brlt2*brgt2))*coshH
                             -8.0*poly*brgt6
                        )/(argH*argH*argH);
     }//endif
   }//endif : intermediate argH
  //-----------------------------------------------------------------------
  //Large argH :
   if(argH>=20.0){
     double eeeH   = exp(-argH);
     gaussH = norm*exp(-argE+argH);
     sinchH = 0.5*(1.0-eeeH*eeeH)/argH;
     coshH  = 0.5*(1.0+eeeH*eeeH);
     pw_erf[0] = (derfP)/rgt + (dderfM_rlt_S) - (gaussH)*(sinchH-1.0*eeeH);

     if(lmax>=1){
       pw_erf[1] = (derfP)*(rlt/rgt2) + (dderfM_rlt_S)*(rgt/rlt)
                 - (gaussH)*( (2.0*argE-1.0)*sinchH + coshH - 2.0*brgt2*eeeH)/argH;
     }//endif

     if(lmax>=2){
       pw_erf[2] = (derfP)*(rlt2/(rgt2*rgt)) + (dderfM_rlt_F)*(rgt2/(rlt2))
                 - (gaussH)*((4.0*combo -2.0*argE + 3.0)*sinchH 
                                       +(2.0*argE - 3.0)*coshH - 4.0*brgt4*poly*eeeH)/(argH*argH);
     }//endif

     if(lmax>=3){
       pw_erf[3] = (derfP)*(rlt2*rlt/rgt4)  + (dderfM_rlt_F)*(rgt2*rgt/(rlt2*rlt))
                 - (gaussH)*( (-15.0 + 6.0*argE - 4.0*(brgt4 + brlt4 + 6.0*brgt2*brlt2)
                                     + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6) )*sinchH 
                             +(+15.0 - 6.0*argE + 4.0*(brgt4 + brlt4 + 1.0*brlt2*brgt2))*coshH
                             -8.0*poly*brgt6*eeeH
                        )/(argH*argH*argH);
    }//endif
  }//endif : large argH
 //-----------------------------------------------------------------------
 // Small argH :
   if(argH<=0.2){
     double argH2   = argH*argH;
     double rats    = argH2/6.0;
     double ratc    = argH2/2.0;
     double sinchH1 = 0.0;
     double coshH1  = 0.0;
     for(int i = 3; i<=13;i+=2){
       sinchH1 += rats;
       coshH1  += ratc;
       double di = (double) i;
       rats *= (argH2/((di+2.0)*(di+1.0)));
       ratc *= (argH2/(di*(di+1.0)));
     }//endfor

     pw_erf[0] = (derfP)/rgt + (dderfM_rlt_S) - (gaussH)*(sinchH1);

     if(lmax>=1){
       pw_erf[1] = (derfP)*(rlt/rgt2) + (dderfM_rlt_S)*(rgt/rlt)
                 - (gaussH)*( (2.0*argE-1.0)*sinchH1 + coshH1 + 2.0*argE - 2.0*brgt2)/argH;
     }//endif

     if(lmax>=2){
       pw_erf[2] = (derfP)*(rlt2/(rgt2*rgt)) + (dderfM_rlt_F)*(rgt2/(rlt2))
                 - (gaussH)*((4.0*combo -2.0*argE + 3.0)*sinchH1  + 4.0*combo
                                       +(2.0*argE - 3.0)*coshH1 - 4.0*brgt4*poly)/(argH*argH);
     }//endif

     if(lmax>=3){
       pw_erf[3] = (derfP)*(rlt2*rlt/rgt4)  + (dderfM_rlt_F)*(rgt2*rgt/(rlt2*rlt))
                 - (gaussH)*( (-15.0 + 6.0*argE - 4.0*(brgt4 + brlt4 + 6.0*brgt2*brlt2)
                                     + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6) )*sinchH1
                                     - 20.0*brgt2*brlt2 + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6)
                             +(+15.0 - 6.0*argE + 4.0*(brgt4 + brlt4 + 1.0*brlt2*brgt2))*coshH1
                             -8.0*poly*brgt6
                            )/(argH*argH*argH);
     }//endif
   }//endif : small argH
//==========================================================================
// nan warning

  for(int il=0;il<=3;il++){ 
    if(isnan(pw_erf[il]) || !isfinite(pw_erf[il])){
      printf("analtyical_pw_%d(%g %g) nan\n",il,rlt,rgt);
    }//endif
  }//endfor

//==========================================================================
 }// end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// readtoendofline: Function to read to end of line in read_coord files     
//==========================================================================
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol){ch=fgetc(fp);}
  if(ch==EOF){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("error: unexpected end of file reached          \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }//endif
}// end routine 
//==========================================================================
