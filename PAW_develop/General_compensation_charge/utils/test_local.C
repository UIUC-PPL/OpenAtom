//===========================================================================
// fast erfc and derivative
//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================
// include files
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <time.h>

//===========================================================================
// defines
#define M_PI_QI 3.14159265358979323846264338327950288419716939937510
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

#define PERFC  (0.3614)
#define CERFC1 (0.2041422096422003)
#define CERFC2 (0.1997535956961481)
#define CERFC3 (0.2213176596405576)
#define CERFC4 (0.03360430734640255)
#define CERFC5 (0.4732592578721755) 
#define CERFC6 (-0.509078520069735)
#define CERFC7 (0.6772631491947646)
#define CERFC8 (-0.369912979092217)
#define CERFC9 (0.06965131976970335)
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
//===========================================================================
// Function prototypes
int main();
inline double gerfc(double , double , double *);

//===========================================================================

//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================
// Generic error function
//===========================================================================
 int main(){
//===========================================================================
// declarations
   int n         = 10000000;
   long int seed = 1458431515;
   double  a     = 0.8142;
   double pre    = 2.0/sqrt(M_PI_QI);

   clock_t start,end;

   double *r   = new double [n];
   double *en  = new double [n];
   double *f   = new double [n];
   double *gen = new double [n];
   double *gf  = new double [n];
   double *ggf = new double [n];

//===========================================================================
// Intialize some distances -- not so slow.
   srand48(seed);
   for (int i=0;i<n;i++){r[i] = 3.0*drand48()+0.1;}
   printf("\n");

//===========================================================================
// Slow erfc(ar)/r and -d/dr (erfc(ar)/r)
   start = clock();
   for (int i=0;i<n;i++){
     double x = a*r[i];
     double erfc_tmp  = erfc(x);
	 double rinv = 1.0/r[i];
     en[i] = erfc_tmp*rinv;
     f[i]  = (a*pre*exp(-x*x) + erfc_tmp*rinv)*rinv;
   }//endfor
   end = clock();
   printf("slow energy and force timing %.10g\n",((double)(end-start))/CLOCKS_PER_SEC);

//===========================================================================
// Fast erfc(ar)/r and -d/dr (erfc(ar)/r)
   start = clock();
   for (int i=0;i<n;i++){
    double x        = a*r[i];
    double eee      = exp(-x*x);
    double tt       = 1.0/(1.0+PERFC*x);
    double gerfc_tmp  = ((((((((CERFC9*tt+CERFC8)*tt+CERFC7)*tt+CERFC6)*tt+CERFC5)*tt+CERFC4)*tt+CERFC3)*tt+CERFC2)*tt+CERFC1)*tt*eee;
    double fgerfc = ((((((((DCERFC9*tt+DCERFC8)*tt+DCERFC7)*tt+DCERFC6)*tt+DCERFC5)*tt+DCERFC4)*tt+DCERFC3)*tt+DCERFC2)*tt+DCERFC1)*tt*tt*eee*PERFC
                           +2.0*gerfc_tmp*x;
	double rinv = 1.0/r[i];
    fgerfc       *= a;   // -d/dr ( erfc(a*r)  )
     gen[i] = gerfc_tmp*rinv;
     gf[i]  = (fgerfc + gerfc_tmp*rinv)*rinv;
   }//endfor
   end = clock();
   printf("fast energy and force timing %.10g\n",((double)(end-start))/CLOCKS_PER_SEC);

//===========================================================================
// Fast erfc(ar)/r and -d/dr (erfc(ar)/r)
   start = clock();
   for (int i=0;i<n;i++){
    double x        = a*r[i];
    double eee      = exp(-x*x);
    double tt       = 1.0/(1.0+PERFC*x);
    double gerfc_tmp  = ((((((((CERFC9*tt+CERFC8)*tt+CERFC7)*tt+CERFC6)*tt+CERFC5)*tt+CERFC4)*tt+CERFC3)*tt+CERFC2)*tt+CERFC1)*tt*eee;
    double fgerfc = a*pre*eee; 
	double rinv = 1.0/r[i];
     gen[i] = gerfc_tmp*rinv;
     ggf[i]  = (fgerfc + gerfc_tmp*rinv)*rinv;
   }//endfor
   end = clock();
   printf("test energy and force timing %.10g\n",((double)(end-start))/CLOCKS_PER_SEC);

//===========================================================================
// Fast erfc(ar)/r and -d/dr (erfc(ar)/r)
   start = clock();
   for (int i=0;i<n;i++){
	double fgerfc;
	double gerfc_tmp = gerfc(r[i], a, &fgerfc);
	double rinv = 1.0/r[i];
     gen[i] = gerfc_tmp*rinv;
     ggf[i]  = (fgerfc + gerfc_tmp*rinv)*rinv;
   }//endfor
   end = clock();
   printf("test energy and force timing %.10g\n",((double)(end-start))/CLOCKS_PER_SEC);

//===========================================================================
// Find the max erro
   double max_en_err = 0.0;
   double max_f_err = 0.0;
   double max_gf_err = 0.0;
   for (int i=0;i<n;i++){
     max_en_err = MAX(fabs(gen[i]-en[i]),max_en_err);
     max_f_err = MAX(fabs(gf[i]-f[i]),max_f_err);
     max_gf_err = MAX(fabs(ggf[i]-gf[i]),max_gf_err);
   }//endfor
   printf("Max error %.10g %.10g %.10g\n",max_en_err,max_f_err, max_gf_err);

//===========================================================================
// Check one exzmple   
   double rn = 1.0;
   double x = a*rn;
   double fgerfc_now;
   double gerfc_now = gerfc(rn,a,&fgerfc_now);
   double ferfc_now = a*pre*exp(-x*x);

   double delta   = 0.00001;
   double rp      = (rn+delta);
   double rm      = (rn-delta);
   double xp      = rp*a;
   double xm      = rm*a;
   double erfc_p  = erfc(xp);
   double erfc_m  = erfc(xm);
   double fgerfc_tmp;
   double gerfc_p = gerfc(rp,a,&fgerfc_tmp);
   double gerfc_m = gerfc(rm,a,&fgerfc_tmp);
   printf("Derivative check: %.10g %.10g %.10g %.10g\n",-0.5*(gerfc_p-gerfc_m)/delta,fgerfc_now,
                                      -0.5*(erfc_p-erfc_m)/delta,ferfc_now);
   printf("\n");

//===========================================================================
   return 1;
}  // end routine main
//===========================================================================

//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================
// Generic error function
//===========================================================================
  inline double gerfc(double r, double a, double *fgerfc){
//===========================================================================
    double x        = a*r;
    double eee      = exp(-x*x);
    double tt       = 1.0/(1.0+PERFC*x);
    double gerfc_x  = ((((((((CERFC9*tt+CERFC8)*tt+CERFC7)*tt+CERFC6)*tt+CERFC5)*tt+CERFC4)*tt+CERFC3)*tt+CERFC2)*tt+CERFC1)*tt*eee;
//    double fgerfc_x = ((((((((DCERFC9*tt+DCERFC8)*tt+DCERFC7)*tt+DCERFC6)*tt+DCERFC5)*tt+DCERFC4)*tt+DCERFC3)*tt+DCERFC2)*tt+DCERFC1)*tt*tt*eee*PERFC
//                         +2.0*gerfc_x*x;
    fgerfc[0]       = a*eee*PRE_ERFC;   // -d/dr ( erfc(a*r)  )
//===========================================================================
    return gerfc_x;
  }//end routine gerf
//===========================================================================
