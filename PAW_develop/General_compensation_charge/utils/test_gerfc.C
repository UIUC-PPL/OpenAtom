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
//   int n         = 10000000;
	int n = 1;
   long int seed = 1458431515;
   double  a     = 0.8142;
   double pre    = 2.0/sqrt(M_PI_QI);

   clock_t start,end;

   double *r   = new double [n];
   double *en  = new double [n];
   double *f   = new double [n];
   double *gen = new double [n];
   double *gf  = new double [n];

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
     en[i] = erfc_tmp/r[i];
     f[i]  = a*pre*exp(-x*x)/r[i] + erfc_tmp/(r[i]*r[i]);
   }//endfor
   end = clock();
   printf("slow energy and force timing %.10g\n",((double)(end-start))/CLOCKS_PER_SEC);

//===========================================================================
// Fast erfc(ar)/r and -d/dr (erfc(ar)/r)
   start = clock();
   for (int i=0;i<n;i++){
     double fgerfc;
     double gerfc_tmp = gerfc(r[i],a,&fgerfc);
     gen[i] = gerfc_tmp/r[i];
     gf[i]  = fgerfc/r[i] + gerfc_tmp/(r[i]*r[i]);       
   }//endfor
   end = clock();
   printf("fast energy and force timing %.10g\n",((double)(end-start))/CLOCKS_PER_SEC);

//===========================================================================
// Find the max erro
   double max_en_err = 0.0;
   double max_f_err = 0.0;
   for (int i=0;i<n;i++){
     max_en_err = MAX(fabs(gen[i]-en[i]),max_en_err);
     max_f_err = MAX(fabs(gf[i]-f[i]),max_f_err);
   }//endfor
   printf("Max error %.10g %.10g\n",max_en_err,max_f_err);

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
    const double p  = 0.3614;
    const double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
    const double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
    const double e5 = 0.4732592578721755, e6 =-0.509078520069735;
    const double e7 = 0.6772631491947646, e8 =-0.369912979092217;
    const double e9 = 0.06965131976970335;
    const double de1 = 1.0*e1;
    const double de2 = 2.0*e2;
    const double de3 = 3.0*e3;
    const double de4 = 4.0*e4;
    const double de5 = 5.0*e5;
    const double de6 = 6.0*e6;
    const double de7 = 7.0*e7;
    const double de8 = 8.0*e8;
    const double de9 = 9.0*e9;
//===========================================================================
    double x        = a*r;
    double eee      = exp(-x*x);
    double tt       = 1.0/(1.0+p*x);
    double gerfc_x  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt+e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
    double fgerfc_x = ((((((((de9*tt+de8)*tt+de7)*tt+de6)*tt+de5)*tt+de4)*tt+de3)*tt+de2)*tt+de1)*tt*tt*eee*p
                           +2.0*gerfc_x*x;
	printf("x: %g, eee: %g, tt: %g, gerfc_x: %g\n",x, eee, tt, gerfc_x);
    fgerfc[0]       = a*fgerfc_x;   // -d/dr ( erfc(a*r)  )
//===========================================================================
    return gerfc_x;
  }//end routine gerf
//===========================================================================
