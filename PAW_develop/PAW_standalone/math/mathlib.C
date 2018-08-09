#include "standard_include.h"
#include "mathlib.h"

/*===============================================================*/
/*  Inverse of a 3x3 */
/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/

void gethinv(double *hmat, double *hmati, double *deth, int iperd)

  /*===============================================================*/
{/*begin routine */
  /*===============================================================*/
  double vol;
  int i;
  /*===============================================================*/
  /* gets inverse, hmati, of the iperd dimensional matrix hmat */
  /* (stored as a 3 x 3) */

  *deth = 0.0;
  for(i=1;i<=9;i++){hmati[i]=0.0;}

  /*===============================================================*/
  /* Perd=3 */

  if (iperd == 3) {
    vol = (hmat[1] * (hmat[5] * hmat[9] - hmat[8] * hmat[6]) +
        hmat[4] * (hmat[8] * hmat[3] - hmat[2] * hmat[9]) +
        hmat[7] * (hmat[2] * hmat[6] - hmat[5] * hmat[3]));
    *deth = vol;
    hmati[1] = (hmat[5] * hmat[9] - hmat[8] * hmat[6]) / vol;
    hmati[5] = (hmat[1] * hmat[9] - hmat[7] * hmat[3]) / vol;
    hmati[9] = (hmat[1] * hmat[5] - hmat[4] * hmat[2]) / vol;
    hmati[4] = (hmat[7] * hmat[6] - hmat[4] * hmat[9]) / vol;
    hmati[2] = (hmat[3] * hmat[8] - hmat[2] * hmat[9]) / vol;
    hmati[7] = (hmat[4] * hmat[8] - hmat[7] * hmat[5]) / vol;
    hmati[3] = (hmat[2] * hmat[6] - hmat[3] * hmat[5]) / vol;
    hmati[8] = (hmat[7] * hmat[2] - hmat[8] * hmat[1]) / vol;
    hmati[6] = (hmat[3] * hmat[4] - hmat[6] * hmat[1]) / vol;
  }/*endif*/

  /*===============================================================*/
  /* Perd=2 */

  if (iperd == 2) {
    vol = hmat[1] * hmat[5] - hmat[4] * hmat[2];
    hmati[1] = hmat[5] / vol;
    hmati[5] = hmat[1] / vol;
    hmati[4] = -hmat[4] / vol;
    hmati[2] = -hmat[2] / vol;
    hmati[9] = 1. / hmat[9];
    *deth = vol * hmat[9];
  }/*endif*/

  /*===============================================================*/
  /* Perd=1,0,cluster_ewald */

  if(iperd <=1 || iperd==4) {
    *deth = hmat[1]*hmat[5]*hmat[9];
    hmati[1] = 1.0/hmat[1];
    hmati[5] = 1.0/hmat[5];
    hmati[9] = 1.0/hmat[9];
  }/*endif*/

  /*===============================================================*/
  /* Errors */

  if((*deth)==0.0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The present volume is zero.                 \n");
    printf("If this is not an error in your input data, \n");
    printf("contact technical support                  \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  } /*endif*/

  /*---------------------------------------------------------------*/
} /* gethinv */
/*===============================================================*/

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the distance between two points
//===================================================================
inline double dist(double dx, double dy, double dz)
//===================================================================
{ // begin routine
//===================================================================
  double result = sqrt(dx*dx + dy*dy + dz*dz);
  return result;
//===================================================================
} // end routine
//===================================================================

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the distance between two points
//===================================================================
inline double erfc_a_r_over_r(double r, double a, double gerfc, double fgerfc)
//===================================================================
{ // begin routine
//===================================================================
//  double result = 2.0*a*r/sqrt(M_PI_QI)*exp(-a*a*r*r) - erf(a*r);
  double result = r*fgerfc + gerfc - 1.0;
  return result;
//===================================================================
} // end routine
//===================================================================

//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================
// Generic error function
//===========================================================================
inline double gerfc(double r, double a, double *fgerfc){
//===========================================================================
#define _STRICT_ERFC_DER_OFF_
  double x        = a*r;
  double eee      = exp(-x*x);
  double tt       = 1.0/(1.0+PERFC*x);
  double gerfc_x  = ((((((((CERFC9*tt+CERFC8)*tt+CERFC7)*tt+CERFC6)*tt+CERFC5)*tt+CERFC4)*tt+CERFC3)*tt+CERFC2)*tt+CERFC1)*tt*eee;
#ifdef _STRICT_ERFC_DER_
    double fgerfc_x = ((((((((DCERFC9*tt+DCERFC8)*tt+DCERFC7)*tt+DCERFC6)*tt+DCERFC5)*tt+DCERFC4)*tt+DCERFC3)*tt+DCERFC2)*tt+DCERFC1)*tt*tt*eee*PERFC
                         +2.0*gerfc_x*x;
#else
    fgerfc[0]       = a*eee*PRE_ERFC;   // -d/dr ( erfc(a*r)  )
#endif
//===========================================================================
    return gerfc_x;
}//end routine gerf
//===========================================================================

