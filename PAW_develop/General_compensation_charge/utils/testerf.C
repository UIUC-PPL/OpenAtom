//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  This program computes the compensation charg energy for a frozen Gaussian
//	core density
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "standard_include.h"
inline double erfc_a_r_over_r(double , double );
inline double gerf(double);
inline double gerfc(double);
inline double gerfc_a_r_over_r(double , double );

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller to manage compensation charge energy
//==========================================================================
int main (int argc, char *argv[]){
	int N = 10000000;
	double * r = new double [N];
	double * y = new double [N];
	long int seed = 2351361;
	srand48(seed);
	for (int i=0; i<N; i++) { 
		r[i] = drand48()*3.0;
	} // endfor

	clock_t start, end;

	start = clock();
	for (int i=0; i<N; i++) { 
		y[i] = 1.0/r[i];
	} // endfor
	end = clock();
	PRINTF("10M 1/r finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

	start = clock();
	for (int i=0; i<N; i++) { 
		y[i] = erfc(r[i])/r[i];
	} // endfor
	end = clock();
	PRINTF("10M erfc(r)/r finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

	start = clock();
	for (int i=0; i<N; i++) { 
		y[i] = erfc_a_r_over_r(r[i], 1.0);
	} // endfor
	end = clock();
	PRINTF("10M my_func finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

	start = clock();
	for (int i=0; i<N; i++) { 
		y[i] = gerf(r[i])/r[i];
	} // endfor
	end = clock();
	PRINTF("10M gerf finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

	start = clock();
	for (int i=0; i<N; i++) { 
		y[i] = gerfc_a_r_over_r(r[i], 1.0);
	} // endfor
	end = clock();
	PRINTF("10M g_func finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);
}//end routine
//==========================================================================

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the distance between two points
//===================================================================
inline double erfc_a_r_over_r(double r, double a)
//===================================================================
{ // begin routine
//===================================================================
    double result = 2.0*a*r/sqrt(M_PI_QI)*exp(-a*a*r*r) - erf(a*r);
    return result;
//===================================================================
} // end routine

/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
inline double gerfc(double x)
{
  /*===============================================================*/
  /*  Local variables */

  const double p = 0.3614;
  const double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  const double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  const double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  const double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  const double e9 = 0.06965131976970335;
  double eee,tt,gerfc;

  /*===============================================================*/
  /* Calculate the error function */

  eee    = exp(-x*x);
  tt     = 1.0/(1.0+p*x);
  gerfc  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
            +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;

  /*===============================================================*/
  return gerfc;
  /*===============================================================*/
}/* end function */
/*===============================================================*/

/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
inline double gerf(double x)
{
  /*===============================================================*/
  /*  Local variables */

  const double p=0.3614;
  const double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  const double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  const double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  const double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  const double e9 = 0.06965131976970335;
  double eee,tt,gerf;

  /*===============================================================*/
  /* Calculate the error function */

  eee    = exp(-x*x);
  tt     = 1.0/(1.0+p*x);
  gerf   = 1.0 - ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
            +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;

  /*===============================================================*/
  return gerf;
  /*===============================================================*/
}/* end function */
/*===============================================================*/

/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/
inline double gerfc_a_r_over_r(double r, double a)
{
  /*===============================================================*/
  /*  Local variables */

  const double p=0.3614;
  const double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  const double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  const double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  const double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  const double e9 = 0.06965131976970335;
  const double pre = 2.0/sqrt(M_PI_QI);
  double eee,tt,gerf;

  /*===============================================================*/
  /* Calculate the error function */
  double x = r*a;
  eee    = exp(-x*x);
  tt     = 1.0/(1.0+p*x);
  gerf   = 1.0 - ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
            +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
  double result = pre*x*eee - gerf;

// double result = 2.0*a*r/sqrt(M_PI_QI)*exp(-a*a*r*r) - erf(a*r);

  /*===============================================================*/
  return result;
  /*===============================================================*/
}/* end function */
/*===============================================================*/
