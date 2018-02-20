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
// include files

#include "standard_include.h"

#include "gen_Gauss_quad_driver.h"
#include "gen_Gauss_quad_entry.h"

#define _LOG_METHOD_

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller for generalized Gaussian quadrature using 
//				   companion matrix 
//==========================================================================
void gen_Gauss_quad_driver (int type, int n, double alpha, int iopt, double * nodes_dbl, double * wghts_dbl){
//==========================================================================
// Local variables
	int ierr_zero; // bad nodes
	int ierr_ortho; // bad orthogonality

//==========================================================================
//	Fetch the integrals of powers of x over the desired weight 
//  Hard coding a for now

	long double * Int_xpow     = new long double [2*n + 1]; 
	long double * Int_xpow_log = new long double [2*n + 1]; 
	long double * Int_xpow_sgn = new long double [2*n + 1]; 
	int * Int_xpow_zero   = new int [2*n + 1]; 
	switch(type) {
		case 0: // w(x) = 1; [-1, 1] 
			fetch_Int_xpow_uniform(n, Int_xpow, Int_xpow_log, Int_xpow_sgn, Int_xpow_zero); break;
		case 1: // w(x) = e^(-x^2); [-infty,infty]
			fetch_Int_xpow_Gaussfull(n, Int_xpow, Int_xpow_log, Int_xpow_sgn, Int_xpow_zero); break;
		case 2: // w(x) = e^(-x^2); [-infty,infty]
			fetch_Int_xpow_Gausshalf(n, Int_xpow, Int_xpow_log, Int_xpow_sgn, Int_xpow_zero); break;
		default:
			PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			PRINTF("Parameters out of range\n");
			PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			FFLUSH(stdout);
			EXIT(1);
	}//end switch

	long double * w = new long double [n]; // Gauss quad weight
	long double * x = new long double [n]; // Gauss quad node

#ifdef _LOG_METHOD_
	for (int i=0; i<2*n+1; i++) {
		if (Int_xpow_zero[i] == 0) {Int_xpow[i] = Int_xpow_sgn[i]*exp(Int_xpow_log[i]);}
	} // end for i
#endif

	gen_Gauss_quad(n, Int_xpow, Int_xpow_log, Int_xpow_sgn, Int_xpow_zero, w, x, iopt, &ierr_zero, &ierr_ortho);

	for (int i=0; i<n; i++) {
		nodes_dbl[i] = ((double) x[i])/alpha;
		wghts_dbl[i] = ((double) w[i])/alpha;
	}// end for i

	delete [] Int_xpow;
	delete [] Int_xpow_log;
	delete [] Int_xpow_sgn;
	delete [] Int_xpow_zero;

	switch(type) {
		case 0: // w(x) = 1; [-1, 1] 
			check_nodes_uniform(n, x); break;
		case 1: // w(x) = e^(-x^2); [-infty,infty]
				break;
		case 2: // w(x) = e^(-x^2); [-infty,infty]
			check_nodes_Gausshalf(n, x); break;
		default:
			PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			PRINTF("Parameters out of range\n");
			PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			FFLUSH(stdout);
			EXIT(1);
	}//end switch

//==========================================================================
//	Tell the user we are done

	PRINTF("\n");
	PRINT_LINE_DASH;
	PRINTF("Completed generalized Gaussian quadrature for weight type = %d and order n = %d\n", type, n);
	PRINT_LINE_STAR;

//==========================================================================
	delete [] x;
	delete [] w;
} // end routine
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over uniform weight on [-1,1]
//==========================================================================
void fetch_Int_xpow_uniform(int n, long double * Int_xpow, long double * Int_xpow_log, long double * Int_xpow_sgn, int * Int_xpow_zero){
	for (int i=0; i<2*n+1; i++) {
		if (i%2 == 0) {
			long double ai = (long double) i;
			Int_xpow[i] = 2.0L/(1.0L + ai);
			Int_xpow_log[i] = log(2.0L) - log(1.0L + ai);
			Int_xpow_sgn[i] = 1.0L;
			Int_xpow_zero[i] = 0;
		} else {
			Int_xpow[i] = 0.0L;
			Int_xpow_sgn[i] = 1.0L;
			Int_xpow_zero[i] = 1;
		} // end if	
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [-inf,inf]
//==========================================================================
void fetch_Int_xpow_Gaussfull(int n, long double * Int_xpow, long double * Int_xpow_log, long double * Int_xpow_sgn, int * Int_xpow_zero){
	Int_xpow[0] = sqrt(M_PI_QI);
	Int_xpow_log[0] = 0.5L*log(M_PI_QI);
	Int_xpow_sgn[0] = 1.0L;
	Int_xpow_zero[0] = 0;
	for (int i=1; i<2*n+1; i++) {
		if (i%2 == 1) {
			Int_xpow[i] = 0.0L;
			Int_xpow_sgn[i] = 1.0L;
			Int_xpow_zero[i] = 1;
		} else {
			long double tmp = 0.5L*(((long double) i)-1.0L);
			Int_xpow[i] = Int_xpow[i-2]*tmp;
			Int_xpow_log[i] = Int_xpow_log[i-2] + log(tmp);
			Int_xpow_sgn[i] = 1.0L;
			Int_xpow_zero[i] = 0;
		} // end if
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [0,inf]
//==========================================================================
void fetch_Int_xpow_Gausshalf(int n, long double * Int_xpow, long double * Int_xpow_log, long double * Int_xpow_sgn, int * Int_xpow_zero){
	Int_xpow[0]     = 0.5L*sqrt(M_PI_QI);
	Int_xpow_log[0] = 0.5L*log(M_PI_QI) + log(0.5L);
	Int_xpow_sgn[0] = 1.0L;
	Int_xpow_zero[0] = 0;

	Int_xpow[1]      = 0.5L;
	Int_xpow_log[1]  = log(0.5L);
	Int_xpow_sgn[1]  = 1.0L;
	Int_xpow_zero[1] = 0;

	for (int i=2; i<2*n+1; i++) {
		long double tmp = 0.5L*(((long double) i)-1.0L);
		Int_xpow[i]      = Int_xpow[i-2]*tmp;
		Int_xpow_log[i]  = Int_xpow_log[i-2] + log(tmp);
		Int_xpow_sgn[i]  = 1.0L;
		Int_xpow_zero[i] = 0;
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Check the Legrendre nodes, although there is a theorem that the nodes lie
//   on [-1,1], make sure they do.
//==========================================================================
void check_nodes_uniform(int n, long double * x){
	int ierr = 0;
	for (int i=0; i<n; i++) {
		if (x[i] > 1.0L || x[i] < -1.0L) {
			ierr++;
			PRINTF("x[%d] = %.10Lg out of range!\n",i,x[i]);
		} // end if
	}// end for
	if (ierr > 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("You have %d nodes out of range!\n", ierr);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if
} // end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Check the Half space Gaussian nodes, although there is a theorem that the
//    nodes lie on [0,\infty], check!
//==========================================================================
void check_nodes_Gausshalf(int n, long double * x){
	int ierr = 0;
	for (int i=0; i<n; i++) {
		if (x[i] < 0.0L) {
			ierr++;
			PRINTF("x[%d] = %.10Lg out of range!\n",i,x[i]);
		} // end if
	}// end for
	if (ierr > 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("You have %d nodes out of range!\n", ierr);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if
} // end routine

