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
extern "C" {
#include <quadmath.h>
}
#include "gen_Gauss_quad_driver_entry.h"
#include "gen_Gauss_quad_driver_local.h"
#include "gen_Gauss_quad_entry.h"

#define _LOG_METHOD_

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller for generalized Gaussian quadrature using 
//				   companion matrix 
//==========================================================================
void gen_Gauss_quad_driver (int type, int n, int iopt, double * nodes_dbl, double * wghts_dbl, int * ierr_zero_out, int * ierr_ortho_out){
//==========================================================================
//  Local variables
	int ierr_zero;
	int ierr_ortho;

//==========================================================================
//  Tell the user what you are doing 

    PRINTF("\n");
    PRINT_LINE_STAR;
    PRINTF("Generating generalized Gaussian quadrature for weight type = %d and order n = %d\n", type, n);
    PRINT_LINE_DASH;
    PRINTF("\n");
	
//==========================================================================
//	Fetch the integrals of powers of x over the desired weight 
//  Hard coding a for now

	__float128 * Int_xpow     = new __float128 [2*n + 1]; 
	__float128 * Int_xpow_log = new __float128 [2*n + 1]; 
	__float128 * Int_xpow_sgn = new __float128 [2*n + 1]; 
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

#ifdef _LOG_METHOD_
	for (int i=0; i<2*n+1; i++) {
		if (Int_xpow_zero[i] == 0) {Int_xpow[i] = Int_xpow_sgn[i]*exp(Int_xpow_log[i]);}
	} // end for i
#endif

//==========================================================================
//	Fetch the weights and nodes 

	__float128 * w = new __float128 [n]; // Gauss quad weight
	__float128 * x = new __float128 [n]; // Gauss quad node

	gen_Gauss_quad(n, Int_xpow, Int_xpow_log, Int_xpow_sgn, Int_xpow_zero, w, x, iopt, &ierr_zero, &ierr_ortho);

	ierr_zero_out[0] = ierr_zero; 
	ierr_ortho_out[0] = ierr_ortho; 

	for (int i=0; i<n; i++) {
		nodes_dbl[i] = ((double) x[i]);
		wghts_dbl[i] = ((double) w[i]);
	}// end for i

	delete [] Int_xpow;
	delete [] Int_xpow_log;
	delete [] Int_xpow_sgn;
	delete [] Int_xpow_zero;

//==========================================================================
//	Check the weights and nodes 

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
void fetch_Int_xpow_uniform(int n, __float128 * Int_xpow, __float128 * Int_xpow_log, __float128 * Int_xpow_sgn, int * Int_xpow_zero){
	for (int i=0; i<2*n+1; i++) {
		if (i%2 == 0) {
			__float128 ai = (__float128) i;
			Int_xpow[i] = 2.q/(1.q + ai);
			Int_xpow_log[i] = log(2.q) - log(1.q + ai);
			Int_xpow_sgn[i] = 1.q;
			Int_xpow_zero[i] = 0;
		} else {
			Int_xpow[i] = 0.q;
			Int_xpow_sgn[i] = 1.q;
			Int_xpow_zero[i] = 1;
		} // end if	
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [-inf,inf]
//==========================================================================
void fetch_Int_xpow_Gaussfull(int n, __float128 * Int_xpow, __float128 * Int_xpow_log, __float128 * Int_xpow_sgn, int * Int_xpow_zero){
	Int_xpow[0] = sqrt(M_PI_QI);
	Int_xpow_log[0] = 0.5q*log(M_PI_QI);
	Int_xpow_sgn[0] = 1.q;
	Int_xpow_zero[0] = 0;
	for (int i=1; i<2*n+1; i++) {
		if (i%2 == 1) {
			Int_xpow[i] = 0.q;
			Int_xpow_sgn[i] = 1.q;
			Int_xpow_zero[i] = 1;
		} else {
			__float128 tmp = 0.5q*(((__float128) i)-1.q);
			Int_xpow[i] = Int_xpow[i-2]*tmp;
			Int_xpow_log[i] = Int_xpow_log[i-2] + log(tmp);
			Int_xpow_sgn[i] = 1.q;
			Int_xpow_zero[i] = 0;
		} // end if
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [0,inf]
//==========================================================================
void fetch_Int_xpow_Gausshalf(int n, __float128 * Int_xpow, __float128 * Int_xpow_log, __float128 * Int_xpow_sgn, int * Int_xpow_zero){
	Int_xpow[0]     = 0.5q*sqrt(M_PI_QI);
	Int_xpow_log[0] = 0.5q*log(M_PI_QI) + log(0.5q);
	Int_xpow_sgn[0] = 1.q;
	Int_xpow_zero[0] = 0;

	Int_xpow[1]      = 0.5q;
	Int_xpow_log[1]  = log(0.5q);
	Int_xpow_sgn[1]  = 1.q;
	Int_xpow_zero[1] = 0;

	for (int i=2; i<2*n+1; i++) {
		__float128 tmp = 0.5q*(((__float128) i)-1.q);
		Int_xpow[i]      = Int_xpow[i-2]*tmp;
		Int_xpow_log[i]  = Int_xpow_log[i-2] + log(tmp);
		Int_xpow_sgn[i]  = 1.q;
		Int_xpow_zero[i] = 0;
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Check the Legrendre nodes, although there is a theorem that the nodes lie
//   on [-1,1], make sure they do.
//==========================================================================
void check_nodes_uniform(int n, __float128 * x){
	int ierr = 0;
	for (int i=0; i<n; i++) {
		if (x[i] > 1.q || x[i] < -1.q) {
			ierr++;
			PRINTF("x[%d] = %.10Lg out of range!\n",i,(long double)x[i]);
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
void check_nodes_Gausshalf(int n, __float128 * x){
	int ierr = 0;
	for (int i=0; i<n; i++) {
		if (x[i] < 0.q) {
			ierr++;
			PRINTF("x[%d] = %.10Lg out of range!\n",i,(long double)x[i]);
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

