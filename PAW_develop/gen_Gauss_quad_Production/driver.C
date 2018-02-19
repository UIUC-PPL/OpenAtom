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

#include "driver.h"
#include "gen_Gauss_quad_entry.h"

int main(int, char **);
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller for generalized Gaussian quadrature using 
//				   companion matrix 
//==========================================================================
int main (int argc, char *argv[]){
//==========================================================================
// Local variables
	int type;    // type of weighting function
	int n;       // power = order of Gauss quadrature
	int iopt;    // output option
	
//==========================================================================
// read in the power and the type of weighting function
// 	 type 0: uniform on -1 to 1
//	 type 1: Gaussian on -infty to infty
//	 type 2: Gaussian on 0 to infty
//	 power n: > 0
//	 iopt = 0 is quiet, and iopt = 1 is verbose
//=========================================================================

	if(argc < 4) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("Insufficient parameters specified\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	}/*endif*/

	sscanf(argv[1],"%d", &type);
	sscanf(argv[2],"%d", &n);
	sscanf(argv[3],"%d", &iopt);

	if (n < 1 || type < 0 || type > 2 || iopt < 0 || iopt > 1) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("Parameters out of range\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	}/*endif*/
		
	PRINTF("\n");
	PRINT_LINE_STAR;
	PRINTF("Generating generalized Gaussian quadrature for weight type = %d and order n = %d\n", type, n);
	PRINT_LINE_DASH;
	PRINTF("\n");

//==========================================================================
//	Fetch the integrals of powers of x over the desired weight 
//  Hard coding a for now

	double * Int_xpow     = new double [2*n + 1]; 
	double * Int_xpow_log = new double [2*n + 1]; 
	double * Int_xpow_sgn = new double [2*n + 1]; 
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

	double * w = new double [n]; // Gauss quad weight
	double * x = new double [n]; // Gauss quad node

	gen_Gauss_quad(n, Int_xpow, Int_xpow_log, Int_xpow_sgn, Int_xpow_zero, w, x, iopt);

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
	return 1;
} // end routine
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over uniform weight on [-1,1]
//==========================================================================
void fetch_Int_xpow_uniform(int n, double * Int_xpow, double * Int_xpow_log, double * Int_xpow_sgn, int * Int_xpow_zero){
	for (int i=0; i<2*n+1; i++) {
		if (i%2 == 0) {
			double ai = (double) i;
			Int_xpow[i] = 2.0/(1.0 + ai);
			Int_xpow_log[i] = log(2.0) - log(1.0 + ai);
			Int_xpow_sgn[i] = 1.0;
			Int_xpow_zero[i] = 0;
		} else {
			Int_xpow[i] = 0.0;
			Int_xpow_sgn[i] = 1.0;
			Int_xpow_zero[i] = 1;
		} // end if	
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [-inf,inf]
//==========================================================================
void fetch_Int_xpow_Gaussfull(int n, double * Int_xpow, double * Int_xpow_log, double * Int_xpow_sgn, int * Int_xpow_zero){
	Int_xpow[0] = sqrt(M_PI);
	Int_xpow_log[0] = 0.5*log(M_PI);
	Int_xpow_sgn[0] = 1.0;
	Int_xpow_zero[0] = 0;
	for (int i=1; i<2*n+1; i++) {
		if (i%2 == 1) {
			Int_xpow[i] = 0.0;
			Int_xpow_sgn[i] = 1.0;
			Int_xpow_zero[i] = 1;
		} else {
			double tmp = 0.5*(((double) i)-1.0);
			Int_xpow[i] = Int_xpow[i-2]*tmp;
			Int_xpow_log[i] += log(tmp);
			Int_xpow_sgn[i] = 1.0;
			Int_xpow_zero[i] = 0;
		} // end if
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [0,inf]
//==========================================================================
void fetch_Int_xpow_Gausshalf(int n, double * Int_xpow, double * Int_xpow_log, double * Int_xpow_sgn, int * Int_xpow_zero){
	Int_xpow[0]     = 0.5*sqrt(M_PI);
	Int_xpow_log[0] = 0.5*log(M_PI) + log(0.5);
	Int_xpow_sgn[0] = 1.0;
	Int_xpow_zero[0] = 0;

	Int_xpow[1]      = 0.5;
	Int_xpow_log[1]  = log(0.5);
	Int_xpow_sgn[1]  = 1.0;
	Int_xpow_zero[1] = 0;

	for (int i=2; i<2*n+1; i++) {
		double tmp = 0.5*(((double) i)-1.0);
		Int_xpow[i]      = Int_xpow[i-2]*tmp;
		Int_xpow_log[i]  = Int_xpow_log[i-2] + log(tmp);
		Int_xpow_sgn[i]  = 1.0;
		Int_xpow_zero[i] = 0;
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Check the Legrendre nodes, although there is a theorem that the nodes lie
//   on [-1,1], make sure they do.
//==========================================================================
void check_nodes_uniform(int n, double * x){
	int ierr = 0;
	for (int i=0; i<n; i++) {
		if (x[i] > 1.0 || x[i] < -1.0) {
			ierr++;
			PRINTF("x[%d] = %.10g out of range!\n",i,x[i]);
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
void check_nodes_Gausshalf(int n, double * x){
	int ierr = 0;
	for (int i=0; i<n; i++) {
		if (x[i] < 0.0) {
			ierr++;
			PRINTF("x[%d] = %.10g out of range!\n",i,x[i]);
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

