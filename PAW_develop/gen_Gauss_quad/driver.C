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
//#define _DEBUG_

//==========================================================================
// POLY structure
typedef struct POLY{
        int 	order;		// order
        double *c;          // c matrix 
        double *a;          // a matrix
        double *O;          // O integral 
		double Norm;		// N_ii
}POLY;

//==========================================================================
// function prototypes 

int main (int , char **);
void gen_Gauss_quad(int, double *, double *, double *, int);

void fetch_Int_xpow_uniform(int, double *); 
void fetch_Int_xpow_Gaussfull(int, double *); 
void fetch_Int_xpow_Gausshalf(int, double *);
void horner(int, double *, double, double *);
void testgrid(int, double *, double *, double *);
extern "C" void dggev( char* JOBVL,  char* JOBVR,  int* N,
                       double* AA,  int* LDA,  double* B,  int* LDB,
                      double* ALPHAR, double* ALPHAI, double* BETA,
                      double* VL,  int* LDVL, double* VR,  int* LDVR,
                      double* WORK,  int* LWORK, int* INFO);

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
	
//==========================================================================
// read in the power and the type of weighting function
// 	 type 0: uniform on -1 to 1
//	 type 1: Gaussian on -infty to infty
//	 type 2: Gaussian on 0 to infty
//	 power n: > 0
//=========================================================================

	if(argc < 2) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("No parameters specified\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	}/*endif*/

	sscanf(argv[1],"%d", &type);
	sscanf(argv[2],"%d", &n);

	if (n < 1 || type < 0 || type > 2) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("Parameters out of range\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	}/*endif*/
		
	PRINTF("\n");
	PRINT_LINE_STAR;
	PRINTF("Generating generalized Gaussian quadrature for weight type = %d and order n = %d\n", type, n);
	PRINT_LINE_STAR;
	PRINTF("\n");

//==========================================================================
//	Fetch the integrals of powers of x over the desired weight 
//  Hard coding a for now

	double * Int_xpow = new double [2*n + 1]; 
	switch(type) {
		case 0: // w(x) = 1; [-1, 1] 
			fetch_Int_xpow_uniform(n, Int_xpow); break;
		case 1: // w(x) = e^(-x^2); [-infty,infty]
			fetch_Int_xpow_Gaussfull(n, Int_xpow); break;
		case 2: // w(x) = e^(-x^2); [-infty,infty]
			fetch_Int_xpow_Gausshalf(n, Int_xpow); break;
		default:
			PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			PRINTF("Parameters out of range\n");
			PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			FFLUSH(stdout);
			EXIT(1);
	}//end switch

	double * w = new double [n]; // Gauss quad weight
	double * x = new double [n]; // Gauss quad node
	int iopt = 1; // test and output option

	gen_Gauss_quad(n, Int_xpow, w, x, iopt);

	delete [] Int_xpow;

//==========================================================================
	return 1;
} // end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Compute the Gaussian quadature weights and nodes of order n for a general weighting 
//	function, given the first 2n moments of weighting function
//  Int_xpow (k) = \int_a^b dx w(x) x^k   k=0 ... 2n 
//==========================================================================
void gen_Gauss_quad(int n, double * Int_xpow, double * wght, double * node, int iopt){
//==========================================================================
//	First: Compute the polynomial coefficients a of the orthogonal polynomials over the weight using Gram Schmidt
//  	We provide a function for each type which computes the integral over the
//		weighting function of x^k (first 2n moments of w(x) stored in Int_xpow)
//		We have a structure Opoly[i].a[k] where a[i+1] is the coefficient of the
//		i+1 th order orthogonal polynomial
//		Opoly[0].a[0] = 1 starts the Gram Schmidt
//		You make Opoly[1]...Opoly[n], you only create i = 2, 3, ..., n
//		for (int i=1; i<=n; i++) { GramSmit (i,Opoly);}
//		Opoly[n].a[0..n] are the desired coefficients
//		GramSmit just requires the evaluation of overlaps over the weighting function
//		for which we have derived a general formula involving the moments which is Int_xpow.

	POLY * poly = new POLY [n + 1];

	for (int i=0; i<n+1; i++) {
		poly[i].order = i;
		poly[i].c	  = new double [i+1];
		poly[i].a	  = new double [i+1];
		poly[i].O	  = new double [i+1];
		poly[i].Norm  = 0.0;
	}// end for i
		
	poly[0].c[0] = 1.0;
	poly[0].a[0] = 1.0;
	poly[0].O[0] = Int_xpow[0];
	poly[0].Norm = Int_xpow[0];
	for (int k=1; k<n+1; k++) {
		poly[k].c[k] = 1.0; // set the last coeff
		for (int j=0; j<k; j++) {
			double Onow = 0.0;
			for (int l=0; l < j; l++) {
				Onow += poly[j].c[l]*poly[k].O[l];
			} // end for l
			poly[k].O[j] = Onow + poly[j].c[j]*Int_xpow[k+j];
			poly[k].c[j] = -poly[k].O[j]/poly[j].Norm;
		} // end for j
		double normtmp =  poly[k].c[k]*Int_xpow[2*k];
		for (int j=0; j<k; j++) {
			normtmp += poly[k].c[j]*poly[k].c[j]*poly[j].Norm + 2.0*poly[k].c[j]*poly[k].O[j];
		} // end for j
		poly[k].Norm = normtmp;
	} // end for k

	for (int k=1; k<n+1; k++) {
		poly[k].a[k] = 1.0;
//		PRINTF("Norm: %d: %g \n",k,poly[k].Norm);
		for (int j=0; j<k; j++) {
			poly[k].a[j] = 0.0;
			for (int l=j; l < k; l++) {
				poly[k].a[j] += poly[l].a[j]*poly[k].c[l];
			} // end for l
//			PRINTF("%d %d : %g \n",k,j,poly[k].a[j]);
		} // end for j
	} // end for k

	double *a = poly[n].a;

//==========================================================================
//	Construct the companion matrix, A, using the polynomial coefficients a

    double* AA      = new double [n*n];
    for (int i=0; i <n*n; i++) {AA[i] 			= 0.0;}
   	for (int i=0; i < n ; i++) {AA[i*n + n - 1] = -a[i];}
    for (int i=1; i < n ; i++) {AA[i*n + i - 1] = 1.0;}

//==========================================================================
//	Find the zeros of the n-th orthogonal polynomial by diagonalizing its companion matrix	

	int N 	 		= n;
	int LDA 		= n; 
	double* B 		= new double [n*n]; 
	int LDB 		= n;
    double* ALPHAR 	= new double [n]; 
	double* ALPHAI	= new double [n]; 
	double* BETA	= new double [n];
    double* VL		= new double [n]; 
	int LDVL		= 1; 
	double* VR		= new double [n]; 
	int LDVR		= 1;
	int LWORK		= -1; 
    double* WORK    = NULL;	
	int INFO		= 1;	
	
	for (int i=0; i<n*n; i++) {B[i] = 0.0;}

	for (int i=0; i<n; i++) {
		int k = i*n + i;
		B[k] = 1.0;
	}// end for i

	double tmp = 1.0;
 	dggev("N", "N", &N, AA, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       		VL, &LDVL, VR, &LDVR, &tmp, &LWORK, &INFO);
	
	LWORK = (int)tmp;
	WORK  = new double [LWORK];
 	dggev("N", "N", &N, AA, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       		VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);

	if (INFO != 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("dggev failed for unknown reason!\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

	int iii = 0;
    for (int i=0; i<n; i++) {
		if (ALPHAI[i]*ALPHAI[i] > 1.0e-16) {
			iii++;
        	PRINTF("Bad root: (%g,%g) ",ALPHAR[i],ALPHAI[i]);
		}// end if
    }// end for i

	if (iii != 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("The roots of the polynomial must be real!\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

//==========================================================================
//	We construct the weights using the zeros and the derivatives of the polynomial

	for (int i=0; i<n; i++) {node[i] = ALPHAR[i];}
	sort(node, node + n);

	double Normn1 = poly[n-1].Norm; 
	double *pnp = new double [n];
	for (int i=0; i<n; i++) {
		pnp[i] = ((double)i+1.0)*a[i+1];
	} // end for
	double *pn1 = poly[n-1].a;

	for (int i=0; i<n; i++) {
		double pn1v = 0, pnpv = 0;
		horner(n, pn1, node[i], &pn1v);
		horner(n, pnp, node[i], &pnpv);
		wght[i] = Normn1/(pn1v*pnpv);
	}// end for
	delete [] pnp;

//==========================================================================
//	Test the weights and nodes	

	if (iopt == 1) {
		PRINTF("=================================================\n");
		PRINTF("Generalized Gaussian quadrature weights and nodes:\n");
		PRINTF("-------------------------------------------------\n");
    	for (int i=0; i<n; i++) {
    	    PRINTF("node[%d] %.10g   wght[%d] %.10g\n",i, node[i], i, wght[i]);
    	}// end for i
		PRINTF("=================================================\n\n");

		PRINTF("=================================================\n");
		PRINTF("Test integrals of 2n-1 moments using the quadrature:\n");
		PRINTF("-------------------------------------------------\n");
		testgrid(n, wght, Int_xpow, node);
		PRINTF("=================================================\n\n");
	} // end if

//==========================================================================
//  Clean up the memory
	
	for (int i=0; i<n+1; i++) {
		delete [] poly[i].c;
		delete [] poly[i].a;
		delete [] poly[i].O;
	}// end for i
	delete [] poly;

	delete [] AA;
	delete [] B; 
    delete [] ALPHAR; 
	delete [] ALPHAI; 
	delete [] BETA;
    delete [] VL; 
	delete [] VR; 
	delete [] WORK;

//==========================================================================
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over uniform weight on [-1,1]
//==========================================================================
void fetch_Int_xpow_uniform(int n, double * Int_xpow){
	for (int i=0; i<2*n+1; i++) {
		int tmp = (i % 2 == 0) ? 1:-1;
		Int_xpow[i] = (1.0 + (double) tmp)/(1.0 + (double) i);
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [-inf,inf]
//==========================================================================
void fetch_Int_xpow_Gaussfull(int n, double * Int_xpow){ 
	Int_xpow[0] = sqrt(M_PI);
	for (int i=1; i<2*n+1; i++) {
		if (i%2 == 1) Int_xpow[i] = 0;
		else Int_xpow[i] = Int_xpow[i-2]*((double)i-1.0)/2.0;
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Integrals of powers of x over Gaussian weight on [0,inf]
//==========================================================================
void fetch_Int_xpow_Gausshalf(int n, double * Int_xpow){
	Int_xpow[0] = sqrt(M_PI)/2.0;
	Int_xpow[1] = 1.0/2.0;
	for (int i=2; i<2*n+1; i++) {
		Int_xpow[i] = Int_xpow[i-2]*((double)i-1.0)/2.0;
	}// end for i
}//end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Valuating the p_n(x_i) using horner's method.
//==========================================================================
void horner(int n, double *a, double x, double *value){
	*value = a[n-1];
	for (int k=n-2; k>=0; k--){
		*value = (*value)*x + a[k];
	} // end for
}// end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Test the correctness
//==========================================================================
void testgrid(int n, double *w, double * Int_xpow, double *x) {
	for (int power=0; power<11; power++) {
		double result = 0;
		for (int i= 0; i<n; i++) {
			result += w[i]*pow(x[i],power);
		} //end for
	PRINTF("%d  %.10g  %.10g\n", power, result, Int_xpow[power]);
	} //end for
} //end routine
