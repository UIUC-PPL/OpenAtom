//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Compute the Gaussian quadature weights and nodes of order n for a general weighting 
//	function, given the first 2n moments of weighting function
//  Int_xpow (k) = \int_a^b dx w(x) x^k   k=0 ... 2n 
//==========================================================================
// include files

#include "standard_include.h"
#include "gen_Gauss_quad_entry.h"
#include "gen_Gauss_quad_local.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Controller for the Gaussian quadature weights and nodes of order n for a general weighting 
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
		for (int j=0; j<k; j++) {
			poly[k].a[j] = 0.0;
			for (int l=j; l < k; l++) {
				poly[k].a[j] += poly[l].a[j]*poly[k].c[l];
			} // end for l
		} // end for j
	} // end for k

	double *a = poly[n].a;
	double *c = poly[n].c;

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
        	PRINTF("Bad root: (%g,%g)\n",ALPHAR[i],ALPHAI[i]);
		}// end if
    }// end for i

	if (iii != 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("The roots of the polynomial must be real!\nAre the 'a' coefficients ill-posed?\n");
		double cmax = fabs(c[0]); double amax = fabs(a[0]);
		double cmin = fabs(c[0]); double amin = fabs(a[0]);
		for (int i=1; i<n+1; i++) {
			cmax = MAX(cmax,fabs(c[i]));
			amax = MAX(amax,fabs(a[i]));
			if (fabs(c[i]) > 0.0) { cmin = MIN(cmin,fabs(c[i]));}
			if (fabs(a[i]) > 0.0) { amin = MIN(amin,fabs(a[i]));}
		} // end for
		PRINTF("a range = [%.10g, %.10g], c range = [%.10g, %.10g],\n",amin, amax, cmin, cmax);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

//==========================================================================
// Test the nodes
	for (int i=0; i<n; i++) {node[i] = ALPHAR[i];}
	sort(node, node + n);
	testnodes(n, node, a);

//==========================================================================
//	We construct the weights using the zeros and the derivatives of the polynomial

	for (int i=0; i<n; i++) {node[i] = ALPHAR[i];}
	sort(node, node + n);

	double Normn1 = poly[n-1].Norm; 
	double *pnp = new double [n];
	for (int i=0; i<n; i++) {
		double tmp = ((double)i)+1.0;
		pnp[i] = tmp*a[i+1];
	} // end for
	double *pn1 = poly[n-1].a;

	for (int i=0; i<n; i++) {
		double pn1v = 0.0, pnpv = 0.0;
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
		testgrid(n, wght, Int_xpow, node, iopt);
		PRINTF("=================================================\n\n");
	} else {
		testgrid(n, wght, Int_xpow, node, iopt);
	}// end if	

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
// Test the nodes
//==========================================================================
void testnodes(int n, double * x, double * a) {
	int ierr = 0;
	double err = 0.0;
	for (int i=0; i<n; i++) {
		double zero = 0.0;
		horner(n+1, a, x[i], &zero);
		err = MAX(err,fabs(zero));
		if (fabs(zero) > 1.0e-8) {ierr++;}
	}//end for
	PRINTF("The maximum err in your zeros is %g\n", err);
	if (ierr > 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("You have %d zeros out of error range! Write a root refinement routine!\n", ierr);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if
} // end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Test the correctness
//==========================================================================
void testgrid(int n, double *w, double * Int_xpow, double *x, int iopt) {
	if (iopt == 1) {
		for (int power=0; power<2*n-1; power++) {
			double result = 0;
			for (int i= 0; i<n; i++) {
				result += w[i]*pow(x[i],power);
			} //end for
			PRINTF("%d  %.10g  %.10g\n", power, result, Int_xpow[power]);
		} //end for
	} // end if

	int ierr = 0;
	double errmax = 0.0;
	for (int power=0; power<2*n-1; power++) {
		double result = 0;
		for (int i= 0; i<n; i++) {
			result += w[i]*pow(x[i],power);
		} //end for
		double err = result - Int_xpow[power];
		if (fabs(Int_xpow[power]) > 1.0e-5) { err /= fabs(Int_xpow[power]);}
		errmax = MAX(err,errmax);
		if (fabs(err) > 1.0e-7) {
			ierr++;
			PRINTF("Error out of range %g for moment %d  %.10g  %.10g\n", err, power, result, Int_xpow[power]);
		} // end if
	} //end for
	PRINTF("The maximum err in your moments is %g\n", errmax);
	if (ierr > 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("You have %d moments out of error range!\n", ierr);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if
} //end routine
