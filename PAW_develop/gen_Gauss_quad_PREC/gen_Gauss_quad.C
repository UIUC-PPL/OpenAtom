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
void gen_Gauss_quad(int n, double * Int_xpow, double * Int_xpow_log, double * Int_xpow_sgn, 
					int * Int_xpow_zero, double * wght, double * node, int iopt){
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
		poly[i].Norm  = 0.0;
		poly[i].c	  = new double [i+1];
		poly[i].a	  = new double [i+1];
		poly[i].O	  = new double [i+1];
		poly[i].p_at_nodes  = new double [n];
	}// end for i
		
	poly[0].c[0] = 1.0;
	poly[0].a[0] = 1.0;
	poly[0].O[0] = Int_xpow[0];
	poly[0].Norm = Int_xpow[0];
	poly[0].c[0] /= sqrt(fabs(poly[0].Norm));
	poly[0].a[0] /= sqrt(fabs(poly[0].Norm));
	poly[0].O[0] /= sqrt(fabs(poly[0].Norm));
	poly[0].Norm = 1.0;
	for (int k=1; k<n+1; k++) {
		poly[k].c[k] = 1.0; // set the last coeff of the naught
		for (int j=0; j<k; j++) {
			double Onow = 0.0;
			for (int l=0; l < j; l++) {
				Onow += poly[j].c[l]*poly[k].O[l];
			} // end for l
			poly[k].O[j] = Onow + poly[j].c[j]*Int_xpow[k+j];
			double tmp   = log(fabs(poly[k].O[j])) - log(fabs(poly[j].Norm));
			double tmpas = ((poly[k].O[j]*poly[j].Norm >= 0.0) ? 1.0 : -1.0);
			poly[k].c[j] = -tmpas*exp(tmp);
//			poly[k].c[j] = -poly[k].O[j]/poly[j].Norm;
		} // end for j
		double normtmp =  poly[k].c[k]*poly[k].c[k]*Int_xpow[2*k];
		for (int j=0; j<k; j++) {
			normtmp += poly[k].c[j]*poly[k].c[j]*poly[j].Norm + 2.0*poly[k].c[j]*poly[k].O[j];
		} // end for j
		double sgn_tmp = ((normtmp >= 0 )? 1: -1);
		double scale   = sgn_tmp/sqrt(fabs(normtmp));
		if (iopt == 1) {PRINTF("k = %d scale =%.10g\n", k, scale);}
		for (int j=0; j<k+1; j++) {
			poly[k].c[j] *= scale;
			poly[k].O[j] *= scale;
		} // end for j
		poly[k].Norm = 1.0;
	} // end for k

	for (int k=1; k<n+1; k++) {
		poly[k].a[k] = poly[k].c[k];
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

	double regular = 0.0;
	for (int i=0; i<n+1; i++) {
		a[i] /= a[n];
		if (fabs(a[i]) > 1.0e-5) { regular += log(fabs(a[i])); }
	} // end for i
	regular /= ((double) n);
	regular = exp(-regular);
	regular = 1.0;
	PRINTF("regular = %.10g\n",regular);

    double* AA      = new double [n*n];
    for (int i=0; i <n*n; i++) {AA[i] 			= 0.0;}
   	for (int i=0; i < n ; i++) {AA[i*n + n - 1] = -a[i]*regular;}
    for (int i=1; i < n ; i++) {AA[i*n + i - 1] = 1.0*regular;}

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
		B[k] = 1.0*regular;
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
		ALPHAR[i] /= regular;
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
//  Test the nodes

	for (int i=0; i<n; i++) {node[i] = ALPHAR[i];}
	sort(node, node + n);
	testnodes(n, node, a);

//==========================================================================
//  Evaluate all the polynomials at the nodes

	for (int j=0; j<n; j++) {poly[0].p_at_nodes[j] = poly[0].a[0];}

	for (int i=1; i<n+1; i++) {
		double * a_tmp = poly[i].a;
		for (int j=0; j<n; j++) {
			double p_tmp;
			horner(i+1, a_tmp, node[j], &p_tmp);
			poly[i].p_at_nodes[j] = p_tmp;
		} // end for j
	} // end for i

//==========================================================================
//	We construct the weights using the zeros and the derivatives of the polynomial

	double an1	= poly[n-1].a[n-1];
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
		double pn1v_log = log(fabs(pn1v));
		double pn1v_sgn = ((pn1v >= 0.0) ? 1.0: -1.0);
		double pnpv_log = log(fabs(pnpv));
		double pnpv_sgn = ((pnpv >= 0.0) ? 1.0: -1.0);
		double Normn1_log = log(fabs(Normn1));
		double Normn1_sgn = ((Normn1 >= 0.0) ? 1.0: -1.0);
		double an1_log = log(fabs(an1));
		double an1_sgn = ((an1 >= 0.0) ? 1.0: -1.0);
//		wght[i] = Normn1/(pn1v*pnpv);
		wght[i] = Normn1_sgn*pn1v_sgn*pnpv_sgn*exp(Normn1_log - pn1v_log - pnpv_log - an1_log)*an1_sgn;
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
		testgrid(n, wght, poly, iopt);
		PRINTF("=================================================\n\n");
	} else {
		testgrid(n, wght, poly, iopt);
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
	double tmp = a[n-1];
	for (int k=n-2; k>=0; k--){
		tmp = tmp*x + a[k];
	} // end for
	value[0] = tmp;
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
void testgrid(int n, double * w, POLY * poly, int iopt) {
//==========================================================================
//	Test the orthogonality of all the polynomials

	int ierr = 0;
	double errmax = 0.0;
	for (int k=0; k<n+1; k++) {
		int iup = ((k != n)? k+1: n); 
		double * p_k = poly[k].p_at_nodes;
		for (int kp=0; kp<iup; kp++) {
			double * p_kp = poly[kp].p_at_nodes;
			double result = 0.0;
			for (int i=0; i<n; i++) {
				result += w[i]*p_k[i]*p_kp[i];
			} //end for i
			double ans = ((k == kp )? 1.0:0.0);
			double err = fabs(result - ans);
			errmax = MAX(err,errmax);
			if (iopt == 1) {
				PRINTF("Result of ortho test for pair (%d, %d) is %.10g with error %.10g\n", k, kp, result, err);
			} // end if
			if (fabs(err) > 1.0e-7 && iopt == 0) {
				ierr++;
				PRINTF("Error out of range %.10g for ortho (%d, %d)\n", err, k, kp);
			} // end if
		}// end for kp
	} //end for k

	PRINTF("The maximum err in your orthogonality test is %g\n", errmax);
	if (ierr > 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("You have %d moments out of error range!\n", ierr);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

} //end routine
