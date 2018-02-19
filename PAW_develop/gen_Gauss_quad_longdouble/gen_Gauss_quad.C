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
void gen_Gauss_quad(int n, long double * Int_xpow, long double * Int_xpow_log, long double * Int_xpow_sgn, 
					int * Int_xpow_zero, long double * wght, long double * node, int iopt, int * ierr_zero, int * ierr_ortho){
//==========================================================================
//	First: Compute the polynomial coefficients "a" of the nth orthogonal polynomials over the weight using Gram Schmidt
//  	        User provides table ofc the first 2n moments of w(x) stored in Int_xpow
//		A datas tructure Ooly[i].a[k] where a[i+1] is the coefficient of thei+1 th order orthogonal polynomial
//		    Opoly[0].a[0] = 1 starts the Gram Schmidt
//		    Create poly[1]...poly[n], you only create k = 1, 2, ..., n by GS
//		    Opoly[n].a[0..n] are the desired coefficients
//              In more detail, Gramchmidt only requires the evaluation of moments over the weighting function stored as Int_xpos.
//              it is peformed in the c-rep -the expansion of the polynomial in lower ortho polynomials with
//              last term c_nx^n.  The a-rep poly = \sum_k a_k x^k is created from the c-rep via an analytic transformation   

//-------------------------------------------------------------------------
// Constructor for poly
	POLY * poly = new POLY [n + 1];
	for (int i=0; i<n+1; i++) {
	    poly[i].order       = i;
	    poly[i].Norm        = 0.0L;
	    poly[i].c	        = new long double [i+1];
	    poly[i].a	        = new long double [i+1];
	    poly[i].O	        = new long double [n+1];  //upper triangular matrix so over allocating - not worth trouble to jazz up
        poly[i].p_at_nodes  = new long double [n];
	}// end for i

//------------------------------------------------------------
// Set up 0th order poly and normalize it and its moments
	poly[0].c[0] = 1.0L;
	poly[0].a[0] = 1.0L;
	for(int i=0;i<n+1;i++){poly[0].O[i] = Int_xpow[i];}

	long double normtmp  = Int_xpow[0];
	long double sgn_tmp  = ((normtmp >= 0.0L )? 1.0L: -1.0L);
	long double scale    = sgn_tmp/sqrt(fabsl(normtmp));

	poly[0].c[0] *= scale;
	poly[0].a[0] *= scale;
	for(int i=0;i<n+1;i++){poly[0].O[i] *=scale;}
	poly[0].Norm = 1.0L;

//------------------------------------------------------------------------------------------
// Set up orders 1 ...n
	for (int k=1; k<n+1; k++) {
	      //++++++++++++++++++++++++++++++++++++++++
          // set the last coeff of kth order poly for stability then create the rest of the c^(0)
    	poly[k].c[k]   = poly[k-1].c[k-1];
		for (int j=0; j<k; j++) {
   	    	poly[k].c[j] = -poly[j].O[k]*poly[k].c[k]/poly[j].Norm;
		}//endfor
	      //++++++++++++++++++++++++++++++++++++++++
      	  // compute the moments O of the kth poly over w(x) starting at j=k (O is upper triable) O=\int x^j w(x) poly_k(x)
		for (int j=k; j<n+1; j++) {
   	    	poly[k].O[j] = poly[k].c[k]*Int_xpow[j+k];
 			for (int i=0; i<k; i++) {
		    	poly[k].O[j] += poly[k].c[i]*poly[i].O[j];
		  	}//endfor
		} //endfor
	      //++++++++++++++++++++++++++++++++++++++++
	      // Compute the current norm
		long double normtmp =  poly[k].c[k]*poly[k].c[k]*Int_xpow[2*k];
		for (int j=0; j<k; j++) {
        	normtmp +=  poly[k].c[j]*(poly[k].c[j]*poly[j].Norm + 2.0L*poly[j].O[k]*poly[k].c[k]);
		} // end for j
		long double sgn_tmp = ((normtmp >= 0.0L )? 1.0L: -1.0L);
		long double scale   = sgn_tmp/sqrt(fabsl(normtmp));
		if (iopt == 1) {PRINTF("  k = %d scale =%.10Lg\n", k, scale);}
	      //++++++++++++++++++++++++++++++++++++++++
	      // Scale the O and the C's to generate unit norm for stability
 		for (int j=0; j<k+1; j++) {poly[k].c[j] *= scale;}
		for (int j=k; j<n+1; j++) {poly[k].O[j] *= scale;}
		poly[k].Norm = 1.0L;
	     //++++++++++++++++++++++++++++++++++++++++
         // Compute the a-representation from the normalized c-rep
 		poly[k].a[k] = poly[k].c[k];
		for (int j=0; j<k; j++) {
		    poly[k].a[j] = 0.0L;
		    for (int l=j; l < k; l++) {
				poly[k].a[j] += poly[l].a[j]*poly[k].c[l];
		    } // end for l
		} // end for j
	} // end for k : iteration of GS to nth order
	if (iopt == 1) {PRINTF("\n"); }

	// local pointer to the nth polynomial whose roots we want
	long double *a = poly[n].a;
	long double *c = poly[n].c;

//==========================================================================
//	Construct the companion matrix, A, using the polynomial coefficients a with ann= 1.0
//         nth order poly(x) = det(A-Bx)  where B=I

	long double atmp = a[n];
	for (int i=0; i <n+1; i++) {a[i] /= atmp;}
	a[n] = 1.0L;

    double* A = new double [n*n];
    for (int i=0; i <n*n; i++) {A[i] 	       = 0.0d;}
    for (int i=0; i < n ; i++) {A[i*n + n - 1] = -a[i];}
    for (int i=1; i < n ; i++) {A[i*n + i - 1] = 1.0d;}

    double* B  = new double [n*n]; 
    for (int i=0; i<n*n; i++) {B[i] = 0.0d;}
    for (int i=0; i<n; i++) {
		int k = i*n + i;
		B[k] = 1.0d;
    }// end for i

//==========================================================================
//    Find the zeros of the n-th orthogonal polynomial by diagonalizing its companion matrix	

	int N 	 	 	 = n;
	int LDA 	 	 = n; 
	int LDB 	     = n;
    double* ALPHAR 	 = new double [n]; 
	double* ALPHAI	 = new double [n]; 
	double* BETA	 = new double [n];
    double* VL	     = new double [n]; 
	int LDVL	     = 1; 
	double* VR	     = new double [n]; 
	int LDVR	     = 1;
	int LWORK	     = -1; 
    double* WORK     = NULL;	
	int INFO	     = 1;
	char NL[5]; 	 strcpy(NL,"N");
	char NR[5];      strcpy(NR,"N");

	double tmp = 1.0d;
 	DGGEV(NL,NR, &N, A, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       		VL, &LDVL, VR, &LDVR, &tmp, &LWORK, &INFO);
	
	LWORK = (int)tmp;
	WORK  = new double [LWORK];
 	DGGEV(NL, NR, &N, A, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       		VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);

	if (INFO != 0) {
		PRINTF("  @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("  dggev failed for unknown reason!\n");
		PRINTF("  @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

	int iii = 0;
    for (int i=0; i<n; i++) {
		if (ALPHAI[i]*ALPHAI[i] > 1.0e-16) {
		   iii++;
           PRINTF("  Bad root: (%g,%g)\n",ALPHAR[i],ALPHAI[i]);
	    }// end if
    }// end for i

	if (iii != 0) {
		PRINTF("  @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("  The roots of the polynomial must be real!\n  Are the 'a' coefficients ill-posed?\n");
		long double cmax = fabsl(c[0]); long double amax = fabsl(a[0]);
		long double cmin = fabsl(c[0]); long double amin = fabsl(a[0]);
		for (int i=1; i<n+1; i++) {
			cmax = MAX(cmax,fabsl(c[i]));
			amax = MAX(amax,fabsl(a[i]));
			if (fabsl(c[i]) > 0.0L) { cmin = MIN(cmin,fabsl(c[i]));}
			if (fabsl(a[i]) > 0.0L) { amin = MIN(amin,fabsl(a[i]));}
		} // end for
		PRINTF("  a range = [%.10Lg, %.10Lg], c range = [%.10Lg, %.10Lg],\n",amin, amax, cmin, cmax);
		PRINTF("  @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

//==========================================================================
//  Copy out, sort the nodes and make sure they are good zeros of the polynomial

	for (int i=0; i<n; i++) {node[i] = ((long double) ALPHAR[i]);}
	sort(node, node + n);

	testnodes(n, node, a, iopt, ierr_zero);

//==========================================================================
//  Evaluate all the polynomials at the nodes for testing later

	for (int j=0; j<n; j++) {poly[0].p_at_nodes[j] = poly[0].a[0];}
	for (int i=1; i<n+1; i++) {
	    long double * a_tmp = poly[i].a;
	    for (int j=0; j<n; j++) {
		long double p_tmp;
		horner(i+1, a_tmp, node[j], &p_tmp);
		poly[i].p_at_nodes[j] = p_tmp;
	    } // end for j
	} // end for i

//==========================================================================
//  Construct the weights using the zeros, derivatinves of the nth poly and the n-1 poly

//  n-1st order poly
	long double *pn1       = poly[n-1].a;
	long double an1	      = poly[n-1].a[n-1];
	long double Normn1     = poly[n-1].Norm; 

//  derivative of nth order poly in the a-representation
	long double *pnp = new long double [n];
	for (int i=0; i<n; i++) {
		long double tmp = ((long double)i)+1.0L;
		pnp[i] = tmp*a[i+1];
	} // end for

//  the wghts are set using the above functions evaluated at the nodes
	for (int i=0; i<n; i++) {
   	    long double pn1v = 0.0L;  horner(n, pn1, node[i], &pn1v);
		long double pnpv = 0.0L;  horner(n, pnp, node[i], &pnpv);
		wght[i] = Normn1/(pn1v*pnpv*an1);
	}// end for
	delete [] pnp;

//==========================================================================
//  Test the weights and nodes  using the overlap of the ortho polynomials

	if (iopt == 1) {
		PRINTF("  =================================================\n");
		PRINTF("  Generalized Gaussian quadrature weights and nodes:\n");
		PRINTF("  -------------------------------------------------\n");
        for (int i=0; i<n; i++) {
			PRINTF("  node[%d] %.22Lg   wght[%d] %.22Lg\n",i, node[i], i, wght[i]);
		}// end for i
		PRINTF("  =================================================\n\n");

		PRINTF("  =================================================\n");
		PRINTF("  Test integrals of 2n-1 moments using the quadrature:\n");
		PRINTF("-  ------------------------------------------------\n");
		testgrid(n, wght, poly, iopt, ierr_ortho);
		PRINTF("  =================================================\n\n");
	} else {
		testgrid(n, wght, poly, iopt, ierr_ortho);
	}// end if
	
	if (ierr_zero[0] > 0 || ierr_ortho[0] > 0) {
		PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		PRINTF("  You have %d nodes out of error range!\n", ierr_zero[0]);
		PRINTF("  You have %d ortho pairs out of error range!\n", ierr_ortho[0]);
		PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		FFLUSH(stdout);
	} // end if
//==========================================================================
//  Clean up the memory
	
// Destructor for poly
	for (int i=0; i<n+1; i++) {
		delete [] poly[i].c;
		delete [] poly[i].a;
		delete [] poly[i].O;
		delete [] poly[i].p_at_nodes;
	}// end for i
	delete [] poly;

// lapack function DDGEV memory
	delete [] A;
	delete [] B; 
    delete [] ALPHAR; 
	delete [] ALPHAI; 
	delete [] BETA;
    delete [] VL; 
	delete [] VR; 
	delete [] WORK;

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Valuating the p_n(x_i) using horner's method.
//==========================================================================
void horner(int n, long double *a, long double x, long double *value){
	long double tmp = a[n-1];
	for (int k=n-2; k>=0; k--){
		tmp = tmp*x + a[k];
	} // end for
	value[0] = tmp;
}// end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Test the nodes
//==========================================================================
void testnodes(int n, long double * x, long double * a, int iopt, int * ierrout) {
	int ierr = 0;
	long double err = 0.0L;
	for (int i=0; i<n; i++) {
		long double zero = 0.0L;
		horner(n+1, a, x[i], &zero);
		if (fabsl(x[i]) > 1.0e-2) { zero /= fabsl(x[i]);}
		err = MAX(err,fabsl(zero));
		if (fabsl(zero) > 1.0e-5) {
		  ierr++;
		  if(iopt==1){PRINTF("  The %dth zero, %.10Lg, has eror %.10Lg\n", i,x[i],err);}
		}//endif
	}//end for
	if (ierr > 0) {
		if (iopt == 1) {PRINTF("\n");}
		PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	  	PRINTF("    Max err in zeros %.10Lg > tolerance\n", err);
		PRINTF("    for %d cases. Try root refinement\n", ierr);
		PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
		FFLUSH(stdout);
		if (err > 1.0L) {EXIT(1);}
	}else{
	  	PRINTF("  The maximum err in your zeros is %Lg\n\n", err);
	} // end if
	ierrout[0] = ierr;
//------------------------------------------------------------------
  } // end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Test the correctness
//==========================================================================
void testgrid(int n, long double * w, POLY * poly, int iopt, int * ierrout) {
//==========================================================================
//  Test the orthogonality of all the polynomials integarate using the weights and nodes
//---------------------------------------------------------------------------

	int ierr = 0;
	long double errmax = 0.0L;
	for (int k=0; k<n+1; k++) {
		int iup = ((k != n)? k+1: n); 
		long double * p_k = poly[k].p_at_nodes;
		for (int kp=0; kp<iup; kp++) {
			long double * p_kp = poly[kp].p_at_nodes;
			long double result = 0.0L;
			for (int i=0; i<n; i++) {
				result += w[i]*p_k[i]*p_kp[i];
			} //end for i
			long double ans = ((k == kp )? 1.0L : 0.0L);
			long double err = fabsl(result - ans);
			errmax = MAX(err,errmax);
			if (iopt == 1) {
				PRINTF("  Ortho test %d, %d) = %.10Lg with error %.10Lg\n", k, kp, result, err);
			} // end if
			if (fabsl(err) > 1.0e-5 && iopt == 0) {
				ierr++;
				PRINTF("  Error out of range %.10Lg for ortho (%d, %d)\n", err, k, kp);
			} // end if
		}// end for kp
	} //end for k

        if (iopt == 1) {PRINTF("\n");}
	PRINTF("  The maximum err in your orthogonality test is %Lg\n", errmax);
        if (iopt == 0) {PRINTF("\n");}

	if (ierr > 0) {
		PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		PRINTF("  You have %d ortho pairs out of error range!\n", ierr);
		PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
		FFLUSH(stdout);
//		EXIT(1);
	} // end if
	ierrout[0] = ierr;
//------------------------------------------------------------------
  } //end routine
//==========================================================================
