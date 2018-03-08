//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Compute the Gaussian quadature weights and nodes of order n for a general 
//  positive semi-definite weighting function, given the first 2n moments of 
//  the weighting function Int_xpow (k) = \int_a^b dx w(x) x^k   k=0 ... 2n
//  Uses the companion matrix method to compute zeros of an orthogonal poly
//  over the weighting function constructed by Gram-Schimdt
//==========================================================================
// include files

#include "standard_include.h"
#include "gen_Gauss_quad_entry.h"
#include "gen_Gauss_quad_local.h"

#define _QGGEV_ON_
//#define _QGGEV_OFF_

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Controller for the Gaussian quadature weights and nodes of order n for
//    a general weighting fnc
//==========================================================================
void gen_Gauss_quad(int n, __float128 * Int_xpow, __float128 * Int_xpow_log, __float128 * Int_xpow_sgn,
		    int * Int_xpow_zero, __float128 * wght, __float128 * node, int iopt, int * ierr_zero,
		    int * ierr_ortho)
//==========================================================================
  {//Begin routine
//==========================================================================
// 1) Compute the polynomial coefficients "a" of the nth orthogonal polynomials over the weight
//    using Gram Schmidt. User provides table of the first 2n moments of w(x) stored in Int_xpow is sufficient
//    A data structure Ooly[i].c[k] where a[i+1] is the coefficient of the i+1st order orthogonal polynomial
//		    Opoly[0].c[0] = 1 starts the Gram Schmidt
//		    Create poly[1]...poly[n], you only create k = 1, 2, ..., n by GS
//		    Opoly[n].a[0..n] are the desired coefficients in the x^k basis.
//   In more detail, Gram-Schmidt only requires the evaluation of moments over the weighting function
//   stored as Int_xpos. GS is peformed in the c-rep -the expansion of the polynomial in lower ortho polynomials with
//   last term c_nx^n.  The a-rep poly = \sum_k a_k x^k is created from the c-rep via an analytic transformation
//
// 2) The zeros or nodes are computed using the companion matrix approach (can probably do better). 
//
// 3) Given Opoly[n], Opoly[n-1] and nodes, compute weights
//
// 4) Test
//
//==========================================================================
// Set the calculation method - a-rep (method=0) or c-rep (method=1)

    int method = 1; // use c or a-rep for overlap computation
//==========================================================================    
//  1) Construct Ortho polys
//-------------------------------------------------------------------------
// Constructor for poly
	POLY * poly = new POLY [n + 1];
	for (int i=0; i<n+1; i++) {
	    poly[i].order       = i;
	    poly[i].Norm        = 0.q;
	    poly[i].c	        = new __float128 [i+1];
	    poly[i].a	        = new __float128 [i+1];
	    poly[i].O	        = new __float128 [n+1]; //upper triangular matrix so over allocate (can jazz up to avoid)
            poly[i].p_at_nodes  = new __float128 [n];
	}// end for i

//------------------------------------------------------------
// Set up 0th order poly and normalize it and its moments
	poly[0].c[0] = 1.q;
	poly[0].a[0] = 1.q;
	for(int i=0;i<n+1;i++){poly[0].O[i] = Int_xpow[i];}

	__float128 normtmp  = Int_xpow[0];
	__float128 sgn_tmp  = ((normtmp >= 0.q )? 1.q: -1.q);
	__float128 scale    = sgn_tmp/sqrt(fabsq(normtmp));

	poly[0].c[0] *= scale;
	poly[0].a[0] *= scale;
	for(int i=0;i<n+1;i++){poly[0].O[i] *=scale;}
	poly[0].Norm = 1.q;

//------------------------------------------------------------------------------------------
// Set up orders 1 ...n
        if(iopt==1){
           PRINTF("  ===========================\n");
           PRINTF("  Reporing scaling factor:\n");
           PRINTF("  ---------------------------\n");
        }//ednfi
	for (int k=1; k<n+1; k++) {
	  //++++++++++++++++++++++++++++++++++++++++
          // set the last coeff of kth order poly for stability then create the rest of the c^(0)
    	        poly[k].c[k]   = poly[k-1].c[k-1];
		for (int j=0; j<k; j++) {
   	    	   poly[k].c[j] = -poly[j].O[k]*poly[k].c[k]/poly[j].Norm;
		}//endfor
	  //++++++++++++++++++++++++++++++++++++++++
      	  // compute the moments O of the kth poly over w(x) starting at j=k (upper triable) O=\int x^j w(x) poly_k(x)
		if(method==1){
		  for (int j=k; j<n+1; j++) {
		    poly[k].O[j]      = poly[k].c[k]*Int_xpow[j+k];
		    for (int i=0; i<k; i++) {
		      poly[k].O[j] += poly[k].c[i]*poly[i].O[j];
		    }//endfor
		  } //endfor
		}//endif
	  //++++++++++++++++++++++++++++++++++++++++
	  // Compute the current norm
		__float128 normtmp = poly[k].c[k]*poly[k].c[k]*Int_xpow[2*k];
		for (int j=0; j<k; j++) {
        	   normtmp +=  poly[k].c[j]*(poly[k].c[j]*poly[j].Norm + 2.q*poly[j].O[k]*poly[k].c[k]);
		} // end for j
		__float128 sgn_tmp = ((normtmp >= 0.q )? 1.q: -1.q);
		__float128 scale   = sgn_tmp/sqrt(fabsq(normtmp));
		if (iopt == 1) {
		  PRINTF("    k = %d scale =%.10Lg\n", k, (long double) scale);
		}//endif
	  //++++++++++++++++++++++++++++++++++++++++
	  // Scale the O and the C's to generate unit norm for stability
 		for (int j=0; j<k+1; j++) {poly[k].c[j] *= scale;}
		if(method==1){
		  for (int j=k; j<n+1; j++) {poly[k].O[j] *= scale;}
		}//enbdif
		poly[k].Norm = 1.q;
	 //++++++++++++++++++++++++++++++++++++++++
         // Compute the a-representation from the normalized c-rep
 		poly[k].a[k] = poly[k].c[k];
		for (int j=0; j<k; j++) {
		    poly[k].a[j] = 0.q;
		    for (int l=j; l < k; l++) {
		      poly[k].a[j] += poly[l].a[j]*poly[k].c[l];
		    } // end for l
		} // end for j
	  //++++++++++++++++++++++++++++++++++++++++
      	  // compute the moments O of the kth poly over w(x) starting at j=k (upper triable) O=\int x^j w(x) poly_k(x)
		if(method==0){
		  for (int j=k; j<n+1; j++) {
    		    poly[k].O[j] = 0.q;
		    for (int i=0; i<=k; i++) {
		      poly[k].O[j] += poly[k].a[i]*Int_xpow[j+i];
		    }//endfor
		  } //endfor
		}//endif
	} // end for k : iteration of GS to nth order
        if(iopt==1){
           PRINTF("  ===========================\n");
	   PRINTF("\n");
        }//ednfi

	// local pointer to the nth polynomial whose roots we want
	__float128 *a = poly[n].a;
	__float128 *c = poly[n].c;

//==========================================================================
//  2) Construct the companion matrix, A, using the polynomial coefficients a with ann= 1.0
//         nth order poly(x) = det(A-Bx)  where B=I and diagonlize to get the zeros
//----------------------------------------------------------------------------
// Companion matrix
     __float128 atmp = a[n];
     for (int i=0; i <n+1; i++) {a[i] /= atmp;}
     a[n] = 1.q;

#ifdef _QGGEV_ON_
    __float128* A = new __float128 [n*n];
    for (int i=0; i <n*n; i++) {A[i] 	       = 0.q;}
    for (int i=0; i < n ; i++) {A[i*n + n - 1] = -a[i];}
    for (int i=1; i < n ; i++) {A[i*n + i - 1] = 1.q;}

    __float128* B  = new __float128 [n*n]; 
    for (int i=0; i<n*n; i++) {B[i] = 0.q;}
    for (int i=0; i<n; i++) {
	int k = i*n + i;
	B[k] = 1.q;
    }// end for i
#else
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
#endif

//----------------------------------------------------------------------------
//    Find the zeros of the n-th orthogonal polynomial by diagonalizing its companion matrix	

	int N 	 	 = n;
	int LDA 	 = n; 
	int LDB 	 = n;
	int LDVL         = 1; 
	int LDVR         = 1;
	int LWORK        = -1; 
	int INFO	 = 27;
	char NL[2]; 	 strcpy(NL,"N");
	char NR[2];      strcpy(NR,"N");

#ifdef _QGGEV_ON_	
        __float128* ALPHAR   = new __float128 [n]; 
	__float128* ALPHAI   = new __float128 [n]; 
	__float128* BETA     = new __float128 [n];
        __float128* VL       = new __float128 [n]; 
	__float128* VR       = new __float128 [n]; 
        __float128* WORK     = NULL;	
	__float128 tmp       = 99.q; // assigned for testing passsing to fortan 

 	QGGEV(NL,NR, &N, A, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
      	      VL, &LDVL, VR, &LDVR, &tmp, &LWORK, &INFO);

	LWORK = (int)tmp;
	WORK  = new __float128 [LWORK];
 	QGGEV(NL, NR, &N, A, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       	      VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
#else
        double* ALPHAR 	 = new double [n]; 
	double* ALPHAI	 = new double [n]; 
	double* BETA	 = new double [n];
        double* VL       = new double [n]; 
	double* VR       = new double [n]; 
        double* WORK     = NULL;	
	double tmp       = 1.0d;
 	DGGEV(NL,NR, &N, A, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       		VL, &LDVL, VR, &LDVR, &tmp, &LWORK, &INFO);
	
	LWORK = (int)tmp;
	WORK  = new double [LWORK];
 	DGGEV(NL, NR, &N, A, &LDA, B, &LDB, ALPHAR, ALPHAI, BETA,
       	      VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
#endif
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
#ifdef _QGGEV_ON_		
                PRINTF("  Bad root: (%Lg,%Lg)\n",(long double)ALPHAR[i],(long double)ALPHAI[i]);
#else
                PRINTF("  Bad root: (%g,%g)\n",ALPHAR[i],ALPHAI[i]);
#endif		
	    }// end if
        }// end for i

	if (iii != 0) {
		PRINTF("  @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("  The roots of the polynomial must be real!\n  Are the 'a' coefficients ill-posed?\n");
		__float128 cmax = fabsq(c[0]); __float128 amax = fabsq(a[0]);
		__float128 cmin = fabsq(c[0]); __float128 amin = fabsq(a[0]);
		for (int i=1; i<n+1; i++) {
			cmax = MAX(cmax,fabsq(c[i]));
			amax = MAX(amax,fabsq(a[i]));
			if (fabsq(c[i]) > 0.q) { cmin = MIN(cmin,fabsq(c[i]));}
			if (fabsq(a[i]) > 0.q) { amin = MIN(amin,fabsq(a[i]));}
		} // end for
		PRINTF("  a range = [%.10Lg, %.10Lg], c range = [%.10Lg, %.10Lg],\n",
		       (long double)amin, (long double)amax, (long double)cmin, (long double)cmax);
		PRINTF("  @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
	} // end if

//==========================================================================
// 2.B Copy out, sort the nodes and make sure they are good zeros of the polynomial

#ifdef _QGGEV_ON_
	for (int i=0; i<n; i++) {node[i] = ALPHAR[i];}	
#else
	for (int i=0; i<n; i++) {node[i] = ((__float128) ALPHAR[i]);}
#endif
	sort(node, node + n);

//==========================================================================
// 2.C Evaluate all the polynomials at the nodes for testing later
//     use horner for a-rep and just sumnming for c-rep							   

        method = 0; // seems a-rep is better here.
        if(iopt==1){
           PRINTF("  ===========================\n");
           PRINTF("  Calcuation method\n");
           PRINTF("  ---------------------------\n");
        }//ednfi
 	for (int j=0; j<n; j++) {poly[0].p_at_nodes[j] = poly[0].a[0];}

	switch (method){
	  //++++++++++++++
	  // a-rep
	   case 0:
	     PRINTF("  Using a-rep ...\n");
	     for (int k=1; k<n+1; k++) { // all polynomials
	       __float128 * a_tmp = poly[k].a;
	       for (int j=0; j<n; j++) { // al nodes
		 __float128 p_tmp;
		 horner(k+1, a_tmp, node[j], &p_tmp);
		 poly[k].p_at_nodes[j] = p_tmp;
	       } // end for j
	     } // end for k
	   break;
	  //++++++++++++++
	  // c-rep
 	   case 1:
     	     PRINTF("   Using c-rep ...\n");
	     for (int k=1; k<n+1; k++) {// all polynomials
	       __float128 * c_k  = poly[k].c;
	       __float128 ak     = ((__float128) k);
	       __float128 ck_sgn = ((c_k[k] >= 0.q ) ? 1.q : -1.q);
	       __float128 ck_log = log(fabsq(c_k[k]));
	       for (int j=0; j<n; j++) {// all nodes
		  __float128 node_sgn = ((node[j] >= 0.q ) ? 1.q : -1.q);
		  __float128 node_log = log(fabsq(node[j]));
		  __float128 p_tmp    = 0.q;
 	          for (int i=0; i<k; i++) { // all lower polynomials
		    p_tmp += c_k[i]*poly[i].p_at_nodes[j];
		  }//endfor
		  if(k%2==0){
		    p_tmp += ck_sgn*exp(ck_log+ak*node_log);
		  }else{
    		    p_tmp += ck_sgn*node_sgn*exp(ck_log+ak*node_log);
		  }//endif
 		  poly[k].p_at_nodes[j] = p_tmp;	  
	       }//endfor
	     }//endfor
          break;
	}//switch : c rep or a-rep

	if(iopt==0){PRINTF("\n");}
	if(iopt==1){PRINTF("  ========================\n\n");}

	testnodes(n, node, poly[n].p_at_nodes, iopt, ierr_zero);

//==========================================================================
// 3) Construct the weights using the zeros, derivatinves of the nth poly and the n-1 poly

//  n-1st order poly
	__float128 an1	     = poly[n-1].a[n-1];
	__float128 Normn1    = poly[n-1].Norm; 

//  derivative of nth order poly in the a-representation
	__float128 *pnp = new __float128 [n];
	for (int i=0; i<n; i++) {
		__float128 tmp = ((__float128)i)+1.q;
		pnp[i] = tmp*a[i+1];
	} // end for

//  the wghts are set using the above functions evaluated at the nodes
	for (int i=0; i<n; i++) {
   	    __float128 pn1v = poly[n-1].p_at_nodes[i];  
	    __float128 pnpv = 0.q;  horner(n, pnp, node[i], &pnpv);
	    wght[i] = Normn1/(pn1v*pnpv*an1);
	}// end for
	delete [] pnp;

//==========================================================================
// 4) Test the weights and nodes  using the overlap of the ortho polynomials

 	if (iopt == 1) {
		PRINTF("  =================================================\n");
		PRINTF("  Test integrals of 2n-1 moments using the quadrature:\n");
		PRINTF("  ------------------------------------------------\n");
		testgrid(n, wght, poly, iopt, ierr_ortho);
		PRINTF("  =================================================\n");
        } else {
		testgrid(n, wght, poly, iopt, ierr_ortho);
	}// end if
	
//==========================================================================
// 5) Output the results if verbose

 	if (iopt == 1) {
	        PRINTF("\n");
		PRINTF("  =================================================\n");
		PRINTF("  Generalized Gaussian quadrature weights and nodes:\n");
		PRINTF("  -------------------------------------------------\n");
                for (int i=0; i<n; i++) {
			PRINTF("     node[%d] %.22Lg   wght[%d] %.22Lg\n",i,
			       (long double)node[i], i, (long double)wght[i]);
		}// end for i
                PRINTF("  =================================================\n\n");
	}//endif

//==========================================================================
// 6) Error Summary

 	if (iopt == 1) {PRINTF("  =================================================\n");}
        PRINTF("   Error summary:\n");
 	if (iopt == 1) {PRINTF("  ------------------------------------------------\n");}
	if (ierr_zero[0] > 0 || ierr_ortho[0] > 0) {
		PRINTF("    $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		PRINTF("      You have %d nodes out of error range!\n", ierr_zero[0]);
		PRINTF("      You have %d ortho pairs out of error range!\n", ierr_ortho[0]);
		PRINTF("    $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		FFLUSH(stdout);
	}else{
	  PRINTF("     No errors detected\n");
	} // end if
 	if (iopt == 1) {PRINTF("  =================================================\n");}

//==========================================================================
// 7) Clean up the memory allowing multiple calls of the routine
	
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
void horner(int n, __float128 *a, __float128 x, __float128 *value){
	__float128 tmp = a[n-1];
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
void testnodes(int n, __float128 * x, __float128 * zero_vec, int iopt, int * ierrout) {
 	if (iopt == 1) {
   	   PRINTF("  =================================================\n");
	   PRINTF("  Test zeros obtained from companion matrix\n");
	   PRINTF("  -------------------------------------------------\n");
	}//endif
        int ierr = 0;
	__float128 err     = 0.q;
	__float128 err_avg = 0.q;
	for (int i=0; i<n; i++) {
  	    __float128 zero = zero_vec[i];
	    if (fabsq(x[i]) > 1.0e-2) { zero /= fabsq(x[i]);}
	    err = MAX(err,fabsq(zero));
	    err_avg += fabsq(zero);
	    if (fabsq(zero) > 1.0e-5) {
	      ierr++;
	      if(iopt==1){PRINTF("  The %dth zero, %.10Lg, has eror %.10Lg\n", i,(long double)x[i],(long double)err);}
	    }//endif
	}//end for
	err_avg /= ((__float128) n);

	if(ierr>0 && iopt==1){PRINTF("\n");}
        PRINTF("    The max and average err in your zeros are (%Lg, %Lg)\n", (long double)err,(long double)err_avg);
	if (ierr > 0) {
		if (iopt == 1) {PRINTF("\n");}
		PRINTF("     $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	  	PRINTF("       Max err in zeros %.10Lg > tolerance\n", (long double)err);
		PRINTF("       for %d cases. Try root refinement\n", ierr);
		PRINTF("     $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		FFLUSH(stdout);
		if (err > 1.q) {EXIT(1);}
	}//endif

 	if (iopt == 1) {PRINTF("  =================================================\n");}
	PRINTF("\n");

	ierrout[0] = ierr;
//------------------------------------------------------------------
  } // end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Test the correctness
//==========================================================================
void testgrid(int n, __float128 * w, POLY * poly, int iopt, int * ierrout) {
//==========================================================================
//  Test the orthogonality of all the polynomials integarate using the weights and nodes
//---------------------------------------------------------------------------

	int ierr = 0;
	__float128 errmax = 0.q;
	__float128 err_avg = 0.q;
	__float128 count   = 0.q;
	for (int k=0; k<n+1; k++) {
		int iup = ((k != n)? k+1: n); 
		__float128 * p_k = poly[k].p_at_nodes;
		for (int kp=0; kp<iup; kp++) {
			__float128 * p_kp = poly[kp].p_at_nodes;
			__float128 result = 0.q;
			for (int i=0; i<n; i++) {
				result += w[i]*p_k[i]*p_kp[i];
			} //end for i
			__float128 ans = ((k == kp )? 1.q : 0.q);
			__float128 err = fabsq(result - ans);
			errmax = MAX(err,errmax);
			err_avg += err;
			count += 1.0L;
			if (iopt == 1) {
				PRINTF("      Ortho test %d, %d) = %.10Lg with error %.10Lg\n", k, kp,
				       (long double)result, (long double)err);
			} // end if
			if (fabsq(err) > 1.0e-5 && iopt == 0) {
				ierr++;
				if(ierr==1){PRINTF("   Ortho test tolerence errors\n");}
				PRINTF("      Ortho pair (%d, %d) error out of range %.10Lg \n", k, kp,(long double)err);
			} // end if
		}// end for kp
	} //end for k
	err_avg /= count;

        if (iopt == 1 || ierr>0) {PRINTF("\n");}

	PRINTF("    The max and avg err in your orthogonality test are (%Lg , %Lg)\n",
	       (long double)errmax,(long double)err_avg);

	if (ierr > 0) {
		PRINTF("    $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		PRINTF("       You have %d ortho pairs out of error range!\n", ierr);
		PRINTF("    $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
		FFLUSH(stdout);
//		EXIT(1);
	} // end if
        if (iopt == 0) {PRINTF("\n");}

	ierrout[0] = ierr;
//------------------------------------------------------------------
  } //end routine
//==========================================================================
