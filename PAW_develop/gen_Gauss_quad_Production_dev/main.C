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

#include "gen_Gauss_quad_driver_entry.h"

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
	int ierr_zero; // bad nodes
	int ierr_ortho; // bad orthogonality
	
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
		
//==========================================================================
//	Fetch the integrals of powers of x over the desired weight 
//  Hard coding a for now

	double * w = new double [n]; // Gauss quad weight
	double * x = new double [n]; // Gauss quad node

	gen_Gauss_quad_driver (type, n, iopt, w, x, &ierr_zero, &ierr_ortho);

	double alpha = 1.0;

    for (int i=0; i<n; i++) {
        w[i] /= alpha;
        x[i] /= alpha;
    }// end for i

	PRINTF("\n\n");
	PRINTF("gen_Gauss_quad_driver completed with errors (%d, %d)\n\n", ierr_zero, ierr_ortho);

	return 1;
} // end routine
//==========================================================================

