//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
// This program generates the f_grid
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "standard_include.h"

//==========================================================================
// Function prototypes

#include "grid.h"
#include "myGaussian.h"

//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Grid generation : Controller
//==========================================================================
void gen_fgrid (int rorder, int thetaorder, int phiorder, int lmax, double expalpha, FGRID *fgrid) {
  //==========================================================================
  // Local variables

	int mmax = lmax;
    double *wr, *xr, *wcostheta, *wphi;  // grid weights and nodes

  //==========================================================================
  //==========================================================================
  // Tell everyone what you are doing 
	PRINT_LINE_STAR;
	PRINTF("Generating the spherical polar weights and nodes, {nr, ntheta, nphi} = {%d, %d, %d}\n",rorder, thetaorder, phiorder);
	PRINT_LINE_DASH;

	wr = new double [rorder];
	xr = new double [rorder];
	wcostheta = new double [thetaorder];

	if (thetaorder < lmax){ 
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("Error! thetaorder too small for lmax!\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		EXIT(1);
  	} //end if
	if (phiorder < 2*mmax) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("Error! phiorder too small for mmax! mmax = %d, phiorder = %d\n",mmax, phiorder);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		EXIT(1);
	} //end if

  //=========================================================================
  // Generate r, theta, phi grids using DVR: Hermite, Legendre, Fourier (equal space)

	int kind = 1; double aaa = -1; double bbb = 1;
//	int kindr = 6; double aar = 0.0;
	int type = 2; int iopt = 0; 
	control_quad_rule(kind, thetaorder, aaa, bbb, wcostheta, fgrid->xcostheta); // Legendre
//	control_quad_rule(kindr, rorderfull, aar, expalpha*expalpha, wrfull, xrfull); // Hermite
	gen_Gauss_quad(type, rorder, expalpha, wr, xr, 0);
	wphi = new double [phiorder];
	genphigrid(phiorder,wphi,fgrid->xphi); // Fourier (equal space)
	
  //=========================================================================
  // Generate the spherical polar grid, Jacobian r^2 absorbed into weight
    double *wf = fgrid->wf;
    double *xf = fgrid->xf;
    double *yf = fgrid->yf;
    double *zf = fgrid->zf;
    double *rf = fgrid->rf;
	if (fgrid->nf != rorder*thetaorder*phiorder) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Bad f grid size! %d != %d\n",fgrid->nf, rorder*thetaorder*phiorder);
        PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
	} // end if

	int f = 0;
	for (int ir=0; ir<rorder; ir++) {
		for (int itheta=0; itheta<thetaorder; itheta++) {
			for (int iphi=0; iphi<phiorder; iphi++) {
				double xsintheta = sqrt(1.0-(fgrid->xcostheta[itheta])*(fgrid->xcostheta[itheta]));
		    	wf[f] = wr[ir]*xr[ir]*xr[ir]*wcostheta[itheta]*wphi[iphi];
				xf[f] = xr[ir]*xsintheta*cos(fgrid->xphi[iphi]);
				yf[f] = xr[ir]*xsintheta*sin(fgrid->xphi[iphi]);
				zf[f] = xr[ir]*fgrid->xcostheta[itheta];
			    rf[f] = xr[ir];
				f++;
			} // end for iphi
		} // end for itheta
	} // end for ir

  //==========================================================================
  // Clean up the memory
	delete [] wcostheta;
  	delete [] wr;
  	delete [] xr;
	delete [] wphi;

  //==========================================================================
  // Tell everyone you are done
	PRINT_LINE_DASH;
	PRINTF("Completed the spherical polar weights and nodes\n");
	PRINT_LINE_STAR;
	PRINTF("\n");

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================
