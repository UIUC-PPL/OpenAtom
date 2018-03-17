//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
// This program generates the Ylm on the f_grid
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "standard_include.h"
#include "ckcomplex.h"
#include "fgrid.h"

//==========================================================================
// Function prototypes

#include "grid.h"

//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Program generates the Ylm on the f-grid: Controller
//==========================================================================
void gen_Ylmf (int rorder, int thetaorder, double *xcostheta, int phiorder, double *xphi, int l, int m, complex *Ylm) {
  //==========================================================================
    complex **Ylmlocal = new complex*[thetaorder];
    for (int i=0; i<thetaorder; i++) Ylmlocal[i] = new complex [phiorder];

    SphericalHarmonicOnGrid(l, m, thetaorder, xcostheta, phiorder, xphi, Ylmlocal);

	int f = 0;

	for (int ir=0; ir<rorder; ir++) {
		for (int itheta=0; itheta<thetaorder; itheta++) {
			for (int iphi=0; iphi<phiorder; iphi++) {
				Ylm[f] = Ylmlocal[itheta][iphi];
				f++;
			}
		}
	}

	delete [] Ylmlocal;
}//end routine
//==========================================================================
