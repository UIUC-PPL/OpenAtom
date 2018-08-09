//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*       Projector Augmented Wave: compute core exchange correlaton contribution       */
//============================================================================
/* Standard Include files and class defs for */
#include "standard_include.h"       // has charmm++.h and pup.h in addition to usual suspects
//#include "allclass_cp.h"
#include "class_PAW_info.h"         // c++ pawinfo readonly class declarations (including pup)
#include "class_cpPAW.h"

//============================================================================
/*  compute PAW core exchange correlation
  	\sum_J \sum_f (Exc(rho_J^total(rf)) - Exc(rho_J^smooth(rf))) 
		This routine get passed rho_J^total(rf) and rho_J^smooth(rf) for a subset of atoms
		from myAtmStart to myAtmEnd
*/
//============================================================================
/*
	This routine will be invoked by PAWrhoChare or whatever Ragu calls it, under LDA
	it will send us the core density for the atoms assigned. Under LSDA, it will send
	us the up and down densities.
*/
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPPAW::PAW_exc_calc(int myAtms, int myAtmStart, int myAtmEnd, PAWCORERHO *PAWcoreRho,
													double *ExcCore_ret)
  //============================================================================
{ /* Begin Function */

  //----------------------------------------------------------------------------
  // Local Variables and local pointers
//  CP           *cp           = CP::get();
  PAWINFO      *pawinfo      = PAWINFO::get();
	int *pawAtmtyp = pawinfo->atmtyp;
	PAWFGRID *pawFgrid = pawinfo->fgrid;
	double ExcCore = 0.0;
  //============================================================================

	for (int j=0; j<myAtms; j++) {
		int jtyp = pawAtmtyp[j+myAtmStart];
		int nf = PAWcoreRho[j].nf;
		double *rho = PAWcoreRho[j].rhotot;
		double *rhoS = PAWcoreRho[j].rhoS;
		double *vxctot = PAWcoreRho[j].vxctot;
		double *vxcS = PAWcoreRho[j].vxcS;
		double *wf = pawFgrid[jtyp].wf;
		if (cp_lda == 1) {
// CPXCFNCTS::PAWexcLDA (nf, wf, rho, rhoS, vxctot, vxcS, &ExcCore);
		} else {
			double *rho_up = PAWcoreRho[j].rhotot_up;
			double *rhoS_up = PAWcoreRho[j].rhoS_up;
			double *vxctot_up = PAWcoreRho[j].vxctot_up;
			double *vxcS_up = PAWcoreRho[j].vxcS_up;
			double *rho_dn = PAWcoreRho[j].rhotot_dn;
			double *rhoS_dn = PAWcoreRho[j].rhoS_dn;
			double *vxctot_dn = PAWcoreRho[j].vxctot_dn;
			double *vxcS_dn = PAWcoreRho[j].vxcS_dn;
// CPXCFNCTS::PAWexcLSDA (nf, wf, rho_up, rhoS_up, vxctot_up, vxcS_up, rho_dn, rhoS_dn, vxctot_dn, vxcS_dn, &ExcCore);
		} // endif
	} // end for j

	ExcCore_ret[0] += ExcCore;

  //==========================================================================
}//end routine
//==========================================================================
