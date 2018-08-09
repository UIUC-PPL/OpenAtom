//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                     Module: fastAtoms.C
//
//
// This subprogram handles the atoms 
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "standard_include.h"
#include "cell.h"
#include "atom_maps.h"
#include "fastAtoms.h"
#include "mathlib.h"
#include "dictionary.h"
#include "proto_handle_entry.h"

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void FASTATOMS::readatoms(CELL *cell, ATOM_MAPS *atom_maps, NAME fnameIn, int NatmChunk)
  //===================================================================================
{//begin routine
  //===================================================================================
  int iperd     = cell->iperd;
  double animg  = cell->animg;
  double Ecut   = cell->Ecut;
  double Rcut   = cell->Rcut;
  double beta_unitless = cell->beta_unitless;

  int natm_typ;          // number of atom types
  int *index_atm_typ;    // index of atom type of each atom
  int *natm_atm_typ;     // the number of atoms of each type
  int *lmax_by_typ;      // lmax by type
	double *alp_by_typ;    // alp by type
	double *beta_by_typ;   // beta by type
	double *Rpc_by_typ;    // Rpc by type
  int **list_atm_by_typ; // list of atoms sorted by atom type
  NAME *atm_typ;         // names of the atom types
  double alpb;           // Ewald alpha               - in "cell" but recomputed below for clarity
  double gcut;           // state g space cutoff      - in "cell" but recomputed below for clarity
//  double Gcut;         // rho g space cutoff      - in "cell" but recomputed below for clarity

  double hmat[10];       // the simulation box
  double hmati[10];      // inverse simulation box
  double volume;         // simulation box volume

//==========================================================================
// Tell the world what you are doing 

  PRINTF("  ");PRINT_LINE_STAR;
  PRINTF("  Reading atom parameters and positions from %s\n",fnameIn);
  PRINTF("  ");PRINT_LINE_DASH;

 //==========================================================================
 // Read in the atoms and parameters
  FILE *fp = fopen(fnameIn,"r");
    fscanf(fp,"%d %d", &natm, &natm_typ); readtoendofline(fp);
    allocate();
    atm_typ      = new NAME [natm_typ];
    natm_atm_typ = new int [natm_typ];
    lmax_by_typ  = new int [natm_typ];
    alp_by_typ   = new double [natm_typ];
    beta_by_typ  = new double [natm_typ];
    Rpc_by_typ   = new double [natm_typ];
    int natm_atm_typ_max = 0;
    int lmax = 0;
    for (int i=0; i<natm_typ; i++) {
      int j;
      fscanf(fp,"%d %d %s %lf %d", &j, &natm_atm_typ[i], atm_typ[i], &Rpc_by_typ[i], &lmax_by_typ[i]); readtoendofline(fp);
      if (j != i) {
        PRINTF("    @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("        atom %d type out of order!\n", i);
        PRINTF("    @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      } //end if 
      lmax = MAX(lmax, lmax_by_typ[i]);
      natm_atm_typ_max = MAX(natm_atm_typ_max,natm_atm_typ[i]);
      PRINTF("     Number of atoms of type %s: %d\n",atm_typ[i], natm_atm_typ[i]);
      alp_by_typ[i]  = 1.8/Rpc_by_typ[i];
      beta_by_typ[i] = alp_by_typ[i]*beta_unitless;
    } // end for i
    int lmax_rho = 2*lmax;
    PRINTF("\n");

    index_atm_typ = new int [natm];
    list_atm_by_typ = new int *[natm_typ];
    for (int i=0; i<natm_typ; i++) { list_atm_by_typ[i] = new int [natm_atm_typ_max]; }
    double Rpc_max = 0.0;
    for (int i=0; i<natm; i++) {
      fscanf(fp,"%lf %lf %lf %lf %lf %lf %d",&x[i],&y[i],&z[i], &q[i], &qt[i], &Rpc[i], &index_atm_typ[i]); readtoendofline(fp);
			if (Rpc[i] != Rpc_by_typ[index_atm_typ[i]]) {
        PRINTF("    @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("        atom %d PAW radius does not match its type! %g != %g\n", i, Rpc[i], Rpc_by_typ[index_atm_typ[i]]);
        PRINTF("    @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      } //end if 
      alp[i] = 1.8/Rpc[i];
      beta[i] = alp[i]*beta_unitless;
      Rpc_max = MAX(Rpc_max, Rpc[i]);
    } //end for i
    fscanf(fp, "%lf %lf %lf", &hmat[1], &hmat[4], &hmat[7]); readtoendofline(fp);
    fscanf(fp, "%lf %lf %lf", &hmat[2], &hmat[5], &hmat[8]); readtoendofline(fp);
    fscanf(fp, "%lf %lf %lf", &hmat[3], &hmat[6], &hmat[9]); readtoendofline(fp);
  fclose(fp);

//==========================================================================
// Consistency checks  (some variables recalculated for clarity of reading)
  
  gcut = sqrt(Ecut);          // state g-space cutoff: gcut = sqrt(2*me*Ecut/hbar^2) = sqrt(2*Ecut) in atomic units
                              //for Ecut in Ryd, gcut = sqrt(Ecut)
  // Gcut = 2.0*gcut          // density g-space cutoff
  double gammasq = gcut*Rcut; // Gcut = 2*gcut, so Gcut^2/4 = gcut^2 = Ecut (in Ryd) = gamma^4/Rcut^2;
                              //gamma^2 = gcut*Rcut
  double gamma_conv = sqrt(gammasq);           // PAW convergence criteria
  alpb = gamma_conv/Rcut;                      // Ewald convergence parameter
  double gcuthat = gcut*hmat[1]/(2.0*M_PI_QI); // this is the state g space
  double Gcuthat = 2.0*gcuthat;                // this is the density g space
  PRINTF("     Pw cutoff E_cut           : %g Ry\n",   Ecut);
  PRINTF("     Real space cutoff R_c     : %g bohr\n", Rcut);
  PRINTF("     alpha_ewald*R_c = gamma   : %g\n", gamma_conv);
  PRINTF("     alpha_ewald*hmat[1]/2.0   : %g\n", alpb*hmat[1]/2.0);
  PRINTF("     Ewald \\hat{g}_cut space   : %g\n", Gcuthat);
  if(NatmChunk>natm || NatmChunk<1){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("     The number of atm chunks must be 1=< %d <= %d",NatmChunk,natm);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif
  if (Gcuthat < alpb*hmat[1]) {
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("     The g space is too small to converge the Ewald sum! gamma = %10g\n", gamma_conv);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  } // end if
  if (gamma_conv < 3.0) {
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("     The gamma(.%10g) = alpha_bar(.%10g)*Rcut(.%10g) is too small!\n", gamma_conv, alpb, Rcut);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  } // end if
  if (Rcut > (animg +1.0)*hmat[1]*0.5) {
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("     The box (%.10g) is too small for your Rcut (%.10g) !\n", hmat[1], Rcut);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  } // end if
  if (Rcut + Rpc_max > (animg + 1.0)*hmat[1]*0.5) {
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("     The box (%.10g) is too small for your Rcut (%.10g) to treat eN without images Rpc_max = %.10g!\n",
	         hmat[1], Rcut, Rpc_max);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
  } // end if
  if (Rcut + 2.0*Rpc_max > (animg+1.0)*hmat[1]*0.5) {
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("     The box (%.10g) is too small for your Rcut (%.10g) to treat Hartree");
    PRINTF("     without images 2Rpc_max = %.10g!\n", hmat[1], Rcut, 2.0*Rpc_max);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
  } // end if

  for (int i=0; i<natm; i++) {
 	fx[i]   = 0; fy[i]   = 0; fz[i]   = 0;
  	fxa[i]  = 0; fya[i]  = 0; fza[i]  = 0;
  } // end for

  //========================================================================
  // compute and store the atom maps

  int *natm_typ_now = new int[natm_typ];
  for (int i=0; i<natm_typ; i++) { natm_typ_now[i] = 0; }
  for (int i=0; i<natm; i++) {
    int ityp = index_atm_typ[i];
    list_atm_by_typ[ityp][natm_typ_now[ityp]] = i;
    natm_typ_now[ityp]++;
  } // end for
  delete [] natm_typ_now;
  atom_maps->natm_typ         = natm_typ;
  atom_maps->natm             = natm;
  atom_maps->lmax             = lmax;
  atom_maps->lmax_rho         = lmax_rho;
  atom_maps->index_atm_typ    = index_atm_typ;
	atom_maps->lmax_by_typ      = lmax_by_typ;
  atom_maps->natm_atm_typ     = natm_atm_typ;
  atom_maps->list_atm_by_typ  = list_atm_by_typ;
  atom_maps->atm_typ          = atm_typ;
  atom_maps->natm_atm_typ_max = natm_atm_typ_max;
	atom_maps->alp_by_typ       = alp_by_typ;
	atom_maps->beta_by_typ      = beta_by_typ;
	atom_maps->Rpc_by_typ       = Rpc_by_typ;

  gethinv(hmat, hmati, &volume, iperd);

  //========================================================================
  // store the new simulation cell information 

  cell->Rpc_max = Rpc_max;
  for (int i=1; i<10; i++) { cell->hmat[i] = hmat[i]; }
  gethinv(hmat, hmati, &volume, iperd);
  for (int i=1; i<10; i++) { cell->hmati[i] = hmati[i]; }  
  cell->volume = volume;

  //========================================================================
  // Tell the world you are done reading the atoms
   PRINTF("  ");PRINT_LINE_DASH;
   PRINTF("  Completed reading atom parameters and positions\n");
   PRINTF("  ");PRINT_LINE_STAR;
   PRINTF("\n");
  //-----------------------------------------------------------------------------
}//end routine
//========================================================================

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void FASTATOMS::output(int myPe, CELL *cell, ATOM_MAPS *atom_maps)
  //===================================================================================
{//begin routine
  //===================================================================================
	// Local variables
	NAME outputFile;
	int natm_typ = atom_maps->natm_typ;
	double *hmat = cell->hmat;
	NAME * atm_typ = atom_maps->atm_typ;
	int * natm_atm_typ = atom_maps->natm_atm_typ;
	int * lmax_by_typ = atom_maps->lmax_by_typ;
	int * index_atm_typ = atom_maps->index_atm_typ;
	double * Rpc_by_typ = atom_maps->Rpc_by_typ;
	double Rcut = cell->Rcut;
	double Ecut = cell->Ecut;

  //===================================================================================
	// output the atoms
	sprintf(outputFile,"%d_atoms.out", myPe);
	FILE * fp = fopen(outputFile,"w");
    fprintf(fp,"%d %d\n", natm, natm_typ);
    for (int i=0; i<natm_typ; i++) {
      fprintf(fp,"%d %d %s %lf %d\n", i, natm_atm_typ[i], atm_typ[i], Rpc_by_typ[i], lmax_by_typ[i]);
    } // end for i

    for (int i=0; i<natm; i++) {
      fprintf(fp,"%lf %lf %lf %lf %lf %lf %d\n",x[i],y[i],z[i],q[i],qt[i],Rpc[i],index_atm_typ[i]);
    } //end for i
    fprintf(fp, "%lf %lf %lf\n", hmat[1], hmat[4], hmat[7]);
    fprintf(fp, "%lf %lf %lf\n", hmat[2], hmat[5], hmat[8]);
    fprintf(fp, "%lf %lf %lf\n", hmat[3], hmat[6], hmat[9]);
	fclose(fp);
  //===================================================================================
}//end routine
//========================================================================
