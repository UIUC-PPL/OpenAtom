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
#include "ckcomplex.h"
#include "fgrid.h"
#include "compchargePAW.h"
#include "grid.h"
#include "gen_Gauss_quad_driver_entry.h"
int main (int, char *[]);

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller to manage compensation charge energy
//==========================================================================
int main (int argc, char *argv[]){
	//==========================================================================
	// Local variables
	ATOM_MAPS atom_maps;
	ATOM_POS  atom_pos;
	ATOM_POS  atom_pos_dummy;
	CELL      cell;

	int iperd;
	int rorder;
	int thetaorder;
	int phiorder;    // grid sizes
	int lmax;									// maximum angular momentum
	int nimg;			// number of images
	double animg;
	double beta_unitless;       // unitless beta for screening, beta[J] = alpha[J]*beta_unitless

	//==========================================================================
	ESTRUCT energy;		// compensation charge energy terms stored nicely
#ifdef _FORCECHECK_
	double delta = 1.0e-5;
	ESTRUCT *energy_plus  = new ESTRUCT [3];
	ESTRUCT *energy_minus = new ESTRUCT [3];
#endif

	fillEstruct(&energy);

	char fnameIn[MAXLINE];

	//==========================================================================
	// Tell everyone what you are doing

	PRINTF("\n");
	PRINT_LINE_STAR
	PRINTF("Test a model PAW energy\n");
	PRINT_LINE_DASH

  //=========================================================================
  //             Check for input file

  if(argc < 9) {
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("No input file specified!\n");
    PRINTF("Run it like: ./model_PAW.x PAW.in rorder thetaorder phiorder lmax beta_unitless iperd nimg model_type\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  //==========================================================================
  // Read the user specified input file: atoms, Ewald parameters, core density parameters and box
  // store the information

	rorder 	 	    = atoi(argv[2]);
	thetaorder    = atoi(argv[3]);
	phiorder      = atoi(argv[4]);
	lmax 	 	      = atoi(argv[5]);
	beta_unitless = atof(argv[6]);
	iperd 		    = atoi(argv[7]);
  nimg		      = atoi(argv[8]);
	animg		      = (double) nimg;

	int model = 0;
	if (strcasecmp(argv[9],"bessel") == 0) {model = 1;}
	if (strcasecmp(argv[9],"gauss") == 0) {model = 2;}
	if (model == 0) {
   	PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
   	PRINTF("unknown model!\n");
   	PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
   	FFLUSH(stdout);
   	EXIT(1);
	} // end if

	//========================================================================
	// reading file input for atom positions and parameters
	strcpy(fnameIn, argv[1]);
	cell.iperd = iperd;
	cell.nimg  = nimg;
	cell.animg = animg;
	fill_atompos_maps(fnameIn, &atom_maps, &atom_pos, &atom_pos_dummy, &cell, beta_unitless);
	int natm = atom_pos.natm;
	int natm_typ = atom_maps.natm_typ;

	//========================================================================
	// compute the grids
	if (model == 1) { rorder++; }
	FGRID *fgrid = new FGRID [natm_typ];      // fgrid structure
	fill_fgrid(fgrid, &atom_maps, &atom_pos, &cell, rorder, thetaorder, phiorder, model);

  //========================================================================
  // Compute the real space energy analytically

	clock_t start, end;
	PRINTF("Computing the PAW energy analytically with a Gaussian basis set for 0D and short range 3D:\n");
	start = clock();
	computePAWreal(&atom_maps, &atom_pos, &cell, &energy);
	end = clock();
	PRINTF("step1 finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

  //========================================================================
  // Compute the long range energy analytically and on grid

	PRINTF("Computing the 3D long range PAW energy with grid and Gaussian basis analytically:\n");
	start = clock();
	computePAWlong(&atom_maps, &atom_pos, &cell, &energy, fgrid);
	end = clock();
	PRINTF("step2 finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

  //========================================================================
  // Compute the real space energy on grid

	PRINTF("Computing the 0D and 3D short range PAW energy with grid:\n");
	start = clock();
	computePAWGrid(lmax, &atom_maps, &atom_pos, &cell, &energy, fgrid);
	end = clock();
	PRINTF("step3 finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

	energy.Etot0D.E           = energy.ENN.E + energy.EeN.E + energy.EHar.E;
	energy.Etot0Dscr.E        = energy.ENN.E + energy.EeN.E + energy.EHarscr.E;
	energy.Etot0D.EGrid       = energy.ENN.EGrid + energy.EeN.EGrid + energy.EHar.EGrid;
	energy.Etot0Dscr.EGrid    = energy.ENN.EGrid + energy.EeN.EGrid + energy.EHarscr.EGrid;
	energy.Etot3D.E           = energy.ENNshort.E + energy.EeNshort.E + energy.EHarshort.E + energy.Elong.E - energy.ENNselflong.E;
	energy.Etot3D.EGrid       = energy.ENNshort.EGrid + energy.EeNshort.EGrid + energy.EHarshort.EGrid + energy.Elong.EGrid - energy.ENNselflong.EGrid;
	energy.Etot3Dscr.E        = energy.ENNshort.E + energy.EeNshort.E + energy.EHarshortscr.E + energy.Elong.E - energy.ENNselflong.E;
	energy.Etot3Dscr.EGrid    = energy.ENNshort.EGrid + energy.EeNshort.EGrid + energy.EHarshortscr.EGrid + energy.Elong.EGrid - energy.ENNselflong.EGrid;

	//========================================================================
	// Print out the energy

	PRINT_LINE_STAR;
	PRINTF(" Printing out results: analytical, f-grid, and error\n");
	PRINT_LINE_DASH;
	PRINTF("Energy                |    Analytic   |   Grid       |  Diff       |  %%Diff\n");
	PRINT_LINE_DASH;
	//-----------------------------------------
	//Printing the 0D energies
	energy.ENN.pres();
	energy.EeN.pres();
	energy.EeNself.pres();
	energy.EHar.pres();
	energy.EHarself.pres();
	energy.EHarscr.pres();
	energy.EHarselfscr.pres();
	//-----------------------------------------
	//Printing the 3D short range energies
	PRINT_LINE_DASH;
	energy.ENNshort.pres();
	energy.EeNshort.pres();
	energy.EeNshortself.pres();
	energy.EHarshort.pres();
	energy.EHarshortself.pres();
	energy.EHarshortscr.pres();
	energy.EHarshortselfscr.pres();
	//-----------------------------------------
	//Printing the 3D long range energies
	PRINT_LINE_DASH;
	energy.ENNselflong.pres();
	energy.Elong.pres();
	//-----------------------------------------
	//Printing the 0D and 3D total energies
	PRINT_LINE_DASH;
	energy.Etot0D.pres();
	energy.Etot3D.pres();
	energy.Etot0Dscr.pres();
	energy.Etot3Dscr.pres();
	PRINT_LINE_DASH;
	PRINTF(" Completed output\n");
	PRINT_LINE_STAR;
	PRINTF("\n");
	//==========================================================================
// } //end for loop over r grid sizes
  //==========================================================================
  //==========================================================================
  // Tell everyone you are done

  PRINT_LINE_DASH
  PRINTF("Finished computing the model PAW energy\n");
  PRINT_LINE_STAR
  PRINTF("\n");

  //==========================================================================
  // Tell everyone you are starting the forces

  PRINT_LINE_STAR
  PRINTF("Checking the forces of the model PAW energy\n");
  PRINT_LINE_DASH
  PRINTF("\n");

#ifdef _FORCECHECK_

	int check_grid = 1;
	double * xd = atom_pos_dummy.x;
	double * yd = atom_pos_dummy.y;
	double * zd = atom_pos_dummy.z;
	double * x = atom_pos.x;
	double * y = atom_pos.y;
	double * z = atom_pos.z;
	double * fxg = atom_pos.fxg;
	double * fyg = atom_pos.fyg;
	double * fzg = atom_pos.fzg;
	double * fx = atom_pos.fx;
	double * fy = atom_pos.fy;
	double * fz = atom_pos.fz;
	double * fx0g = atom_pos.fx0g;
	double * fy0g = atom_pos.fy0g;
	double * fz0g = atom_pos.fz0g;
	double * fx0 = atom_pos.fx0;
	double * fy0 = atom_pos.fy0;
	double * fz0 = atom_pos.fz0;
  PRINTF("Enter the atom you want me to check:\n");
  int iii;
  scanf("%d", &iii);
  double deltau;
  ESTRUCT * edummy;

  for (int i=0; i<2; i++) {
	for (int j=0; j<3; j++) {
		switch (i) {
			case 0: deltau = delta;  edummy = &energy_plus[j]; break;
			case 1: deltau = -delta; edummy = &energy_minus[j]; break;
		} // end switch i
		switch (j) {
			case 0: xd[iii] = x[iii] + deltau; break;
			case 1: yd[iii] = y[iii] + deltau; break;
			case 2: zd[iii] = z[iii] + deltau; break;
		} // end switch j
   		computePAWreal(&atom_maps, &atom_pos_dummy, &cell, edummy);
   		computePAWlong(&atom_maps, &atom_pos_dummy, &cell, edummy, fgrid);
   		computePAWGrid(lmax, &atom_maps, &atom_pos_dummy, &cell, edummy, fgrid);
  		edummy->Etot0D.E     = edummy->ENN.E + edummy->EeN.E + edummy->EHar.E;
  		edummy->Etot0D.EGrid = edummy->ENN.EGrid + edummy->EeN.EGrid + edummy->EHar.EGrid;
  		edummy->Etot3D.E     = edummy->ENNshort.E + edummy->EeNshort.E + edummy->EHarshort.E + edummy->Elong.E - edummy->ENNselflong.E;
  		edummy->Etot3D.EGrid = edummy->ENNshort.EGrid + edummy->EeNshort.EGrid + edummy->EHarshort.EGrid
								+ edummy->Elong.EGrid - edummy->ENNselflong.EGrid;
		switch (j) {
			case 0: xd[iii] = x[iii]; break;
			case 1: yd[iii] = y[iii]; break;
			case 2: zd[iii] = z[iii]; break;
		} // end switch j
	} // end for j
  }// end for i

  double * fdummy = new double [3];
	for (int j=0; j<3; j++) {
		if (check_grid == 1) {
			if (iperd == 3) {
				fdummy[j] = (energy_minus[j].Etot3D.EGrid -  energy_plus[j].Etot3D.EGrid)/(2.0*delta);
			} else {
				fdummy[j] = (energy_minus[j].Etot0D.EGrid -  energy_plus[j].Etot0D.EGrid)/(2.0*delta);
			}// end if
		} else {
			if (iperd == 3) {
				fdummy[j] = (energy_minus[j].Etot3D.E -  energy_plus[j].Etot3D.E)/(2.0*delta);
			} else {
				fdummy[j] = (energy_minus[j].Etot0D.E -  energy_plus[j].Etot0D.E)/(2.0*delta);
			}// end if
		} // end if
	} // end for j
	if (iperd == 3) {
		if (check_grid != 1) {
			PRINTF("finite diff 3D: %.10g %.10g %.10g \n", fdummy[0], fdummy[1], fdummy[2]);
			PRINTF(" analytical 3D: %.10g %.10g %.10g\n", fx[iii], fy[iii], fz[iii]);
		} else {
			PRINTF("finite diff 3D: %.10g %.10g %.10g \n", fdummy[0], fdummy[1], fdummy[2]);
			PRINTF("       grid 3D: %.10g %.10g %.10g\n", fxg[iii], fyg[iii], fzg[iii]);
		} // end if
	} else {
		if (check_grid != 1) {
			PRINTF("finite diff 0D: %.10g %.10g %.10g \n", fdummy[0], fdummy[1], fdummy[2]);
			PRINTF(" analytical 0D: %.10g %.10g %.10g\n", fx0[iii], fy0[iii], fz0[iii]);
		} else {
			PRINTF("finite diff 0D: %.10g %.10g %.10g \n", fdummy[0], fdummy[1], fdummy[2]);
			PRINTF("       grid 0D: %.10g %.10g %.10g\n", fx0g[iii], fy0g[iii], fz0g[iii]);
		} // end if
	}// end if iperd
#endif // _FORCECHECK_

  //========================================================================
  // Print out the forces

	PRINT_LINE_STAR;
	PRINTF(" Printing out forces\n");
	PRINT_LINE_DASH;
	for (int i=0; i<natm; i++) {
	  PRINTF("fx0  [%d]: % 14.12f  |  % 14.12f  |  % 14.12f\n",i,atom_pos.fx0[i],atom_pos.fy0[i],atom_pos.fz0[i]);
	  PRINTF("fx0g [%d]: % 14.12f  |  % 14.12f  |  % 14.12f\n",i,atom_pos.fx0g[i],atom_pos.fy0g[i],atom_pos.fz0g[i]);
	} // end for
	PRINT_LINE_DASH;
	for (int i=0; i<natm; i++) {
	  PRINTF("fx   [%d]: % 14.12f  |  % 14.12f  |  % 14.12f\n",i,atom_pos.fx[i],atom_pos.fy[i],atom_pos.fz[i]);
	  PRINTF("fxg  [%d]: % 14.12f  |  % 14.12f  |  % 14.12f\n",i,atom_pos.fxg[i],atom_pos.fyg[i],atom_pos.fzg[i]);
	} // end for

  //==========================================================================
  // Tell everyone you are done
  PRINT_LINE_DASH
  PRINTF("Finished computing the model PAW forces\n");
  PRINT_LINE_STAR
  PRINTF("\n");
  return 0;

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
void fillEstruct(ESTRUCT * energy) {
//==========================================================================
  strcpy(energy->ENN.name                  , "E_NN_0D                  ");
  strcpy(energy->EeN.name                  , "E_eN_0D                  ");
  strcpy(energy->EeNself.name              , "E_eNself_0D              ");
  strcpy(energy->EHar.name                 , "E_Har_0D                 ");
  strcpy(energy->EHarself.name             , "E_Har_self_0D            ");
  strcpy(energy->ENNshort.name             , "E_NNshort_3D             ");
  strcpy(energy->EeNshort.name             , "E_eN_short_3D            ");
  strcpy(energy->EeNshortself.name         , "E_eN_short_self_3D       ");
  strcpy(energy->EHarshort.name            , "E_Har_short_3D           ");
  strcpy(energy->EHarshortself.name        , "E_Har_short_self_3D      ");
  strcpy(energy->Elong.name                , "E_long_3D                ");
  strcpy(energy->ENNselflong.name          , "E_NNselflong_3D          ");
  strcpy(energy->Etot0D.name               , "E_tot_0D                 ");
  strcpy(energy->Etot3D.name               , "E_tot_3D                 ");
  strcpy(energy->EHarselfscr.name          , "E_Har_self_0D_scr        ");
  strcpy(energy->EHarshortselfscr.name     , "E_Har_short_self_3D_scr  ");
  strcpy(energy->Etot0Dscr.name            , "E_tot_0D_scr             ");
  strcpy(energy->Etot3Dscr.name            , "E_tot_3D_scr             ");
  strcpy(energy->EHarscr.name              , "E_Har_0D_scr             ");
  strcpy(energy->EHarshortscr.name         , "E_Har_short_3D_scr       ");
} // end subroutine

//==========================================================================
void fill_atompos_maps(char *fnameIn, ATOM_MAPS *atom_maps, ATOM_POS *atom_pos, ATOM_POS *atom_pos_dummy, CELL *cell, double beta_unitless){
//==========================================================================
	int iperd = cell->iperd;
	int nimg  = cell->nimg;
	double animg  = cell->animg;
	int natm_typ;       	// number of atom types
	int natm;				// number of atoms
	int *index_atm_typ;	// index of atom type of each atom
	int *natm_atm_typ;	// the number of atoms of each type
	int **list_atm_by_typ;// list of atoms sorted by atom type
	NAME *atm_typ;	    // names of the atom types
	double alpb;			// Ewald alpha
	double gcut;		    // state g space cutoff
	double Gcut;		    // density g space cutoff = 2*gcut
	double Ecut, Rcut;

	double *x,*y,*z,*q,*qt,*alp,*beta, *Rpc;    // atom positions and core Gaussian parameters
	double hmat[10];				  // the simulation box
	double hmati[10];				  // inverse simulation box
	double volume;				  // simulation box volume

//==========================================================================
	PRINTF("\nReading input parameters and atom positions from %s\n\n",fnameIn);
	FILE *fp = fopen(fnameIn,"r");
		fscanf(fp,"%d %d", &natm, &natm_typ); readtoendofline(fp);
		atm_typ = new NAME [natm_typ];
		natm_atm_typ = new int [natm_typ];
		int natm_atm_typ_max = 0;
		for (int i=0; i<natm_typ; i++) {
			int j;
			fscanf(fp,"%d %d %s", &j, &natm_atm_typ[i], atm_typ[i]); readtoendofline(fp);
			if (j != i) {
				PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
				PRINTF("atom type out of order!\n");
				PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
				FFLUSH(stdout);
				EXIT(1);
			} //end if
			natm_atm_typ_max = MAX(natm_atm_typ_max,natm_atm_typ[i]);
			PRINT_LINE_DASH;
			PRINTF("Number of atom %s: %d\n",atm_typ[i], natm_atm_typ[i]);
		} // end for i
		PRINT_LINE_STAR;
		PRINTF("\n");
		fscanf(fp,"%lf %lf",&Rcut, &Ecut); readtoendofline(fp);

		x = new double [natm]; y = new double [natm]; z = new double [natm];
		q = new double [natm]; qt = new double [natm];
		Rpc = new double [natm]; alp = new double [natm]; beta = new double [natm];

		double *vx   = new double [natm]; double *vy   = new double [natm]; double *vz   = new double [natm];
		double *fx0  = new double [natm]; double *fy0  = new double [natm]; double *fz0  = new double [natm];
		double *fx   = new double [natm]; double *fy   = new double [natm]; double *fz   = new double [natm];
		double *fx0g = new double [natm]; double *fy0g = new double [natm]; double *fz0g = new double [natm];
		double *fxg  = new double [natm]; double *fyg  = new double [natm]; double *fzg  = new double [natm];

		index_atm_typ = new int [natm];
		list_atm_by_typ = new int *[natm_typ];
		for (int i=0; i<natm_typ; i++) { list_atm_by_typ[i] = new int [natm_atm_typ_max]; }
		double Rpc_max = 0.0;
		for (int i=0; i<natm; i++) {
			fscanf(fp,"%lf %lf %lf %lf %lf %lf %d",&x[i],&y[i],&z[i], &q[i], &qt[i], &Rpc[i], &index_atm_typ[i]); readtoendofline(fp);
			alp[i] = 1.8/Rpc[i];
			beta[i] = alp[i]*beta_unitless;
			Rpc_max = MAX(Rpc_max, Rpc[i]);
		} //end for i
		fscanf(fp, "%lf %lf %lf", &hmat[1], &hmat[4], &hmat[7]); readtoendofline(fp);
		fscanf(fp, "%lf %lf %lf", &hmat[2], &hmat[5], &hmat[8]); readtoendofline(fp);
		fscanf(fp, "%lf %lf %lf", &hmat[3], &hmat[6], &hmat[9]); readtoendofline(fp);
  fclose(fp);

	gcut = sqrt(Ecut); // gcut = sqrt(2*me*Ecut/hbar^2) = sqrt(2*Ecut) in atomic units, for Ecut in Ryd, gcut = sqrt(Ecut)
	Gcut = 2.0*gcut; // the density g cutoff is twice the pw g cutoff
	double gammasq = gcut*Rcut; // Gcut = 2*gcut, so Gcut^2/4 = gcut^2 = Ecut (in Ryd) = gamma^4/Rcut^2; gamma^2 = gcut*Rcut
	double gamma_conv = sqrt(gammasq);
	alpb = gamma_conv/Rcut;
	double gcuthat = gcut*hmat[1]/(2.0*M_PI_QI); // this is the state g space
	double Gcuthat = 2.0*gcuthat; // this is the density g space
	PRINTF("Pw cutoff E_cut is %lf Ry\n", Ecut);
	PRINTF("Real space cutoff R_c %lf bohr\n", Rcut);
	PRINTF("Gamma = alpb*R_c is %lf\n", gamma_conv);
	PRINTF("Ewald g space = %lf\n",Gcuthat);
	PRINTF("Ewald alpb*hmat[1]/2.0 = %lf\n",alpb*hmat[1]/2.0);
	if (Gcuthat < alpb*hmat[1]) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The g space is too small to converge the Ewald sum! gamma = %10g\n", gamma_conv);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
	} // end if
	if (gamma_conv < 3.0) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The gamma(.%10g) = alpha_bar(.%10g)*Rcut(.%10g) is too small!\n", gamma_conv, alpb, Rcut);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
	} // end if
	if (Rcut > (animg +1.0)*hmat[1]*0.5) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The box (%.10g) is too small for your Rcut (%.10g) !\n", hmat[1], Rcut);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
	} // end if
	if (Rcut + Rpc_max > (animg + 1.0)*hmat[1]*0.5) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The box (%.10g) is too small for your Rcut (%.10g) to treat eN without images Rpc_max = %.10g!\n", hmat[1], Rcut, Rpc_max);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
	} // end if
	if (Rcut + 2.0*Rpc_max > (animg+1.0)*hmat[1]*0.5) {
      PRINTF("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("The box (%.10g) is too small for your Rcut (%.10g) to tream Hartree without images 2Rpc_max = %.10g!\n", hmat[1], Rcut, 2.0*Rpc_max);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
	} // end if

  for (int i=0; i<natm; i++) {
	vx[i]   = 0; vy[i]   = 0; vz[i]   = 0;
	fx0[i]  = 0; fy0[i]  = 0; fz0[i]  = 0;
	fx[i]   = 0; fy[i]   = 0; fz[i]   = 0;
	fx0g[i] = 0; fy0g[i] = 0; fz0g[i] = 0;
	fxg[i]  = 0; fyg[i]  = 0; fzg[i]  = 0;
  } // end for

#ifdef _FORCECHECK_
//========================================================================
// make a copy of the atoms: give the dummy structure its own memory!!!!!!

	double *xd = new double [natm]; double *yd = new double [natm]; double *zd = new double [natm];
	double *qd = new double [natm]; double *qtd = new double [natm];
	double *alpd = new double [natm]; double *betad = new double [natm]; double *Rpcd = new double [natm];

	double *vxd   = new double [natm]; double *vyd   = new double [natm]; double *vzd   = new double [natm];
	double *fx0d  = new double [natm]; double *fy0d  = new double [natm]; double *fz0d  = new double [natm];
	double *fxd   = new double [natm]; double *fyd   = new double [natm]; double *fzd   = new double [natm];
	double *fx0gd = new double [natm]; double *fy0gd = new double [natm]; double *fz0gd = new double [natm];
	double *fxgd  = new double [natm]; double *fygd  = new double [natm]; double *fzgd  = new double [natm];


  for (int i=0; i<natm; i++) {
	vxd[i]   = 0; vyd[i]   = 0; vzd[i]   = 0;
	fx0d[i]  = 0; fy0d[i]  = 0; fz0d[i]  = 0;
	fxd[i]   = 0; fyd[i]   = 0; fzd[i]   = 0;
	fx0gd[i] = 0; fy0gd[i] = 0; fz0gd[i] = 0;
	fxgd[i]  = 0; fygd[i]  = 0; fzgd[i]  = 0;

	xd[i]   = x[i]; yd[i]   = y[i]; zd[i]   = z[i];
	qd[i]   = q[i]; qtd[i]  = qt[i];
    alpd[i] = alp[i];
	betad[i] = beta[i];
	Rpcd[i] = Rpc[i];
  } // end for
#endif // _FORCECHECK_

  //========================================================================
  // store the atom information: this is safe because the dummy has its own memory

	atom_pos->natm = natm;
	atom_pos->x    = x;   atom_pos->y 	 = y;    atom_pos->z   = z;
	atom_pos->q    = q;
	atom_pos->qt   = qt;
	atom_pos->alp  = alp; atom_pos->beta = beta; atom_pos->Rpc = Rpc;

	atom_pos->vx	 = vx;   atom_pos->vy   = vy;   atom_pos->vz   = vz;
	atom_pos->fx0	 = fx0;  atom_pos->fy0  = fy0;  atom_pos->fz0  = fz0;
	atom_pos->fx	 = fx;   atom_pos->fy   = fy;   atom_pos->fz   = fz;
	atom_pos->fx0g = fx0g; atom_pos->fy0g = fy0g; atom_pos->fz0g = fz0g;
	atom_pos->fxg	 = fxg;  atom_pos->fyg  = fyg;  atom_pos->fzg  = fzg;

#ifdef _FORCECHECK_
	atom_pos_dummy->natm = natm;
	atom_pos_dummy->x 	 = xd;   atom_pos_dummy->y 	  = yd;    atom_pos_dummy->z   = zd;
	atom_pos_dummy->q 	 = qd;
	atom_pos_dummy->qt 	 = qtd;
	atom_pos_dummy->alp  = alpd; atom_pos_dummy->beta = betad; atom_pos_dummy->Rpc = Rpcd;

	atom_pos_dummy->vx	 = vxd;   atom_pos_dummy->vy   = vyd;   atom_pos_dummy->vz   = vzd;
	atom_pos_dummy->fx0	 = fx0d;  atom_pos_dummy->fy0  = fy0d;  atom_pos_dummy->fz0  = fz0d;
	atom_pos_dummy->fx	 = fxd;   atom_pos_dummy->fy   = fyd;   atom_pos_dummy->fz   = fzd;
	atom_pos_dummy->fx0g = fx0gd; atom_pos_dummy->fy0g = fy0gd; atom_pos_dummy->fz0g = fz0gd;
	atom_pos_dummy->fxg	 = fxgd;  atom_pos_dummy->fyg  = fygd;  atom_pos_dummy->fzg  = fzgd;
#endif // _FORCECHECK_

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
	atom_maps->natm_typ 			  = natm_typ;
	atom_maps->natm				      = natm;
	atom_maps->index_atm_typ    = index_atm_typ;
	atom_maps->natm_atm_typ     = natm_atm_typ;
	atom_maps->list_atm_by_typ  = list_atm_by_typ;
	atom_maps->atm_typ			    = atm_typ;
	atom_maps->natm_atm_typ_max = natm_atm_typ_max;

	gethinv(hmat, hmati, &volume, iperd);

	//========================================================================
	// store the simulation cell information
	cell->iperd = iperd;
	cell->alpb	= alpb;
	cell->gcut	= gcut;
	cell->Gcut	= Gcut;
	cell->Rcut	= Rcut;
	cell->nimg  = nimg;
	cell->animg = animg;

	for (int i=1; i<10; i++) {
		cell->hmat[i]	= hmat[i];
	}
	gethinv(hmat, hmati, &volume, iperd);
	for (int i=1; i<10; i++) {
		cell->hmati[i] = hmati[i];
	}
	cell->volume = volume;
} // end subroutine


//========================================================================
void fill_fgrid(FGRID *fgrid, ATOM_MAPS *atom_maps, ATOM_POS *atom_pos, CELL *cell, int rorder, int thetaorder, int phiorder, int model) {
//========================================================================
	int nf = rorder*thetaorder*phiorder;
	int natm_typ = atom_maps->natm_typ;
	int **list_atm_by_typ = atom_maps->list_atm_by_typ;
	double *alp = atom_pos->alp;
	double *beta = atom_pos->beta;
	double *Rpc = atom_pos->Rpc;
	double alpb = cell->alpb;

	double * xr_master = new double [rorder]; double * wr_master = new double [rorder];
	double * wr_bare = new double [rorder];
	double * xphi_master = new double [phiorder]; double * wphi_master = new double [phiorder];
	double * xtheta_master = new double [thetaorder]; double * wtheta_master = new double [thetaorder];
	double * psi0_2 = new double [rorder];

	int kind = 1; double aaa = -1; double bbb = 1;
	control_quad_rule(kind, thetaorder, aaa, bbb, wtheta_master, xtheta_master); // Legendre

	if (model == 2) {
		int type = 2; int iopt = 0;
		int ierr_zero; int ierr_ortho;
		gen_Gauss_quad_driver(type, rorder, iopt, xr_master, wr_master, &ierr_zero, &ierr_ortho);
		double norm = 4.0/sqrt(M_PI_QI);
		double check = 0.0;
		for (int i=0; i < rorder; i++) {
			double x2 = xr_master[i]*xr_master[i];
			psi0_2[i] = x2*norm;  // the weight function of the quadrature has the Gaussian
			wr_bare[i] = wr_master[i];
			wr_master[i] *= psi0_2[i]; // unitless
			check += wr_master[i];
		} // end for
#ifdef _INTEGRATION_CHECK_
		PRINTF("check = %g\n", check);
#endif
	} // end if

	if (model == 1) {
		double Rpc_master = 1.0;
		double delta_master = Rpc_master/((double) (rorder - 1));
		double check = 0.0;
		for (int i=0; i < rorder; i++) {
			double r = ((double) i)*delta_master;
			xr_master[i] = r; // depends on Rpc - only guy that needs to change if Rpc != 1.0
			double ak1 = M_PI_QI/Rpc_master;
			psi0_2[i] = sin(ak1*r)*sin(ak1*r)*2.0/Rpc_master; // ak1*r is unitless
			wr_bare[i] = delta_master;
			wr_master[i] = delta_master*psi0_2[i]; // unitless
			check += wr_master[i];
		} // end for
		wr_master[0] *= 0.5; wr_master[rorder-1] *= 0.5;
		wr_bare[0] *= 0.5; wr_bare[rorder-1] *= 0.5;
#ifdef _INTEGRATION_CHECK_
		PRINTF("check = %g\n", check);
#endif
	} // end if


  genphigrid(phiorder,wphi_master,xphi_master); // Fourier (equal space)


  //---------------------------------------------------------------------
  // use the master quadratures to create the full fgrid

  for (int i=0; i<natm_typ; i++) {
		int nang = thetaorder*phiorder;
		int nr = rorder;
  	fgrid[i].nf = nf;
  	fgrid[i].nr = rorder;
  	fgrid[i].ntheta = thetaorder;
  	fgrid[i].nphi = phiorder;
		fgrid[i].nang = nang;
  	fgrid[i].wf = new double [nf];
  	fgrid[i].xf = new double [nf];
  	fgrid[i].yf = new double [nf];
  	fgrid[i].zf = new double [nf];
  	fgrid[i].rf = new double [nf];
		fgrid[i].rho = new double [nf];
  	fgrid[i].xcostheta = new double [thetaorder];
  	fgrid[i].xphi = new double [phiorder];
		fgrid[i].Ylmf = new complex [nf];
  	fgrid[i].xr = new double [rorder];
  	fgrid[i].wr = new double [rorder];
  	fgrid[i].wr_bare = new double [rorder];
		fgrid[i].ylm = new double[nang];
		fgrid[i].wang = new double[nang];
		fgrid[i].rho_lm = new double[rorder];

		double **rho_scr    = new double *[nr];
		double *contig     = new double [nr*nang];
		for(int ir=0;ir<nr;ir++){rho_scr[ir] = &contig[(ir*nang)];}
		fgrid[i].rho_scr = rho_scr;

    double **pw_erfA = new double *[nr];
    double *contigrA  = new double [nr*nr];
    for(int ir=0;ir<nr;ir++){pw_erfA[ir] = &contigrA[(ir*nr)];}
		fgrid[i].pw_erfA = pw_erfA;

    double **pw_erfB = new double *[nr];
    double *contigrB  = new double [nr*nr];
    for(int ir=0;ir<nr;ir++){pw_erfB[ir] = &contigrB[(ir*nr)];}
		fgrid[i].pw_erfB = pw_erfB;

    double **pw_coul = new double *[nr];
    double *contigc  = new double [nr*nr];
    for(int ir=0;ir<nr;ir++){pw_coul[ir] = &contigc[(ir*nr)];}
		fgrid[i].pw_coul = pw_coul;

		int J = list_atm_by_typ[i][0];
		double alp_tmp = alp[J];
		double beta_tmp = beta[J];
		double Rpc_tmp = Rpc[J];
		fgrid[i].alp = alp_tmp;
		fgrid[i].beta = beta_tmp;
		fgrid[i].Rpc = Rpc_tmp;

    double *wf = fgrid[i].wf;
    double *xf = fgrid[i].xf;
    double *yf = fgrid[i].yf;
    double *zf = fgrid[i].zf;
    double *rf = fgrid[i].rf;
		double *rho = fgrid[i].rho;
    double *xr = fgrid[i].xr;
    double *wr = fgrid[i].wr;
    double *wr_bare_tmp = fgrid[i].wr_bare;
		double *xcostheta = fgrid[i].xcostheta;
		double *xphi = fgrid[i].xphi;
		double *ylm = fgrid[i].ylm;
		double *wang = fgrid[i].wang;

		double scale = 1.0;
		if (model == 1) {scale = Rpc_tmp;}
		if (model == 2) {scale = alp_tmp;}
		for (int ir=0; ir<rorder; ir++) {
			xr[ir] = xr_master[ir]*scale;
			wr[ir] = wr_master[ir];
			wr_bare_tmp[ir] = wr_bare[ir];
		} // end for ir

    double result4 = 0.0, resultA = 0.0;
    double result3 = 0.0;
    double beta2 = beta_tmp*beta_tmp;
		double alpb2 = alpb*alpb;
    for (int ir=0; ir<rorder; ir++) {
        for (int jr=0; jr<rorder; jr++) {
            double rgt = MAX(xr[ir],xr[jr]);
            double rlt = MIN(xr[ir],xr[jr]);
            double rd = rgt - rlt;
            double rs = rgt + rlt;
            double complicated, complicatedA;
            double simple;
            if (ir != 0 && jr != 0) {
                double part1 = (exp(-beta2*rs*rs) - exp(-beta2*rd*rd))/(2.0*beta_tmp*sqrt(M_PI_QI)*rgt*rlt);
                double part2 = (rd*erfc(beta_tmp*rd) - rs*erfc(beta_tmp*rs))/(2.0*rgt*rlt);
                double part3 = 1.0/rgt;
                complicated = part1 + part2 + part3;
                double part1A = (exp(-alpb2*rs*rs) - exp(-alpb2*rd*rd))/(2.0*alpb*sqrt(M_PI_QI)*rgt*rlt);
                double part2A = (rd*erfc(alpb*rd) - rs*erfc(alpb*rs))/(2.0*rgt*rlt);
                double part3A = 1.0/rgt;
                complicatedA = part1A + part2A + part3A;
                simple = 1.0/rgt;
            } else {
                if (ir !=0 || jr != 0) {
                    complicated = erf(beta_tmp*rgt)/rgt;
                    complicatedA = erf(alpb*rgt)/rgt;
                    simple = 1.0/rgt;
                } else {
                    complicated = 2.0*beta_tmp/sqrt(M_PI_QI);
                    complicatedA = 2.0*alpb/sqrt(M_PI_QI);
                    simple = 0.0;
                } // end if
            } // end if
            result3 += wr[ir]*wr[jr]*simple;
            result4 += wr[ir]*wr[jr]*complicated;
            resultA += wr[ir]*wr[jr]*complicatedA;
						pw_erfA[ir][jr] = complicatedA;
						pw_erfB[ir][jr] = complicated;
						pw_coul[ir][jr] = simple;
        } // end for jr
    } // end for ir
#ifdef _INTEGRATION_CHECK_
    PRINTF("2*HartreeselfGrid_pw:        %.14g\n", result3);
		PRINTF("2*HartreeselfGrid_pw_screen: %.14g\n", result4);
		PRINTF("2*HartreeselfGrid_pw_Ewald:  %.14g\n", resultA);
#endif

		double gamma2_inv = 0.5/(alp_tmp*alp_tmp) + 0.25/(beta_tmp*beta_tmp);
		double gamma = 1.0/sqrt(gamma2_inv);
#ifdef _INTEGRATION_CHECK_
		PRINTF("2*Hartreeself_screen: %.10g %.10g beta_tmp = %g\n", gamma/sqrt(M_PI_QI), result4, beta_tmp);
#endif
		for (int itheta=0; itheta<thetaorder; itheta++) { xcostheta[itheta] = xtheta_master[itheta];}
		for (int iphi=0; iphi<phiorder; iphi++) { xphi[iphi] = xphi_master[iphi];}

    int f = 0; int iii = 0;
		double result5 = 0.0, result6 = 0.0;
		for (int itheta=0; itheta<thetaorder; itheta++) {
       for (int iphi=0; iphi<phiorder; iphi++) {
					wang[iii] = wtheta_master[itheta]*wphi_master[iphi];
					ylm[iii]  = 1.0/sqrt(4.0*M_PI_QI);
    			for (int ir=0; ir<rorder; ir++) {
                double xsintheta = sqrt(1.0-(xtheta_master[itheta])*(xtheta_master[itheta]));
                wf[f] = wr[ir]*wtheta_master[itheta]*wphi_master[iphi];
                xf[f] = xr[ir]*xsintheta*cos(xphi_master[iphi]);
                yf[f] = xr[ir]*xsintheta*sin(xphi_master[iphi]);
                zf[f] = xr[ir]*xtheta_master[itheta];
                rf[f] = xr[ir];
								rho[f] = psi0_2[ir]/(4.0*M_PI_QI);
								result6 += rho[f]*wr_bare[ir]*wang[iii];
								result5 += wf[f];
                f++;
          } // end for iphi
					iii++;
        } // end for itheta
    } // end for ir
#ifdef _INTEGRATION_CHECK_
		PRINTF("result5: %.10g\n", result5/(4.0*M_PI_QI));
		PRINTF("result6: %.10g\n", result6);
#endif

  } // end for i: atom type

	for (int i=0; i<natm_typ; i++) {
		int lmax_tmp = 0; int nr = fgrid[i].nr; int nang = fgrid[i].nang;
		double *ylm = fgrid[i].ylm; double *wang = fgrid[i].wang;
		double *wr_bare = fgrid[i].wr_bare;
		double **pw_coul = fgrid[i].pw_coul;
		double **pw_erfA = fgrid[i].pw_erfA;
		double **pw_erfB = fgrid[i].pw_erfB;
		double *rho_in = fgrid[i].rho; double **rho_scr = fgrid[i].rho_scr;
		double *rho_lm = fgrid[i].rho_lm; double vself_coul = 0.0, vselfA = 0.0, vselfB = 0.0;
		screen_self_Hartree(lmax_tmp, nr, nang, ylm, wang, wr_bare, pw_coul, rho_in, rho_scr, rho_lm, &vself_coul);
		screen_self_Hartree(lmax_tmp, nr, nang, ylm, wang, wr_bare, pw_erfA, rho_in, rho_scr, rho_lm, &vselfA);
		screen_self_Hartree(lmax_tmp, nr, nang, ylm, wang, wr_bare, pw_erfB, rho_in, rho_scr, rho_lm, &vselfB);
#ifdef _INTEGRATION_CHECK_
		PRINTF("vself_coul, vself_screenA, vself_screenB: %.14g %.14g %.14g\n", 2.0*vself_coul, 2.0*vselfA, 2.0*vselfB);
#endif
  } // end for i: atom type
} // end routine
