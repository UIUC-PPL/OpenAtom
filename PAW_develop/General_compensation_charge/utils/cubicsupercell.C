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

void readtoendofline(FILE *);
void gethinv(double *, double *, double *, int);

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
	CELL      cell;

	int natm_typ;       	// number of atom types
	int natm;				// number of atoms
	int natmnew;		 	// new number of atoms
	int *index_atm_typ;		// index of atom type of each atom
	int *natm_atm_typ;		// the number of atoms of each type
	int **list_atm_by_typ;	// list of atoms sorted by atom type
	NAME *atm_typ;	    	// names of the atom types

	int iperd = 3;			// periodicity

	double *x,*y,*z,*q,*qt,*alp,*beta, *Rpc;    						 // atom positions and core Gaussian parameters
	double hmat[10];				  // the simulation box
	double hmatnew[10];				  // the new simulation box
	double hmati[10];				  // inverse simulation box
	double hmatinew[10];			  // newinverse simulation box
	double volume;				  	  // simulation box volume
	double volumenew;				  // new simulation box volume

	char fnameIn[MAXLINE];
	char fnameOut[MAXLINE];

	FILE *fp;

	//==========================================================================
	// Tell everyone what you are doing

	PRINTF("\n");
	PRINT_LINE_STAR
	PRINTF("replicate the unit cell\n");
	PRINT_LINE_DASH

    //=========================================================================
    //             Check for input file                                 

    if(argc < 6) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("No input file specified!\n");
		PRINTF("Run it like: ./cubicsupercell.x atom.in nx ny nz atom.out\n");
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
    }/*endif*/

	//==========================================================================
	// Read the user specified input file: atoms, Ewald parameters, core density parameters and box
	// store the information

	int nx = atoi(argv[2]);
	int ny = atoi(argv[3]);
	int nz = atoi(argv[4]);
	if (nx <= 0 || ny <= 0 || nz <= 0) {
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("Number of replicas (%d,%d,%d) must be > 0!\n", nx, ny, nz);
		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
    }/*endif*/

	strcpy(fnameIn, argv[1]);
	PRINTF("\nReading input parameters and atom positions from %s\n\n",fnameIn);
	fp = fopen(fnameIn,"r");
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
			} 
			natm_atm_typ_max = MAX(natm_atm_typ_max,natm_atm_typ[i]);
			PRINT_LINE_DASH;
			PRINTF("Number of atom %s: %d\n",atm_typ[i], natm_atm_typ[i]);
		} // end for i
		PRINT_LINE_STAR;
		PRINTF("\n");
		double a,b;
		fscanf(fp,"%lf %lf",&a, &b); readtoendofline(fp);

		x = new double [natm]; y = new double [natm]; z = new double [natm];
		q = new double [natm]; qt = new double [natm];
		Rpc = new double [natm]; alp = new double [natm]; beta = new double [natm];
		double beta_unitless = 1.0;

		index_atm_typ = new int [natm];
		list_atm_by_typ = new int *[natm_typ];
		for (int i=0; i<natm_typ; i++) { list_atm_by_typ[i] = new int [natm_atm_typ_max]; }
		for (int i=0; i<natm; i++) {
			fscanf(fp,"%lf %lf %lf %lf %lf %lf %d",&x[i],&y[i],&z[i], &q[i], &qt[i], &Rpc[i], &index_atm_typ[i]); readtoendofline(fp);
			alp[i] = 1.8/Rpc[i];
			beta[i] = alp[i]*beta_unitless;
		} //end for i
		fscanf(fp, "%lf %lf %lf", &hmat[1], &hmat[4], &hmat[7]); readtoendofline(fp);
		fscanf(fp, "%lf %lf %lf", &hmat[2], &hmat[5], &hmat[8]); readtoendofline(fp);
		fscanf(fp, "%lf %lf %lf", &hmat[3], &hmat[6], &hmat[9]); readtoendofline(fp);
	fclose(fp);

//=======================================================================================================
//replicate the old atoms

	natmnew = natm*nx*ny*nz;
	hmatnew[1] = hmat[1]*((double) nx);
	hmatnew[5] = hmat[5]*((double) ny);
	hmatnew[9] = hmat[9]*((double) nz);

	strcpy(fnameOut, argv[5]);
	PRINTF("\nWriting input parameters and new atom positions to %s\n\n",fnameOut);
	fp = fopen(fnameOut,"w");
		fprintf(fp,"%d %d\n", natmnew, natm_typ);
		for (int i=0; i<natm_typ; i++) { 
			fprintf(fp,"%d %d %s\n", i, natm_atm_typ[i]*nx*ny*nz, atm_typ[i]);
			PRINTF("Number of atom %s: %d\n",atm_typ[i], natm_atm_typ[i]*nx*ny*nz);
		} // end for i
		PRINT_LINE_STAR;
		PRINTF("\n");
		fprintf(fp,"%lf %lf\n",a, b);

		for (int i=0; i<natm; i++) {
			for (int ix=0; ix<nx; ix++) {
			for (int iy=0; iy<ny; iy++) {
			for (int iz=0; iz<nz; iz++) {
				double xxx = x[i] + ((double) ix)*hmat[1];
				double yyy = y[i] + ((double) iy)*hmat[5];
				double zzz = z[i] + ((double) iz)*hmat[9];
				fprintf(fp,"%lf %lf %lf %lf %lf %lf %d\n",xxx, yyy, zzz, q[i], qt[i], Rpc[i], index_atm_typ[i]);
			}}}
		} //end for i
		fprintf(fp, "%lf %lf %lf\n", hmatnew[1], hmatnew[4], hmatnew[7]);
		fprintf(fp, "%lf %lf %lf\n", hmatnew[2], hmatnew[5], hmatnew[8]);
		fprintf(fp, "%lf %lf %lf\n", hmatnew[3], hmatnew[6], hmatnew[9]);
	fclose(fp);

	return 0;

  //--------------------------------------------------------------------------
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// readtoendofline: Function to read to end of line in read_coord files     
//==========================================================================
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol){ch=fgetc(fp);}
  if(ch==EOF){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("ERROR: Unexpected end of file reached          \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif
}// end routine 
//==========================================================================


void gethinv(double *hmat, double *hmati, double *deth, int iperd)

  /*===============================================================*/
{/*begin routine */
  /*===============================================================*/
  double vol;
  int i;
  /*===============================================================*/
  /* gets inverse, hmati, of the iperd dimensional matrix hmat */
  /* (stored as a 3 x 3) */

  *deth = 0.0;
  for(i=1;i<=9;i++){hmati[i]=0.0;}

  /*===============================================================*/
  /* Perd=3 */

  if (iperd == 3) {
    vol = (hmat[1] * (hmat[5] * hmat[9] - hmat[8] * hmat[6]) +
        hmat[4] * (hmat[8] * hmat[3] - hmat[2] * hmat[9]) +
        hmat[7] * (hmat[2] * hmat[6] - hmat[5] * hmat[3]));
    *deth = vol;
    hmati[1] = (hmat[5] * hmat[9] - hmat[8] * hmat[6]) / vol;
    hmati[5] = (hmat[1] * hmat[9] - hmat[7] * hmat[3]) / vol;
    hmati[9] = (hmat[1] * hmat[5] - hmat[4] * hmat[2]) / vol;
    hmati[4] = (hmat[7] * hmat[6] - hmat[4] * hmat[9]) / vol;
    hmati[2] = (hmat[3] * hmat[8] - hmat[2] * hmat[9]) / vol;
    hmati[7] = (hmat[4] * hmat[8] - hmat[7] * hmat[5]) / vol;
    hmati[3] = (hmat[2] * hmat[6] - hmat[3] * hmat[5]) / vol;
    hmati[8] = (hmat[7] * hmat[2] - hmat[8] * hmat[1]) / vol;
    hmati[6] = (hmat[3] * hmat[4] - hmat[6] * hmat[1]) / vol;
  }/*endif*/

  /*===============================================================*/
  /* Perd=2 */

  if (iperd == 2) {
    vol = hmat[1] * hmat[5] - hmat[4] * hmat[2];
    hmati[1] = hmat[5] / vol;
    hmati[5] = hmat[1] / vol;
    hmati[4] = -hmat[4] / vol;
    hmati[2] = -hmat[2] / vol;
    hmati[9] = 1. / hmat[9];
    *deth = vol * hmat[9];
  }/*endif*/

  /*===============================================================*/
  /* Perd=1,0,cluster_ewald */

  if(iperd <=1 || iperd==4) {
    *deth = hmat[1]*hmat[5]*hmat[9];
    hmati[1] = 1.0/hmat[1];
    hmati[5] = 1.0/hmat[5];
    hmati[9] = 1.0/hmat[9];
  }/*endif*/

  /*===============================================================*/
  /* Errors */

  if((*deth)==0.0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The present volume is zero.                 \n");
    printf("If this is not an error in your input data, \n");
    printf("contact technical support                  \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  } /*endif*/

  /*---------------------------------------------------------------*/
} /* gethinv */
/*===============================================================*/

