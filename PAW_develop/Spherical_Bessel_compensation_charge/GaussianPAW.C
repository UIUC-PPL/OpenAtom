//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  This program takes classical input and make the PIMD ring polymers
//  To run this program : executable inputfile
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "standard_include.h"
#include "ckcomplex.h"
#include "fgrid.h"
#include "GaussianPAW.h"

//=================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================
//  Generate Uniform Random Numbers
//=================================================================
#define four_pi_inv    (1.0/(4.0*M_PI_QI))
#define MODULUS_R    2147483647 // DON'T CHANGE THIS VALUE       
#define MULTIPLIER_R 48271      // DON'T CHANGE THIS VALUE       
//=================================================================
// Random returns a pseudo-random real number uniformly distributed 
// between 0.0 and 1.0. 
//=================================================================
double altRandom(long *seed){
  long t;
  const long Q = MODULUS_R / MULTIPLIER_R;
  const long R = MODULUS_R % MULTIPLIER_R;

  t = MULTIPLIER_R * (seed[0] % Q) - R * (seed[0] / Q);
  if(t > 0){
    seed[0] = t;
  }else {
    seed[0] = t + MODULUS_R;
  }//endif
  return ((double) seed[0] / MODULUS_R);
}
//=================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Generate Gaussian Random numbers 
//===============================================================
void gaussran(int nran, long *seed, double *gauss)
  //========================================================================
{//begin routine
  //========================================================================
  //             Local variable declarations                                
  int i,loop;
  double twopi,rad2,al,r,phi,arg;
  //========================================================================
  // I) Constants 

  twopi = 2.0*M_PI_QI;
  rad2  = sqrt(2.0);
  loop  = nran/2;

  //========================================================================
  // II) Make nran (or nran-1 if nran odd) Gaussian random numbers          
  for(i=1;i<=loop;i++){
    //------------------------------------------------------------------------
    // A) uniform random numbers in r and phi 
    r   = altRandom(seed); 
    phi = altRandom(seed); 
    r   = MAX(r,1e-30);
    r   = MIN(r,1.0);
    //------------------------------------------------------------------------
    // B) Gaussify in x and y
    al  = sqrt(-log(r))*rad2;
    arg = twopi*phi;
    gauss[2*i-1] = al*cos(arg);
    gauss[2*i]   = al*sin(arg);
  }//endfor
  //========================================================================
  // III) Make one more if nran is odd 

  if((nran % 2)!=0){
    r   = altRandom(seed); 
    phi = altRandom(seed); 
    r   = MAX(r,1e-30);
    r   = MIN(r,1.0);
    arg = twopi*phi;
    al  = sqrt(-log(r))*rad2;
    gauss[nran] = al*cos(arg);
  }//endif
  //------------------------------------------------------------------------
}//end routine
//========================================================================



/*===============================================================*/
/*  Inverse of a 3x3 */
/*===============================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===============================================================*/

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



//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the PAW energies in real space analytically
//===================================================================
void computePAWreal(ATOM_MAPS *atom_maps, ATOM_POS *atom_pos, CELL *cell, ESTRUCT *energy) 
//===================================================================
{ // begin routine
//===================================================================
// read in parameters from the structures
	int natm 	= atom_maps->natm;
	double *x 	= atom_pos->x;   double *y   = atom_pos->y;   double *z   = atom_pos->z; 
	double *fx0	= atom_pos->fx0; double *fy0 = atom_pos->fy0; double *fz0 = atom_pos->fz0;
	double *fx	= atom_pos->fx;  double *fy	 = atom_pos->fy;  double *fz  = atom_pos->fz;
	double *q 	= atom_pos->q;   double *qt  = atom_pos->qt;
	double *alp = atom_pos->alp; double *beta = atom_pos->beta; 
	double *Rpc = atom_pos->Rpc;
	double alpb = cell->alpb;    double *hmat = cell->hmat;    double *hmati = cell->hmati; double Rcut = cell->Rcut;
	int iperd 	= cell->iperd;

//===================================================================
// zero all the energies 
  double ENN = 0.0, EeNself = 0.0, EeN = 0.0, EHarself = 0.0, EHar = 0.0, EHarselfscr = 0.0, EHarscr = 0.0;

//===================================================================
// energy contributions for clusters
//===================================================================
 if (iperd == 0) { 
	double prec = 1.0/sqrt(2.0*M_PI_QI);
	double precH = 0.5/sqrt(M_PI_QI);
	double preceN = -2.0/sqrt(M_PI_QI);
	for (int i=0; i<natm; i++) {
		double gamma2_inv   = 0.5/(alp[i]*alp[i]) + 0.25/(beta[i]*beta[i]);
		double gamma   = 1.0/sqrt(gamma2_inv);
		EHarself      += qt[i]*qt[i]*alp[i];
		EHarselfscr   += qt[i]*qt[i]*gamma; 
		EeNself       += qt[i]*q[i]*alp[i];
	} // end for
	EHarself    *= prec;
	EHarselfscr *= precH;
	EeNself     *= preceN;

  for (int i=0; i<natm; i++) {
      for (int j=0; j<natm; j++) {
          if (j != i) {
              double alp_ij  = alp[i]*alp[j]/sqrt(alp[i]*alp[i] + alp[j]*alp[j]);
              double R_ij    = dist((x[i]-x[j]), (y[i]-y[j]), (z[i]-z[j]));
			  double fgerfc_i, fgerfc_j, fgerfc_ij;
			  double gerfc_tmp_i =  gerfc(R_ij, alp[i], &fgerfc_i);
			  double gerfc_tmp_j =  gerfc(R_ij, alp[j], &fgerfc_j);
			  double gerfc_tmp_ij =  gerfc(R_ij, alp_ij, &fgerfc_ij);
			  double gerf_tmp_i =  1.0 - gerfc_tmp_i;
			  double gerf_tmp_ij = 1.0 - gerfc_tmp_ij;
			  double R_ij_inv = 1.0/R_ij;
              ENN           += 0.5*q[i]*q[j]*R_ij_inv;
              EeN           -= 1.0*qt[i]*q[j]*gerf_tmp_i*R_ij_inv;
              EHar          += 0.5*qt[i]*qt[j]*gerf_tmp_ij*R_ij_inv;
			  double coeff1  = -q[i]*q[j];
			  double coeff2  = -(qt[i]*q[j]*erfc_a_r_over_r(R_ij,alp[i],gerfc_tmp_i, fgerfc_i) 
								+ q[i]*qt[j]*erfc_a_r_over_r(R_ij,alp[j],gerfc_tmp_j, fgerfc_j));
			  double coeff3  = qt[i]*qt[j]*erfc_a_r_over_r(R_ij,alp_ij, gerfc_tmp_ij, fgerfc_ij);
			  double ftmp = (coeff1 + coeff2 + coeff3)*R_ij_inv*R_ij_inv*R_ij_inv;
			  fx0[i]		  += ftmp*(x[j] - x[i]);  // 0D, no periodic images
			  fy0[i]		  += ftmp*(y[j] - y[i]); 
			  fz0[i]		  += ftmp*(z[j] - z[i]); 
          } // end if
      } // end for 
  } // end for
	EeN  += EeNself;
	EHarscr = EHar + EHarselfscr;
	EHar += EHarself;
 } // end if iperd = 0

//===================================================================
// zero all the energies 
	double EeNshortself = 0.0, EeNshort = 0.0, EHarshortself = 0.0, EHarshort = 0.0, EHarshortscr = 0.0;
	double EHarshortselfscr = 0.0;
	double ENNshort = 0.0;
	double ENNselflong = 0.0;

//===================================================================
// energy contributions for 3D PDC
//===================================================================

 if (iperd == 3) { 
  for (int i=0; i<natm; i++) {
	double alpb_i     = alpb*alp[i]/sqrt(alp[i]*alp[i] + alpb*alpb);
	double alpb_ii 	  = alpb*alp[i]*alp[i]/sqrt(alp[i]*alp[i]*alp[i]*alp[i] + 2*alpb*alpb*alp[i]*alp[i]);
	double betab_ii   = beta[i]*alp[i]*alp[i]/sqrt(alp[i]*alp[i]*alp[i]*alp[i] + 2*beta[i]*beta[i]*alp[i]*alp[i]);
	EeNshortself  	 += 1.0/sqrt(M_PI_QI)*(-2.0*(alp[i] - alpb_i))*qt[i]*q[i];
	EHarshortself 	 += 1.0/sqrt(M_PI_QI)*(alp[i]/sqrt(2.0) - alpb_ii)*qt[i]*qt[i];
	EHarshortselfscr += 1.0/sqrt(M_PI_QI)*(betab_ii - alp[i]/sqrt(2.0))*qt[i]*qt[i];
	ENNselflong   	 += 1.0*alpb/sqrt(M_PI_QI)*q[i]*q[i];
  } // end for
	EHarshortselfscr += EHarshortself;

  for (int i=0; i<natm; i++) {
    double alpb_i = alpb*alp[i]/sqrt(alp[i]*alp[i] + alpb*alpb);
    for (int j=0; j<natm; j++) {
        if (j != i) {
    		double alpb_j = alpb*alp[j]/sqrt(alp[j]*alp[j] + alpb*alpb);
            double alp_ij  = alp[i]*alp[j]/sqrt(alp[i]*alp[i] + alp[j]*alp[j]);
            double alpb_ij = alpb*alp[i]*alp[j]/sqrt(alp[i]*alp[i]*alp[j]*alp[j] + alpb*alpb*alp[i]*alp[i] + alpb*alpb*alp[j]*alp[j]);
			double dx = x[j]-x[i];
			double dy = y[j]-y[i];
			double dz = z[j]-z[i];
			dx -= hmat[1]*NINT(dx*hmati[1]);
			dy -= hmat[5]*NINT(dy*hmati[5]);
			dz -= hmat[9]*NINT(dz*hmati[9]);
			double r2 = dx*dx + dy*dy + dz*dz;
			double R_ij = sqrt(r2);
			if (R_ij < Rcut + Rpc[i] + Rpc[j]) { // Hartree allowed
                double fgerfc_bar_ij, fgerfc_ij;
                double gerfc_tmp_ij =  gerfc(R_ij, alp_ij, &fgerfc_ij);
                double gerfc_tmp_bar_ij =  gerfc(R_ij, alpb_ij, &fgerfc_bar_ij);
                double R_ij_inv = 1.0/R_ij;
            	EHarshort     += 0.5*qt[i]*qt[j]*(gerfc_tmp_bar_ij - gerfc_tmp_ij)*R_ij_inv;
				double coeff1 = 0.0;
				double coeff2 = 0.0;
				double coeff3  = qt[i]*qt[j]*(erfc_a_r_over_r(R_ij,alp_ij,gerfc_tmp_ij,fgerfc_ij) - 
											  erfc_a_r_over_r(R_ij,alpb_ij, gerfc_tmp_bar_ij, fgerfc_bar_ij));
				if (R_ij < Rcut + Rpc[i]) { // eN allowed
				    double fgerfc_i, fgerfc_j, fgerfc_bar_i, fgerfc_bar_j;
				    double gerfc_tmp_i =  gerfc(R_ij, alp[i], &fgerfc_i);
				    double gerfc_tmp_j =  gerfc(R_ij, alp[j], &fgerfc_j);
				    double gerfc_tmp_bar_i =  gerfc(R_ij, alpb_i, &fgerfc_bar_i);
				    double gerfc_tmp_bar_j =  gerfc(R_ij, alpb_j, &fgerfc_bar_j);
					coeff2  = -(qt[i]*q[j]*(erfc_a_r_over_r(R_ij,alp[i],gerfc_tmp_i,fgerfc_i) - erfc_a_r_over_r(R_ij,alpb_i,gerfc_tmp_bar_i,fgerfc_bar_i)) 
								+q[i]*qt[j]*(erfc_a_r_over_r(R_ij,alp[j],gerfc_tmp_j,fgerfc_j) - erfc_a_r_over_r(R_ij,alpb_j,gerfc_tmp_bar_j,fgerfc_bar_j)));
            		EeNshort      -= 1.0*qt[i]*q[j]*(gerfc_tmp_bar_i - gerfc_tmp_i)*R_ij_inv;
					if (R_ij < Rcut) { // NN allowed
					   double fgerfc_bar;
					   double gerfc_tmp_bar =  gerfc(R_ij, alpb, &fgerfc_bar);
            			ENNshort      += 0.5*q[i]*q[j]*gerfc_tmp_bar*R_ij_inv;
						coeff1  = -q[i]*q[j]*(1.0 + erfc_a_r_over_r(R_ij,alpb,gerfc_tmp_bar,fgerfc_bar));
					} // end if NN allowed 
				} // end if eN allowed
				double ftmp = (coeff1 + coeff2 + coeff3)*R_ij_inv*R_ij_inv*R_ij_inv;
				fx[i]		  += ftmp*dx; 
				fy[i]		  += ftmp*dy; 
				fz[i]		  += ftmp*dz; 
			} // end if Hartree allowed
     	} // end if j!= i
    } // end for  j
  } // end for i
    EeNshort  += EeNshortself;
	EHarshortscr = EHarshort + EHarshortselfscr;
    EHarshort += EHarshortself;
 } // end if iperd = 3


//===================================================================
// Put the energies in the structure

	energy->ENN.E                = ENN; 
	energy->ENNshort.E           = ENNshort; 
	energy->EeNself.E            = EeNself; 
	energy->EeN.E                = EeN; 
	energy->EHarself.E           = EHarself; 
	energy->EHar.E               = EHar; 
	energy->EHarscr.E            = EHarscr; 
	energy->EeNshortself.E       = EeNshortself; 
	energy->EeNshort.E           = EeNshort; 
	energy->EHarshortself.E      = EHarshortself; 
	energy->EHarshort.E          = EHarshort; 
	energy->EHarshortscr.E       = EHarshortscr; 
	energy->ENNselflong.E        = ENNselflong;
	energy->EHarselfscr.E        = EHarselfscr;
	energy->EHarshortselfscr.E   = EHarshortselfscr;


//===================================================================
} // end routine
//===================================================================

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the PAW real space energies on Grid
//===================================================================
void computePAWGrid(int lmax, ATOM_MAPS *atom_maps, ATOM_POS *atom_pos, CELL *cell, ESTRUCT *energy, FGRID *fgrid) 
//===================================================================
{ // begin routine
//===================================================================
// read in parameters from the structures
	int natm			  = atom_maps->natm;
	int natm_typ 		  = atom_maps->natm_typ;
	int *natm_atm_typ 	  = atom_maps->natm_atm_typ;
	int **list_atm_by_typ = atom_maps->list_atm_by_typ;
	double *x 			= atom_pos->x;    double *y    = atom_pos->y;    double *z    = atom_pos->z; 
	double *fx0g		= atom_pos->fx0g; double *fy0g = atom_pos->fy0g; double *fz0g = atom_pos->fz0g;
	double *fxg			= atom_pos->fxg;  double *fyg  = atom_pos->fyg;  double *fzg  = atom_pos->fzg;
	double *q 			= atom_pos->q;    double *qt   = atom_pos->qt;
	double *alp  = atom_pos->alp;
	double *beta = atom_pos->beta;
	double *Rpc  = atom_pos->Rpc;
	double alpb   		= cell->alpb;     double *hmat = cell->hmat;	 double *hmati = cell->hmati;
	int iperd 			= cell->iperd;
	double Rcut = cell->Rcut;

//===================================================================
// local variables
	int rorder     = fgrid[0].nr;
	int thetaorder = fgrid[0].ntheta;
	int phiorder   = fgrid[0].nphi;
	int nf         = fgrid[0].nf;
//-------------------------------------------------------------------
// nuclear-nuclear energy is the same
    double ENNGrid         = energy->ENN.E;
    double ENNshortGrid    = energy->ENNshort.E;
    double ENNselflongGrid = energy->ENNselflong.E;

//===================================================================
// zero all the energies 
    double EeNGrid      = 0.0, EeNselfGrid = 0.0, EHarGrid = 0.0, EHarscrGrid = 0.0, EHarselfGrid = 0.0;
    double EeNshortGrid = 0.0, EeNshortselfGrid = 0.0, EHarshortGrid = 0.0, EHarshortscrGrid = 0.0, EHarshortselfGrid = 0.0;
	double EHarshortselfscrGrid = 0.0;
	double EHarselfscrGrid = 0.0;

//===================================================================
// energy contributions for clusters
//===================================================================
 if (iperd == 0) { 

clock_t start, end;
//===================================================================
// Adding the ENN part force 
	start = clock();
  for (int i=0; i<natm; i++) {
      for (int j=0; j<natm; j++) {
          if (j != i) {
              double R_ij    = dist((x[i]-x[j]), (y[i]-y[j]), (z[i]-z[j]));
			  double coeff1  = -q[i]*q[j];
			  fx0g[i]		+= (coeff1)/(R_ij*R_ij*R_ij)*(x[j] - x[i]); 
			  fy0g[i]		+= (coeff1)/(R_ij*R_ij*R_ij)*(y[j] - y[i]); 
			  fz0g[i]		+= (coeff1)/(R_ij*R_ij*R_ij)*(z[j] - z[i]); 
          } // end if
      } // end for 
  } // end for
  end = clock();
  PRINTF("0D ENN finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);
	

// 0D eN and Har no self
	int pair = 0;
	start = clock();
	for (int jtyp=0; jtyp<natm_typ; jtyp++) {
		double *xf_jtyp     = fgrid[jtyp].xf;
		double *yf_jtyp     = fgrid[jtyp].yf;
		double *zf_jtyp     = fgrid[jtyp].zf;
		double *wf_jtyp     = fgrid[jtyp].wf;
		//for(int ktyp=jtyp; ktyp<natm_typ; ktyp++){
		for(int ktyp=0; ktyp<natm_typ; ktyp++){
			double *xf_ktyp     = fgrid[ktyp].xf;
			double *yf_ktyp     = fgrid[ktyp].yf;
			double *zf_ktyp     = fgrid[ktyp].zf;
			double *wf_ktyp     = fgrid[ktyp].wf;
			for(int j=0; j<natm_atm_typ[jtyp];j++){
				int J = list_atm_by_typ[jtyp][j];
				double Ncoeff_J = four_pi_inv;
				//int kstart = (ktyp == jtyp ? (j+1):0);
				//for(int k=kstart;k<natm_atm_typ[ktyp];k++){
				for(int k=0;k<natm_atm_typ[ktyp];k++){
					int K = list_atm_by_typ[ktyp][k];
					if (K != J) {
						pair++;
						double Ncoeff_K = four_pi_inv;
						for (int f1=0; f1 < nf; f1++) {
							double dx_J   = xf_jtyp[f1]-x[K]+x[J];
							double dy_J   = yf_jtyp[f1]-y[K]+y[J];
							double dz_J   = zf_jtyp[f1]-z[K]+z[J];
							double r2_J   = dx_J*dx_J + dy_J*dy_J + dz_J*dz_J;
							double r_J    = sqrt(r2_J);
							double tmp_J  = -Ncoeff_J*qt[J]*q[K]*wf_jtyp[f1]/r_J;
							EeNGrid 	 += tmp_J;
							double dx_K   = xf_ktyp[f1]-x[J]+x[K];
							double dy_K   = yf_ktyp[f1]-y[J]+y[K];
							double dz_K   = zf_ktyp[f1]-z[J]+z[K];
							double r2_K   = dx_K*dx_K + dy_K*dy_K + dz_K*dz_K;
							double r_K    = sqrt(r2_K);
							double tmp_K  = -Ncoeff_K*qt[K]*q[J]*wf_ktyp[f1]/r_K;
							double coeffJ = tmp_J/(r_J*r_J);
							double coeffK = tmp_K/(r_K*r_K);
							fx0g[J]		 += coeffJ*dx_J - coeffK*dx_K;
							fy0g[J]		 += coeffJ*dy_J - coeffK*dy_K;
							fz0g[J]		 += coeffJ*dz_J - coeffK*dz_K;
							for (int f2=0; f2 < nf; f2++) {
								double Hdx_J 	= xf_jtyp[f1]-xf_ktyp[f2]-x[K]+x[J];
								double Hdy_J 	= yf_jtyp[f1]-yf_ktyp[f2]-y[K]+y[J];
								double Hdz_J 	= zf_jtyp[f1]-zf_ktyp[f2]-z[K]+z[J];
								double Hr2_J 	= Hdx_J*Hdx_J + Hdy_J*Hdy_J + Hdz_J*Hdz_J;
								double Hr_J 	= sqrt(Hr2_J);
								double Htmp_J  	= 0.5*Ncoeff_J*Ncoeff_K*qt[J]*qt[K]*wf_jtyp[f1]*wf_ktyp[f2]/Hr_J;
								EHarGrid 	   += Htmp_J; 
								double Hdx_K 	= xf_ktyp[f1]-xf_jtyp[f2]+x[K]-x[J];
								double Hdy_K 	= yf_ktyp[f1]-yf_jtyp[f2]+y[K]-y[J];
								double Hdz_K 	= zf_ktyp[f1]-zf_jtyp[f2]+z[K]-z[J];
								double Hr2_K 	= Hdx_K*Hdx_K + Hdy_K*Hdy_K + Hdz_K*Hdz_K;
								double Hr_K 	= sqrt(Hr2_K);
								double Htmp_K 	= 0.5*Ncoeff_K*Ncoeff_J*qt[J]*qt[K]*wf_ktyp[f1]*wf_jtyp[f2]/Hr_K;
								double HcoeffJ  = Htmp_J/(Hr_J*Hr_J);
								double HcoeffK  = Htmp_K/(Hr_K*Hr_K);
								fx0g[J]		   += HcoeffJ*Hdx_J - HcoeffK*Hdx_K; 
								fy0g[J]		   += HcoeffJ*Hdy_J - HcoeffK*Hdy_K; 
								fz0g[J]		   += HcoeffJ*Hdz_J - HcoeffK*Hdz_K; 
							} // end for f2
						} // end if
					} // end for f1
				} // end for k
			} //end for j
		} // end for ktyp
 	} // end for jtyp
  end = clock();
  PRINTF("(%d) pairs\n",pair);
  PRINTF("0D EeN and EHar noself finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);
//==============================================================================
//eN and Har self term
	start = clock();
	for (int jtyp=0; jtyp<natm_typ; jtyp++) {
        double *wf        = fgrid[jtyp].wf;
        double *rf        = fgrid[jtyp].rf;
        double *xf        = fgrid[jtyp].xf;
        double *yf        = fgrid[jtyp].yf;
        double *zf        = fgrid[jtyp].zf;
		double *xcostheta = fgrid[jtyp].xcostheta;
		double *xphi      = fgrid[jtyp].xphi;
		complex *Ylmf     = fgrid[jtyp].Ylmf;
		double alpJ		  = fgrid[jtyp].alp;
		double betaJ	  = fgrid[jtyp].beta;
		int J = list_atm_by_typ[jtyp][0];
		double Ncoeff = four_pi_inv;
		double wght = ((double) natm_atm_typ[jtyp]);

		for (int f1=0; f1 < nf; f1++) {
			double r = rf[f1];
			if (r > 0) {
				EeNselfGrid += -1.0*wght*Ncoeff*qt[J]*q[J]*wf[f1]/r;
			} // end if
		} // end for f1

#ifdef _FGRIDTEST_		
		complex *Ylmf_test       = new complex [nf];
		complex *Ylpmpf_test     = new complex [nf];
		double atmp			 = fgrid[jtyp].alp;
		for (int l=0; l <= lmax; l++) {
			for (int lp = 0; lp <= lmax; lp++) {
				for (int m=-l; m<= l; m++) {
					for (int mp = -lp; mp <= lp ; mp++) {
						complex result  	= complex(0.0,0.0);
						double  result2x  	= 0.0;
						double  result2y  	= 0.0;
						double  result2z  	= 0.0;
						gen_Ylmf (rorder, thetaorder, xcostheta, phiorder, xphi, l, m, Ylmf_test);
						gen_Ylmf (rorder, thetaorder, xcostheta, phiorder, xphi, lp, mp, Ylpmpf_test);
						for (int f=0; f < nf; f++) {
							result   += wf[f]*Ylmf_test[f]*Ylpmpf_test[f].conj()*(4.0/sqrt(M_PI))*atmp*atmp*atmp;
							result2x += wf[f]*xf[f]*xf[f];
							result2y += wf[f]*yf[f]*yf[f];
							result2z += wf[f]*zf[f]*zf[f];
						} // end for f
						PRINTF("\nfgrid test result:\n");
						PRINTF("int Y(%d %d), Y(%d %d) exp(-a^2*r^2)r^2 = (%lf, %lf) !\n",l,m,lp,mp,result.re, result.im);
						PRINTF("int x^2 r^2 exp(-a^2^r2)= (%lf) !\n",result2x);
						PRINTF("int y^2 r^2 exp(-a^2^r2)= (%lf) !\n",result2y);
						PRINTF("int z^2 r^2 exp(-a^2^r2)= (%lf) !\n",result2z);
					} // end for mp
				}// end for m
			} // end for lp
		} // end for l
#endif
		complex tmpsum = (0.0, 0.0);

		for (int l=0; l <= lmax; l++) {
			double l_re = (double) l;
			for (int m=-l; m<=l; m++) {
				gen_Ylmf (rorder, thetaorder, xcostheta, phiorder, xphi, l, m, Ylmf);
				for (int f1=0; f1 < nf; f1++) {
					for (int f2=0; f2 < nf; f2++) {
						double rmin = MIN(rf[f1],rf[f2]);
						double rmax = MAX(rf[f1],rf[f2]);
						double rrat = rmin/rmax;
						double pref = Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
						double pref_harm = 4.0*M_PI_QI/(2.0*l_re + 1.0);
						if (rmax > 0) {
							double pref_r = pow(rrat,l_re)/rmax;
							complex temp = 0.5*pref*pref_harm*pref_r*(Ylmf[f1]*Ylmf[f2].conj());
							tmpsum += temp*wght;
						} // end if
					} //end for f2
				} //end for f1
			} //end for m
		} //end for l
		EHarselfGrid += tmpsum.re;
		
		for (int f1=0; f1 < nf; f1++) {
			for (int f2=0; f2 < nf; f2++) {
				double pref = 0.5*wght*Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
				if (f1 != f2) {
					double dx   = xf[f1] - xf[f2];
					double dy   = yf[f1] - yf[f2];
					double dz   = zf[f1] - zf[f2];
					double r2   = dx*dx + dy*dy + dz*dz;
					double r    = sqrt(r2);
					if (r > 0) { 
						double fgerfc;
						double gerfc_tmp = gerfc(r, betaJ, &fgerfc);
						double gerf_tmp =  1.0 - gerfc_tmp;
						double r_inv = 1.0/r;
						EHarselfscrGrid += pref*gerf_tmp*r_inv;
					} // end if r>0
				} else {
					EHarselfscrGrid += pref*PRE_ERFC*betaJ;
				} // end if
			} //end for f2
		} //end for f1
	} //end for jtyp
	end = clock();
	PRINTF("0D EeN and EHar self finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

// register the energies to the struct
	EeNGrid      += EeNselfGrid;
	EHarscrGrid   = EHarGrid + EHarselfscrGrid;
	EHarGrid     += EHarselfGrid;
} // end if iperd = 0

//===================================================================
// energy contributions for 3D PDC
//===================================================================
 if (iperd == 3) { 

clock_t start, end;
//===================================================================
// Adding the ENN part force 
  start = clock();
  for (int J=0; J<natm; J++) {
      for (int K=0; K<natm; K++) {
          if (J != K) {
			  double dx = x[K]-x[J];
			  double dy = y[K]-y[J];
			  double dz = z[K]-z[J];
			  dx -= hmat[1]*NINT(dx*hmati[1]);
			  dy -= hmat[5]*NINT(dy*hmati[5]);
			  dz -= hmat[9]*NINT(dz*hmati[9]);
			  double r2 = dx*dx + dy*dy + dz*dz;
			  double R_ij = sqrt(r2);
			  double coeff1  = -q[J]*q[K];
              double fgerfc_bar;
              double gerfc_tmp_bar =  gerfc(R_ij, alpb, &fgerfc_bar);
              double R_ij_inv = 1.0/R_ij;
			  double tmp    = coeff1*R_ij_inv*R_ij_inv*R_ij_inv*(1.0 + erfc_a_r_over_r(R_ij, alpb,gerfc_tmp_bar,fgerfc_bar));
			  fxg[J]		+= tmp*dx; 
			  fyg[J]		+= tmp*dy; 
			  fzg[J]		+= tmp*dz; 
          } // end if
      } // end for 
  } // end for
  end = clock();
  PRINTF("3D ENN finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);
	

	int pair = 0;
	// short range eN and Hartree without self
	start = clock();
	for (int jtyp=0; jtyp<natm_typ; jtyp++) {
		double *xf_jtyp     = fgrid[jtyp].xf;
		double *yf_jtyp     = fgrid[jtyp].yf;
		double *zf_jtyp     = fgrid[jtyp].zf;
		double *wf_jtyp     = fgrid[jtyp].wf;
		//for(int ktyp=typ; ktyp<natm_typ; ktyp++){
		for(int ktyp=0; ktyp<natm_typ; ktyp++){
			double *xf_ktyp     = fgrid[ktyp].xf;
			double *yf_ktyp     = fgrid[ktyp].yf;
			double *zf_ktyp     = fgrid[ktyp].zf;
			double *wf_ktyp     = fgrid[ktyp].wf;
			for(int j=0; j<natm_atm_typ[jtyp];j++){
				int J = list_atm_by_typ[jtyp][j];
				double Ncoeff_J = four_pi_inv;
				for(int k=0;k<natm_atm_typ[ktyp];k++){
					int K = list_atm_by_typ[ktyp][k];
					double Ncoeff_K = four_pi_inv;
					if (K != J) {
                        double dx = x[K]-x[J];
                        double dy = y[K]-y[J];
                        double dz = z[K]-z[J];
                        dx -= hmat[1]*NINT(dx*hmati[1]);
                        dy -= hmat[5]*NINT(dy*hmati[5]);
                        dz -= hmat[9]*NINT(dz*hmati[9]);
                        double r2 = dx*dx + dy*dy + dz*dz;
                        double R_ij = sqrt(r2);
						if (R_ij < Rcut + Rpc[J] + Rpc[K]) { // allows Hartree and e-N to be calculated
							pair++;
							for (int f1=0; f1 < nf; f1++) {
								if (R_ij < Rcut + Rpc[J]) {  // alows e-N
									double dx_J = xf_jtyp[f1] - dx;
									double dy_J = yf_jtyp[f1] - dy;
									double dz_J = zf_jtyp[f1] - dz;
//									dx_J -= hmat[1]*NINT(dx_J*hmati[1]);
//									dy_J -= hmat[5]*NINT(dy_J*hmati[5]);
//									dz_J -= hmat[9]*NINT(dz_J*hmati[9]);
									double r2_J = dx_J*dx_J + dy_J*dy_J + dz_J*dz_J;
									double r_J = sqrt(r2_J);

									double fgerfc_bar_J;
									double gerfc_tmp_bar_J =  gerfc(r_J, alpb, &fgerfc_bar_J);
									double r_J_inv = 1.0/r_J;

									double tmp_eJ = -Ncoeff_J*qt[J]*q[K]*wf_jtyp[f1]*gerfc_tmp_bar_J*r_J_inv;
									double tmp_J = -Ncoeff_J*qt[J]*q[K]*wf_jtyp[f1]*r_J_inv;
									EeNshortGrid += tmp_eJ;

									double dx_K   = xf_ktyp[f1] + dx;
									double dy_K   = yf_ktyp[f1] + dy;
									double dz_K   = zf_ktyp[f1] + dz;
//									dx_K -= hmat[1]*NINT(dx_K*hmati[1]);
//									dy_K -= hmat[5]*NINT(dy_K*hmati[5]);
//									dz_K -= hmat[9]*NINT(dz_K*hmati[9]);
									double r2_K   = dx_K*dx_K + dy_K*dy_K + dz_K*dz_K;
									double r_K    = sqrt(r2_K);

									double fgerfc_bar_K;
									double gerfc_tmp_bar_K =  gerfc(r_K, alpb, &fgerfc_bar_K);
									double r_K_inv = 1.0/r_K;

									double tmp_K  = -Ncoeff_K*qt[K]*q[J]*wf_ktyp[f1]*r_K_inv;
									double coeffJ = tmp_J*r_J_inv*r_J_inv;
									double coeffK = tmp_K*r_K_inv*r_K_inv;

									double eNtmp_J = coeffJ*(1.0 + erfc_a_r_over_r(r_J, alpb,gerfc_tmp_bar_J,fgerfc_bar_J));
									double eNtmp_K = coeffK*(1.0 + erfc_a_r_over_r(r_K, alpb,gerfc_tmp_bar_K,fgerfc_bar_K));
									fxg[J]		 += (eNtmp_J*dx_J - eNtmp_K*dx_K);
									fyg[J]		 += (eNtmp_J*dy_J - eNtmp_K*dy_K);
									fzg[J]		 += (eNtmp_J*dz_J - eNtmp_K*dz_K);
								} // ends e-N
								for (int f2=0; f2 < nf; f2++) {
									double Hdx_J 	= xf_jtyp[f1]-xf_ktyp[f2] - dx;
									double Hdy_J 	= yf_jtyp[f1]-yf_ktyp[f2] - dy;
									double Hdz_J 	= zf_jtyp[f1]-zf_ktyp[f2] - dz;
//									Hdx_J -= hmat[1]*NINT(Hdx_J*hmati[1]);
//									Hdy_J -= hmat[5]*NINT(Hdy_J*hmati[5]);
//									Hdz_J -= hmat[9]*NINT(Hdz_J*hmati[9]);
									double Hr2_J 	= Hdx_J*Hdx_J + Hdy_J*Hdy_J + Hdz_J*Hdz_J;
									double Hr_J 	= sqrt(Hr2_J);

									double fgerfc_bar_HJ;
									double gerfc_tmp_bar_HJ =  gerfc(Hr_J, alpb, &fgerfc_bar_HJ);
									double Hr_J_inv = 1.0/Hr_J;

									double Htmp_eJ  = 0.5*Ncoeff_J*Ncoeff_K*qt[J]*qt[K]*wf_jtyp[f1]*wf_ktyp[f2]*gerfc_tmp_bar_HJ*Hr_J_inv;
									double Htmp_J   = 0.5*Ncoeff_J*Ncoeff_K*qt[J]*qt[K]*wf_jtyp[f1]*wf_ktyp[f2]*Hr_J_inv;
									EHarshortGrid  += Htmp_eJ;

									double Hdx_K 	= xf_ktyp[f1]-xf_jtyp[f2] + dx;
									double Hdy_K 	= yf_ktyp[f1]-yf_jtyp[f2] + dy;
									double Hdz_K 	= zf_ktyp[f1]-zf_jtyp[f2] + dz;
//									Hdx_K -= hmat[1]*NINT(Hdx_K*hmati[1]);
//									Hdy_K -= hmat[5]*NINT(Hdy_K*hmati[5]);
//									Hdz_K -= hmat[9]*NINT(Hdz_K*hmati[9]);
									double Hr2_K 	= Hdx_K*Hdx_K + Hdy_K*Hdy_K + Hdz_K*Hdz_K;
									double Hr_K 	= sqrt(Hr2_K);

									double fgerfc_bar_HK;
									double gerfc_tmp_bar_HK =  gerfc(Hr_K, alpb, &fgerfc_bar_HK);
									double Hr_K_inv = 1.0/Hr_K;

									double Htmp_K 	= 0.5*Ncoeff_K*Ncoeff_J*qt[J]*qt[K]*wf_ktyp[f1]*wf_jtyp[f2]*Hr_K_inv;
									double HcoeffJ  = Htmp_J*Hr_J_inv*Hr_J_inv;
									double HcoeffK  = Htmp_K*Hr_K_inv*Hr_K_inv;
									double Hartmp_J   = HcoeffJ*(1.0 + erfc_a_r_over_r(Hr_J, alpb, gerfc_tmp_bar_HJ, fgerfc_bar_HJ));
									double Hartmp_K   = HcoeffK*(1.0 + erfc_a_r_over_r(Hr_K, alpb, gerfc_tmp_bar_HK, fgerfc_bar_HK));
									fxg[J]		   += (Hartmp_J*Hdx_J - Hartmp_K*Hdx_K); 
									fyg[J]		   += (Hartmp_J*Hdy_J - Hartmp_K*Hdy_K); 
									fzg[J]		   += (Hartmp_J*Hdz_J - Hartmp_K*Hdz_K); 
								} // end for f2
							} // end for f1
						} // end allows Hartree and eN
					} // end if K != J
				} // end for k
			} //end for j
		} // end for ktyp
 	} // end for jtyp
	end = clock();
	PRINTF("(%d) pairs\n",pair);
	PRINTF("3D short EeN and Har noself finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);
//==============================================================================
//eN and Har self term (short)
	start = clock();
	for (int jtyp=0; jtyp<natm_typ; jtyp++) {
		double *xf		  = fgrid[jtyp].xf;
		double *yf		  = fgrid[jtyp].yf;
		double *zf		  = fgrid[jtyp].zf;
        double *wf        = fgrid[jtyp].wf;
        double *rf        = fgrid[jtyp].rf;
		double *xcostheta = fgrid[jtyp].xcostheta;
		double *xphi      = fgrid[jtyp].xphi;
		complex *Ylmf     = fgrid[jtyp].Ylmf;
		double alpJ			= fgrid[jtyp].alp;
		double betaJ	  	= fgrid[jtyp].beta;
		int J = list_atm_by_typ[jtyp][0];
		double Ncoeff = four_pi_inv;
		double wght = ((double) natm_atm_typ[jtyp]);

		for (int f1=0; f1 < nf; f1++) {
			double r = rf[f1];
			if (r > 0) { EeNshortselfGrid += -1.0*wght*Ncoeff*qt[J]*q[J]*wf[f1]*erfc(alpb*r)/r; }
		} // end for f1

		double erf_lim = alpb/sqrt(M_PI_QI);
		for (int f1=0; f1 < nf; f1++) {
			for (int f2=0; f2 < nf; f2++) {
				double pref = wght*Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
				double dx = xf[f1] - xf[f2];
                double dy = yf[f1] - yf[f2];
                double dz = zf[f1] - zf[f2];
                double r2 = dx*dx + dy*dy + dz*dz;
                double r = sqrt(r2);
				if (f1 != f2) {
					if (r > 0) { 
						double fgerfc;
						double gerfc_tmp = gerfc(r, alpb, &fgerfc);
						double gerf_tmp =  1.0 - gerfc_tmp;
						double r_inv = 1.0/r;
						EHarshortselfGrid -= 0.5*pref*gerf_tmp*r_inv; 
					}// end if r > 0
				} else {
					EHarshortselfGrid -= erf_lim*pref;
				} // end if
			} //end for f2
		} //end for f1

		for (int f1=0; f1 < nf; f1++) {
			for (int f2=0; f2 < nf; f2++) {
				double pref = 0.5*wght*Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
				if (f1 != f2) {
					double dx   = xf[f1] - xf[f2];
					double dy   = yf[f1] - yf[f2];
					double dz   = zf[f1] - zf[f2];
					double r2   = dx*dx + dy*dy + dz*dz;
					double r    = sqrt(r2);
					if (r > 0) { 
						double fgerfc_bar, fgerfc_beta;
						double gerfc_tmp_bar  = gerfc(r, alpb, &fgerfc_bar);
						double gerfc_tmp_beta = gerfc(r, betaJ, &fgerfc_beta);
						double r_inv = 1.0/r;
						EHarshortselfscrGrid += pref*(gerfc_tmp_bar - gerfc_tmp_beta)*r_inv; 
					}// end if r > 0
				} else {
					EHarshortselfscrGrid += pref*2.0*(betaJ - alpb)/sqrt(M_PI_QI);
				} // end if
			} //end for f2
		} //end for f1

		for (int l=0; l <= lmax; l++) {
			for (int m=-l; m<=l; m++) {
				gen_Ylmf (rorder, thetaorder, xcostheta, phiorder, xphi, l, m, Ylmf);
				for (int f1=0; f1 < nf; f1++) {
					for (int f2=0; f2 < nf; f2++) {
						double rmin = MIN(rf[f1],rf[f2]);
						double rmax = MAX(rf[f1],rf[f2]);
						double pref = Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
						complex temp = pref*(2.0*M_PI_QI)/(2*l+1)*pow(rmin,l)/pow(rmax,l+1)*(Ylmf[f1]*Ylmf[f2].conj());
						if (rmax > 0) {EHarshortselfGrid += temp.re*wght;}
					} //end for f2
				} //end for f1
			} //end for m
		} //end for l
	} //end for jtyp
	end = clock();
	PRINTF("3D eN and Har self finishes in (%lf seconds) \n",((double) end-start)/CLOCKS_PER_SEC);

	// register the energies in the energy struct
	EeNshortGrid      += EeNshortselfGrid;
	EHarshortscrGrid   = EHarshortGrid + EHarshortselfscrGrid;
	EHarshortGrid     += EHarshortselfGrid;
 } // end if iperd = 3D

//===================================================================
// Put the energies in the structure

	energy->ENN.EGrid            = ENNGrid;
	energy->ENNshort.EGrid       = ENNshortGrid;
	energy->ENNselflong.EGrid    = ENNselflongGrid;
	energy->EeN.EGrid            = EeNGrid; 
	energy->EeNself.EGrid        = EeNselfGrid; 
	energy->EHar.EGrid           = EHarGrid; 
	energy->EHarscr.EGrid        = EHarscrGrid; 
	energy->EHarself.EGrid       = EHarselfGrid; 
	energy->EeNshortself.EGrid   = EeNshortselfGrid; 
	energy->EeNshort.EGrid       = EeNshortGrid; 
	energy->EHarshortself.EGrid  = EHarshortselfGrid; 
	energy->EHarshort.EGrid      = EHarshortGrid; 
	energy->EHarshortscr.EGrid   = EHarshortscrGrid; 
	energy->ENNselflong.EGrid    = ENNselflongGrid; 
	energy->EHarselfscr.EGrid    = EHarselfscrGrid; 
	energy->EHarshortselfscr.EGrid    = EHarshortselfscrGrid; 

//===================================================================
} // end routine
//===================================================================

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the long range part energy (with and without f grid)
//===================================================================
void computePAWlong(ATOM_MAPS *atom_maps, ATOM_POS *atom_pos, CELL *cell, ESTRUCT *energy, FGRID *fgrid) 
//===================================================================
{ // begin routine
int iperd = cell->iperd;
double Elong = 0.0, ElongGrid = 0.0, Elongzero = 0.0;
if (iperd == 3) {
//===================================================================
// read in parameters from the structures
	int natm			  = atom_maps->natm;
	double hmati[10];
    for (int i=1; i<10; i++) hmati[i] = cell->hmati[i];
	double Gcut			  = cell->Gcut;
    double vol			  = cell->volume;
    int natm_typ          = atom_maps->natm_typ;
    int *natm_atm_typ     = atom_maps->natm_atm_typ;
    int **list_atm_by_typ = atom_maps->list_atm_by_typ;
  	double *x 	= atom_pos->x;   double *y   = atom_pos->y;   double *z   = atom_pos->z; 
  	double *fxg	= atom_pos->fxg; double *fyg = atom_pos->fyg; double *fzg = atom_pos->fzg;
  	double *fx	= atom_pos->fx;  double *fy	 = atom_pos->fy;  double *fz  = atom_pos->fz;
  	double *q 	= atom_pos->q;   double *qt  = atom_pos->qt;  double *alp = atom_pos->alp; 
    double alpb           = cell->alpb;
 
//===================================================================
// local variables
    int nf         = fgrid[0].nf;

//===================================================================
// zero the energies and compute the g-vectors

	int Ixgmax = (int) (Gcut/(2.0*M_PI_QI*hmati[1]));
	int Iygmax = (int) (Gcut/(2.0*M_PI_QI*hmati[5]));
	int Izgmax = (int) (Gcut/(2.0*M_PI_QI*hmati[9]));
	double xgmax = (double) Ixgmax;
	double ygmax = (double) Iygmax;
	double zgmax = (double) Izgmax;
	
//	PRINTF("xgmax = %d, ygmax = %d, zgmax = %d\n", xgmax, ygmax, zgmax);

//===================================================================
// compute the g=0 term 

	double delq = 0.0, qtalp = 0.0;
	for (int i=0; i<natm; i++) {
		delq += (q[i] - qt[i]);
		qtalp += qt[i]/(alp[i]*alp[i]);
	} // end for
	Elongzero = M_PI_QI*delq/vol*(1.0*qtalp - delq/(2.0*alpb*alpb)); // fix me for non-Gaussian

//===================================================================
// compute the g-space sum (g != 0)
	
	complex * ng_form_fact  = new complex [natm_typ];
	complex * ngf_form_fact = new complex [natm_typ];
	for (double gx = -xgmax; gx <= xgmax; gx++) {
		double ggx = gx*(2.0*M_PI_QI)*hmati[1];
		for (double gy = -ygmax; gy <= ygmax; gy++) {
			double ggy = gy*(2.0*M_PI_QI)*hmati[5];
			for (double gz = -zgmax; gz <= zgmax; gz++) {
				double ggz = gz*(2.0*M_PI_QI)*hmati[9];
				double g2 = ggx*ggx + ggy*ggy + ggz*ggz;
				complex ng = complex (0.0,0.0);
				complex ngf = complex (0.0,0.0);
				if (g2 <= Gcut*Gcut && g2 > 0.0) { 
					double kernel = 2.0*M_PI_QI/(g2*vol)*exp(-g2/(4.0*alpb*alpb));
					// compute the form factors, which only depends on atom type and g
					for (int jtyp=0; jtyp<natm_typ; jtyp++) {
						double *xf		     = fgrid[jtyp].xf;
						double *yf		     = fgrid[jtyp].yf;
						double *zf		     = fgrid[jtyp].zf;
				        double *wf           = fgrid[jtyp].wf;
 			            int K 				 = list_atm_by_typ[jtyp][0];
						double Ncoeff 		 = four_pi_inv;
						ng_form_fact[jtyp]   = complex(q[K] - qt[K]*exp(-g2/(4.0*alp[K]*alp[K])),0.0);
						ngf_form_fact[jtyp]	 = complex(q[K],0.0);
						double ngf_pre	     = qt[K]*Ncoeff;
						for (int f=0; f<nf; f++) {
							double gdotRf = -(ggx*xf[f] + ggy*yf[f] + ggz*zf[f]);
							ngf_form_fact[jtyp] -= (ngf_pre*wf[f])*CkExpIm(gdotRf);
						}// end for f
					} // end for jtyp
					// compute the atom type resolved density and add it to the total density with the correct form factor
					for (int jtyp=0; jtyp<natm_typ; jtyp++) {
						complex ng_atm_typ = complex(0.0,0.0);
						for(int j=0; j<natm_atm_typ[jtyp];j++){
							int J = list_atm_by_typ[jtyp][j];
							double gdotR = -(ggx*x[J] + ggy*y[J] + ggz*z[J]);
							ng_atm_typ   += CkExpIm(gdotR);
						} // end for j
						ng  += ng_form_fact[jtyp]*ng_atm_typ;
						ngf += ngf_form_fact[jtyp]*ng_atm_typ;
					} // end for jtyp
					// compute the energy using the kernel and the density
					double ng_mag_square = ng.getMagSqr();
					double ngf_mag_square = ngf.getMagSqr();
					Elong += kernel*ng_mag_square;
					ElongGrid += kernel*ngf_mag_square;
					// get the forces on the atoms
					for (int jtyp=0; jtyp<natm_typ; jtyp++) {
						for(int j=0; j<natm_atm_typ[jtyp]; j++) {
 			        		int J = list_atm_by_typ[jtyp][j];
							double  gdotR 		= -(ggx*x[J] + ggy*y[J] + ggz*z[J]);
							complex expgdotR	= CkExpIm(gdotR);
							complex div_ng      = ng_form_fact[jtyp]*expgdotR*complex(0.0,-1.0);
							complex div_ngf     = ngf_form_fact[jtyp]*expgdotR*complex(0.0,-1.0);
							complex fg_tmp_c	= ng*div_ng.conj();
							complex fgf_tmp_c	= ngf*div_ngf.conj();
  							double  fg_tmp		= -2.0*fg_tmp_c.re*kernel;
  							double  fgf_tmp		= -2.0*fgf_tmp_c.re*kernel;
							fx[J]			   += fg_tmp*ggx;
							fy[J]			   += fg_tmp*ggy;
							fz[J]			   += fg_tmp*ggz;
							fxg[J]			   += fgf_tmp*ggx;
							fyg[J]			   += fgf_tmp*ggy;
							fzg[J]			   += fgf_tmp*ggz;
						} // end for j
					}// end for jtyp
				} // end if
			} // end for gz
     	} // end for gy
	} // end for gx	
} // end if iperd = 3
//===================================================================
// put the energies in the structure 
	Elong      			   += Elongzero;
	ElongGrid 			   += Elongzero;
	energy->Elong.E         = Elong;
	energy->Elong.EGrid     = ElongGrid;
//===================================================================
} // end routine
//===================================================================

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the distance between two points
//===================================================================
inline double dist(double dx, double dy, double dz) 
//===================================================================
{ // begin routine
//===================================================================
	double result = sqrt(dx*dx + dy*dy + dz*dz);
	return result;
//===================================================================
} // end routine
//===================================================================

//===================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// Compute the distance between two points
//===================================================================
inline double erfc_a_r_over_r(double r, double a, double gerfc, double fgerfc) 
//===================================================================
{ // begin routine
//===================================================================
//	double result = 2.0*a*r/sqrt(M_PI_QI)*exp(-a*a*r*r) - erf(a*r);
	double result = r*fgerfc + gerfc - 1.0;
	return result;
//===================================================================
} // end routine
//===================================================================

//===========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===========================================================================
// Generic error function
//===========================================================================
inline double gerfc(double r, double a, double *fgerfc){
//===========================================================================
#define _STRICT_ERFC_DER_OFF_
	double x        = a*r;
	double eee      = exp(-x*x);
	double tt       = 1.0/(1.0+PERFC*x);
	double gerfc_x  = ((((((((CERFC9*tt+CERFC8)*tt+CERFC7)*tt+CERFC6)*tt+CERFC5)*tt+CERFC4)*tt+CERFC3)*tt+CERFC2)*tt+CERFC1)*tt*eee;
#ifdef _STRICT_ERFC_DER_
    double fgerfc_x = ((((((((DCERFC9*tt+DCERFC8)*tt+DCERFC7)*tt+DCERFC6)*tt+DCERFC5)*tt+DCERFC4)*tt+DCERFC3)*tt+DCERFC2)*tt+DCERFC1)*tt*tt*eee*PERFC
                         +2.0*gerfc_x*x;
#else
    fgerfc[0]       = a*eee*PRE_ERFC;   // -d/dr ( erfc(a*r)  )
#endif
//===========================================================================
    return gerfc_x;
}//end routine gerf
//===========================================================================
