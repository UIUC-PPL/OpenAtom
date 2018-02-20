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
	#include "GaussianPAW.h"

//=================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================
//  Generate Uniform Random Numbers
//=================================================================
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
  int natm 		= atom_maps->natm;
  double *x 	= atom_pos->x;   double *y   = atom_pos->y;   double *z   = atom_pos->z; 
  double *fx0	= atom_pos->fx0; double *fy0 = atom_pos->fy0; double *fz0 = atom_pos->fz0;
  double *fx	= atom_pos->fx;  double *fy	 = atom_pos->fy;  double *fz  = atom_pos->fz;
  double *q 	= atom_pos->q;   double *qt  = atom_pos->qt;  double *alp = atom_pos->alp; 
  double alpb   = cell->alpb;
  int iperd 	= cell->iperd;

//===================================================================
// zero all the energies 
  double ENN = 0.0, EeNself = 0.0, EeN = 0.0, EHarself = 0.0, EHar = 0.0;

//===================================================================
// energy contributions for clusters
//===================================================================
// if (iperd == 0) { 
  for (int i=0; i<natm; i++) {
	EHarself      += 0.5*alp[i]*qt[i]*(sqrt(2.0)*qt[i])/sqrt(M_PI_QI);
	EeNself       += 0.5*alp[i]*qt[i]*(-4.0*q[i])/sqrt(M_PI_QI);
  } // end for

  for (int i=0; i<natm; i++) {
      for (int j=0; j<natm; j++) {
          if (j != i) {
              double alp_ij  = alp[i]*alp[j]/sqrt(alp[i]*alp[i] + alp[j]*alp[j]);
              double R_ij    = dist((x[i]-x[j]), (y[i]-y[j]), (z[i]-z[j]));
              ENN           += 0.5*q[i]*q[j]/R_ij;
              EeN           -= 1.0*qt[i]*q[j]*erf(alp[i]*R_ij)/R_ij;
              EHar          += 0.5*qt[i]*qt[j]*erf(alp_ij*R_ij)/R_ij;
			  double coeff1  = -q[i]*q[j];
			  double coeff2  = -(qt[i]*q[j]*erfc_a_r_over_r(R_ij,alp[i])+ q[i]*qt[j]*erfc_a_r_over_r(R_ij,alp[i]));
			  double coeff3  = qt[i]*qt[j]*erfc_a_r_over_r(R_ij,alp_ij);
			  fx0[i]		  += (coeff1 + coeff2 + coeff3)/(R_ij*R_ij*R_ij)*(x[j] - x[i]); 
			  fy0[i]		  += (coeff1 + coeff2 + coeff3)/(R_ij*R_ij*R_ij)*(y[j] - y[i]); 
			  fz0[i]		  += (coeff1 + coeff2 + coeff3)/(R_ij*R_ij*R_ij)*(z[j] - z[i]); 
//			  fx0[i]		  += (coeff3)/(R_ij*R_ij*R_ij)*(x[j] - x[i]); 
//			  fy0[i]		  += (coeff3)/(R_ij*R_ij*R_ij)*(y[j] - y[i]); 
//			  fz0[i]		  += (coeff3)/(R_ij*R_ij*R_ij)*(z[j] - z[i]); 
          } // end if
      } // end for 
  } // end for
	EeN  += EeNself;
	EHar += EHarself;
// } // end if

//===================================================================
// zero all the energies 
  double EeNshortself = 0.0, EeNshort = 0.0, EHarshortself = 0.0, EHarshort = 0.0;
  double ENNshort = 0.0;
  double ENNselflong = 0.0;

//===================================================================
// energy contributions for 3D PDC
//===================================================================

// if (iperd == 3) { 
  for (int i=0; i<natm; i++) {
	double alpb_i  = alpb*alp[i]/sqrt(alp[i]*alp[i] + alpb*alpb);
	double alpb_ii = alpb*alp[i]*alp[i]/sqrt(alp[i]*alp[i]*alp[i]*alp[i]+ 2*alpb*alpb*alp[i]*alp[i]);
	EeNshortself  += 1.0/sqrt(M_PI_QI)*(-2.0*(alp[i] - alpb_i))*qt[i]*q[i];
	EHarshortself += 1.0/sqrt(M_PI_QI)*(alp[i]/sqrt(2.0) - alpb_ii)*qt[i]*qt[i];
	ENNselflong   += 1.0*alpb/sqrt(M_PI_QI)*q[i]*q[i];
  } // end for

  for (int i=0; i<natm; i++) {
    double alpb_i = alpb*alp[i]/sqrt(alp[i]*alp[i] + alpb*alpb);
    for (int j=0; j<natm; j++) {
        if (j != i) {
    		double alpb_j = alpb*alp[j]/sqrt(alp[j]*alp[j] + alpb*alpb);
            double alp_ij  = alp[i]*alp[j]/sqrt(alp[i]*alp[i] + alp[j]*alp[j]);
            double alpb_ij = alpb*alp[i]*alp[j]/sqrt(alp[i]*alp[i]*alp[j]*alp[j] + alpb*alpb*alp[i]*alp[i] + alpb*alpb*alp[j]*alp[j]);
            double R_ij    = dist((x[i]-x[j]), (y[i]-y[j]), (z[i]-z[j]));
            ENNshort      += 0.5*q[i]*q[j]*erfc(alpb*R_ij)/R_ij;
            EeNshort      -= 1.0*qt[i]*q[j]*(erfc(alpb_i*R_ij)-erfc(alp[i]*R_ij))/R_ij;
            EHarshort     += 0.5*qt[i]*qt[j]*(erfc(alpb_ij*R_ij)-erfc(alp_ij*R_ij))/R_ij;
			double coeff1  = -q[i]*q[j]*(1 + erfc_a_r_over_r(R_ij,alpb));
			double coeff2  = -(qt[i]*q[j]*(erfc_a_r_over_r(R_ij,alp[i]) - erfc_a_r_over_r(R_ij,alpb_i))+ q[i]*qt[j]*(erfc_a_r_over_r(R_ij,alp[j]) - erfc_a_r_over_r(R_ij,alpb_j)));
			double coeff3  = qt[i]*qt[j]*(erfc_a_r_over_r(R_ij,alp_ij) - erfc_a_r_over_r(R_ij,alpb_ij));
			fx[i]		  += (coeff1 + coeff2 + coeff3)/(R_ij*R_ij*R_ij)*(x[j] - x[i]); 
			fy[i]		  += (coeff1 + coeff2 + coeff3)/(R_ij*R_ij*R_ij)*(y[j] - y[i]); 
			fz[i]		  += (coeff1 + coeff2 + coeff3)/(R_ij*R_ij*R_ij)*(z[j] - z[i]); 
//			fx[i]		  += (coeff1)/(R_ij*R_ij*R_ij)*(x[j] - x[i]); 
//			fy[i]		  += (coeff1)/(R_ij*R_ij*R_ij)*(y[j] - y[i]); 
//			fz[i]		  += (coeff1)/(R_ij*R_ij*R_ij)*(z[j] - z[i]); 
     	} // end if
    } // end for 
  } // end for
    EeNshort  += EeNshortself;
    EHarshort += EHarshortself;
// } // end if


//===================================================================
// Put the energies in the structure

  energy->ENN.E                = ENN; 
  energy->ENNshort.E           = ENNshort; 
  energy->EeNself.E            = EeNself; 
  energy->EeN.E                = EeN; 
  energy->EHarself.E           = EHarself; 
  energy->EHar.E               = EHar; 
  energy->EeNshortself.E       = EeNshortself; 
  energy->EeNshort.E           = EeNshort; 
  energy->EHarshortself.E      = EHarshortself; 
  energy->EHarshort.E          = EHarshort; 
  energy->ENNselflong.E        = ENNselflong;

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
  int natm				= atom_maps->natm;
  int natm_typ 			= atom_maps->natm_typ;
  int *natm_atm_typ 	= atom_maps->natm_atm_typ;
  int **list_atm_by_typ = atom_maps->list_atm_by_typ;
  double *x 			= atom_pos->x;    double *y    = atom_pos->y;    double *z    = atom_pos->z; 
  double *fx0g			= atom_pos->fx0g; double *fy0g = atom_pos->fy0g; double *fz0g = atom_pos->fz0g;
  double *fxg			= atom_pos->fxg;  double *fyg  = atom_pos->fyg;  double *fzg  = atom_pos->fzg;
  double *q 			= atom_pos->q;    double *qt   = atom_pos->qt;   double *alp  = atom_pos->alp; 
  double alpb   		= cell->alpb;
  int iperd 			= cell->iperd;

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
    double EeNGrid      = 0.0, EeNselfGrid = 0.0, EHarGrid = 0.0, EHarselfGrid = 0.0;
    double EeNshortGrid = 0.0, EeNshortselfGrid = 0.0, EHarshortGrid = 0.0, EHarshortselfGrid = 0.0;

//===================================================================
// energy contributions for clusters
//===================================================================
// if (iperd == 0) { 

//===================================================================
// Adding the ENN part force 
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
// } // end if
	
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
				double Ncoeff_J = pow(alp[J]*alp[J]/(M_PI_QI),1.5);
				//int kstart = (ktyp == jtyp ? (j+1):0);
				//for(int k=kstart;k<natm_atm_typ[ktyp];k++){
				for(int k=0;k<natm_atm_typ[ktyp];k++){
					int K = list_atm_by_typ[ktyp][k];
					if (K != J) {
						double Ncoeff_K = pow(alp[K]*alp[K]/(M_PI_QI),1.5);
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
//==============================================================================
//eN and Har self term
	for (int jtyp=0; jtyp<natm_typ; jtyp++) {
        double *wf        = fgrid[jtyp].wf;
        double *rf        = fgrid[jtyp].rf;
		double *xcostheta = fgrid[jtyp].xcostheta;
		double *xphi      = fgrid[jtyp].xphi;
		complex *Ylmf     = fgrid[jtyp].Ylmf;
		int J = list_atm_by_typ[jtyp][0];
		double Ncoeff = pow(alp[J]*alp[J]/(M_PI_QI),1.5);
		double wght = ((double) natm_atm_typ[jtyp]);

		for (int f1=0; f1 < nf; f1++) {
			double r = rf[f1];
			EeNselfGrid += -1.0*wght*Ncoeff*qt[J]*q[J]*wf[f1]/r;
		} // end for f1

		for (int l=0; l <= lmax; l++) {
			for (int m=-l; m<=l; m++) {
				gen_Ylmf (rorder, thetaorder, xcostheta, phiorder, xphi, l, m, Ylmf);
				for (int f1=0; f1 < nf; f1++) {
					for (int f2=0; f2 < nf; f2++) {
						double rmin = MIN(rf[f1],rf[f2]);
						double rmax = MAX(rf[f1],rf[f2]);
						double pref = Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
						complex temp = pref*(2.0*M_PI_QI)/(2*l+1)*pow(rmin,l)/pow(rmax,l+1)*(Ylmf[f1]*Ylmf[f2].conj());
						EHarselfGrid += temp.re*wght;
					} //end for f2
				} //end for f1
			} //end for m
		} //end for l
	} //end for jtyp
	EeNGrid      += EeNselfGrid;
	EHarGrid     += EHarselfGrid;
// } // end if

//===================================================================
// energy contributions for 3D PDC
//===================================================================
// if (iperd == 3) { 

//===================================================================
// Adding the ENN part force 
  for (int i=0; i<natm; i++) {
      for (int j=0; j<natm; j++) {
          if (j != i) {
              double R_ij    = dist((x[i]-x[j]), (y[i]-y[j]), (z[i]-z[j]));
			  double coeff1  = -q[i]*q[j];
			  fxg[i]		+= (coeff1)/(R_ij*R_ij*R_ij)*(x[j] - x[i])*(1 + erfc_a_r_over_r(R_ij, alpb)); 
			  fyg[i]		+= (coeff1)/(R_ij*R_ij*R_ij)*(y[j] - y[i])*(1 + erfc_a_r_over_r(R_ij, alpb)); 
			  fzg[i]		+= (coeff1)/(R_ij*R_ij*R_ij)*(z[j] - z[i])*(1 + erfc_a_r_over_r(R_ij, alpb)); 
          } // end if
      } // end for 
  } // end for
// } // end if
	
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
				double Ncoeff_J = pow(alp[J]*alp[J]/(M_PI_QI),1.5);
				//int kstart = (ktyp == jtyp ? (j+1):0);
				//for(int k=kstart;k<natm_atm_typ[ktyp];k++){
				for(int k=0;k<natm_atm_typ[ktyp];k++){
					int K = list_atm_by_typ[ktyp][k];
					double Ncoeff_K = pow(alp[K]*alp[K]/(M_PI_QI),1.5);
					if (K != J) {
						for (int f1=0; f1 < nf; f1++) {
							double dx_J = xf_jtyp[f1]-x[K]+x[J];
							double dy_J = yf_jtyp[f1]-y[K]+y[J];
							double dz_J = zf_jtyp[f1]-z[K]+z[J];
							double r2_J = dx_J*dx_J + dy_J*dy_J + dz_J*dz_J;
							double r_J = sqrt(r2_J);
							double tmp_eJ = -Ncoeff_J*qt[J]*q[K]*wf_jtyp[f1]*erfc(alpb*r_J)/r_J;
							double tmp_J = -Ncoeff_J*qt[J]*q[K]*wf_jtyp[f1]/r_J;
							EeNshortGrid += tmp_eJ;
							double dx_K   = xf_ktyp[f1]-x[J]+x[K];
							double dy_K   = yf_ktyp[f1]-y[J]+y[K];
							double dz_K   = zf_ktyp[f1]-z[J]+z[K];
							double r2_K   = dx_K*dx_K + dy_K*dy_K + dz_K*dz_K;
							double r_K    = sqrt(r2_K);
							double tmp_K  = -Ncoeff_K*qt[K]*q[J]*wf_ktyp[f1]/r_K;
							double coeffJ = tmp_J/(r_J*r_J);
							double coeffK = tmp_K/(r_K*r_K);
							fxg[J]		 += coeffJ*dx_J*(1 + erfc_a_r_over_r(r_J, alpb)) - coeffK*dx_K*(1 + erfc_a_r_over_r(r_K, alpb));
							fyg[J]		 += coeffJ*dy_J*(1 + erfc_a_r_over_r(r_J, alpb)) - coeffK*dy_K*(1 + erfc_a_r_over_r(r_K, alpb));
							fzg[J]		 += coeffJ*dz_J*(1 + erfc_a_r_over_r(r_J, alpb)) - coeffK*dz_K*(1 + erfc_a_r_over_r(r_K, alpb));
							for (int f2=0; f2 < nf; f2++) {
								double Hdx_J 	= xf_jtyp[f1]-xf_ktyp[f2]-x[K]+x[J];
								double Hdy_J 	= yf_jtyp[f1]-yf_ktyp[f2]-y[K]+y[J];
								double Hdz_J 	= zf_jtyp[f1]-zf_ktyp[f2]-z[K]+z[J];
								double Hr2_J 	= Hdx_J*Hdx_J + Hdy_J*Hdy_J + Hdz_J*Hdz_J;
								double Hr_J 	= sqrt(Hr2_J);
								double Htmp_eJ  = 0.5*Ncoeff_J*Ncoeff_K*qt[J]*qt[K]*wf_jtyp[f1]*wf_ktyp[f2]*erfc(alpb*Hr_J)/Hr_J;
								double Htmp_J   = 0.5*Ncoeff_J*Ncoeff_K*qt[J]*qt[K]*wf_jtyp[f1]*wf_ktyp[f2]/Hr_J;
								EHarshortGrid  += Htmp_eJ;
								double Hdx_K 	= xf_ktyp[f1]-xf_jtyp[f2]+x[K]-x[J];
								double Hdy_K 	= yf_ktyp[f1]-yf_jtyp[f2]+y[K]-y[J];
								double Hdz_K 	= zf_ktyp[f1]-zf_jtyp[f2]+z[K]-z[J];
								double Hr2_K 	= Hdx_K*Hdx_K + Hdy_K*Hdy_K + Hdz_K*Hdz_K;
								double Hr_K 	= sqrt(Hr2_K);
								double Htmp_K 	= 0.5*Ncoeff_K*Ncoeff_J*qt[J]*qt[K]*wf_ktyp[f1]*wf_jtyp[f2]/Hr_K;
								double HcoeffJ  = Htmp_J/(Hr_J*Hr_J);
								double HcoeffK  = Htmp_K/(Hr_K*Hr_K);
								fxg[J]		   += HcoeffJ*Hdx_J*(1 + erfc_a_r_over_r(Hr_J, alpb)) - HcoeffK*Hdx_K*(1 + erfc_a_r_over_r(Hr_K, alpb)); 
								fyg[J]		   += HcoeffJ*Hdy_J*(1 + erfc_a_r_over_r(Hr_J, alpb)) - HcoeffK*Hdy_K*(1 + erfc_a_r_over_r(Hr_K, alpb)); 
								fzg[J]		   += HcoeffJ*Hdz_J*(1 + erfc_a_r_over_r(Hr_J, alpb)) - HcoeffK*Hdz_K*(1 + erfc_a_r_over_r(Hr_K, alpb)); 
							} // end for f2
						} // end for f1
					} // end if
				} // end for k
			} //end for j
		} // end for ktyp
 	} // end for jtyp
//==============================================================================
//eN and Har self term (short)
	for (int jtyp=0; jtyp<natm_typ; jtyp++) {
		double *xf		  = fgrid[jtyp].xf;
		double *yf		  = fgrid[jtyp].yf;
		double *zf		  = fgrid[jtyp].zf;
        double *wf        = fgrid[jtyp].wf;
        double *rf        = fgrid[jtyp].rf;
		double *xcostheta = fgrid[jtyp].xcostheta;
		double *xphi      = fgrid[jtyp].xphi;
		complex *Ylmf     = fgrid[jtyp].Ylmf;
		int J = list_atm_by_typ[jtyp][0];
		double Ncoeff = pow(alp[J]*alp[J]/(M_PI_QI),1.5);
		double wght = ((double) natm_atm_typ[jtyp]);

		for (int f1=0; f1 < nf; f1++) {
			double r = rf[f1];
			EeNshortselfGrid += -1.0*wght*Ncoeff*qt[J]*q[J]*wf[f1]*erfc(alpb*r)/r;
		} // end for f1

		double erf_lim = alpb/sqrt(M_PI_QI);
		for (int f1=0; f1 < nf; f1++) {
			for (int f2=0; f2 < nf; f2++) {
				double pref = Ncoeff*Ncoeff*qt[J]*qt[J]*wf[f1]*wf[f2];
				double dx = xf[f1] - xf[f2];
                double dy = yf[f1] - yf[f2];
                double dz = zf[f1] - zf[f2];
                double r2 = dx*dx + dy*dy + dz*dz;
                double r = sqrt(r2);
				if (f1 != f2) {
					EHarshortselfGrid -= 0.5*pref*erf(alpb*r)/r;
				} else {
					EHarshortselfGrid -= erf_lim*pref;
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
//						double dx = xf[f1] - xf[f2];
//                      double dy = yf[f1] - yf[f2];
//                      double dz = zf[f1] - zf[f2];
//                      double r2 = dx*dx + dy*dy + dz*dz;
						EHarshortselfGrid += temp.re*wght;
					} //end for f2
				} //end for f1
			} //end for m
		} //end for l
	} //end for jtyp
	EeNshortGrid      += EeNshortselfGrid;
	EHarshortGrid     += EHarshortselfGrid;
// } // end if

//===================================================================
// Put the energies in the structure

  energy->ENN.EGrid            = ENNGrid;
  energy->ENNshort.EGrid       = ENNshortGrid;
  energy->ENNselflong.EGrid    = ENNselflongGrid;
  energy->EeN.EGrid            = EeNGrid; 
  energy->EeNself.EGrid        = EeNselfGrid; 
  energy->EHar.EGrid           = EHarGrid; 
  energy->EHarself.EGrid       = EHarselfGrid; 
  energy->EeNshortself.EGrid   = EeNshortselfGrid; 
  energy->EeNshort.EGrid       = EeNshortGrid; 
  energy->EHarshortself.EGrid  = EHarshortselfGrid; 
  energy->EHarshort.EGrid      = EHarshortGrid; 
  energy->ENNselflong.EGrid    = ENNselflongGrid; 

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
//===================================================================
// read in parameters from the structures
	int natm			  = atom_maps->natm;
	double hmati[10];
    for (int i=1; i<10; i++) hmati[i] = cell->hmati[i];
	double gcut			  = cell->gcut;
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

	double Elong = 0.0, ElongGrid = 0.0;
	int Ixgmax = (int) (gcut/(2.0*M_PI_QI*hmati[1]));
	int Iygmax = (int) (gcut/(2.0*M_PI_QI*hmati[5]));
	int Izgmax = (int) (gcut/(2.0*M_PI_QI*hmati[9]));
	double xgmax = (double) Ixgmax;
	double ygmax = (double) Iygmax;
	double zgmax = (double) Izgmax;
	
//	PRINTF("xgmax = %d, ygmax = %d, zgmax = %d\n", xgmax, ygmax, zgmax);

//===================================================================
// compute the g=0 term 

	double Elongzero, delq = 0.0, qtalp = 0.0;
	for (int i=0; i<natm; i++) {
		delq += (q[i] - qt[i]);
		qtalp += qt[i]/(alp[i]*alp[i]);
	} // end for
	Elongzero = M_PI_QI*delq/vol*(1.0*qtalp - delq/(2.0*alpb*alpb));

//===================================================================
// compute the g-space sum (g != 0)

	for (double gx = -xgmax; gx <= xgmax; gx++) {
		double ggx = gx*(2.0*M_PI_QI)*hmati[1];
		for (double gy = -ygmax; gy <= ygmax; gy++) {
			double ggy = gy*(2.0*M_PI_QI)*hmati[5];
			for (double gz = -zgmax; gz <= zgmax; gz++) {
				double ggz = gz*(2.0*M_PI_QI)*hmati[9];
				double g2 = ggx*ggx + ggy*ggy + ggz*ggz;
				complex ng = complex (0.0,0.0);
				complex ngf = complex (0.0,0.0);
				if (g2 <= gcut*gcut && g2 > 0.0) { 
					double prefact = 2.0*M_PI_QI/(g2*vol)*exp(-g2/(4.0*alpb*alpb));
					for (int jtyp=0; jtyp<natm_typ; jtyp++) {
						double *xf		    = fgrid[jtyp].xf;
						double *yf		    = fgrid[jtyp].yf;
						double *zf		    = fgrid[jtyp].zf;
				        double *wf          = fgrid[jtyp].wf;
						complex ng_atm_typ  = complex(0.0,0.0);
						complex ngf_atm_typ1 = complex(0.0,0.0);
						complex ngf_atm_typ2 = complex(0.0,0.0);
 			            int K 				= list_atm_by_typ[jtyp][0];
						double Ncoeff 		= pow(alp[K]*alp[K]/(M_PI_QI),1.5);
						double ng_pref     = q[K] - qt[K]*exp(-g2/(4.0*alp[K]*alp[K]));
						double ngf_pref1	= q[K];
						double ngf_pref2	= qt[K]*Ncoeff;
						for(int j=0; j<natm_atm_typ[jtyp];j++){
 			                int J = list_atm_by_typ[jtyp][j];
							double gdotR = -(ggx*x[J] + ggy*y[J] + ggz*z[J]);
     						ng_atm_typ   += CkExpIm(gdotR);
     						ngf_atm_typ1 += CkExpIm(gdotR);
							for (int f=0; f<nf; f++) {
								double gdotRf = -(ggx*(x[J] + xf[f]) + ggy*(y[J] + yf[f]) + ggz*(z[J] + zf[f]));
								ngf_atm_typ2 -= wf[f]*CkExpIm(gdotRf);
							}// end for f
						} // end for j
						ng  += ng_pref*ng_atm_typ;
						ngf += ngf_pref1*ngf_atm_typ1 + ngf_pref2*ngf_atm_typ2;
					} // end for jtyp
					double normng = ng.getMagSqr();
					double normngf = ngf.getMagSqr();
					Elong += prefact*normng;
					ElongGrid += prefact*normngf;
					for (int jtyp=0; jtyp<natm_typ; jtyp++) {
						double *xf		    = fgrid[jtyp].xf;
						double *yf		    = fgrid[jtyp].yf;
						double *zf		    = fgrid[jtyp].zf;
				        double *wf          = fgrid[jtyp].wf;
 			            int K 				= list_atm_by_typ[jtyp][0];
						double Ncoeff 		= pow(alp[K]*alp[K]/(M_PI_QI),1.5);
						for(int j=0; j<natm_atm_typ[jtyp];j++){
 			        		int J = list_atm_by_typ[jtyp][j];
							double  gdotR 		= -(ggx*x[J] + ggy*y[J] + ggz*z[J]);
							double  ng_pref     = q[J] - qt[J]*exp(-g2/(4.0*alp[J]*alp[J]));
							complex ngf_pref    = complex (q[J],0.0);
							double  fg_tmp		= 2.0*(ng*CkExpIm(-gdotR)).im;
							fx[J]			   += fg_tmp*ng_pref*ggx*prefact;
							fy[J]			   += fg_tmp*ng_pref*ggy*prefact;
							fz[J]			   += fg_tmp*ng_pref*ggz*prefact;
							for (int f=0; f<nf; f++) {
								double gdotf    = ggx*xf[f] + ggy*yf[f] + ggz*zf[f];
								ngf_pref	   -= qt[J]*wf[f]*CkExpIm(gdotf)*Ncoeff;
							}// end for f
							double  fgf_tmp		= 2.0*(ng*CkExpIm(-gdotR)*ngf_pref).im;
							fxg[J]			   += fgf_tmp*ggx*prefact;
							fyg[J]			   += fgf_tmp*ggy*prefact;
							fzg[J]			   += fgf_tmp*ggz*prefact;
						} // end for j
					}// end for jtyp
				} // end if
			} // end for gz
     	} // end for gy
	} // end for gx
	

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
inline double erfc_a_r_over_r(double r, double a) 
//===================================================================
{ // begin routine
//===================================================================
	double result = 2.0*a*r/sqrt(M_PI_QI)*exp(-a*a*r*r) - erf(a*r);
	return result;
//===================================================================
} // end routine
//===================================================================
