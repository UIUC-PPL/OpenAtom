//==========================================================================
// FGRID structure: spherical polar with DVR
typedef struct FGRID{
	int nf;             // total number of f points
	int nr;             // number of r grid points
	int ntheta;         // number of theta grid points
	int nphi;           // number of phi grid points
	int nrfull;         // number of r grid points*2
	int nang;					  // number of angular grid points = ntheta*nphi
	double alp;         // Gaussian parameter
	double beta;        // screen parameter
	double Rpc;					// Rpc
	double *wf;         // weight including the Jacobian - size[nf]
	double *xf;         // x coordinate - size[nf]
	double *yf;         // y - size[nf]
	double *zf;         // z - size[nf]
	double *rf;         // r - size[nf]
	double *rho;				// density w/o r or angular weights - size[nf]
	complex *Ylmf;      // memory to store 1 spehrical harmonic on the f grid - size[nf]
	double *ylm;        // memory to store 1 spehrical harmonic on the angular grid - size[nang]
	double *wang;				// weight for the angular integral - size[nang]
	double *rho_lm;			// the Ylm component of the density - size[nr]
	double *xcostheta;  // xcostheta - used to generate Ylmf - size[ntheta]
	double *xphi;       // xphi - used to generate Ylmf - size[nphi]
	double *xr;         // this is the r grid only - size[nr]
	double *wr;         // this is the r grid only - size[nr]
	double *wr_bare;
	double **rho_scr;	  // the scratch storage for the density - size[nr][nang]
	double **pw_erfB;		// the s-wave partial wave of erf(B|r-r'|)/|r-r'| - size[nr][nr]
	double **pw_erfA;		// the s-wave partial wave of erf(A|r-r'|)/|r-r'| - size[nr][nr]
	double **pw_coul;		// the s-wave partial wave of r<^l/r>^(l+1) Coulumb - size[nr][nr]
}FGRID;
