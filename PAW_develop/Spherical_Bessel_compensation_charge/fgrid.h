//==========================================================================
// FGRID structure: spherical polar with DVR
typedef struct FGRID{
    int nf;             // total number of f points
    int nr;             // number of r grid points
    int ntheta;         // number of theta grid points
    int nphi;           // number of phi grid points
    int nrfull;         // number of r grid points*2
    double alp;         // Gaussian parameter
    double beta;        // screen parameter
    double *wf;         // weight including the Jacobian - size[nf]
    double *xf;         // x coordinate - size[nf]
    double *yf;         // y - size[nf]
    double *zf;         // z - size[nf]
    double *rf;         // r - size[nf]
    complex *Ylmf;      // memory to store 1 spehrical harmonic on the f grid - size[nf]
    double *xcostheta;  // xcostheta - used to generate Ylmf - size[ntheta]
    double *xphi;       // xphi - used to generate Ylmf - size[nphi]
    double *xr;         // this is the r grid only - size[nr]
    double *wr;         // this is the r grid only - size[nr]
}FGRID;
