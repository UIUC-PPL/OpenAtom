//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  The function prototypes for the grid.C
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

// Function prototypes
void readtoendofline(FILE *);
void SphericalHarmonicOnGrid(int, int, int, double *, int, double *, complex **);
void testHermite(int, double *, double *);
void testLegendre(int, double *, double *);
void genphigrid(int, double *, double *);
void testphi(int, double *, double *);
void calcHermiteHalfSpace(int, int, double *, double *, double *, double *);
void control_quad_rule(int, int, double, double, double *, double *);
void gen_fgrid (int, int, int, int, double, FGRID *);
//==========================================================================
