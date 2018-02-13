
//==========================================================================
// POLY structure
typedef struct POLY{
        int 	order;		// order
        double *c;          // c matrix 
        double *a;          // a matrix
        double *O;          // O integral 
		double Norm;		// N_ii
}POLY;

//==========================================================================
// function prototypes 
void gen_Gauss_quad(int, int, double, double *, double *, int);

void fetch_Int_xpow_uniform(int, double *); 
void fetch_Int_xpow_Gaussfull(int, double *); 
void fetch_Int_xpow_Gausshalf(int, double *);
void horner(int, double *, double, double *);
void testgrid(int, double *, double *, double *);
extern "C" void dggev( char* JOBVL,  char* JOBVR,  int* N,
                       double* AA,  int* LDA,  double* B,  int* LDB,
                      double* ALPHAR, double* ALPHAI, double* BETA,
                      double* VL,  int* LDVL, double* VR,  int* LDVR,
                      double* WORK,  int* LWORK, int* INFO);
