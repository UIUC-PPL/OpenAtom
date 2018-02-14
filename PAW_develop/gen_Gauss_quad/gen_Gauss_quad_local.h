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

void horner(int, double *, double, double *);
void testgrid(int, double *, double *, double *, int);
void testnodes(int, double * , double *);
extern "C" void dggev( char* JOBVL,  char* JOBVR,  int* N,
                       double* AA,  int* LDA,  double* B,  int* LDB,
                      double* ALPHAR, double* ALPHAI, double* BETA,
                      double* VL,  int* LDVL, double* VR,  int* LDVR,
                      double* WORK,  int* LWORK, int* INFO);

