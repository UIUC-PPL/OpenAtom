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
        int 	 order;		 // order
		double   Norm;		 // N_ii
        double * c;          // c matrix (order+1)
        double * a;          // a matrix (order+1)
        double * O;          // O integral (order+1)
		double * p_at_nodes; // size n
}POLY;

//==========================================================================
// function prototypes 

void horner(int, double *, double, double *);
void testgrid(int, double *, POLY *, int);
void testnodes(int, double * , double *,int);

#ifdef FORT_UNDER
#define DGGEV dggev_
#else
#define DGGEV dggev
#endif

extern "C" void DGGEV( char* JOBVL,  char* JOBVR,  int* N,
                       double* AA,  int* LDA,  double* B,  int* LDB,
                      double* ALPHAR, double* ALPHAI, double* BETA,
                      double* VL,  int* LDVL, double* VR,  int* LDVR,
                      double* WORK,  int* LWORK, int* INFO);

