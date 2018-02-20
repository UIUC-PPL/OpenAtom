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
		long double   Norm;		 // N_ii
        long double * c;          // c matrix (order+1)
        long double * a;          // a matrix (order+1)
        long double * O;          // O integral (n+1), its upper triangle, too lazy to do the pointers
		long double * p_at_nodes; // size n
}POLY;

//==========================================================================
// function prototypes 

void horner(int, long double *, long double, long double *);
void testgrid(int, long double *, POLY *, int, int *);
void testnodes(int, long double * , long double *,int, int *);

#ifdef FORT_UNDER
#define DGGEV dggev_
#else
#define DGGEV dggev
#endif

extern "C" void DGGEV( char* JOBVL,  char* JOBVR,  int* N,
                       double* A,  int* LDA,  double* B,  int* LDB,
                      double* ALPHAR, double* ALPHAI, double* BETA,
                      double* VL,  int* LDVL, double* VR,  int* LDVR,
                      double* WORK,  int* LWORK, int* INFO);

