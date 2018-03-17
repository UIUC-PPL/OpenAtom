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
	__float128   Norm;		 // N_ii
        __float128 * c;          // c matrix (order+1)
        __float128 * a;          // a matrix (order+1)
        __float128 * O;          // O integral (n+1), its upper triangle, too lazy to do the pointers
	__float128 * p_at_nodes; // size n
}POLY;

//==========================================================================
// function prototypes 

void horner(int, __float128 *, __float128, __float128 *);
void testgrid(int, __float128 *, POLY *, int, int *);
void testnodes(int, __float128 * , __float128 *,int, int *);

#ifdef FORT_UNDER
#define DGGEV dggev_
#else
#define DGGEV dggev
#endif

// this guys is different because of use_iso_c and bind(c) in the fortan source
#define QGGEV qggev

extern "C" {
void DGGEV(char* JOBVL,  char* JOBVR,  int* N,
           double* A,  int* LDA,  double* B,  int* LDB,
           double* ALPHAR, double* ALPHAI, double* BETA,
           double* VL,  int* LDVL, double* VR,  int* LDVR,
	   double* WORK,  int* LWORK, int* INFO);
}

extern "C" {
void QGGEV(char* JOBVL,  char* JOBVR,  int* N,
           __float128* A,  int* LDA,  __float128* B,  int* LDB,
           __float128* ALPHAR, __float128* ALPHAI, __float128* BETA,
           __float128* VL,  int* LDVL, __float128* VR,  int* LDVR,
	   __float128* WORK,  int* LWORK, int* INFO);
}
