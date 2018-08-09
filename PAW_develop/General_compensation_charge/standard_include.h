//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
// This is the standard include files
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "gperftools/profiler.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <time.h>
//extern "C" {
//#include <quadmath.h>
//}
//#include "ckcomplex.h"

using namespace std;
//==========================================================================
// PINY constants, inline functions and typedefs

#define BOLTZ 315777.0       // 1/k_b in atomic units
#define BOHR  0.529177       // bohr to angstrom
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI_QI 3.14159265358979323846264338327950288419716939937510

// constants for approx to erfc
#define PERFC  (0.3614)
#define CERFC1 (0.2041422096422003)
#define CERFC2 (0.1997535956961481)
#define CERFC3 (0.2213176596405576)
#define CERFC4 (0.03360430734640255)
#define CERFC5 (0.4732592578721755) 
#define CERFC6 (-0.509078520069735)
#define CERFC7 (0.6772631491947646)
#define CERFC8 (-0.369912979092217)
#define CERFC9 (0.06965131976970335)
#define DCERFC1 (1.0*CERFC1)
#define DCERFC2 (2.0*CERFC2)
#define DCERFC3 (3.0*CERFC3)
#define DCERFC4 (4.0*CERFC4)
#define DCERFC5 (5.0*CERFC5)
#define DCERFC6 (6.0*CERFC6)
#define DCERFC7 (7.0*CERFC7)
#define DCERFC8 (8.0*CERFC8)
#define DCERFC9 (9.0*CERFC9)
#define PRE_ERFC (2.0/sqrt(M_PI_QI))

#define MAXWORD   80         // length of a string in OA
#define MAXLINE  100         // length of a line in OA
typedef char NAME[MAXWORD];  

// some generic definitionss we overwrite to charm++
#define PRINTF printf        
#define FFLUSH fflush
#define EXIT(N) {exit(N);}

// some nice formatting
#define PRINT_LINE_STAR {PRINTF("==============================================================================\n");}
#define PRINT_LINE_DASH {PRINTF("------------------------------------------------------------------------------\n");}

#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))
#define NINT(X) ( (int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)) )

//==========================================================================
// ATOM_MAPS structure
typedef struct ATOM_MAPS{
	int natm_typ;         // number of atom types
	int natm;             // number of atoms
	int natm_atm_typ_max; // max number of atoms of any type
	int *index_atm_typ;   // index of atom type of each atom
	int *natm_atm_typ;    // the number of atoms of each type
	int **list_atm_by_typ;// list of atoms sorted by atom type
	NAME *atm_typ;        // names of the atom types
}ATOM_MAPS;

//==========================================================================
// ATOM_POS structure
typedef struct ATOM_POS{
	int natm;
	double *x,*y,*z;			      // atom coordinates
	double *vx,*vy,*vz;			    // atom velocities
	double *fx0,*fy0,*fz0;		  // atom forces iperd = 0 from analytical model
	double *fx,*fy,*fz;			    // atom forces iperd = 3 from analytical model
	double *fx0g,*fy0g,*fz0g;   // atom forces iperd = 0 from the grid
	double *fxg,*fyg,*fzg;      // atom forces iperd = 3 from the grid
	double *q,*qt;				      // charges
	double *alp;  		    	    // Gaussian parameter for the core 
	double *beta;		            // beta screener beta[J]
	double *Rpc;				        // Rpc
}ATOM_POS;

//==========================================================================
// CELL structure
typedef struct CELL{
	int iperd;  			  // periodicity
	int nimg;						// number of images
	double animg;				// (double) nimg
	double alpb;        // Ewald alpha
	double Rcut;				// Ewald cutoff real space
	double gcut;        // state g space
	double Gcut;				// density g cutoff, twice gcut
	double hmat[10];    // the simulation box
	double hmati[10];   // inverse simulation box
	double volume;      // simulation box volume
}CELL;

