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

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include "ckcomplex.h"

using namespace std;
//==========================================================================
// PINY constants, inline functions and typedefs

#define BOLTZ 315777.0       // 1/k_b in atomic units
#define BOHR  0.529177       // bohr to angstrom
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
	double *x,*y,*z;			// atom coordinates
	double *vx,*vy,*vz;			// atom velocities
	double *fx0,*fy0,*fz0;		// atom forces iperd = 0 from analytical model
	double *fx,*fy,*fz;			// atom forces iperd = 3 from analytical model
	double *fx0g,*fy0g,*fz0g;   // atom forces iperd = 0 from the grid
	double *fxg,*fyg,*fzg;      // atom forces iperd = 3 from the grid
	double *q,*qt;				// charges
	double *alp;  		    	// Gaussian parameter for the core 
}ATOM_POS;

//==========================================================================
// CELL structure
typedef struct CELL{
    int iperd;  			        // periodicity
    double alpb;          			// Ewald alpha
    double gcut;          			// Ewald reciprocal space cutoff 
    double hmat[10];                   // the simulation box
    double hmati[10];                  // inverse simulation box
    double volume;                  // simulation box volume
}CELL;

//==========================================================================
// FGRID structure: spherical polar with DVR
typedef struct FGRID{
	int nf;				// total number of f points
	int nr;				// number of r grid points
	int ntheta;			// number of theta grid points
	int nphi;			// number of phi grid points
	int nrfull;			// number of r grid points*2
	double *wf;			// weight including the Jacobian
	double *xf;			// x coordinate
	double *yf;			// y
	double *zf;			// z
	double *rf;			// r
	double *xcostheta;  // xcostheta - used to generate Ylmf
	double *xphi;		// xphi - used to generate Ylmf
	complex *Ylmf;		// memory to store 1 spehrical harmonic on the f grid
}FGRID;

