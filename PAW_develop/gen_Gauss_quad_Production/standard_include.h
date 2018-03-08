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
extern "C" {
#include <quadmath.h>
}

using namespace std;
//==========================================================================
// PINY constants, inline functions and typedefs

#define BOLTZ 315777.0       // 1/k_b in atomic units
#define BOHR  0.529177       // bohr to angstrom

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif
#define M_PI_QI 3.14159265358979323846264338327950288419716939937510

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

