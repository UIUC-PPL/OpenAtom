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
#ifndef _STD_INCLUDE_
#define _STD_INCLUDE_

#ifndef CHARM_OFF
#include "charm++.h"
#include <pup.h>
#endif
//#include "gperftools/profiler.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <unistd.h>
#include <sstream>
#include <string>
#include <assert.h>
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
#define MAXWORD   80         // length of a string in OA
#define MAXLINE  100         // length of a line in OA
typedef char NAME[MAXWORD];
//==========================================================
// some generic definitionss we overwrite to charm++

#ifndef CHARM_OFF
#define PRINTF CkPrintf
#ifdef PUP_PRINTING_ON
#define PUP_PRINTF CkPrintf
#else
#define PUP_PRINTF
#endif
#define SCANF  CkScanf
#define FFLUSH fflush
#define EXIT(N) {CkExit();}
#endif

#ifdef CHARM_OFF
#define PUP_PRINTF printf
#define PRINTF printf
#define FFLUSH fflush
#define SCANF  scanf
#define EXIT(N) {exit(N);}
#endif

//==========================================================
// some nice formatting
#define PRINT_LINE_STAR {PRINTF("==============================================================================\n");}
#define PRINT_LINE_DASH {PRINTF("------------------------------------------------------------------------------\n");}

//==========================================================
// Useful functions
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))
#define NINT(X) ( (int) ((X)>=0.0 ? ((X)+0.5):((X)-0.5)) )

#endif
