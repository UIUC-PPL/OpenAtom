//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <quadmath.h>

using namespace std;
int main(int, char **);
extern "C" void temp_(__float128 *);

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller for generalized Gaussian quadrature using 
//				   companion matrix 
//==========================================================================
int main (int argc, char *argv[]){
//==========================================================================

  __float128 junk = 23.q;
  __float128 junk2 = 0.5q;
  temp_(&junk);
  printf("%d %.30Lg %.30Lg\n",(int) sizeof(junk),(long double) junk,(long double) junk2);

//--------------------------------------------------------------------------
	return 1;
} // end routine
//==========================================================================
