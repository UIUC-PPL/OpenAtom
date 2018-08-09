/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

#ifndef _PAWATMRHO_
#define _PAWATMRHO_

#include "atom_maps.h"
#include "cell.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class PAWATMRHO {
//============================================================================
 public:

  int natm;
  double *x,*y,*z;            // no masses or velocities cause these are working vectors
  double *fx,*fy,*fz;					// atom forces grid
  double *fxa,*fya,*fza;      // atom forces analytical
  double *q,*qt;              // charges
  double *alp;                // Gaussian parameter for the core 
  double *beta;               // beta screener beta[J]
  double *Rpc;                // Rpc

  PAWATMRHO(){natm=0;x=y=z=fx=fy=fz=NULL;}
 ~PAWATMRHO(){}

	void allocate(){
 	     x = new double [natm];    y = new double [natm];   z = new double [natm];
	    fx = new double [natm];   fy = new double [natm];  fz = new double [natm]; // not sent in message
	   fxa = new double [natm];  fya = new double [natm]; fza = new double [natm]; // not sent in message
	     q = new double [natm];   qt = new double [natm];
	   alp = new double [natm]; beta = new double [natm]; Rpc = new double [natm]; 
	} // end routine

	void destroy(){
 	  delete []  x;  delete []  y;   delete []  z;
	  delete [] fx;  delete [] fy;   delete [] fz;
	  delete [] fxa; delete [] fya;  delete [] fza;
	  delete []  q;  delete [] qt;
	  delete [] alp; delete [] beta; delete [] Rpc;
	} //end routine

	void readatoms(CELL *, ATOM_MAPS *, NAME, int); // done once on Proc 0 at start up
	void output(int, CELL *, ATOM_MAPS *); //write output into a file
  
};//end PAWATMRHO
//============================================================================
#endif
