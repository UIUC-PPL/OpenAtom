//==========================================================================
//                PAW info                          
//				PAW atom info, corresponds to ATOM_POS in the old code
//                                                                          
//==========================================================================

#ifndef _ATMSGRP_
#define _ATMSGRP_
#include "fastAtoms.h"        // Fastatom class

class ATMSGRPMSG : public CMessage_ATMSGRPMSG
{ // size of the msg is (8*8 + 1) sizeof(int)
 public:
	int natm;
	double *x,*y,*z;            // no masses or velocities cause these are working vectors
	double *q,*qt;              // charges
	double *alp;                // Gaussian parameter for the core 
	double *beta;               // beta screener beta[J]
	double *Rpc;                // Rpc
};

class ATMSGRP : public CBase_ATMSGRP {
public:
  //======================================
  // integer member variables
  FASTATOMS fastatoms; // instance of fastatoms
  bool atmsReady;
  int nElements;
  int me;
  int natm;
  //==========================
  // Member functions
	ATMSGRP(ATMSGRPMSG *);
	void AtmsGrpConstructed(CkReductionMsg *);
	void AtmsGrpConstructed0();
};
#endif
