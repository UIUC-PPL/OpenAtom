//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*       Projector Augmented Wave Standalone Code test 1d Chare array       */
//============================================================================
/* Standard Include files and class defs for */
#include "standard_include.h"       // has charmm++.h and pup.h in addition to usual suspects
#include "class_PAW_info.h"         // c++ pawinfo readonly class declarations (including pup)
#include "cleanexit.decl.h"         // charm++ declarations for this module - 1D chare array
#include "PAWrho.decl.h"            // charm++ declarations for this module - 1D chare array
#include "ATMSGRP.decl.h"
#include "cleanexit.h"              // c++ cleanexit declarations
#include "PAWrho.h"                 // c++ pawrho declarations 
#include "atmsGrp.h"

//============================================================================
/* charm++ global variable : instantiate in pawmain.C */
extern PAWINFO  readonly_pawinfo;
extern CProxy_ATMSGRP ATMSGRP_Proxy;  // proxy to atmsGrp chare array for launch
extern CProxy_cleanexit cleanexit_Proxy;  // proxy to cleanexit chare array for launch
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*           PAwRho constructor */
//============================================================================
PAWrho::PAWrho(void) {
  PAWINFO * pawinfo = PAWINFO::get();
	ATMSGRP * atmsGrpLoc = ATMSGRP_Proxy.ckLocalBranch();
	FASTATOMS * fastatoms = &(atmsGrpLoc->fastatoms);
  me        = thisIndex;
  nElements = pawinfo->NatmChunk;
  Natm      = pawinfo->Natm;
  numMyAtm  = Natm/nElements;
  int rem   = (Natm % nElements);
  if (me < rem) {numMyAtm++;}
  if (me < rem) {
    myAtmStart = numMyAtm*me;
  } else {
    myAtmStart = numMyAtm*me + rem;
  }// end if
  myAtmEnd = myAtmStart + numMyAtm;
  CkPrintf("PAWrho chare array index %d of 0-%d reporting to constructor has %d atoms index %d-%d\n",
					me, nElements-1, numMyAtm, myAtmStart, myAtmEnd-1);
  CkPrintf("  Starting atom: %g %g %g\n", fastatoms->x[myAtmStart], fastatoms->y[myAtmStart], fastatoms->z[myAtmStart]);
  CkPrintf("    Ending atom: %g %g %g\n", fastatoms->x[myAtmEnd-1], fastatoms->y[myAtmEnd-1], fastatoms->z[myAtmEnd-1]);
  if (nElements > 1) {
    CkCallback cb(CkIndex_PAWrho::PAWrhoConstructed(NULL),thisProxy); // send to everyone
    int i = 1;
    contribute(sizeof(int), &i, CkReduction::sum_int, cb);
  } else {
    PAWrhoConstructed0();
  } // end if
//----------------------------------------------------------------------------------------------------
 }//end constuctor
//=====================================================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for PAWrho terminating round robin: invoke scalar exit after output  */
//=====================================================================================================
void PAWrho::PAWrhoConstructed(CkReductionMsg *m){
  int value = ((int *) m->getData())[0];
  PRINTF("PAWrho chare array index %d has arrived in the reduction client RhoConstructed: result %d correct answer %d\n",
	 me,value,nElements);
  delete m;  // reductions don't have a no keep and you must delete messages for suffer memory leaks
  cleanexit_Proxy[0].recvPAWrhoMsg();
}// end rotuine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for PAWrho terminating round robin: invoke scalar exit after output  */
//=====================================================================================================
void PAWrho::PAWrhoConstructed0(){
  int value = 1;
  PRINTF("PAWrho chare array index %d has arrived in the reduction client RhoConstructed0: result %d correct answer %d\n",
      me,value,nElements);
  cleanexit_Proxy[0].recvPAWrhoMsg();
}//end routine
//============================================================================

//============================================================================
/* Strange charm++ input that is needed at the bottom of files */
#include "PAWrho.def.h"
//============================================================================
