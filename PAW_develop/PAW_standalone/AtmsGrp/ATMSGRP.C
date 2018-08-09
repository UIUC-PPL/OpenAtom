//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*       Projector Augmented Wave Standalone Code test 1d Chare array       */
//============================================================================
/* Standard Include files and class defs for */
#include "standard_include.h"        // has charmm++.h and pup.h
#include "class_PAW_info.h"
#include "ATMSGRP.decl.h"             // charm++ declarations for this module - 1D chare array
#include "cleanexit.decl.h"
#include "atmsGrp.h"                  // c++ declarations for this chare array
#include "cleanexit.h"

//============================================================================
/* charm++ global variable : instantiate in pawmain.C */
extern PAWINFO  readonly_pawinfo;
extern CProxy_cleanexit cleanexit_Proxy;  // proxy to cleanexit chare array for 

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  ATMSGRP constructor */
//=====================================================================================================
ATMSGRP::ATMSGRP(ATMSGRPMSG * atmsGrpMsg)  {
	//=================================================================================
	// get your identity and output it
  me             = thisIndex;
  nElements      = CkNumPes();
  CkPrintf("atmsGrp chare array index %d of 0-%d reporting to constructor \n",me, nElements-1);
	
	//=================================================================================
	// prepare the fastatoms
  atmsReady      = false;
  natm           = atmsGrpMsg->natm;
  fastatoms.natm = natm;
  fastatoms.allocate();  // asllocate space for the atoms
  for(int i=0; i<natm; i++) {
    fastatoms.x[i]    = atmsGrpMsg->x[i];
    fastatoms.y[i]    = atmsGrpMsg->y[i];
    fastatoms.z[i]    = atmsGrpMsg->z[i];
    fastatoms.q[i]    = atmsGrpMsg->q[i];
    fastatoms.qt[i]   = atmsGrpMsg->qt[i];
    fastatoms.alp[i]  = atmsGrpMsg->alp[i];
    fastatoms.beta[i] = atmsGrpMsg->beta[i];
    fastatoms.Rpc[i]  = atmsGrpMsg->Rpc[i];
  } // end for    
	PAWINFO * pawinfo = PAWINFO::get();
  CELL *cell = &(pawinfo->cell);
  ATOM_MAPS *atom_maps = &(pawinfo->atom_maps);
  fastatoms.output(me, cell, atom_maps);
	
	//=================================================================================
	// create an PAW f-grid for all the atom types
	 

	//=================================================================================
	// perform a reduction so everyone knows we are done
  if (nElements > 1) {
    CkCallback cb(CkIndex_ATMSGRP::AtmsGrpConstructed(NULL),thisProxy); // send to everyone
    int i = 1;
    contribute(sizeof(int), &i, CkReduction::sum_int, cb);
  } else {
    AtmsGrpConstructed0();
  } // end if
}//end constuctor
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for AtmsGrp terminating round robin: invoke scalar exit after output  */
//=====================================================================================================
void ATMSGRP::AtmsGrpConstructed(CkReductionMsg *m){
  int value = ((int *) m->getData())[0];
  PRINTF("AtmsGrp chare array index %d has arrived in the reduction client AtmsGrpConstructed: result %d correct answer %d\n",
   me,value,nElements);
  delete m;  // reductions don't have a no keep and you must delete messages for suffer memory leaks
  cleanexit_Proxy[0].recvAtmsGrpMsg();
}// end rotuine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for AtmsGrp terminating round robin: invoke scalar exit after output  */
//=====================================================================================================
void ATMSGRP::AtmsGrpConstructed0(){
  int value = 1;
  PRINTF("AtmsGrp chare array index %d has arrived in the reduction client AtmsGrpConstructed0: result %d correct answer %d\n",
      me,value,nElements);
  cleanexit_Proxy[0].recvAtmsGrpMsg();
}//end routine
//============================================================================

//============================================================================
/* Strange charm++ input that is needed at the bottom of files */
#include "ATMSGRP.def.h"
//============================================================================
