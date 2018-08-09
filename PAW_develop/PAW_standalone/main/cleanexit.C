//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*       Projector Augmented Wave Standalone Code test 1d Chare array       */
//============================================================================
/* Standard Include files and class defs for */
#include "standard_include.h"        // has charmm++.h and pup.h
#include "class_PAW_info.h"    // c++ pawinfo readonly class declarations (including pup)
#include "cleanexit.decl.h"             // charm++ declarations for this module - 1D chare array
#include "cleanexit.h"                  // c++ declarations for this chare array
#include "scalar_cleanexit.h"        // scalar cleanexit routine

//============================================================================
/* charm++ global variable : instantiate in pawmain.C */
extern CProxy_cleanexit cleanexit_Proxy;
extern PAWINFO  readonly_pawinfo;

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* Constuctor invokes this rouine : array index 0 starts a round robin        */
//============================================================================
void cleanexit::startRoundRobin(){
   int next = (me+1 < nElements ? me+1 : 0);
   if(me==0){thisProxy[next].SayHi(me);sends++;}
}//end routine
//============================================================================
  
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*   Entry method : Receive a message from a index-1 element and send to index+1 */
//============================================================================
void cleanexit::SayHi(int indexm1)  { //
   recvs++;
   int next = (me+1 < nElements ? me+1 : 0);
   if(sends<nrounds){thisProxy[next].SayHi(me);sends++;}
   PRINTF("cleanexit chare array index %d received greetings from array index %d for the %d time\n",me,indexm1,recvs);
   fflush(stdout);
   if(recvs==nrounds && sends==nrounds){
     if (nElements > 1) {
       int i=1;
       PRINTF("cleanexit chare array index %d contributes to exit reduction : %d %d\n",me,recvs,sends);
       fflush(stdout);
       CkCallback cb(CkIndex_cleanexit::PAWcharmCleanExitSelf(NULL),thisProxy[0]);// send to array element 0 only
//        CkCallback cb(CkIndex_cleanexit::PAWcharmCleanExitSelf(NULL),thisProxy);// send to all array elements
        contribute(sizeof(int), &i, CkReduction::sum_int, cb);
     } else {
        PAWcharmCleanExitSelf0();
     } // end if
   }//endif
//----------------------------------------------------------------------------
 }// end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for cleanexit terminating round robin: invoke scalar exit after output  */
//=====================================================================================================
void cleanexit::PAWcharmCleanExitSelf(CkReductionMsg *m){
  doneRoundRobin = true;
  int value = ((int *) m->getData())[0];
  PRINTF("cleanexit chare array index %d has arrived in the reduction client Self: result %d correct answer %d: %d %d %d\n",
	 me,value,nElements,donePAWrho,doneRoundRobin,doneAtmsGrp);
  fflush(stdout);
  delete m;  // reductions don't have a no keep and you must delete messages for suffer memory leaks
  if (donePAWrho && doneRoundRobin && doneAtmsGrp) {PAWscalarCleanExit();}
}//end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for cleanexit terminating round robin: invoke scalar exit after output  */
//=====================================================================================================
void cleanexit::PAWcharmCleanExitSelf0(){
  doneRoundRobin = true;
  int value = 1;
  PRINTF("cleanexit chare array index %d has arrived in the reduction client Self0: result %d correct answer %d: %d %d ^%d\n",
	 me,value,nElements,donePAWrho,doneRoundRobin,doneAtmsGrp);
    fflush(stdout);
  if (donePAWrho && doneRoundRobin && doneAtmsGrp) {PAWscalarCleanExit();}
}// end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for cleanexit terminating PAWrho: invoke scalar exit after output  */
//=====================================================================================================
void cleanexit::PAWcharmCleanExitRho(CkReductionMsg *m){
  donePAWrho = true;
  int value = ((int *) m->getData())[0];
  PRINTF("cleanexit chare array index %d has arrived in the reduction client Rho: result %d correct answer %d: %d %d %d\n",
	 me,value,nElements,donePAWrho,doneRoundRobin,doneAtmsGrp);
  fflush(stdout);
  delete m;  // reductions don't have a no keep and you must delete messages for suffer memory leaks
  if (donePAWrho && doneRoundRobin && doneAtmsGrp) {PAWscalarCleanExit();}
}//endif
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for cleanexit terminating AtmsGrp: invoke scalar exit after output  */
//=====================================================================================================
void cleanexit::PAWcharmCleanExitAtmsGrp(CkReductionMsg *m){
  doneAtmsGrp = true;
  int value = ((int *) m->getData())[0];
  PRINTF("cleanexit chare array index %d has arrived in the reduction client AtmsGrp: result %d correct answer %d: %d %d %d\n",
	 me,value,nElements, donePAWrho,doneRoundRobin,doneAtmsGrp);
  fflush(stdout);
  delete m;  // reductions don't have a no keep and you must delete messages for suffer memory leaks
  if (donePAWrho && doneRoundRobin && doneAtmsGrp) {PAWscalarCleanExit();}
}//endif
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for cleanexit terminating Rho on 1 proc: invoke scalar exit after output  */
//=====================================================================================================
void cleanexit::PAWcharmCleanExitRho0(){
  donePAWrho = true;
  int value = 1;
  PRINTF("cleanexit chare array index %d has arrived in the reduction client Rho0: result %d correct answer %d: %d %d %d\n",
	 me,value,nElements,donePAWrho,doneRoundRobin,doneAtmsGrp);
  fflush(stdout);
  if (donePAWrho && doneRoundRobin && doneAtmsGrp) {PAWscalarCleanExit();}
}//end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Reduction callback function for cleanexit terminating AtmsGrp on 1 proc: invoke scalar exit after output  */
//=====================================================================================================
void cleanexit::PAWcharmCleanExitAtmsGrp0(){
  doneAtmsGrp = true;
  int value = 1;
  PRINTF("cleanexit chare array index %d has arrived in the reduction client AtmsGrp0: result %d correct answer %d: %d %d %d\n",
	 me,value,nElements,donePAWrho,doneRoundRobin,doneAtmsGrp);
  fflush(stdout);
  if (donePAWrho && doneRoundRobin && doneAtmsGrp) {PAWscalarCleanExit();}
}//end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Each PAWrho chare sends a msg here when it's done, when everyone reports, bcast PAWrhoDone()*/
//=====================================================================================================
void cleanexit::recvPAWrhoMsg(){
   PAWINFO * pawinfo = PAWINFO::get();
   int NatmChunk = pawinfo->NatmChunk;
   numPAWrhoMsg++;
   PRINTF("cleanexit chare array index %d has arrived in recvPAWrhoMsg for the %d time\n", me,numPAWrhoMsg);
   fflush(stdout);
   if (numPAWrhoMsg == NatmChunk) {
     thisProxy.bcastPAWrhoDone();
   }//endif
}//end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Each ATMSGRP chare sends a msg here when it's done, when everyone reports, bcast bcastAtmsGrpDone()*/
//=====================================================================================================
void cleanexit::recvAtmsGrpMsg(){
   int numPes = CkNumPes();
   numAtmsGrpMsg++;
   PRINTF("cleanexit chare array index %d has arrived in recvAtmsGrpMsg for the %d time\n", me,numAtmsGrpMsg);
   fflush(stdout);
   if (numAtmsGrpMsg == numPes) {
     thisProxy.bcastAtmsGrpDone();
   }//endif
}//end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Since PAWrho is done, contribute to a reduction will cause the code to exit  */
//=====================================================================================================
void cleanexit::bcastPAWrhoDone(){
  ncontrib++;
  PRINTF("cleanexit chare array index %d has arrived in bcastPAWrhoDone for the %d time\n", me,ncontrib);
  fflush(stdout);
  if (nElements > 1) {
    int i = 1;
    CkCallback cb(CkIndex_cleanexit::PAWcharmCleanExitRho(NULL),thisProxy[0]); // send to array element 0 only
    contribute(sizeof(int), &i, CkReduction::sum_int, cb);
  } else {
    PAWcharmCleanExitRho0();
  } // end if
}//end routine
//============================================================================

//=====================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=====================================================================================================
/*  Since AtmsGrp is done, contribute to a reduction will cause the code to exit  */
//=====================================================================================================
void cleanexit::bcastAtmsGrpDone(){
  PRINTF("cleanexit chare array index %d has arrived in bcastAtmsGrpDone\n", me);
  fflush(stdout);
  if (nElements > 1) {
    int i = 1;
    CkCallback cb(CkIndex_cleanexit::PAWcharmCleanExitAtmsGrp(NULL),thisProxy[0]); // send to array element 0 only
    contribute(sizeof(int), &i, CkReduction::sum_int, cb);
  } else {
    PAWcharmCleanExitAtmsGrp0();
  } // end if
}//end routine
//============================================================================

//============================================================================
/* Strange charm++ input that is needed at the bottom of files */
#include "cleanexit.def.h"
//============================================================================
