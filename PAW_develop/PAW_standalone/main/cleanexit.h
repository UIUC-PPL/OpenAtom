//==========================================================================
//                PAW info                          
//             {Variables needed for mem allocation:                        
//                 vxctyplen,nsplin_g,nsplin_g_tot,n_ang_max,               
//                 num_nl_lst       }                                       
//                                                                          
//==========================================================================

#ifndef _CLEANEXIT_
#define _CLEANEXIT_

class cleanexit : public CBase_cleanexit {
public:
  //======================================
  // integer member variables
    bool doneRoundRobin;
    bool donePAWrho;
    bool doneAtmsGrp;
    int recvs;
    int sends;
    int nrounds;
    int nElements;
    int me;
    int numPAWrhoMsg;
    int numAtmsGrpMsg;
    int ncontrib;
  //==========================
  // Constructor
  cleanexit(void)  {
    doneRoundRobin = false;
    donePAWrho     = false;
    doneAtmsGrp    = false;
    me             = thisIndex;
    nElements      = CkNumPes();
    recvs          = 0;
    sends          = 0;
    nrounds        = 3;
    numPAWrhoMsg   = 0;
    numAtmsGrpMsg  = 0;
    ncontrib       = 0;
    CkPrintf("cleanexit chare array index %d of 0-%d reporting to constructor \n",me, nElements-1);
    startRoundRobin();
  }//end constuctor
  void startRoundRobin();
  void SayHi(int );
  void recvPAWrhoMsg();
  void recvAtmsGrpMsg();
  void PAWcharmCleanExitSelf(CkReductionMsg *);
  void PAWcharmCleanExitRho(CkReductionMsg *);
  void PAWcharmCleanExitAtmsGrp(CkReductionMsg *);
  void PAWcharmCleanExitSelf0();
  void PAWcharmCleanExitRho0();
  void PAWcharmCleanExitAtmsGrp0();
  void bcastPAWrhoDone();
  void bcastAtmsGrpDone();
};
#endif
