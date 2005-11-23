/** \file ckPairCalculator.h
 *
 */

#ifndef _ckPairCalculator_h_
#define _ckPairCalculator_h_
#define _PC_COMMLIB_MULTI_ 0
#include "pairutil.h"
#include "ckmulticast.h"
#include "ckhashtable.h"
#include "cksparsecontiguousreducer.h"
#include "PipeBroadcastStrategy.h"
#include "BroadcastStrategy.h"
#include "DirectMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#include "MultiRingMulticast.h"
#include "NodeMulticast.h"

// Flag to use sparse reduction or regular reduction

// Debugging flag for Verbose output
//#define _PAIRCALC_DEBUG_

// Optimize flags: 
#define _PAIRCALC_FIRSTPHASE_STREAM_
#define _PAIRCALC_USE_ZGEMM_
#define _PAIRCALC_USE_DGEMM_

#define NORMALPC   0
#define KEEPORTHO  1
#define PSIV       2


//#define _PAIRCALC_SECONDPHASE_LOADBAL_
enum redtypes {section=0, machine=1, sparsecontiguous=2};
PUPbytes(redtypes);

#ifdef FORTRANUNDERSCORE
#define ZGEMM zgemm_ 
#define DGEMM dgemm_ 
#define DCOPY dcopy_
#define ZTODO ztodo_
#else
#define ZGEMM zgemm
#define DGEMM dgemm
#define DCOPY dcopy
#define ZTODO ztodo
#endif
extern ComlibInstanceHandle mcastInstanceCP;


#ifdef _PAIRCALC_USE_BLAS_
extern "C" complex ZTODO( const int *N,  complex *X, const int *incX, complex *Y, const int *incY);

#endif

#ifdef _PAIRCALC_USE_DGEMM_

//extern "C" void DGEMM(char *,char *, int *,int *, int *,double *,double *,int *, double *,int *,double *,double *,int *);
extern "C" {void DGEMM (char *, char *, int *, int *, int *,double *,double *,
                        int *, double *, int *, double *, double *, int * );}


#endif

#ifdef _PAIRCALC_USE_ZGEMM_
extern "C" void ZGEMM(char *,char *, int *,int *, int *,complex *,complex *,int *,
                       const complex *,int *,complex *,complex *,int *);

extern "C" void DCOPY(int*,double *,int*, double *,int *);

#endif

typedef void (*FuncType) (complex a, complex b);
PUPmarshallBytes(FuncType);

#include "ckPairCalculator.decl.h"

class IndexAndID {
 public:
  CkArrayIndex4D idx;
  CkArrayID id;
  IndexAndID(CkArrayIndex4D *_idx, CkArrayID *_id) 
    {
      idx=*_idx;
      id=*_id;
    }
  IndexAndID(CkArrayIndex4D _idx, CkArrayID _id): idx(_idx),id(_id)
    {
    }
  void dump()
    {
      CkPrintf("w %d x %d y %d z %d\n",idx.index[0],idx.index[1],idx.index[2],idx.index[3]);
    }
  void sdump(char *string, int size)
    {
      snprintf(string,size,"w %d x %d y %d z %d\n",idx.index[0],idx.index[1],idx.index[2],idx.index[3]);
    }
  IndexAndID()
    {
    }
};

PUPbytes(IndexAndID);

class PairCalcReducer : public Group {
 public:
  PairCalcReducer(CkMigrateMessage *m) { }
  PairCalcReducer(){ 
    acceptCount = 0; numRegistered[0] = 0; numRegistered[1] = 0;
    reduction_elementCount = 0;
    tmp_matrix = NULL;
  }
  ~PairCalcReducer() {if (tmp_matrix !=NULL) delete [] tmp_matrix;}
  void clearRegister()
    {  
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d] clearing register\n",CkMyPe());
#endif
      acceptCount=0;
      reduction_elementCount=0;
      localElements[0].removeAll();
      localElements[1].removeAll();
      numRegistered[0]=0;
      numRegistered[1]=0;
      /*      if(tmp_matrix!=NULL)
	      {
	      delete [] tmp_matrix;
	      tmp_matrix=NULL;
	      }
      */
    }
  void broadcastEntireResult(entireResultMsg *msg);
  void broadcastEntireResult(entireResultMsg2 *msg);
  void broadcastEntireResult(int size, double* matrix, bool symmtype);
  void broadcastEntireResult(int size, double* matrix1, double* matrix2, bool symmtype);
  void doRegister(IndexAndID elem, bool symmtype){

    numRegistered[symmtype]++;
    localElements[symmtype].push_back(elem);
#ifdef _PAIRCALC_DEBUG_
    char string[80];
    localElements[symmtype][localElements[symmtype].length()-1].sdump(string,80);
    CkPrintf("[%d] registered %d %s",CkMyPe(),symmtype,string);
#endif
  }

  void acceptContribute(int size, double* matrix, CkCallback cb, bool isAllReduce, bool symmtype, int offx, int offy, int inputsize);
  
  void startMachineReduction();

 private:
  CkVec<IndexAndID> localElements[2];
  int numRegistered[2];
  int acceptCount;
  int reduction_elementCount;
  double *tmp_matrix;
  bool isAllReduce;
  int size;
  bool symmtype;
  CkCallback cb;

  void pup(PUP::er &p){
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("PairCalcReducer on %d pupping\n",CkMyPe());
#endif
    p|localElements[0];
    p|localElements[1];
    p(numRegistered,2);
    p|acceptCount;
    p|reduction_elementCount;
    p|isAllReduce;
    p|size;
    p|symmtype;
    p|cb;
    if(p.isUnpacking())
      {
	tmp_matrix=NULL;
      }
  }

}; 


class initGRedMsg : public CkMcastBaseMsg, public CMessage_initGRedMsg {
 public:
  CkCallback cb;
  CkGroupID mCastGrpId;
  CkCallback synccb;
  bool lbsync;
  friend class CMessage_initGRedMsg;
};

class initResultMsg : public CkMcastBaseMsg, public CMessage_initResultMsg {
 public:
  int offset;
  int dest;
  CkGroupID mCastGrpId;
  CkCallback synccb;
  bool lbsync;
  friend class CMessage_initGRedMsg;
};

class sendBWsignalMsg : public CMessage_sendBWsignalMsg{
};

class mySendMsg : public CMessage_mySendMsg {
 public:
  int N;
  complex *data;
  friend class CMessage_mySendMsg;
};

class partialResultMsg : public CMessage_partialResultMsg {
 public:
  complex *result;
  int N;
  int myoffset;
  int priority;
  CkCallback cb;

  friend class CMessage_partialResultMsg;
};

class priorSumMsg : public CMessage_priorSumMsg {
 public:
  complex *result;
  int N;
  int priority;
  CkCallback cb;

  friend class CMessage_priorSumMsg;

};

class calculatePairsMsg : public CkMcastBaseMsg, public CMessage_calculatePairsMsg {
 public:
  int size;
  int sender;
  bool fromRow;
  bool flag_dp;
  complex *points;
  bool doPsiV;
  void init(int _size, int _sender, bool _fromRow, bool _flag_dp, complex *_points , bool _doPsiV)
    {
      size=_size;
      sender=_sender;
      fromRow=_fromRow;
      flag_dp=_flag_dp;
      doPsiV=_doPsiV;
      memcpy(points,_points,size*sizeof(complex));
    }
  friend class CMessage_calculatePairsMsg;

};

class multiplyResultMsg : public CkMcastBaseMsg, public CMessage_multiplyResultMsg {
 public:
  double *matrix1;
  double *matrix2;
  int size;
  int size2;
  int actionType;
  void init(int _size, int _size2, double *_points1, double *_points2, bool _actionType)
    {
      size=_size;
      size2=_size2;
      memcpy(matrix1,_points1,size*sizeof(double));
      memcpy(matrix2,_points2,size2*sizeof(double));
      actionType=_actionType;
    }
  void init1(int _size, double *_points1, int _actionType)
    {
      size=_size;
      size2=0;
      memcpy(matrix1,_points1,size*sizeof(double));
      actionType=_actionType;
      // this field does nothing in minimization
      matrix2=NULL;
    }
  friend class CMessage_multiplyResultMsg;
};

class entireResultMsg : public CMessage_entireResultMsg {
 public:
  double *matrix;
  int size;
  bool symmetric;
  void init(int _size, double *_points, bool _symmetric)
    {
      size=_size;
      symmetric=_symmetric;
      memcpy(matrix,_points,size*sizeof(double));
    }
  friend class CMessage_entireResultMsg;
};

class entireResultMsg2 : public CMessage_entireResultMsg2 {
 public:
  double *matrix1;
  double *matrix2;
  int size;
  bool symmetric;
  void init(int _size, double *_points1, double *_points2, bool _symmetric)
    {
      size=_size;
      symmetric=_symmetric;
      memcpy(matrix1,_points1,size*sizeof(double));
      memcpy(matrix2,_points2,size*sizeof(double));
    }
  friend class CMessage_entireResultMsg2;
};




class PairCalculator: public CBase_PairCalculator {
 public:
  PairCalculator(bool sym, int grainSize, int s, int blkSize, CkCallback cb,  CkArrayID final_callbackid, int final_callback_ep, int callback_ep_tol, bool conserveMemory, bool lbpaircalc, redtypes reduce);
    
  PairCalculator(CkMigrateMessage *);
  ~PairCalculator();
  void lbsync() {
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d] atsyncs\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z);
#endif
    resumed=false;
    rck=0;
    AtSync();
  };
  void multiplyForward(bool);
  void ResumeFromSync();
  void initGRed(initGRedMsg *msg);
  void acceptPairData(calculatePairsMsg *msg);
  void sendBWResult(sendBWsignalMsg *msg);
  void multiplyResult(multiplyResultMsg *msg);
  void multiplyPsiV();
  void multiplyResultI(multiplyResultMsg *msg);
  void initResultSection(initResultMsg *msg);
  void pup(PUP::er &);
  // this is never called we leave it here as an historical artifact
  inline double compute_entry(int n, complex *psi1, complex *psi2, int op) 
    {

      int i;
      register double sum = 0;
      for (i = 0; i < n; i++) {
	sum += psi1[i].re * psi2[i].re +  psi1[i].im * psi2[i].im;
      }
        
      return sum;
    }

 private:
  int numRecd, numExpected, grainSize, S, blkSize, N;
  bool symmetric;
  bool conserveMemory;
  bool lbpaircalc;
  redtypes cpreduce;
  CkCallback cb;
  CkArrayID cb_aid;
  int cb_ep;
  int cb_ep_tol;
  bool existsLeft;
  bool existsRight;
  bool existsOut;
  bool existsNew;
  bool resumed;
  CkSectionInfo cookie; 
  complex *mynewData, *othernewData;
  double *inDataLeft, *inDataRight;
  double *outData;
  int actionType;
  
  /* to support the simpler section reduction*/
  int rck;
  CkGroupID mCastGrpId;

  CkSectionInfo *resultCookies;
  CkSectionInfo *otherResultCookies;
};

//forward declaration
CkReductionMsg *sumMatrixDouble(int nMsg, CkReductionMsg **msgs);

#endif
