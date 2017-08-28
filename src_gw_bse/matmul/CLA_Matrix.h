#ifndef CLA_Matrix_H
#define CLA_Matrix_H

#include "CLA_Matrix.decl.h"
#include "ckmulticast.h"
#include "ckcomplex.h"

/* Class definitions */

/* The CLA_Matrix class should not be used directy by the user. See comments
 * below regarding CLA_Matrix_interface. */
class CLA_Matrix : public CBase_CLA_Matrix{
  friend class CLA_Matrix_interface;

  public:
    CLA_Matrix(){}
    CLA_Matrix(CkMigrateMessage *m){}
    ~CLA_Matrix();
    virtual void pup(PUP::er &p);
    virtual void ResumeFromSync();
    void synchronize(void);//{AtSync();}

    /* For 2D algorihtm */
    CLA_Matrix(int M, int K, int N, int m, int k, int n, int strideM,
     int strideN, int strideK, int part,
     CProxy_CLA_Matrix other1, CProxy_CLA_Matrix other2, CkCallback ready);
    void receiveA(CLA_Matrix_msg *m);
    void receiveB(CLA_Matrix_msg *m);

    /* For 3D algorithm */
    CLA_Matrix(CProxy_CLA_MM3D_multiplier p, int M, int K, int N, int m, int k,
     int n, int strideM, int strideK, int strideN, int part, CkCallback cb);//,
//     CkGroupID gid);
    void ready(CkReductionMsg *m);
    void readyC(CkReductionMsg *m);
    void mult_done(CkReductionMsg *m);
  private:
    void multiply(double alpha, double beta, complex *data,
     void (*fptr) (void *), void *usr_data);
    void multiply();

    /* shared */
    int M, K, N; // dimensions of global arrays
    int m, k, n; // (actual) dimensions of local arrays
    int um, uk, un; // dimensions given by users (above diff. if !(m | M), etc)
    int M_chunks, K_chunks, N_chunks;
    int M_stride, K_stride, N_stride;
    int part; // A, B, or C in C += A * B
    int algorithm;
    void (*fcb) (void *obj);
    void *user_data;
    double alpha, beta;

    /* For 2D algorithm */
    CProxySection_CLA_Matrix commGroup2D; // used by A and B
    complex *tmpA, *tmpB;
    complex    *dest; // used by C
    int row_count, col_count; // used by C
    CProxy_CLA_Matrix other1; // For A, B. For B, A. For C, A.
    CProxy_CLA_Matrix other2; // For A, C. For B, C. For C, B.

    /* For 3D algorithm */
    CProxySection_CLA_MM3D_multiplier commGroup3D; // used by all
    /* below used only by C */
    CkCallback init_cb;
    CkReductionMsg *res_msg;
    bool got_start, got_data;
};

/* This class below and the make_multiplier function below are the only way in
 * which a user of the library should interact with the libary. Users should
 * never explicitly create char CLA_Matrix object or call their entry methods.
 */

class CLA_Matrix_interface {
  public:
    CLA_Matrix_interface(){}
    inline void multiply(double alpha, double beta, complex *data,
     void (*fptr) (void *), void *usr_data, int x, int y){
      p(x, y).ckLocal()->multiply(alpha, beta, data, fptr, usr_data);
    }
    inline void sync(int x, int y){
      p(x, y).ckLocal()->synchronize();
    }
    inline void pup(PUP::er &per){ per | p;}
  private:
    CProxy_CLA_Matrix p;
    inline void init(CProxy_CLA_Matrix p){ this->p = p; }

  friend int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
   CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
   CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
   int M, int K, int N, int m, int k, int n, int strideM, int strideK,
   int strideN, CkCallback cbA, CkCallback cbB,
   CkCallback cbC, CkGroupID gid, int algorithm);
};

class CLA_Matrix_msg : public CkMcastBaseMsg, public CMessage_CLA_Matrix_msg {
  public:
    CLA_Matrix_msg(complex *data, int d1, int d2, int fromX, int fromY);
//    ~CLA_Matrix_msg(){delete [] data;}
    complex *data;
    int d1, d2;
    int fromX, fromY;
};

class CLA_MM3D_mult_init_msg : public CkMcastBaseMsg,
 public CMessage_CLA_MM3D_mult_init_msg {
  public:
    CLA_MM3D_mult_init_msg(CkCallback ready, CkCallback reduce){
      this->ready = ready;
      this->reduce = reduce;
    }
    CkCallback ready;
    CkCallback reduce;
};

class CLA_MM3D_mult_reinit_msg : public CkMcastBaseMsg,
 public CMessage_CLA_MM3D_mult_reinit_msg {
  public:
    CLA_MM3D_mult_reinit_msg(CkCallback ready){
      this->ready = ready;
    }
    CkCallback ready;
};

class CLA_MM3D_migrate_msg : public CkMcastBaseMsg,
 public CMessage_CLA_MM3D_migrate_msg {
};

class CLA_MM3D_Map : public CkArrayMap {
  public:
    CLA_MM3D_Map(int mc, int kc, int nc){
      M_chunks = mc;
      K_chunks = kc;
      N_chunks = nc;
      pes = CkNumPes();
    }
    virtual int procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex3D idx3d = *(CkArrayIndex3D *) &idx;
      return (N_chunks * idx3d.index[0] + N_chunks * M_chunks * idx3d.index[2] +
       idx3d.index[1]) % pes;
    }
  private:
    int M_chunks, K_chunks, N_chunks, pes;
};

class CLA_MM3D_multiplier : public CBase_CLA_MM3D_multiplier{
  public:
    CLA_MM3D_multiplier(){};
    CLA_MM3D_multiplier(CkMigrateMessage *m){};
    CLA_MM3D_multiplier(int m, int k, int n);
    ~CLA_MM3D_multiplier(){};
    virtual void pup(PUP::er &p);
    void initialize_reduction(CLA_MM3D_mult_init_msg *m);
    void reinitialize_reduction(CLA_MM3D_mult_reinit_msg *m);
    void receiveA(CLA_Matrix_msg *msg);
    void receiveB(CLA_Matrix_msg *msg);
    void multiply(complex *A, complex *B);
    void synchronize(CLA_MM3D_migrate_msg *msg);
    virtual void ResumeFromSync();
  private:
    int m, k, n;
    bool gotA, gotB;
    CLA_Matrix_msg *data_msg;
    CkSectionInfo sectionCookie;
    CkCallback reduce_CB;
    CkMulticastMgr *redGrp;
/*
    double *A, *B, *C;
*/
};

/* Function below creates the necessary interfaces so that
 * C = beta * C + alpha * A * B
 * can be computed. X will be bound to bindX for
 * X in {A, B, C}. A is M x K, B is K x N, C is M x N. M, K, and N are
 * decomposed into chunks of size m, k, and n, respectively. Along a given
 * dimension, the array must be of size ceil(1.0 * Y / y) (Y in {M, K, N},
 * y in {m, k, n}). If y does not divide Y, the last element must have only
 * Y % y elements along the given dimension. The variables strideY indicate
 * the stride at which elements are to be placed along that dimension.
 * When matrix X is ready, it does a
 * callback to cbX. A CKGroupID gid to be used by the library must be passed
 * in (this can be a simple CProxy_CkMulticastMgr::ckNew()). The algorithm
 * used to multiply the matrices is determined by the value of 'algorithm',
 * which can take the values defined below.
 *
 * Return value: Zero is returned upon success. A negative value is returned
 * if an error occurs.
 *  ERR_INVALID_ALG: an invalid algorithm was selected
 *  ERR_INVALID_DIM: invalid dimensions were given
 */

/* Valid values for 'algorithm' are given below. As new ones are added,
 * they should be given incremental numbers. MM_ALG_MIN should not be changed,
 * and MM_ALG_MAX should have the value of the greatest defined algorithm.
 * MM_ALG_MIN MM_ALG_MAX are used to validate user input, so they should always
 * be updated as new algorithms are added.
 */
#define MM_ALG_MIN  1
#define MM_ALG_2D   1
#define MM_ALG_3D   2
#define MM_ALG_MAX  2

int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
 CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
 CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
 int M, int K, int N, int m, int k, int n, int strideM, int strideK,
 int strideZ, CkCallback cbA, CkCallback cbB,
 CkCallback cbC, CkGroupID gid, int algorithm);

/* Error codes */
#define SUCCESS         0
#define ERR_INVALID_ALG     -1
#define ERR_INVALID_DIM     -2

#endif
