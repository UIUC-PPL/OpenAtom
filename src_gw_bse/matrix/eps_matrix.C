#include "standard_include.h"
#include "allclass_gwbse.h"
#include "eps_matrix.h"
#include "messages.h"
#include "pmatrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"

#include <cstring> // for memcpy
using std::memcpy;

#define eps_rows 20
#define eps_cols 20
#define IDX_eps(r,c) ((r)*eps_cols + (c))

void example_dgemm(int M, int N, int K, double alpha,
                   complex *A, complex *B, complex *C) {
  /* multiply */
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      complex sum = 0.0;
      for (int k = 0; k < K; ++k) {
        sum += A[i*K + k] * B[k*N + j];
      }
      C[N*i + j] = C[N*i + j] + alpha*sum;
    }
  }
}

EpsMatrix::EpsMatrix() {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  data_received = 0;
  total_time = 0.0;
}

EpsMatrix::EpsMatrix(MatrixConfig config) : CBase_EpsMatrix(config) {
  GWBSE* gwbse = GWBSE::get();
  blockSize = config.tile_rows;
  numBlocks = config.chareRows();

  // Set some constants
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  total_time = 0.0;
  data_received = 0;
}

void EpsMatrix::setI(CLA_Matrix_interface mat, bool clean){
  matrix = mat;
  if (clean) {
    delete [] data;
    initialize();
  }
}

void EpsMatrix::receiveFs(Phase3Message* msg) {
  int n = 0;
  // TODO: memcpy
  for(int i=msg->start_i;i<=msg->end_i;i++)
    for(int j=msg->start_j;j<=msg->end_j;j++)
        data[IDX_eps(i,j)] = msg->data[n++];

  data_received+=n;
  if(data_received == eps_cols*eps_rows) {
    CkCallback cb(CkReductionTarget(Controller, epsilon_created), controller_proxy);
    contribute(cb);
  }
}

void EpsMatrix::multiply(double alpha, double beta) {
  matrix.multiply(alpha, beta, data, EpsMatrix::done_cb, (void*) this,
       thisIndex.x, thisIndex.y);
}

void EpsMatrix::add_compl_two() {
  int i = 0;
  complex compl_two(2.0, 0.0);
  if(thisIndex.x==thisIndex.y)
  for(int i=0;i<eps_rows;i++)
    data[IDX_eps(i,i)] += compl_two;

  CkCallback cb(CkReductionTarget(Controller, complement_multiplied), controller_proxy);
  contribute(cb);
}

void inline EpsMatrix::round_done(void) {
  CmiMemoryCheck();
  CkCallback cb(CkReductionTarget(Controller, m_multiplied), controller_proxy);
  contribute(cb);
}

void EpsMatrix::scalar_multiply(double alpha) {
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = alpha*data[IDX_eps(i,j)]; 

  CkCallback cb(CkReductionTarget(Controller, scalar_multiplied), controller_proxy);
  contribute(cb);
}

void EpsMatrix::screenedExchange() {


  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex total_contribution(0.0,0.0);
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;

  for (int k = 0; k < K; k++) {
    for (int i = 0; i < f_cache->getNSize(); i++) {
        complex contribution(0.0,0.0);
  //      for (int j = 0; j < f_cache->getNSize(); j++) { //Performs only <n|Sigma|n> as does the fortran code
  // Uncommenting above loop will perform <n|Sigma|n’>
        for (int l = 0; l < L; l++) {
          complex* fi = f_cache->getFVec(k, i, l, thisIndex.x, eps_rows);
          complex* fj = f_cache->getFVec(k, i, l, thisIndex.y, eps_cols);
          for (int r = 0; r < config.tile_rows; r++) {
            for (int c = 0; c < config.tile_cols; c++) {
              contribution += fi[r]*fj[c].conj()*data[IDX_eps(r,c)];
            }
          }
        }
        contrib_data[ik] = contribution;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution;
    }
  }

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::screenedExchangeComplete(NULL), controller_proxy));
  contribute(msg);
}

void EpsMatrix::bareExchange() {
  complex total_contribution = (0.0,0.0);
  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();

  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;
  std::vector<double> vcoulb = psi_cache->getVCoulb();
  vcoulb[0] = 0.38343474/2.0;

  if(thisIndex.x == thisIndex.y) {
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < f_cache->getNSize(); i++) {//ib = 5 to 8 actually, map to a number from 0
        complex contribution = (0.0,0.0);
        for (int l = 0; l < L; l++) {
          complex* f = f_cache->getFVec(k, i, l, thisIndex.x, eps_rows);
          for(int ii=0; ii < config.tile_rows; ii++) {
            int g = thisIndex.x*eps_rows+ii;
            if(g < vcoulb.size())
              contribution += f[ii]*f[ii].conj()*vcoulb[g];
          }
        }
        contrib_data[ik] = contribution;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution;
      }
    }
  }
  else{
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < f_cache->getNSize(); i++) {
        complex contribution = (0.0,0.0);
        tuple_reduction[ik++] =  CkReduction::tupleElement(sizeof(complex), &contribution, CkReduction::sum_double);
      }
    }
  }

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::bareExchangeComplete(NULL), controller_proxy));
  contribute(msg);
}

void EpsMatrix::coh(){

  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();

  complex* states = psi_cache->getStates();
  int psi_size = nfft[0]*nfft[1]*nfft[2];
  std::vector<int> accept_v = f_cache->getAcceptVector();
  std::vector<int> geps_x = f_cache->getGepsXVector();
  std::vector<int> geps_y = f_cache->getGepsYVector();
  std::vector<int> geps_z = f_cache->getGepsZVector();

  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;
  complex total_contribution = (0.0,0.0);

  for (int k = 0; k < K; k++) {
    complex *f;
    int epsilon_size = 0;

    for(int g=0;g<psi_size;g++)
      if(accept_v[g]) epsilon_size++;

    f = new complex[epsilon_size];

//This could probably be done once per node and cached
    int counter = 0;
    for (int i = 0; i < f_cache->getNSize(); i++){
      for (int j = 0; j < f_cache->getNSize(); j++){
        counter = 0;
        for(int g=0; g < psi_size; g++){
          if(accept_v[g]) {
            f[counter] += states[i*psi_size + g] * states[j*psi_size + g];
            counter++;
          }
        }
        if(counter != epsilon_size) {
          CkPrintf("\nWarning!!! %d != %d\n", counter, epsilon_size);
          CkAbort("\nsize is not equal expected");
        }

      }
    }

    int end_x = config.tile_rows;
    int end_y = config.tile_cols;

    int last_index = epsilon_size/eps_rows;
    if(thisIndex.x == last_index-1) end_x = epsilon_size%eps_rows;
    if(thisIndex.y == last_index-1) end_y = epsilon_size%eps_cols;

    for (int i = 0; i < f_cache->getNSize(); i++) {
      complex contribution = (0.0,0.0);
      for (int j = 0; j < f_cache->getNSize(); j++) {
        for (int r = 0; r < end_x; r++) {
          for (int c = 0; c < end_y; c++) {
            for(int g=0; g<epsilon_size; g++) {
              //apply some filters on which g's should be applied
              //if( (gppvec(1,igp) .eq. gvec(1,gidx(g)))
              //.AND.(gppvec(2,igp) .eq. gvec(2,gidx(g))) &
              //.AND. (gppvec(3,igp) .eq. gvec(3,gidx(g)) )) then
              if(geps_x[thisIndex.y*eps_cols+c]-geps_x[thisIndex.x*eps_rows+r] == geps_x[g] &&
                  geps_y[thisIndex.y*eps_cols+c]-geps_y[thisIndex.x*eps_rows+r] == geps_y[g] &&
                  geps_z[thisIndex.y*eps_cols+c]-geps_z[thisIndex.x*eps_rows+r] == geps_z[g])
              {
                contribution += f[g]*data[IDX_eps(r,c)];
              }
            }
          }
        }
      }
      contrib_data[ik] = contribution;
      tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
      ik++;
      total_contribution += contribution;
    }
  } //end of K loop

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::cohComplete(NULL), controller_proxy));
  contribute(msg);
}

void EpsMatrix::findAlpha() {
  if (config.chareCols() != 1) {
    CkAbort("findAlpha() only implemented for 1D decompositions\n");
  }
  double R = 0;
  for(int i = 0; i < config.tile_cols; i++) {
    R += abs(data[i]);
  }

  CkCallback cb(CkReductionTarget(Controller, found_alpha), controller_proxy);
  contribute(sizeof(long double), &R, CkReduction::max_double, cb);
}

void EpsMatrix::convergence_check(CProxy_EpsMatrix cproxy){
    // TODO: memcpy
    std::vector<complex> data_out(total_data);
    for(int i=0;i<total_data;i++)
      data_out[i] = data[i];

    cproxy(thisIndex.x, thisIndex.y).receiveConvCheck(data_out);
}

void EpsMatrix::receiveConvCheck(std::vector<complex> data_in) {
  double Rmax=0;  // the largest element
  double tmp;
  for(int i=0; i<total_data; i++) {
    tmp = abs(data[i] - data_in[i]);
    if( tmp > Rmax ){ Rmax = tmp; }
  }
  contribute(sizeof(complex), &Rmax, CkReduction::max_double,
      CkCallback(CkReductionTarget(Controller, converge_results), controller_proxy));
}

void EpsMatrix::createTranspose(CProxy_EpsMatrix other, bool todo) {
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
      if (todo) {
        complex tranpose = data[IDX_eps(i,j)];
        tranpose.im *= -1;
        incoming.push_back(tranpose);
      } else {
        incoming.push_back(data[IDX_eps(i,j)]);
      }
    }
  }
  if(todo) {
    other(thisIndex.y, thisIndex.x).receiveTranspose(incoming);
  } else {
    other(thisIndex.x, thisIndex.y).receiveTranspose(incoming);
  }
}

void EpsMatrix::receiveTranspose(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];

  contribute(CkCallback(CkReductionTarget(Controller, transpose_complete), controller_proxy));
}

void EpsMatrix::createCopy(CProxy_EpsMatrix other, bool todo) {
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++) {
    for(int j=0; j < config.tile_cols; j++) {
        incoming.push_back(data[IDX_eps(i,j)]);
    }
  }
    other(thisIndex.x, thisIndex.y).recvCopy(incoming);
}

void EpsMatrix::recvCopy(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];

  contribute(CkCallback(CkReductionTarget(Controller, copy_complete), controller_proxy));
}

void EpsMatrix::createConjugate(CProxy_EpsMatrix other){
  std::vector<complex> incoming;
  for(int i=0; i < config.tile_rows; i++)
    for(int j=0; j < config.tile_cols; j++){
      complex conj = data[IDX_eps(i,j)].conj();
      incoming.push_back(conj);
    }
  other(thisIndex.x, thisIndex.y).receiveConjugate(incoming);
}

void EpsMatrix::receiveConjugate(std::vector<complex> new_data) {
  unsigned n = 0;
  for(int i=0;i<config.tile_rows;i++)
    for(int j=0;j<config.tile_cols;j++)
      data[IDX_eps(i,j)] = new_data[n++];
    
  contribute(CkCallback(CkReductionTarget(Controller, conjugateComplete), controller_proxy));
}

void EpsMatrix::multiply_coulb(){
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  std::vector<double> coulb = psi_cache->getVCoulb();

  coulb[0] = 0.38343474/2.0;

  for(int i=0;i<config.tile_rows;i++){
    for(int j=0;j<config.tile_cols;j++){
      int g = thisIndex.x*config.tile_rows+i;
      int gp = thisIndex.y*config.tile_cols+j;
      if(g==gp && g<coulb.size())
        data[IDX_eps(i,j)] -= 1.0;
      if(g<coulb.size() && gp<coulb.size())
        data[IDX_eps(i,j)] *= sqrt(coulb[g])*sqrt(coulb[gp]);
    }
  }

  contribute(CkCallback(CkReductionTarget(Controller, s_ready), controller_proxy));
}

#include "eps_matrix.def.h"
