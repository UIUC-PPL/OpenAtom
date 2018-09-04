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
  qindex = 0; // The controller will set this
  CkPrintf("\nWarning!! The controller should set this q-point!\n");

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
  qindex = config.qindex; // The controller sets this

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
        contrib_data[ik] = contribution * -1.0;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution * -1.0;
    }
  }

  tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &total_contribution, CkReduction::sum_double);

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::screenedExchangeComplete(NULL), controller_proxy));
  contribute(msg);
  delete[] contrib_data;
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

  if(qindex==0)
    vcoulb[0] = psi_cache->getVCoulb0();

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
        contrib_data[ik] = contribution * -1.0;
        tuple_reduction[ik] =  CkReduction::tupleElement(sizeof(complex), &(contrib_data[ik]), CkReduction::sum_double);
        ik++;
        total_contribution += contribution * -1.0;
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
  delete[] contrib_data;
}

void EpsMatrix::coh(){

  FVectorCache* f_cache = fvector_cache_proxy.ckLocalBranch();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();

  complex* states = psi_cache->getStates();
  int psi_size = nfft[0]*nfft[1]*nfft[2];
  std::vector<int> accept_v = f_cache->getAcceptVector();

  int ga[psi_size];
  int gb[psi_size];
  int gc[psi_size];
  fftidx_to_gidx(ga,gb,gc,nfft);

  int n = f_cache->getNSize();
  int tuple_size = K*n;
  tuple_size += 1;
  CkReduction::tupleElement *tuple_reduction;
  tuple_reduction = new CkReduction::tupleElement[tuple_size];
  complex *contrib_data;
  contrib_data = new complex[tuple_size];
  int ik = 0;
  complex total_contribution = (0.0,0.0);

  GWBSE *gwbse = GWBSE::get();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  GW_SIGMA *gw_sigma = &(gwbse->gw_sigma);
  int n_np = gw_sigma->num_sig_matels;
  int *n_list = gw_sigma->n_list_sig_matels;
  int *np_list = gw_sigma->np_list_sig_matels;

  complex *f = new complex[n_np*psi_size];
  std::vector<int> map(psi_size);
  for (int k = 0; k < K; k++) {
    int epsilon_size = 0;

    for(int g=0;g<psi_size;g++){
      if(accept_v[g]){
        map[epsilon_size] = g;
        epsilon_size++;
      }
    }
    map.resize(epsilon_size);

    int base_index = k*2*n*psi_size;


//This could probably be done once per node and cached
    for (int i = 0; i < f_cache->getNSize(); i++){
      int i_index = n_list[i]-1;
      int j_index = np_list[i]-1;

      int state_index = i*2*psi_size;
      int f_base = i*psi_size;

      for(int g=0; g < psi_size; g++){
        f[f_base+g] = states[base_index + state_index + g].conj() * states[base_index + state_index + g];
      }


      fft_controller->setup_fftw_3d(nfft, -1);
      fftw_complex* in_pointer = fft_controller->get_in_pointer();
      fftw_complex* out_pointer = fft_controller->get_out_pointer();
      // Pack our data, do the fft, then get the output
      put_into_fftbox(nfft, &f[f_base], in_pointer);
      fft_controller->do_fftw();
      fftbox_to_array(psi_size, out_pointer, &f[f_base], 1);
    }

    for (int i = 0; i < f_cache->getNSize(); i++) {
      int f_base = i*psi_size;
      complex contribution = (0.0,0.0);
      for (int r = 0; r < config.tile_rows; r++) {
        int g1 = thisIndex.x*eps_rows+r;
        for (int c = 0; c < config.tile_cols; c++) {
          int g2 = thisIndex.y*eps_cols+c;
          if(g1>=epsilon_size || g2>=epsilon_size) continue;

          int gdiff[3];
          gdiff[0] = ga[map[g1]]-ga[map[g2]];
          gdiff[1] = gb[map[g1]]-gb[map[g2]];
          gdiff[2] = gc[map[g1]]-gc[map[g2]];
          // flip the value and
          // set back to gdiff values

          for (int ii=0; ii<3; ii++){
            if (gdiff[ii] < -nfft[ii]/2){
              gdiff[ii] += nfft[ii];
            }
            if (gdiff[ii] >= nfft[ii]/2){
              gdiff[ii] -= nfft[ii]/2;
            }
          }

          int gdiffIndex = -1;
          for (int ii=0; ii<psi_size; ii++){
            if (gdiff[0]==ga[ii] && gdiff[1]==gb[ii] && gdiff[2]==gc[ii]){
              gdiffIndex = ii;
              break;
            }
          }

          contribution += f[f_base+gdiffIndex]*data[IDX_eps(r,c)];
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
  delete[] contrib_data;
  delete[] f;
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
  if(qindex==0)
    coulb[0] = psi_cache->getVCoulb0();

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
