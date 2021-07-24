#include "controller.h"

#include "standard_include.h"
#include "allclass_gwbse.h"
#include "messages.h"
#include "eps_matrix.h"
#include "pmatrix.h"
#include "mat_mul.h"
#include "main.h"
#include "states.h"
#include "fft_controller.h"
#include "fft_routines.h"
#include "CkLoopAPI.h"
#include "limits.h"

#define eps_rows 20
#define eps_cols 20

#define DEBUG4 1

void init_plan_lock();

Controller::Controller() {
  GWBSE *gwbse = GWBSE::get();
  GW_SIGMA *gw_sigma = &(gwbse->gw_sigma);
  // Set our class variables
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;

  qpts = gwbse->gw_parallel.n_qpt;
  all_qpts = false;
  if(gwbse->gw_parallel.n_qpt >= K)
    all_qpts = true;
  else
    Q = gwbse->gw_parallel.Q;

  Bands = gw_sigma->num_sig_matels;

  bare_x_final = new complex*[K];
  screen_x_final = new complex*[K];
  coh_final = new complex*[K];

  for ( int i=0; i<K; i++) {
    bare_x_final[i] = new complex[Bands];
    screen_x_final[i] = new complex[Bands];
    coh_final[i] = new complex[Bands];
    for( int j=0; j<Bands; j++) {
      bare_x_final[i][j] = (0.0,0.0);
      screen_x_final[i][j] = (0.0,0.0);
      coh_final[i][j] = (0.0,0.0);
    }
  }
  n_list = gw_sigma->n_list_sig_matels;
  np_list = gw_sigma->np_list_sig_matels;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;

  next_K = next_state = total_sent = total_complete = next_report_threshold = 0;

  dimension = gwbse->gw_parallel.n_elems;
  rows = gwbse->gw_parallel.rows_per_chare;

  epsCut = gwbse->gw_epsilon.Ecuteps;
  tol_iter_mtxinv = gwbse->gw_epsilon.tol_iter;
  alat = 10.261200; 
  shift[0] = 0;
  shift[1] = 0;
  shift[2] = 0.001;
  maxiter = gwbse->gw_epsilon.max_iter?gwbse->gw_epsilon.max_iter:MAX_ITERATIONS;
  // TODO: Make these config options
  do_output = true;
  max_sends = M*K;  // For debugging this can be changed to a smaller number
  //maxiter = 1;
  msg_received = 0;
  global_inew = 0;
  max_local_inew = global_inew;
  global_jnew = 0;

  copyCB = CkCallback(CkReductionTarget(Controller, copyComplete), thisProxy);
  readCB = CkCallback(CkReductionTarget(Controller, readComplete), thisProxy);
  writeCB = CkCallback(CkReductionTarget(Controller, writeComplete), thisProxy);
  verifyCB = CkCallback(CkReductionTarget(Controller, verifyComplete), thisProxy);

  qindexCB = CkCallback(CkReductionTarget(Controller, setQIndex), thisProxy);

  p_config = gwbse->gw_io.p_matrix;
  eps_config = gwbse->gw_io.epsilon;
  eps_inv_config = gwbse->gw_io.epsilon_inv;
}


void Controller::computeEpsDimensions() {
  GWBSE *gwbse = GWBSE::get();

  double *this_q, *b1, *b2, *b3;
  b1 = gwbse->gwbseopts.b1;
  b2 = gwbse->gwbseopts.b2;
  b3 = gwbse->gwbseopts.b3;

  this_q = gwbse->gwbseopts.qvec[qindex];

  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;

  int ndata = nfft[0]*nfft[1]*nfft[2];
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();

  fft_controller->get_geps(epsCut, this_q, b1, b2, b3, alat, nfft);
}

void Controller::calc_Geps() {
  
  GWBSE *gwbse = GWBSE::get();

  double *this_q, *b1, *b2, *b3;
  b1 = gwbse->gwbseopts.b1;
  b2 = gwbse->gwbseopts.b2;
  b3 = gwbse->gwbseopts.b3;

  double *a1, *a2, *a3;
  a1 = gwbse->gwbseopts.a1;
  a2 = gwbse->gwbseopts.a2;
  a3 = gwbse->gwbseopts.a3;

  this_q = gwbse->gwbseopts.qvec[qindex];

  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;

  int ndata = nfft[0]*nfft[1]*nfft[2];
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();

  //output - vcoulb
  fft_controller->calc_vcoulb(this_q, a1, a2, a3, b1, b2, b3, shift, alat, gwbse->gwbseopts.nkpt, qindex, gwbse->gwbseopts.nk);
}

void Controller::got_vcoulb(std::vector<double> vcoulb_in, double vcoulb0){
  vcoulb = vcoulb_in;
  psi_cache_proxy.setVCoulb(vcoulb_in, vcoulb0);
}

PsiCache::PsiCache() {
  states_received = 0;
  GWBSE *gwbse = GWBSE::get();
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  int M = gwbse->gw_parallel.M;
  GW_SIGMA *gw_sigma = &(gwbse->gw_sigma);
  n_np = gw_sigma->num_sig_matels;
  n_list = gw_sigma->n_list_sig_matels;
  np_list = gw_sigma->np_list_sig_matels;
  qindex = 0;
  elements = 0;
  int dim = gwbse->gw_parallel.n_elems/gwbse->gw_parallel.rows_per_chare;
  total_elements = dim*dim;
  if(gwbse->gw_parallel.n_qpt < K){
    qindex = gwbse->gw_parallel.Q[0];
  }
  psi_size = gwbse->gw_parallel.n_elems;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;
  received_psis = 0;
  received_chunks = 0;
  psis = new complex**[K];
  for (int k = 0; k < K; k++) {
    psis[k] = new complex*[L+M];
    for (int l = 0; l < L+M; l++) {
      psis[k][l] = new complex[psi_size];
    }
  }
  // shifted k grid psis. Need this for qindex=0
  psis_shifted = new complex**[K];
  for (int k = 0; k < K; k++) {
    psis_shifted[k] = new complex*[L];
    for (int l = 0; l < L; l++) {
      psis_shifted[k][l] = new complex[psi_size];
    }
  }

  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;

  int ndata = nfft[0]*nfft[1]*nfft[2];

  fs = new complex[L*psi_size*pipeline_stages];
  fsave = new complex[L*psi_size];
  f_nop = new complex[L*psi_size];
  states = new complex[K*2*n_np*psi_size];
  P_m = new complex[ndata*ndata];

  umklapp_factor = new complex[psi_size];

  // Variables for chare region registration
  min_row = INT_MAX;
  min_col = INT_MAX;
  max_row = INT_MIN;
  max_col = INT_MIN;
  tile_lock = CmiCreateLock();

  lp = new LAPLACE(gwbse->gw_epsilon.Eocc, gwbse->gw_epsilon.Eunocc);

  char fromFile[200];
  Occ_occ = new double**[1];
  Occ_unocc = new double**[1];
  for (int s = 0; s < 1; s++) {
    Occ_occ[s] = new double*[K];
    Occ_unocc[s] = new double*[K];
    for (int k = 0; k < K; k++) {
      Occ_occ[s][k] = new double[L];
      Occ_unocc[s][k] = new double[M];
    }
  }
for (int s = 0; s < 1; s++) {
    for (int k = 0; k < K; k++) {
  sprintf(fromFile, "./STATES/Spin.%d_Kpt.%d_Bead.0_Temper.0/%s",s,k,"occupations.in");
  FILE* fp_occ = fopen(fromFile, "r");
  if (fp_occ == NULL) {
    PRINTF("Cannot open Occ Value File: %s\n", fromFile);
    EXIT(1);
  }
  for (int i = 0; i < L; i++) {
    fscanf(fp_occ,"%lg",&Occ_occ[s][k][i]);
  }
  for (int i = 0; i < M; i++) {
    fscanf(fp_occ,"%lg",&Occ_unocc[s][k][i]);
  }
 }
} 

  total_time = 0.0;
  contribute(CkCallback(CkReductionTarget(Controller,psiCacheReady), controller_proxy));
}
double PsiCache::get_OccOcc(int k, int iv) {
  return Occ_occ[0][k][iv];//tmp;
}

LAPLACE *PsiCache::getLP() {
  return lp;
}

void PsiCache::setQIndex(int q_index){

  qindex = q_index;
  received_psis = 0;
  received_chunks = 0;
  for (int k = 0; k < K; k++) {
    for (int l = 0; l < L; l++) {
      std::fill(psis[k][l], psis[k][l]+psi_size, 0.0);
    }
  }
  // shifted k grid psis. Need this for qindex=0
  for (int k = 0; k < K; k++) {
    for (int l = 0; l < L; l++) {
      std::fill(psis_shifted[k][l], psis_shifted[k][l]+psi_size, 0.0);
    }
  }

  std::fill(fs, fs+L*psi_size*pipeline_stages, 0.0);
  std::fill(fsave, fsave+L*psi_size, 0.0);
  std::fill(f_nop, f_nop+L*psi_size, 0.0);
  std::fill(states, states+K*2*n_np*psi_size, 0.0);

  std::fill(umklapp_factor, umklapp_factor+psi_size, 0.0);

  // Variables for chare region registration
  min_row = INT_MAX;
  min_col = INT_MAX;
  max_row = INT_MIN;
  max_col = INT_MIN;
  tile_lock = CmiCreateLock();

  total_time = 0.0;
  contribute(CkCallback(CkReductionTarget(Controller,psiCacheReady), controller_proxy));
}

void PsiCache::reportFTime() {
  CkReduction::statisticsElement stats(total_time);
  int tuple_size = 2;
  CkReduction::tupleElement tuple_reduction[] = {
    CkReduction::tupleElement(sizeof(double), &total_time, CkReduction::max_double),
    CkReduction::tupleElement(sizeof(CkReduction::statisticsElement), &stats, CkReduction::statistics) };

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::reportFTime(NULL), controller_proxy));
  contribute(msg);
}

void PsiCache::receivePsi(PsiMessage* msg) {

  if (msg->spin_index != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(msg->k_index < K);
//  CkAssert(msg->state_index < L);
  CkAssert(msg->size == psi_size);
  if(msg->shifted==false){std::copy(msg->psi, msg->psi+psi_size, psis[msg->k_index][msg->state_index]);}
  if(msg->shifted==true){std::copy(msg->psi, msg->psi+psi_size, psis[msg->k_index][msg->state_index]);}
  fflush(stdout);
  states_received++;
#if 1
  if(states_received == (L+M)*K)
    contribute(CkCallback(CkReductionTarget(Controller,cachesFilled), controller_proxy));
#endif

  delete msg;
}

/**
*
* Authors: Robert Kaucic, Eric Mikida
*
* Update chare registration,  optimize fvector computation
* Pmatrix chares now register the region of the pmatrix they correspond to
* with their local psiCache. This data is stored in an ordered, coalesced
* vector of disjoint regions. computeF now calculates fvectors only in
* the regions specified by the vector. psiCache now maintains a vector
* of pointers to chares that have registered with it.
**/

void PsiCache::setRegionData(PMatrix* matrix_chare, int start_row, int start_col, int tile_nrows, int tile_ncols) {
  CmiLock(tile_lock);

  matrix_chares.push_back(matrix_chare);

  std::pair<int, int> rowSpan = std::pair<int, int>(start_row, start_row + tile_nrows);
  std::pair<int, int> colSpan = std::pair<int, int>(start_col, start_col + tile_ncols);

  // Should only ever trigger once per cache.
  if (regions.size() == 0) {
    regions.push_back(rowSpan);
  } else { // Insert rowSpan in order, or coalesce it.
    bool rowInserted = false;
    for (int i = 0; i < regions.size(); i++) {
      std::pair<int, int> curSpan = regions[i];
      if ((rowSpan.first >= curSpan.first && rowSpan.first <= curSpan.second)
         || (rowSpan.second >= curSpan.first && rowSpan.second <= curSpan.second)
         || (rowSpan.first <= curSpan.first && rowSpan.second >= curSpan.second)) { // Reasons to coalesce.
        regions[i].first = std::min(rowSpan.first, curSpan.first);
        regions[i].second = std::max(rowSpan.second, curSpan.second);
        int j;
        for (j = i + 1; j < regions.size(); j++) {
          if (regions[i].second >= regions[j].first) { // Merge the overlapping spans into regions[i].
            regions[i].first = std::min(regions[i].first, regions[j].first);
            regions[i].second = std::max(regions[i].second, regions[j].second);
          } else {
            break;
          }
        }
        if (j > i + 1) { // If at least one additional span was coalesced.
          regions.erase(regions.begin() + i + 1, regions.begin() + j);
        }
        rowInserted = true;
        break;
      } else if (rowSpan.second < curSpan.first) { // rowSpan is disjoint with other spans. Insert.
        regions.insert(regions.begin() + i, rowSpan);
        rowInserted = true;
        break;
      }
    }
    if (!rowInserted) { // rowSpan occurs entirely after all other spans.
      regions.push_back(rowSpan);
    }
  }
  // Insert colSpan in order, or coalesce it.
  bool colInserted = false;
  for (int i = 0; i < regions.size(); i++) {
    std::pair<int, int> curSpan = regions[i];
    if ((colSpan.first >= curSpan.first && colSpan.first <= curSpan.second)
       || (colSpan.second >= curSpan.first && colSpan.second <= curSpan.second)
       || (colSpan.first <= curSpan.first && colSpan.second >= curSpan.second)) { // Reasons to coalesce.
      regions[i].first = std::min(colSpan.first, curSpan.first);
      regions[i].second = std::max(colSpan.second, curSpan.second);
      int j;
      for (j = i + 1; j < regions.size(); j++) {
        if (regions[i].second >= regions[j].first) { // Merge the overlapping spans into regions[i].
          regions[i].first = std::min(regions[i].first, regions[j].first);
          regions[i].second = std::max(regions[i].second, regions[j].second);
        } else {
          break;
        }
      }
      if(j > i + 1) { // If at least one additional span was coalesced.
        regions.erase(regions.begin() + i + 1, regions.begin() + j);
      }
      colInserted = true;
      break;
    } else if (colSpan.second < curSpan.first) { // colSpan is disjoint with other spans. Insert.
      regions.insert(regions.begin() + i, colSpan);
      colInserted = true;
      break;
    }
  }
  if (!colInserted) { // colSpan occurs entirely after all other spans.
    regions.push_back(colSpan);
  }

  CmiUnlock(tile_lock);
}

// Called by CkLoop to spread the computation of f vectors across the node
void computeF(int first, int last, void* result, int count, void* params) {
  FComputePacket* f_packet = (FComputePacket*)params;
  unsigned psi_size = f_packet->size;
  complex* psi_unocc = f_packet->unocc_psi;
  complex* umklapp_factor = f_packet->umklapp_factor;
  double* e_occ = f_packet->e_occ;
  double e_unocc = f_packet->e_unocc;
  complex* fs = f_packet->fs;
  std::vector<std::pair<int, int>>* regions = f_packet->regions;

  for (int l = first; l <= last; l++) {
    complex* f = &(fs[l*psi_size]);
    complex* psi_occ = f_packet->occ_psis[l];
    double scaling_factor = 2/sqrt(e_unocc - e_occ[l]);

    for (int i = 0; i < (*regions).size(); i++) {
      for (int j = (*regions)[i].first; j < (*regions)[i].second; j++) {
        f[j] = psi_occ[j] * psi_unocc[j].conj();
        if (umklapp_factor) {
          f[j] *= umklapp_factor[j];
        }
        #ifdef USE_LAPACK
        // BLAS calls compute the complex conjugate of P, which is hermitian. This
        // change to f corrects that so we get the correct P.
        f[j] = f[j].conj();
        #endif
        f[j] *= scaling_factor;
      }
    }
  }
}

// Called by CkLoop to spread the computation of f vectors across the node
void computeF_all(int first, int last, void* result, int count, void* params) {
  FComputePacket* f_packet = (FComputePacket*)params;
  unsigned psi_size = f_packet->size;
  complex* psi_unocc = f_packet->unocc_psi;
  complex* umklapp_factor = f_packet->umklapp_factor;
  double* e_occ = f_packet->e_occ;
  double e_unocc = f_packet->e_unocc;
  complex* fs = f_packet->fs;
  complex* fsave = f_packet->fsave;

  for (int l = first; l <= last; l++) {
    complex* f = &(fs[l*psi_size]);
    complex* fsave_subset = &(fsave[l*psi_size]);
    complex* psi_occ = f_packet->occ_psis[l];
    double scaling_factor = 2/sqrt(e_unocc - e_occ[l]);

    for (int i = 0; i < psi_size; i++) {
      f[i] = psi_occ[i] * psi_unocc[i].conj();
      if (umklapp_factor) {
        f[i] *= umklapp_factor[i];
      }
#ifdef USE_LAPACK
      // BLAS calls compute the complex conjugate of P, which is hermitian. This
      // change to f corrects that so we get the correct P.
      f[i] = f[i].conj();
#endif
      fsave_subset[i] = f[i];
      f[i] *= scaling_factor;
    }
  }
}

// Receive an unoccupied psi, and split off the computation of all associated f
// vectors across the node using CkLoop.

bool PsiCache::in_np_list(int n_index){
  for(int i=0;i<n_np;i++)
    if(n_index+1==np_list[i]) return true;
  return false;
}

int PsiCache::get_index(int n_index){
  for(int i=0;i<n_np;i++){
    if(n_index+1==np_list[i]) return i;
  }
  return -1;
}

void PsiCache::computeFs(PsiMessage* msg) {
  double start = CmiWallTimer();

  if (msg->spin_index != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(msg->size == psi_size);

  // Compute ikq index and the associated umklapp factor
  // TODO: This should just be a table lookup
  unsigned ikq;
  int umklapp[3];
  kqIndex(msg->k_index, ikq, umklapp);

  bool uproc = false;
  if (umklapp[0] != 0 || umklapp[1] != 0 || umklapp[2] != 0) {
    uproc = true;
    computeUmklappFactor(umklapp);
  }

  GWBSE* gwbse = GWBSE::get();
  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_occ_shifted = gwbse->gw_epsilon.Eocc_shifted;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;

#define TESTING 1
#ifdef TESTING
  if(in_np_list(msg->state_index) && msg->shifted==false){
//Cache this
    int state_index = get_index(msg->state_index)*2*psi_size;
    complex *store_x = &states[(msg->k_index*2*n_np*psi_size)+ state_index];
    complex *load_x = msg->psi;
    memcpy(store_x, load_x, psi_size*sizeof(complex));
  }
#endif

  // Create the FComputePacket for this set of f vectors and start CkLoop
  f_packet.size = psi_size;
  f_packet.unocc_psi = msg->psi;
  f_packet.regions = &regions;
  f_packet.fsave = fsave;

  if(msg->state_index >= L){
    if ( qindex == 0 ) {
      f_packet.occ_psis = psis_shifted[ikq];
      f_packet.e_occ = e_occ_shifted[msg->spin_index][ikq];
    }
    else {
      f_packet.occ_psis = psis[ikq];
      f_packet.e_occ = e_occ[msg->spin_index][ikq];
    }

    if(!msg->sigma) {
      f_packet.e_unocc = e_unocc[msg->spin_index][msg->k_index][msg->state_index-L];
      f_packet.fs = fs + (L*psi_size*(received_chunks%pipeline_stages));

      if (uproc) { f_packet.umklapp_factor = umklapp_factor; }
      else { f_packet.umklapp_factor = NULL; }

  #ifdef USE_CKLOOP
      CkLoop_Parallelize(computeF, 1, &f_packet, L, 0, L - 1);
  #else
      for (int l = 0; l < L; l++) {
        computeF(l,l,NULL,1,&f_packet);
      }
  #endif
      received_chunks++;
    }
  }

#ifdef TESTING
  if(in_np_list(msg->state_index) && msg->sigma)
  {
    if(msg->state_index < L)
      f_packet.e_unocc = f_packet.e_occ[msg->state_index];
    else
      f_packet.e_unocc = e_unocc[msg->spin_index][msg->k_index][msg->state_index-L];

  // Ignore shifted states(q=0) for fvectors, when caching for sigma calc
    f_packet.occ_psis = psis[ikq];
    f_packet.e_occ = e_occ[msg->spin_index][ikq];
    f_packet.fs = f_nop;

    if (uproc) { f_packet.umklapp_factor = umklapp_factor; }
    else { f_packet.umklapp_factor = NULL; }

#ifdef USE_CKLOOP
    CkLoop_Parallelize(computeF_all, 1, &f_packet, L, 0, L - 1);
#else
    for (int l = 0; l < L; l++) {
      computeF_all(l,l,NULL,1,&f_packet);
    }
#endif

    FVectorCache *fvec_cache = fvector_cache_proxy.ckLocalBranch();
    fvec_cache->computeFTilde(fsave);
  //compute ftilde first - similar to ckloop above for all L's
    fvec_cache->applyCutoff(fsave);
    fvec_cache->putFVec(msg->k_index, get_index(msg->state_index), fsave);
  }
#endif

  // Let the matrix chares know that the f vectors are ready
  CkCallback cb;
  if(!msg->sigma)
    cb = CkCallback(CkReductionTarget(PMatrix, applyFs), pmatrix2D_proxy);
  else
    cb = CkCallback(CkReductionTarget(Controller,prepare_epsilon), controller_proxy);

  contribute(cb);

  // Cleanup
  delete msg;
  total_time += CmiWallTimer() - start;
}

complex* PsiCache::getPsi(unsigned ispin, unsigned ikpt, unsigned istate) const {
  if (ispin != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(ikpt >= 0 && ikpt < K);
  CkAssert(istate >= 0 && istate < L);
  return psis[ikpt][istate];
}

complex* PsiCache::getF(unsigned idx, unsigned req_no) const {
  CkAssert(idx >= 0 && idx < L);
  CkAssert(req_no < received_chunks && req_no >= received_chunks - pipeline_stages);
  return &(fs[idx*psi_size+(L*psi_size*(req_no%pipeline_stages))]);
}

complex* PsiCache::getStates() {
  return states;
}

void PsiCache::setVCoulb(std::vector<double> vcoulb_in, double vcoulb0){
  vcoulb = vcoulb_in;
  vcoulb_0 = vcoulb0;
  contribute(CkCallback(CkReductionTarget(Controller,prepare_epsilon), controller_proxy));
}

std::vector<double> PsiCache::getVCoulb() {
  return vcoulb;
}

double PsiCache::getVCoulb0() {
  return vcoulb_0;
}

// TODO: improve this to only be called when REGISTER_REGIONS is active
// at the moment "#ifdef REGISTER REGIONS" doesn't work in controller.ci
void PsiCache::reportInfo() {
#ifdef DEBUG
  if(min_row != -1 && max_row != -1 && min_col != -1 && max_col != -1) {
    CkPrintf("[PSICACHE %i]: minRow: %d, maxRow: %d, minCol: %d, maxCol: %d\n", CkMyNode(), min_row, max_row, min_col, max_col);
  }
#endif
}

void PsiCache::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
  GWBSE* gwbse = GWBSE::get();

  // temporary space to save k/q/k+q vectors
  double *this_k, *this_q;
  double k_plus_q[3], k_plus_q_orig[3];
  this_k = gwbse->gwbseopts.kvec[ikpt];
  if(qindex >= K) CkAbort("Q Index is greater than K, please provide larger K or smaller Q Index");
  this_q = gwbse->gwbseopts.qvec[qindex];

  for (int i=0; i<3; i++) {
    // calculate k+q vector 
    k_plus_q[i] = this_k[i] + this_q[i]; // k+q vector
    k_plus_q_orig[i] = k_plus_q[i]; // save it for Umklapp process
    // if not 0 =< k+q [i] <1, adjust k+q so that k+q[i] is in the Brillouine zone 
    if ( k_plus_q[i] >= 1 ) {
      k_plus_q[i] -= 1;
    }
    else if( k_plus_q[i] < 0 ){
      k_plus_q[i] += 1;
    }
  }
    
  // find k+q vector index
  for (int kk=0; kk < gwbse->gwbseopts.nkpt; kk++) {
    bool match = true;
    this_k = gwbse->gwbseopts.kvec[kk];
    //this_k is now a difference between k and k+q
    for (int i=0; i<3; i++) {
      if (this_k[i] != k_plus_q[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      ikq = kk;
      break;
    }
  }
  // save umklapp scattering information
  for (int i=0; i<3; i++) {
    uklapp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
  }

}


void PsiCache::computeUmklappFactor(int uklpp[3]){

  if (uklpp[0]==0 && uklpp[1]==0 && uklpp[2]==0){
    // do nothing
  }
  else{
    GWBSE *gwbse = GWBSE::get();
    int* nfft;
    nfft = gwbse->gw_parallel.fft_nelems;
    double *a1, *a2, *a3, *b1, *b2, *b3;
    a1 = gwbse->gwbseopts.a1;
    a2 = gwbse->gwbseopts.a2;
    a3 = gwbse->gwbseopts.a3;
    b1 = gwbse->gwbseopts.b1;
    b2 = gwbse->gwbseopts.b2;
    b3 = gwbse->gwbseopts.b3;
    double lattconst = gwbse->gwbseopts.latt;

    double rijk, G0, phase;
    unsigned counter = 0;
    for(int i=0; i<nfft[0]; i++){
      for(int j=0; j<nfft[1]; j++){
        for(int k=0; k<nfft[2]; k++){
          phase = 0;
          for (int l=0; l<3; l++){
            rijk = a1[l]*i/nfft[0] + a2[l]*j/nfft[1] + a3[l]*k/nfft[2];
            G0 = b1[l]*uklpp[0] + b2[l]*uklpp[1] + b3[l]*uklpp[2];
            G0 *= -2*M_PI/lattconst;
            phase += rijk*G0;
          }
          umklapp_factor[counter].re = cos(phase);
          umklapp_factor[counter].im = sin(phase);
          counter += 1;
        }// end k loop
      }// end j loop
    }// end i loop
  }//end if-else statement

}//end function

void FVectorCache::findIndices(){

  int count = 0;
  for(int i=0;i<eps_chares_x;i++){
    for(int j=0;j<eps_chares_y;j++){
      count++;
      if(count == my_chare_start+1){
        eps_start_chare_x = i;
        eps_start_chare_y = j;
      }
      if(count == my_chare_start+my_chare_count){
        eps_end_chare_x = i;
        eps_end_chare_y = j;
        return;
      }
    }
  }

  return;
}

int FVectorCache::getNSize(){
  return n_list_size;
}

void FVectorCache::setDim(int dim, std::vector<int> accept,
  std::vector<int> geps_X, std::vector<int> geps_Y, std::vector<int> geps_Z) {
  eps_chares_x = eps_chares_y = dim;
  accept_vector = accept;
  geps_x = geps_X;
  geps_y = geps_Y;
  geps_z = geps_Z;
  epsilon_size = 0;
  for(int i=0;i<accept_vector.size();i++)
    if(accept_vector[i])
      epsilon_size++;
  PAD(epsilon_size);

  if(CkMyNode() >= eps_chares_x*eps_chares_y){
    storing = false;
    contribute(CkCallback(CkReductionTarget(Controller,fCacheReady), controller_proxy));
    return;
  }
  totalSize = 0;
  GWBSE *gwbse = GWBSE::get();
  L = gwbse->gw_parallel.L;
  K = gwbse->gw_parallel.K;
  n_list_size = gwbse->gw_sigma.num_sig_matels;
  int total_eps_chares = eps_chares_x*eps_chares_y;
  int node_count = CkNumNodes();
  if(CkNumNodes() > eps_chares_x*eps_chares_y) node_count = eps_chares_x*eps_chares_y;
  my_chare_count = total_eps_chares/node_count;

  my_chare_start = CkMyNode()*my_chare_count;
  int remaining = total_eps_chares%node_count;

  if(CkMyNode()>0)
    my_chare_start += remaining;

  if(CkMyNode()==0)
    my_chare_count += remaining;

  my_eps_chare_indices_x = new int[my_chare_count];
  my_eps_chare_indices_y = new int[my_chare_count];

  findIndices();
  int count = 0;
  for(int i=eps_start_chare_x;i<=eps_end_chare_x;i++){
    int j = 0;
    if(i==eps_start_chare_x)
      j = eps_start_chare_y;
    int j_end = eps_chares_y-1;
    if(i==eps_end_chare_x)
      j_end = eps_end_chare_y;
    while(j<=j_end){
      my_eps_chare_indices_x[count] = i;
      my_eps_chare_indices_y[count++] = j;
      j++;
    }
  }

  ndata = padded_epsilon_size;
  psi_size = gwbse->gw_parallel.n_elems;
  data_size_x = ndata/eps_chares_x;
  if(ndata%eps_chares_x > 0)
    data_size_x += 2;
  data_size_y = ndata/eps_chares_y;
    if(ndata%eps_chares_y > 0)
      data_size_y += 2;
  data_offset_x = new int[my_chare_count];
  data_offset_y = new int[my_chare_count];

  for(int i=0;i<my_chare_count;i++){
    data_offset_x[i] = my_eps_chare_indices_x[i]*data_size_x;
    data_offset_y[i] = my_eps_chare_indices_y[i]*data_size_y;
  }

  int size_x = data_size_x;
  int size_y = data_size_y;
  local_offset =  new int[my_chare_count*2];
  global_offset = new int[my_chare_count*2];
  for(int i=0;i<my_chare_count;i++){
    global_offset[2*i] = data_offset_x[i];
    local_offset[2*i] = totalSize;
    totalSize += size_x;

    global_offset[2*i+1] = data_offset_y[i];
    local_offset[2*i+1] = totalSize;
    totalSize += size_y;
  }

  fs = new complex[K*n_list_size*L*totalSize];

  contribute(CkCallback(CkReductionTarget(Controller,fCacheReady), controller_proxy));
}

void FVectorCache::putFVec(int kpt, int n, complex* fs_input){ //fs_input has all L's corresponding to n
 if(!storing) return;
 for(int i=0;i<my_chare_count;i++){
    for(int l=0;l<L;l++){
      int global_x = global_offset[2*i];
      global_x += l*ndata;
      int local_x = local_offset[2*i];
      local_x += kpt*n_list_size*L*totalSize + n*L*totalSize + l*totalSize;

      complex *store_x = &fs[local_x];
      complex *load_x = &fs_input[global_x];
      memcpy(store_x, load_x, data_size_x*sizeof(complex));

      int global_y = global_offset[2*i+1];
      global_y += l*ndata;
      int local_y = local_offset[2*i+1];
      local_y += kpt*n_list_size*L*totalSize + n*L*totalSize + l*totalSize;

      complex *store_y = &fs[local_y];
      complex *load_y = &fs_input[global_y];
      memcpy(store_y, load_y, data_size_y*sizeof(complex));

#ifdef DEBUG_ph4
      CkPrintf("\n[Node-%d], (L=%d) Storing global(%d,%d) at local(%d,%d)\n", CkMyNode(), l, global_x, global_y, local_x, local_y);
#endif
    }
  }
}

complex* FVectorCache::getFVec(int kpt, int n, int l, int chare_start_index, int size){
  if(!storing) return NULL;
  for(int i=0;i<my_chare_count;i++){
    if(chare_start_index == my_eps_chare_indices_x[i]){
      int local_x = local_offset[2*i];
      local_x += kpt*n_list_size*L*totalSize + n*L*totalSize + l*totalSize;
      complex *f = &fs[local_x];
      return f;
    }
    if(chare_start_index == my_eps_chare_indices_y[i]){
      int local_y = local_offset[2*i+1];
      local_y += kpt*n_list_size*L*totalSize + n*L*totalSize + l*totalSize;
      complex *f = &fs[local_y];
      return f;
    }
  }
  return NULL;
}


// Called by CkLoop to spread the computation of f vectors across the node
void fTildeWorkUnit(int first, int last, void* result, int count, void* params) {

  FComputePacket* f_packet = (FComputePacket*)params;
  complex* fs = f_packet->fs;
  GWBSE *gwbse = GWBSE::get();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  int psi_size = nfft[0]*nfft[1]*nfft[2];
  int direction = 1;
  int L = gwbse->gw_parallel.L;
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();


  for (int i=0; i < L; i++){ //for all the L*n_list_size, L are computed in this node
    // First set up the data structures in the FFTController
    fft_controller->setup_fftw_3d(nfft, direction);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();
    // Pack our data, do the fft, then get the output
    put_into_fftbox(nfft, &fs[i*psi_size], in_pointer);
    fft_controller->do_fftw();
    fftbox_to_array(psi_size, out_pointer, &fs[i*psi_size], -1.0); //Now cached on the same partitions
    // replace f_vector to f_tilde_vector
  }
}

//Each node calculates its own ftilde
void FVectorCache::computeFTilde(complex *fs_in){

  // Create the FComputePacket for this set of f vectors and start CkLoop
  f_packet.size = ndata;
  f_packet.fs = fs_in;
  
  
#if 0//ifdef USE_CKLOOP
  CkLoop_Parallelize(fTildeWorkUnit, 1, &f_packet, n_list_size, 0, n_list_size - 1);
#else
    fTildeWorkUnit(0,0,NULL,1,&f_packet);
#endif
}

void FVectorCache::applyCutoff(complex* fs){
  if(!storing) return;
  int count = 0;

  for(int l=0;l<L;l++){
    complex* fk = &(fs[l*psi_size]);

    for(int i=0;i<psi_size;i++)
      if(accept_vector[i]) fs[count++] = fk[i];

    for(int i=epsilon_size;i<padded_epsilon_size;i++)
      fs[count++] = complex(0.0,0.0);
  }
}

#include "psi_cache.def.h"
#include "fvector_cache.def.h"
#include "fft_controller.def.h"
#include "controller.def.h"
