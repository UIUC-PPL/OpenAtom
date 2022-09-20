#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include <cstdlib>
#include "ckcomplex.h"

#include "fft_controller.decl.h"
#include "psi_cache.decl.h"
#include "fvector_cache.decl.h"
#include "class_gw_io.h"
#include "laplace.h"

#include <unordered_set>
#include <unordered_map>

#include "controller.decl.h"
#define MAX_ITERATIONS 100

#define PAD(size)\
  int pad = eps_rows - (epsilon_size%eps_rows);\
  padded_epsilon_size = epsilon_size + pad;

class Stopwatch {
  private:
    struct Time {
      double start, stop;
      bool started, stopped;
      Time() : start(0.0), stop(0.0), started(false), stopped(false) {}
    };
    std::unordered_map<std::string, Time> times;
  public:
    void startTimer(std::string key) {
      if (times[key].started) {
        CkPrintf("Warning: Timer already started for %s\n",key.c_str());
      }
      times[key].start = CmiWallTimer();
      times[key].started = true;
    }
    void stopTimer(std::string key) {
      if (!times[key].started) {
        CkPrintf("Warning: Timer not started for %s\n", key.c_str());
      }
      if (times[key].stopped) {
        CkPrintf("Warning: Timer already stopped for %s\n", key.c_str());
      }
      times[key].stop = CmiWallTimer();
      times[key].stopped = true;
    }
    void resetTimer(std::string key) {
      times[key].started = false;
      times[key].stopped = false;
    }
    double getEnd(std::string key) {
      if (!times[key].stopped) {
        CkPrintf("Warning: Timer not stopped for %s\n", key.c_str());
      }
      return times[key].stop;
    }
    double getElapsed(std::string key) {
      if (times[key].stopped) {
        return times[key].stop - times[key].start;
      } else if (times[key].started) {
        return CmiWallTimer() - times[key].start;
      } else {
        CkPrintf("Warning: Timer not started or stopped for %s\n", key.c_str());
        return 0.0;
      }
    }
    bool isStopped(std::string key) {
      return times[key].stopped;
    }
    bool isStarted(std::string key) {
      return times[key].stopped;
    }
};

class Controller : public CBase_Controller {
  Controller_SDAG_CODE
  public:
    Controller();
    void computeEpsDimensions();
    void calc_Geps();
    void got_geps(std::vector<int> accept, int epsilon_size);
    void got_vcoulb(std::vector<double> vcoulb_in, double vcoulb0);
  private:
    bool do_output;
    int msg_received;
    int iter, maxiter, qpts;
    int iteration;
    int qindex;
    unsigned index;
    unsigned dimension, rows;
    bool resultInsert;
    bool all_qpts;
    double epsCut;
    double tol_iter_mtxinv;
    double alat;
    double vol;
    complex bare_x, screen_x, coh;
    complex **bare_x_final, **screen_x_final, **coh_final;
    std::vector<double> vcoulb;
    double shift[3];
    unsigned K, L, M, *Q, Bands, pipeline_stages;
    unsigned next_K, next_state, total_sent, total_complete;
    unsigned max_sends, next_report_threshold;
    unsigned p_matrix_dimension, num_p_rows;
    int global_inew, global_jnew;
    int max_local_inew;
    int padded_epsilon_size;
    double prev_max;
    int *n_list, *np_list;
    std::vector<int> accept_result;
    Stopwatch stopwatch;

    CkCallback copyCB, readCB, writeCB, verifyCB, qindexCB;

    IOConfig p_config, eps_config, eps_inv_config;

    // Epsilon proxies and matrices
    CLA_Matrix_interface matA, matB, matC;
    CLA_Matrix_interface matA2, matB2, matC2;
    CLA_Matrix_interface matA3, matB3, matC3;
    CProxy_EpsMatrix eps_matrix1D_proxy;
    CProxy_EpsMatrix eps_matrix2D_m_proxy, eps_matrix2D_mT_proxy,
                     eps_matrix2D_X_proxy, eps_matrix2D_A_proxy,
                     eps_matrix2D_M1_proxy, eps_matrix2D_X1_proxy,
                     s_matrix2D_proxy;
};

// A struct containing the required info for computing a set of f vectors for a
// single unoccupied state. For each f the equation is:
// f[i] = psi_occ[i] * psi_unocc[i].conj() * scaling_factor * umklapp_factor[i]
// f, psi vectors, and umklapp_factor all have 'size' elements
// e_unocc and umklapp_factor are the same for every f
// e_occ has an entry for each f to be computed
// The set of occupied psis are all occupied psis needed for the given unocc psi
struct FComputePacket {
  unsigned int size;
  double e_unocc;
  double* e_occ;
  complex* unocc_psi;
  complex** occ_psis;
  complex* umklapp_factor;
  complex* fs;
  complex* fsave;
  std::vector<std::pair<int, int>>* regions;
};

class PsiCache : public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    void computeFs(PsiMessage*);
    void reportFTime();
    complex* getPsi(unsigned, unsigned, unsigned) const;
    complex* getF(unsigned,unsigned) const;
    void setVCoulb(std::vector<double> vcoulb_in, double vcoulb0);
    std::vector<double> getVCoulb();
    double getVCoulb0();
    complex* getStates();
    bool in_np_list(int n_index);
    int get_index(int n_index);
    void setQIndex(int q_index);
    void setRegionData(PMatrix* matrix_chare, int start_row, int start_col,
                       int tile_nrows, int tile_ncols);
    void reportInfo();
    LAPLACE* getLP();

    int elements;
    int total_elements;
    int debug;
    complex*** psis;
    std::vector<std::pair<int, int>> regions;
    double get_OccOcc(int k, int iv);
  private:
    void kqIndex(unsigned, unsigned&, int*);
    void computeUmklappFactor(int*);

    // Used for CkLoop parameters
    FComputePacket f_packet;

    unsigned K, L, M, psi_size, received_psis, qindex, pipeline_stages, received_chunks;
    // TODO: Flatten arrays?
//    complex*** psis;
    complex*** psis_shifted;
    complex* fs;
    complex *fsave;
    complex *f_nop;
    complex *states;
    complex *P_m;
    std::vector<double> vcoulb;
    complex* umklapp_factor;
    int n_np;
    int states_received;
    int *n_list, *np_list;
    double vcoulb_0;

    double total_time;

    // Used for registering fvector regions
    std::vector<PMatrix*> matrix_chares;
    int min_row, min_col, max_row, max_col;
    CmiNodeLock tile_lock;
    LAPLACE *lp;
    double*** Occ_occ;
    double*** Occ_unocc;
};

class FVectorCache : public CBase_FVectorCache {
  FVectorCache_SDAG_CODE
  public:
    FVectorCache(){ storing = true;}
    void setDim(int dim, std::vector<int> accept,
    std::vector<int> geps_X, std::vector<int> geps_Y, std::vector<int> geps_Z);
    void computeFTilde(complex* fs_in);
    void putFVec(int kpt, int n, complex* fs_input);
    complex* getFVec(int kpt, int n, int l, int start, int size);
    void applyCutoff(complex* fs_in);
    void findIndices();
    int getNSize();
    std::vector<int> getAcceptVector() { return accept_vector;}
    std::vector<int> getGepsXVector() { return geps_x;}
    std::vector<int> getGepsYVector() { return geps_y;}
    std::vector<int> getGepsZVector() { return geps_z;}
  private:
    FComputePacket f_packet;
    unsigned K, L, psi_size, fcount, n_list_size, node_count;
    int ndata, totalSize, data_size_x, data_size_y;
    int eps_chares_x, eps_chares_y, my_chare_count, my_chare_start;
    int eps_start_chare_x, eps_start_chare_y, eps_end_chare_x, eps_end_chare_y;
    int epsilon_size, padded_epsilon_size;
    int *my_eps_chare_indices_x;
    int *my_eps_chare_indices_y;
    int *data_offset_x;
    int *data_offset_y;
    int count;
    complex sum;
    complex* fs;
    std::vector<int> accept_vector;
    std::vector<int> geps_x;
    std::vector<int> geps_y;
    std::vector<int> geps_z;
    int num_chares, num_chares_x, num_chares_y, chare_factor;
    int *charesX;
    int *charesY;
    int *local_offset;
    int *global_offset;
    int num_rows, num_cols;
    int *offsets;//[num_chares*2];//already calculated
    bool storing;
};


extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;
extern /* readonly */ CProxy_FVectorCache fvector_cache_proxy;

#endif
