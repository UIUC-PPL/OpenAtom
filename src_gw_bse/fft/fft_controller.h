#ifndef __FFT_CONTROLLER_H__
#define __FFT_CONTROLLER_H__

#include <cstdlib>
#include "ckcomplex.h"
#include "fftw3.h"
#include "gspace.h"
#include "constant.h"

#include "fft_controller.decl.h"

class FFTController : public CBase_FFTController {
  public:
    FFTController();
    void setup_fftw_3d(int nfft[3], int direction);
    void do_fftw();

    fftw_complex* get_in_pointer() const { return in_pointer; }
    fftw_complex* get_out_pointer() const { return out_pointer; }

    void get_geps(double epsCut, double* qvec, double* b1, double* b2, double * b3, 
                              double alat, int nfft[3]);
    void calc_vcoulb(double* qvec, double* a1, double* a2, double* a3,
                      double* b1, double* b2, double * b3,
                      double shift[3],double alat, int nkpt, int iq, int *nk);
    double average_mbz(double pa,double pb,double pc,double *hmatik,
                   int Na,int Nb,int Nc,double tol);
    double compute_average_mbz(const double pa,const double pb,const double pc,
                   double *hmatik,
                   const int Na,const int Nb,const int Nc,const int n);
  private:
    void destroy_fftw_stuff();
    bool first_time;
    int old_direction;
    int old_nfft[3];
    fftw_plan plan;
    fftw_complex* in_pointer;
    fftw_complex* out_pointer;
    GSPACE *geps;
};

extern /* readonly */ CProxy_FFTController fft_controller_proxy;

#endif
