#include <cstdlib>
#include "ckcomplex.h"
#include "fftw3.h"

void set_radix(int nrad_in,int *nrad_ret, int *krad);
void set_fftsize(int, int, int*, int*, int*, int*);

void gidx_to_fftidx(int, int**, int [3], int**);
void put_into_fftbox(int, complex*, int**, int [3], fftw_complex*, bool);
void put_into_fftbox(int [3], complex*, fftw_complex*);
void fftbox_to_array(int, fftw_complex*, complex*, double);
void fftidx_to_gidx(int*, int*, int*, int[3]);
