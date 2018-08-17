#include "fft_controller.h"
#include "fft_routines.h"
#include "controller.h"

#define DEBUG(x) /*CkPrintf x*/
#define FLUSHDEBUG() /*fflush(stdout);*/

static CmiNodeLock fft_plan_lock;
void init_plan_lock() {
  fft_plan_lock = CmiCreateLock();
}

FFTController::FFTController() {
  first_time = true;
  in_pointer = out_pointer = NULL;

  geps = new GSPACE();
  

  // TODO: A group dependency could probably solve this better
  contribute(CkCallback(CkReductionTarget(Controller, fftControllerReady), controller_proxy));
}

void FFTController::do_fftw() {
  fftw_execute(plan);
}

double calc_vol(double* a1, double* a2, double* a3){
  double a[3][3];
  double m[3][3];
  double vol;

  for (int i=0; i<3; i++){
      m[0][i] = a1[i];
      m[1][i] = a2[i];
      m[2][i] = a3[i];
  }

  /* compute matrix of cofactors */
  a[0][0] =  m[1][1]*m[2][2] - m[1][2]*m[2][1];
  a[1][0] = -m[1][0]*m[2][2] + m[1][2]*m[2][0];
  a[2][0] =  m[1][0]*m[2][1] - m[1][1]*m[2][0];
  a[0][1] = -m[0][1]*m[2][2] + m[0][2]*m[2][1];
  a[1][1] =  m[0][0]*m[2][2] - m[0][2]*m[2][0];
  a[2][1] = -m[0][0]*m[2][1] + m[0][1]*m[2][0];
  a[0][2] =  m[0][1]*m[1][2] - m[0][2]*m[1][1];
  a[1][2] = -m[0][0]*m[1][2] + m[0][2]*m[1][0];
  a[2][2] =  m[0][0]*m[1][1] - m[0][1]*m[1][0];

  vol = m[0][0]*a[0][0] + m[0][1]*a[1][0] + m[0][2]*a[2][0];
  return vol;
}

double FFTController::compute_average_mbz(const double pa,const double pb,const double pc,
                   double *hmatik,
                   const int Na,const int Nb,const int Nc,const int n)
{
  // make sure n is even
  if (n%2 != 0 or n < 2 ) {
    CkPrintf("\naverage_integral() : n=%d is not even or less than 2!!\n\n",n);
    CkExit(1);
  }

  // do the integral of size n x n x n
  const double twopi = 2.0*M_PI;
  const double fourpi = 4.0*M_PI;
  double avg = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        // how much grid point (i,j,k) deviates from (pa,pb,pc)
        double dai = (-0.5 + (2.0*(double)i+1.0)/(2.0*(double)n))/(double)Na;
        double dbj = (-0.5 + (2.0*(double)j+1.0)/(2.0*(double)n))/(double)Nb;
        double dck = (-0.5 + (2.0*(double)k+1.0)/(2.0*(double)n))/(double)Nc;
        // add deviation to (pa,pb,pc) to get current point
        double pai = pa + dai;
        double pbj = pb + dbj; // PB&J!
        double pck = pc + dck;
        // convert lattice (pa,pb,pc) to physical (px,py,pz) vectorx
        double px = (pai*hmatik[1] + pbj*hmatik[2] + pck*hmatik[3]) * twopi;
        double py = (pai*hmatik[4] + pbj*hmatik[5] + pck*hmatik[6]) * twopi;
        double pz = (pai*hmatik[7] + pbj*hmatik[8] + pck*hmatik[9]) * twopi;
        // squared length of (px,py,pz)
        double plen  = sqrt( px*px + py*py + pz*pz );
        // accumulate avearge of 4pi/|p|^2
        avg += fourpi / ( plen*plen );
      } // k loop
    } // jloop
  } // i loop

  // divide by number of grid points to get final average
  avg /= (double)(n*n*n);
  return avg;
}

double FFTController::average_mbz(double pa,double pb,double pc,double *hmatik,
                   int Na,int Nb,int Nc,double tol)
{
  DEBUG(("pa=%g pb=%g pc=%g tol=%g : n avg\n",pa,pb,pc,tol));
  // we are going to make a table of computed results as we go
  // actually, only the last two entries are useful...
  const int ntable=10;
  double table[ntable][2];
  int tabrow = 0;
  bool gotbelowtol = false;
  // current average and old average from previous cycle
  double avg = 0.0 , avgold = 0.0;
  for (int n = 50 ; n <= 50 *ntable; n += 50 ) {
    // compute current average
    avgold = avg;
    avg = compute_average_mbz(pa,pb,pc,hmatik,Na,Nb,Nc,n);
    // stuff latest result in table
    table[tabrow][0] = (double)n;
    table[tabrow][1] = avg;
    tabrow++;
    DEBUG(("%3d %.10f\n",n,avg));
    FLUSHDEBUG();
    // are we below tolerance?  Then we are done!
    if ( fabs( (avg-avgold)/avg ) < tol ) {
      gotbelowtol = true;
      DEBUG(("Got below tolerance!  Answer is %.10f\n",avg));
      return avg;
    }
  }

  // if we got here, we did not get below tolerance
  // this is generally the case when pa==pb==pc==0 and the integral
  // is *very* hard to do brute force.  But then it behaves very
  // nicely versus 1/n in that it goes as const1 + const2/n for large n.
  // we will be using this fact to fit it to such a form.
  DEBUG(("Did not get below tolerance...\n"));
  DEBUG(("Last 2 tabulated entries are\n"));
  for (int tabrow = ntable-2; tabrow <= ntable-1; tabrow++)
    DEBUG(("%3d  %.10f\n",(int)table[tabrow][0],table[tabrow][1]));

  // do a linear fit to the last 2 points in the table assuming
  // the asymptotic form is f(n) = a/n + b.  We only care about b
  // which is the extrapolated integral at n=infinity.
  // if f(n1) and f(n2) are known then b=(f(n1)*n1-f(n2)*n2)/(n1-n2).
  double n1 = table[ntable-2][0]; double f1 = table[ntable-2][1];
  double n2 = table[ntable-1][0]; double f2 = table[ntable-1][1];
  double intercept = (f1*n1-f2*n2)/(n1-n2);
  DEBUG(("Linear fit to a/n + b gives b=%.10f\n",intercept));

  return intercept;
}

void FFTController:: calc_vcoulb(double* qvec, double* a1, double* a2, double* a3,
                                 double* b1, double* b2, double * b3, double shift[3],
                                 double alat, int nkpt, int iq, int *nk){

  double* vcoulb;
  vcoulb = new double [geps->ng];
  double gx, gy, gz;
  double gq[3];
  double vol = calc_vol(a1, a2, a3);
  double fact = 4*PI/vol/nkpt;
  double vcoulb0 = 0;

  for (int i=0; i<geps->ng; i++) {
      if (iq==0) {
          if(i==0){
            double hmatik[10];
            hmatik[1] = b1[0]; hmatik[2] = b2[0]; hmatik[3] = b3[0];
            hmatik[4] = b1[1]; hmatik[5] = b2[1]; hmatik[6] = b3[1];
            hmatik[7] = b1[2]; hmatik[8] = b2[2]; hmatik[9] = b3[2];
            for (int k=1; k<10; k++){ hmatik[k] /= alat; }

            double tol = 1.0e-6;
            double avg = average_mbz(geps->ig[i], geps->jg[i], geps->kg[i], hmatik, nk[0], nk[1], nk[2], tol);
            avg /= (vol*nkpt);
            vcoulb0 = avg;
          }

          gx = geps->ig[i] + shift[0];
          gy = geps->jg[i] + shift[1];
          gz = geps->kg[i] + shift[2];
      }
      else {
          gx = geps->ig[i] + qvec[0];
          gy = geps->jg[i] + qvec[1];
          gz = geps->kg[i] + qvec[2];
      }

      vcoulb[i] = 0;
      for (int j=0; j<3; j++) {
          gq[j] =  gx*b1[j] + gy*b2[j] + gz*b3[j];
          gq[j] *= 2*PI/alat;

          vcoulb[i] += gq[j]*gq[j];
      }
      vcoulb[i] = 1/vcoulb[i];
      vcoulb[i] *= fact;
  }

  std::vector<double> vcoulb_v;
  vcoulb_v.resize(geps->ng);
  for(int i=0;i<geps->ng;i++)
    vcoulb_v[i] = vcoulb[i];
  controller_proxy.got_vcoulb(vcoulb_v, vcoulb0);
}

void FFTController::get_geps(double epsCut, double* qvec, double* b1, double* b2, double * b3, 
                              double alat, int nfft[3]){

  int ndata = nfft[0]*nfft[1]*nfft[2];
  bool accept[ndata];
  int *gx, *gy, *gz;
    
  gx = new int [ndata];
  gy = new int [ndata];
  gz = new int [ndata];

  fftidx_to_gidx(gx, gy, gz, nfft);

//Values would need to be sent to Pmatrix geps

   double gxtmp, gytmp, gztmp;
    double vtmp[3];
    double Ekin;
    int eps_size = 0;
      
    for (int i=0; i<ndata; i++) { //can't be 0?
        gxtmp = gx[i] + qvec[0];
        gytmp = gy[i] + qvec[1];
        gztmp = gz[i] + qvec[2];//iq was removed assuming 0 index - might be wrong, since we have only one node
        /* transfer to cartesian unit to calculate energy */
        Ekin = 0;
        for (int j=0; j<3; j++) {
            vtmp[j] = gxtmp*b1[j] + gytmp*b2[j] + gztmp*b3[j];
            vtmp[j] *= 2*PI/alat;
            Ekin += 0.5 * vtmp[j] * vtmp[j];
        }
        
        if (Ekin <= epsCut) {
            accept[i] = true;
            eps_size += 1;
        }
        else{
            accept[i] = false;
        }

    }

    CkPrintf("[FFT CONTROLLER] Dimension of epsilon matrix = %d\n", eps_size);
    // set values
    geps->ng = eps_size;
    geps->ig = new int [eps_size];
    geps->jg = new int [eps_size];
    geps->kg = new int [eps_size];
   
    int j=0;

    for (int i=0; i<ndata; i++) {
        if (accept[i]) {
            geps->ig[j] = gx[i];
            geps->jg[j] = gy[i];
            geps->kg[j] = gz[i];
            j += 1;
        }
    }
   
    if ( j!= eps_size ) {
        CkPrintf(" Oops. Error when reducing gspace!!!");
    }

    geps->ig_diff = new int [eps_size*eps_size];
    geps->jg_diff = new int [eps_size*eps_size];
    geps->kg_diff = new int [eps_size*eps_size];

    int ngdata = eps_size;
    int k=0;

#ifdef DEBUG_2
    for (int i=0; i<ngdata; i++)
      CkPrintf("\ninit_geps->ig[%d] = %d\n", i+1, geps->ig[i]);
#endif
    for (int i=0; i<ngdata; i++) {
      for (int j=0; j<ngdata; j++) {
          geps->ig_diff[k] = geps->ig[j] - geps->ig[i];
          geps->jg_diff[k] = geps->jg[j] - geps->jg[i];
          geps->kg_diff[k] = geps->kg[j] - geps->kg[i];
//          CkPrintf("\ngeps->ig_diff[%d] = %d\n", k+1, geps->ig_diff[k]);
          k += 1;
      }
    }
   
    delete[] gx;
    delete[] gy;
    delete[] gz;

    std::vector<int> accept_v;
    std::vector<int> geps_x(eps_size);
    std::vector<int> geps_y(eps_size);
    std::vector<int> geps_z(eps_size);
    accept_v.resize(ndata);

    int counter = 0;
    for(int i=0;i<ndata;i++){
      if(accept[i]){
        accept_v[i] = 1;

        geps_x[counter] = geps->ig[counter];
        geps_y[counter] = geps->jg[counter];
        geps_z[counter] = geps->kg[counter];
        counter++;
      }
      else{
        accept_v[i] = 0;
      }
    }

    controller_proxy.receiveEpsDimensions(accept_v, geps_x, geps_y, geps_z, eps_size);
}

void FFTController::destroy_fftw_stuff() {
  fftw_destroy_plan(plan);
  fftw_free(in_pointer);
  fftw_free(out_pointer);
  in_pointer = out_pointer = NULL;
}

void FFTController::setup_fftw_3d(int nfft[3], int direction) {
  // check for some dumb input values
  if (nfft[0]<=0 || nfft[1] <=0 || nfft[2] <=0) {
    CkPrintf("setup_fftw_3d routine received illegal value for nfft. \
              nfft should be positive number.");
    CkExit();
  }
  if (!(direction==-1 || direction==1)) {
    CkPrintf("setup_fftw_3d routine received illegal value for direction. \
              FFTW direction must either 1 or -1");
    CkExit();
  }

  // if not the first time and any parameters (size, direction) mismatch
  // we need to destroy old plans to set up new ones
  if (!first_time && (direction != old_direction ||
      nfft[0] != old_nfft[0] ||
      nfft[1] != old_nfft[1] ||
      nfft[2] != old_nfft[2])) {
    destroy_fftw_stuff();
    first_time = true;
  }

  // if first time, we need to set up
  if(first_time) {
    const int ndata = nfft[0]*nfft[1]*nfft[2];
    in_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
    out_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);

    CmiLock(fft_plan_lock);
    plan = fftw_plan_dft_3d(nfft[0], nfft[1], nfft[2],
        in_pointer, out_pointer, direction, FFTW_ESTIMATE);
    CmiUnlock(fft_plan_lock);

  }

  // now the old value changes to new value
  first_time = false;
  old_nfft[0] = nfft[0];
  old_nfft[1] = nfft[1];
  old_nfft[2] = nfft[2];
  old_direction = direction;
}
