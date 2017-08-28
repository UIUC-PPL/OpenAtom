#ifndef __STATES_H__
#define __STATES_H__

#include "states.decl.h"
#include "ckcomplex.h"

class States : public CBase_States {

 public:
  complex *stateCoeff;   // state coefficient in reciprocal space (G space)
  complex *stateCoeffR;  // state in R space 
  complex *stateCoeff_shifted;  // shifted k grid state coefficient in reciprocal space (G space)
  complex *stateCoeffR_shifted; // shifted k grid state in R space 
  int *ga, *gb, *gc;
  int numCoeff;
  int ikpt;        // index for k point
  int ispin;       // index for spin
  int istate;      // index for state
  bool shifted;    // if states are shifted or not
  
  int ibinary_opt; // binary file option to read state file
  bool doublePack; // if only have gamma point, then true (=1). Otherwise, false(=0)
  char fileName[1000];  // file name for state

  int nfft[3]; // number of fft grid in each direction
  
  int nocc; // number of occupied states
  
  /// Constructors ///
  States();
  States(CkMigrateMessage *msg);

  /// Entry Methods ///                                                                         
  void fftGtoR();
  void sendToCache();
  void sendToComputeF();

  /// fftw routines ///
  // void fft_G_to_R(); -> doesn't seem to exist anymore

  /// scalar routines ///
  void readState(char *);
  void readStateShifted(char *);  
};

extern /* readonly */ CProxy_States states_proxy;

#endif //__STATES_H__
