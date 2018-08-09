//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP PAW energies (however calculated)
//==========================================================================

#ifndef _CPPAW_
#define _CPPAW_

#include "Atoms.h"
class CPPAW{

  //---------------------------------------------------------------------------
  public:

    //---------------------------------------------------------------------------
    //con-destruct:
    CPPAW(){};
    ~CPPAW(){};

    //---------------------------------------------------------------------------
    // functions defined in PAWexchangeCorrelation.C
    static void PAW_exc_calc();

    //---------------------------------------------------------------------------
    // functions defined in PAWhartElectronIon.C 
    static void PAW_eNhart_calc();

    //---------------------------------------------------------------------------
    // functions defined in PAWcomputeNN.C 
    static void PAW_compute_NNlist();


    //---------------------------------------------------------------------------
}; //CPPAW
//==========================================================================

#endif 
