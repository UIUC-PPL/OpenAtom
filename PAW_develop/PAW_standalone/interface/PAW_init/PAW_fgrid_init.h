/** \file PAW_fgrid_init.h
 *
 */

#ifndef _PAW_FGRID_INIT_
#define _PAW_FGRID_INIT_

#include "fgrid.h"
#include "atom_maps.h"

//===================================================================================
// constants needed to create accurate partial waves of erf(a*|r-r'|)/|r-r'| 
typedef struct PW_ERF_CONS {
  int lmax;                      // int : max number of erf pw channels
  int ng;                        // int : number of GL integation points
  double *i_n;                   // list[lmax+1] : mod spher bessel fnc 1st kind small arg
  double *i_n_0;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_1;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_2;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_3;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_1r;                // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind large arg
  double *i_n_2r;                // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind large arg
  double *i_n_3r;                // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind large arg
  double *pref;                  // list[lmax+1] : partial wave pre facotr
  double *pw_erf_num;            // list[lmax+1] : partial wave for one given {r>,r<}
  double *pw_erf_math;           // list[lmax+1] : partial wave for one given {r>,r<}
  double *pw_erf_big;            // list[lmax+1] : partial wave for one given {r>,r<}
  double *gx,*wx;                // list[ng] : GL integation pts and nodes
}PW_ERF_CONS ;

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
class PAW_FGRID_INIT {
//===================================================================================
public:
  int natm_typ;
  int nr;
  int ntheta;
  int nphi;
  int nf;
  int nang;
  int lmax;
  double alpb;
  PW_ERF_CONS pw_erf_cons;

  //==================================
  // Class Functions
  //----------------------------------
   PAW_FGRID_INIT(){};
  ~PAW_FGRID_INIT(){};
   void fill_fgrid(FGRID_CONTAINER *, ATOM_MAPS *);
  //==================================

//-----------------------------------------------------------------------------------
}; // end class
//===================================================================================

//===================================================================================

#endif
