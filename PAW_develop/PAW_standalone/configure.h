/** \file configure.h
 *
 */

#ifndef _Configure_
#define _Configure_

#include "dictionary.h"

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
class Config {
//===================================================================================
 public:

  //==================================
  // Class Functions
  //----------------------------------
   Config(){};
  ~Config(){};
   void Finale(int, int, int, int, int, int);

   void set_config_dict_fun    (int *, DICT_WORD **);
   void set_config_dict_pawinfo    (int *, DICT_WORD **);

   void set_config_params_pawinfo  (DICT_WORD *, char *, char *);
  //==================================

//-----------------------------------------------------------------------------------
   }; // end class
//===================================================================================

//===================================================================================

#endif
