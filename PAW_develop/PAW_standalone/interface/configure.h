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
   void Finale(PAWINFO*);
   void FinaleCharm(PAWINFO*);

   void set_config_dict_fun    (int *, DICT_WORD **);
   void set_config_dict_paw    (int *, DICT_WORD **);
	 void set_config_dict_PAWRhoChare(int *,DICT_WORD **);
   void set_config_params_paw  (DICT_WORD *, char *, char *, PAWINFO *);
   void set_config_params_PAWRhoChare  (DICT_WORD *, char *, char *, PAWINFO *);
	 void write_cpaimd_config(FILE *,DICT_WORD *, int, char *);

	 void readConfig(char*, PAWINFO*);
	 void readConfigCharm(char *, PAWINFO*);
	 void simpleRangeCheck();
	 void rangeExit(int, const char *, int);
  //==================================

//-----------------------------------------------------------------------------------
   }; // end class
//===================================================================================

//===================================================================================

#endif
