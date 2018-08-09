//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                         PI_MD:
//             The future of simulation technology
//             ------------------------------------
//                     Module: control_sim_parms.c
//
//
// This subprogram reads in user simulation params and echoes them
// and the default parameters to a file
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "standard_include.h"
#include "class_PAW_info.h"
#include "configure.h"
//#include "proto_friend_lib_entry.h"
#include "proto_handle_entry.h"
#if CMK_PROJECTIONS_USE_ZLIB
#include "zlib.h"
#endif

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfig(char* input_name, PAWINFO* pawinfo)
  //===================================================================================
{//begin routine
  //===================================================================================

  int num_dict_fun;
	int num_dict_paw;
  DICT_WORD *dict_fun;
  DICT_WORD *dict_paw;
  DICT_WORD word; 
	
  int nline;
  int nkey,nfun_key;
  char *fun_key,*fname;
  int ind_key;
  FILE *fp;

  fun_key = new char [PINY_MAXWORD];
  fname   = new char [1024]; 

  //===================================================================================
  // Tell everyone you are busy

  PRINTF("  =============================================================\n");
  PRINTF("  Reading Piny PAW physics input from file : %s\n",input_name);
  PRINTF("  -------------------------------------------------------------\n\n");

  if(PINY_MAXWORD!=MAXWORD || PINY_MAXLINE != MAXLINE){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect word and line sizes\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //===================================================================================
  // Set up the dictionaries

  set_config_dict_fun  (&num_dict_fun  ,&dict_fun);
  set_config_dict_paw  (&num_dict_paw  ,&dict_paw);

  //===================================================================================
  // Read the input file and fill the dictionaries with user input

  fp = fopen((const char *) input_name,"r");
  if(fp==NULL){
     PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("       Input File %s not found\n",input_name);
     PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  nline = 1;  nfun_key=0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,input_name)){
    get_fun_key_index(fun_key,num_dict_fun,dict_fun,nline,nfun_key,
        input_name,&ind_key);
    nkey  = 0;
    while(get_word(fp,&word,&nline,&nkey,nfun_key,input_name)){
      switch(ind_key){
        case 1 : put_word_dict(&word,dict_paw,num_dict_paw,fun_key,nline,
                     nkey,nfun_key,input_name);break;
      }//end switch
    }// end while
  }//end while

  fclose(fp);

  //===================================================================================
  // Take the information out of the dictionary and put it in the class

  set_config_params_paw  (dict_paw,  dict_fun[1].keyword, input_name, pawinfo);

  simpleRangeCheck(); // redundant checking
	Finale(pawinfo);

  //===================================================================================
  // Output your parameter choices to the screen

  // Tokenize the input cpaimd config file name to remove all directory paths
  char *cfgInputName  = new char[strlen(input_name)+1];
  char *cfgOutputName = NULL;
  char dirSeparator[] = "/"; ///< @warning: Will we ever run on Windows and get screwed?
  char *tokenized     = strtok( strcpy(cfgInputName,input_name), dirSeparator);
  while (tokenized != NULL)
  {
    cfgOutputName = tokenized;
    tokenized = strtok(NULL,dirSeparator);
  }
  // The cpaimd config output file is written in the current directory (and not where the input file is located)
  sprintf(fname,"%s.out",cfgOutputName);
  fp = fopen((const char*) fname,"w");
  write_cpaimd_config(fp, dict_paw,  num_dict_paw,  dict_fun[1].keyword);
  fclose(fp);
  delete [] cfgInputName;
  //===================================================================================
  // Free memory :

	delete [] fun_key;
	delete [] fname;
	delete [] dict_fun;
	delete [] dict_paw;

  //===================================================================================
  // Tell Everyone you are done

  PRINTF("  -------------------------------------------------------------\n");
  PRINTF("  Completed reading Piny PAW physics input from file : %s\n",input_name);
  PRINTF("  =============================================================\n\n");

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfigCharm(char *input_name, PAWINFO* pawinfo)
  //===================================================================================
{//begin routine
  //===================================================================================

  int num_dict_fun;
	int num_dict_PAWrho_chare;
  DICT_WORD *dict_fun;
  DICT_WORD *dict_PAWrho_chare;
  DICT_WORD word; 
	
  int nline;
  int nkey,nfun_key;
  char *fun_key,*fname;
  int ind_key;
  FILE *fp;

  fun_key = new char [PINY_MAXWORD];
  fname   = new char [1024]; 

  //===================================================================================
  // Tell everyone you are busy

  PRINTF("  =============================================================\n");
  PRINTF("  Reading PAW charm input from file : %s\n",input_name);
  PRINTF("  -------------------------------------------------------------\n\n");

  if(PINY_MAXWORD!=MAXWORD || PINY_MAXLINE != MAXLINE){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect word and line sizes\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //===================================================================================
  // Set up the dictionaries

  set_config_dict_fun  (&num_dict_fun  ,&dict_fun);
	set_config_dict_PAWRhoChare (&num_dict_PAWrho_chare  ,&dict_PAWrho_chare);

  //===================================================================================
  // Read the input file and fill the dictionaries with user input

  fp = fopen((const char *) input_name,"r");
  if(fp==NULL){
     PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("       Input File %s not found\n",input_name);
     PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  nline = 1;  nfun_key=0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,input_name)){
    get_fun_key_index(fun_key,num_dict_fun,dict_fun,nline,nfun_key,
        input_name,&ind_key);
    nkey  = 0;
    while(get_word(fp,&word,&nline,&nkey,nfun_key,input_name)){
      switch(ind_key){
        case 2 : put_word_dict(&word,dict_PAWrho_chare,num_dict_PAWrho_chare,fun_key,nline,
                     nkey,nfun_key,input_name);break;
      }//end switch
    }// end while
  }//end while

  fclose(fp);

  //===================================================================================
  // Take the information out of the dictionary and put it in the class

  set_config_params_PAWRhoChare (dict_PAWrho_chare, dict_fun[2].keyword, input_name, pawinfo);

  simpleRangeCheck(); // redundant checking
//	FinalePAW(pawinfo);

  //===================================================================================
  // Output your parameter choices to the screen

  // Tokenize the input cpaimd config file name to remove all directory paths
  char *cfgInputName  = new char[strlen(input_name)+1];
  char *cfgOutputName = NULL;
  char dirSeparator[] = "/"; ///< @warning: Will we ever run on Windows and get screwed?
  char *tokenized     = strtok( strcpy(cfgInputName,input_name), dirSeparator);
  while (tokenized != NULL)
  {
    cfgOutputName = tokenized;
    tokenized = strtok(NULL,dirSeparator);
  }
  // The cpaimd config output file is written in the current directory (and not where the input file is located)
  sprintf(fname,"%s.out",cfgOutputName);
  fp = fopen((const char*) fname,"w");
  write_cpaimd_config(fp, dict_PAWrho_chare,  num_dict_PAWrho_chare,  dict_fun[2].keyword);
  fclose(fp);
  delete [] cfgInputName;
  //===================================================================================
  // Free memory :

	delete [] fun_key;
	delete [] fname;
	delete [] dict_fun;
	delete [] dict_PAWrho_chare;

  //===================================================================================
  // Tell Everyone you are done

  PRINTF("  -------------------------------------------------------------\n");
  PRINTF("  Completed reading Piny PAW physics input from file : %s\n",input_name);
  PRINTF("  =============================================================\n\n");

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_fun  (int *num_dict  ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary

  num_dict[0] = 2;
	*dict = new DICT_WORD [num_dict[0] + 1];

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind = 0;
  //------------------------------------------------------------------------------
  //  1)~ProjectorAugmentedWave_def[ ]
  ind++;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"ProjectorAugmentedWave_def");
  strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  2)~PAWRhoChare_def[ ]
  ind++;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"PAWRhoChare_def");
  strcpy((*dict)[ind].keyarg," ");
  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_paw  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary

  num_dict[0] = 11;
	*dict = new DICT_WORD [num_dict[0] + 1];

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind = 0;
  //-----------------------------------------------------------------------------
  //  1)\pawOn{}
  ind++;
  strcpy((*dict)[ind].keyword,"pawOn");
  strcpy((*dict)[ind].keyarg,"Off");
  strcpy((*dict)[ind].error_mes,"Off or On");
  //-----------------------------------------------------------------------------
  //  2)\numRadialPoints{}
  ind++;
  strcpy((*dict)[ind].keyword,"numRadialPoints");
  strcpy((*dict)[ind].keyarg,"10");
  strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  3)\numCosthetaPoints{}
  ind++;
  strcpy((*dict)[ind].keyword,"numCosthetaPoints");
  strcpy((*dict)[ind].keyarg,"6");
  strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  4)\UnitlessScreeningLength{}
  ind++;
  strcpy((*dict)[ind].keyword,"UnitlessScreeningLength");
  strcpy((*dict)[ind].keyarg,"3.0");
  strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  5)\CompensationModelType{}
  ind++;
  strcpy((*dict)[ind].keyword,"CompensationModelType");
  strcpy((*dict)[ind].keyarg,"Gaussian");
  strcpy((*dict)[ind].error_mes,"Gaussian or sphericalBessel");
  //-----------------------------------------------------------------------------
  //  6)\AtomInputFileName{}
  ind++;
  strcpy((*dict)[ind].keyword,"AtomInputFileName");
  strcpy((*dict)[ind].keyarg,"atoms.in");
  strcpy((*dict)[ind].error_mes,"file containing atom positions and types");
  //-----------------------------------------------------------------------------
  //  7)\RealSpaceRadialCutoff{}
  ind++;
  strcpy((*dict)[ind].keyword,"RealSpaceRadialCutoff");
  strcpy((*dict)[ind].keyarg,"4.0");
  strcpy((*dict)[ind].error_mes,"a number > 0 in bohr");
  //-----------------------------------------------------------------------------
  //  8)\PlaneWaveEnergyCutoff{}
  ind++;
  strcpy((*dict)[ind].keyword,"PlaneWaveEnergyCutoff");
  strcpy((*dict)[ind].keyarg,"20.0");
  strcpy((*dict)[ind].error_mes,"a number > 0 in Rydberg");
  //-----------------------------------------------------------------------------
  //  9)\Periodicity{}
  ind++;
  strcpy((*dict)[ind].keyword,"Periodicity");
  strcpy((*dict)[ind].keyarg,"3");
  strcpy((*dict)[ind].error_mes,"periodicity of the atomic supercell, 0/1/2/3");
  //-----------------------------------------------------------------------------
  //  10)\RealSpaceImages{}
  ind++;
  strcpy((*dict)[ind].keyword,"RealSpaceImages");
  strcpy((*dict)[ind].keyarg,"0");
  strcpy((*dict)[ind].error_mes,"a number >=0");
  //-----------------------------------------------------------------------------
  //  11)\LmaxPW{}
  ind++;
  strcpy((*dict)[ind].keyword,"LmaxPW");
  strcpy((*dict)[ind].keyarg,"0");
  strcpy((*dict)[ind].error_mes,"a number >=0");
  //-----------------------------------------------------------------------------
}//end routine
//===================================================================================

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_PAWRhoChare  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary

  num_dict[0] = 1;
	*dict = new DICT_WORD [num_dict[0] + 1];

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind = 0;
  //-----------------------------------------------------------------------------
  //  1)\NatmChunk{}
  ind++;
  strcpy((*dict)[ind].keyword,"NatmChunk");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"1 <= number <= Natm");
  //-----------------------------------------------------------------------------
}//end routine
//===================================================================================

//===================================================================================
// Some simple range checking : needs some love
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::simpleRangeCheck(){
  //===================================================================================

  //  rangeExit(launchNLeesFromRho,"launchNLeesFromRho",0);

  //  rangeExit(invsqr_tolerance,"invsqr_tolerance;",0);

  //---------------------------------------------------------------------------------
}//end routine
//===================================================================================

//============================================================================
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Config::rangeExit(int param, const char *name, int iopt){
  //============================================================================

  switch(iopt){
    case 0:
      if(param<1){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("     The parameter %s must be >0 not %d \n",name,param);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      break;
    case 1:
      if(param<0 || param>1){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("     The parameter %s must be 1(on) or 0 (off) \n",name);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      break;
  }//end switch

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

//===================================================================================
// Consistency Checks on the input
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::Finale(PAWINFO* pawinfo){
  //===================================================================================
	CELL *cell = &(pawinfo->cell);
	double Rcut = cell->Rcut;
	double Ecut = cell->Ecut;
	double gcut = cell->gcut;
	double alpb = cell->alpb;

  double gammasq = gcut*Rcut; // Gcut = 2*gcut, so Gcut^2/4 = gcut^2 = Ecut (in Ryd) = gamma^4/Rcut^2; gamma^2 = gcut*Rcut
  double gamma_conv = sqrt(gammasq);

  //===================================================================================
  // Code consistency checks
	if (pawinfo->nr > 20) {
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("     The PAW radial integration points are very large %d\n", pawinfo->nr);
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	}//end if
	if (pawinfo->ntheta > 12) {
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("     The PAW Costheta integration points are very large %d\n", pawinfo->ntheta);
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	}//end if
	if ((pawinfo->beta_unitless > 5.0) || (pawinfo->beta_unitless < 1.0)) {
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("     The PAW Coulomb screening parameter %g should be >= 1.0 and <= 5.0\n", pawinfo->beta_unitless);
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	}//end if
  if (gamma_conv < 3.0) {
		PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("     The gamma(.%10g) = alpha_bar(.%10g)*Rcut(.%10g) is too small!\n", gamma_conv, alpb, Rcut);
		PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		FFLUSH(stdout);
		EXIT(1);
  } // end if
	if ((Ecut > 100.0) || (Ecut < 5.0)) {
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("     The PW energy cutoff %g should generally be >= 5.0 and <= 100.0\n", Ecut);
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	}//end if
	if ((Rcut > 5.0) || (Rcut < 0.5)) {
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("     The PAW real space cutoff %g should generally be >= 0.5 and <= 5.0\n", Rcut);
    PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
	}//end if

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

//===================================================================================
// Consistency Checks on the input
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::FinaleCharm(PAWINFO* pawinfo){
  //===================================================================================
  // Code consistency checks
  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void Config::write_cpaimd_config(FILE *fp,DICT_WORD *dict, int num_dict, char *fun_key){
  //========================================================================

  int iuset,itype;

  //========================================================================
  //     I) Write out meta key word

  fprintf(fp,"================================================\n");
  fprintf(fp,"cccccccccccccccccccccccccccccccccccccccccccccccc\n");
  fprintf(fp,"================================================\n");
  fprintf(fp,"~%s[\n",fun_key);

  //================================================================================
  //    II) User defined parameters

  iuset = 1;itype = 1;
  fprintf(fp,"------------------------------------------------\n");
  fprintf(fp,"  User defined parameters with software overides\n");
  fprintf(fp,"------------------------------------------------\n");
  dict_print(fp,num_dict,dict,itype,iuset);

  //================================================================================
  //   III) Default parameters

  iuset = 0;itype = 1;
  fprintf(fp,"------------------------------------------------\n");
  fprintf(fp,"  Default parameters with software overides\n");
  fprintf(fp,"------------------------------------------------\n");
  dict_print(fp,num_dict,dict,itype,iuset);

  //================================================================================
  //   IV) End Meta key word  output

  fprintf(fp,"------------------------------------------------\n]\n");
  fprintf(fp,"================================================\n\n\n");

  //========================================================================
}//end routine
//========================================================================

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_paw  (DICT_WORD *dict, char *fun_key, char *input_name, PAWINFO * pawinfo){
  //===================================================================================
	CELL* cell = &(pawinfo->cell);
  FGRID_CONTAINER *fgrid_container = &(pawinfo->fgrid_container);

  int ind = 0;

  //===================================================================================
  //-----------------------------------------------------------------------------
  //  1)\pawOn{}
  ind++;
	pawinfo->paw_on = -1;
	if (strcasecmp(dict[ind].keyarg, "On") == 0) {pawinfo->paw_on = 1;}
	if (strcasecmp(dict[ind].keyarg, "Off") == 0) {pawinfo->paw_on = 0;}
	if (pawinfo->paw_on == -1) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  2)\numRadialPoints{}
	ind++;
	sscanf(dict[ind].keyarg, "%d", &(pawinfo->nr));
	if (pawinfo->nr < 1) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  3)\numCosthetaPoints{}
	ind++;
	sscanf(dict[ind].keyarg, "%d", &(pawinfo->ntheta));
	if (pawinfo->ntheta < 1) {keyarg_barf(dict,input_name,fun_key,ind);}
	pawinfo->nphi = 2*(pawinfo->ntheta);
  //-----------------------------------------------------------------------------
  //  4)\UnitlessScreeningLength{}
	ind++;
	sscanf(dict[ind].keyarg, "%lg", &(pawinfo->beta_unitless));
	sscanf(dict[ind].keyarg, "%lg", &(cell->beta_unitless));
	if (pawinfo->beta_unitless <= 0.0) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  5)\CompensationModelType{}
  ind++;
	pawinfo->model = 0;
	if (strcasecmp(dict[ind].keyarg, "sphericalBessel") == 0) {pawinfo->model = 1;}
	if (strcasecmp(dict[ind].keyarg, "Gaussian") == 0) {pawinfo->model = 2;}
	if (pawinfo->model == 0) {keyarg_barf(dict,input_name,fun_key,ind);}
	if (pawinfo->model == 1) {pawinfo->nr++;}
  //-----------------------------------------------------------------------------
  //  6)\AtomInputFileName{}
  ind++;
  strcpy(pawinfo->atmFnameIn, dict[ind].keyarg);
  FILE * fp = fopen(pawinfo->atmFnameIn, "r");
  if (fp == NULL) {keyarg_barf(dict,input_name,fun_key,ind);}
  fclose(fp);
  //-----------------------------------------------------------------------------
  //  7)\RealSpaceRadialCutoff{}
  ind++;
	sscanf(dict[ind].keyarg, "%lg", &(cell->Rcut));
	if (cell->Rcut < 0) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  8)\PlaneWaveEnergyCutoff{}
  ind++;
	sscanf(dict[ind].keyarg, "%lg", &(cell->Ecut));
	if (cell->Ecut < 0) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  9)\Periodicity{}
  ind++;
	sscanf(dict[ind].keyarg, "%d", &(cell->iperd));
	if (cell->iperd < 0 || cell->iperd > 3) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  10)\RealSpaceImages{}
  ind++;
	sscanf(dict[ind].keyarg, "%d", &(cell->nimg));
	if (cell->nimg < 0) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  11)\LmaxPW{}
  ind++;
	sscanf(dict[ind].keyarg, "%d", &(fgrid_container->lmax_pw));
	if (fgrid_container->lmax_pw < 0) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
	//	assign all the derived sizes
	int nr     = pawinfo->nr;
	int ntheta = pawinfo->ntheta;
	int nphi   = pawinfo->nphi;
  int nf     = nr*ntheta*nphi;
  int nang   = ntheta*nphi;
	pawinfo->nf   = nf;
	pawinfo->nang = nang;
  fgrid_container->nf     = nf;
  fgrid_container->nr     = nr;
  fgrid_container->ntheta = ntheta;
  fgrid_container->nphi   = nphi;
  fgrid_container->nang   = nang;
  
	double Rcut = cell->Rcut;
	double Ecut = cell->Ecut;
	double gcut = sqrt(Ecut); // gcut = sqrt(2*me*Ecut/hbar^2) = sqrt(2*Ecut) in atomic units, for Ecut in Ryd, gcut = sqrt(Ecut)
	double Gcut = 2.0*gcut; // the density g cutoff is twice the pw g cutoff
	double gammasq = gcut*Rcut; // Gcut = 2*gcut, so Gcut^2/4 = gcut^2 = Ecut (in Ryd) = gamma^4/Rcut^2; gamma^2 = gcut*Rcut
	double gamma_conv = sqrt(gammasq);
	double alpb = gamma_conv/Rcut;
	cell->gcut = gcut;
	cell->Gcut = Gcut;
	cell->alpb = alpb;
	cell->animg = (double) (cell->nimg);
	fgrid_container->alpb = alpb;
  //========================================================================
}//end routine
//========================================================================

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_PAWRhoChare  (DICT_WORD *dict, char *fun_key, char *input_name, PAWINFO * pawinfo){
  //===================================================================================

  int ind = 0, my_int;

  //===================================================================================
  //-----------------------------------------------------------------------------
  //  1)\NatmChunk{}
  ind++;
	sscanf(dict[ind].keyarg, "%d", &my_int);
	pawinfo->NatmChunk = my_int;
	if (pawinfo->NatmChunk < 1 ) {keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
}//end routine
//========================================================================
