//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*       Projector Augmented Wave Standalone Code  main module              */
//============================================================================
/* Standard Include files filed by class defs for this module */
#include "standard_include.h"  // has charmm++.h and pup.h
#ifndef CHARM_OFF
#include "pawmain.decl.h"      // charm++ declarations for main   module including extern declaration of cleanexit chare
#include "cleanexit.decl.h"       // charm++ declarations for cleanexit module
#include "ATMSGRP.decl.h"       // charm++ declarations for atmsGrp module
#include "PAWrho.decl.h"       // charm++ declarations for PAWrho module
#endif
#include "pawmain.h"           // c++ declarations for main class / module
#include "class_PAW_info.h"    // c++ pawinfo readonly class declarations (including pup)
#include "configure.h"         // c++ configure input class declarations (non-parallel class)
#include "fastAtoms.h"				// Fastatom class
#ifndef CHARM_OFF
#include "cleanexit.h"            // c++ declarations for cleanexit chare array module
#include "atmsGrp.h"            // c++ declarations for atmsGrp chare array module
#include "PAWrho.h"            // c++ declarations for PAWrho chare array module
#include "PAW_"
void copyFastAtomsToMsg(ATMSGRPMSG *, FASTATOMS *);
#endif
#include "scalar_cleanexit.h"  // scalar clean exit proto type

//============================================================================
/* charm++ globels used in this file */
#ifndef CHARM_OFF
PAWINFO  readonly_pawinfo;   // declaration of global readonly class - PAWINFO : everywhere else is extern.
CProxy_cleanexit cleanexit_Proxy;  // proxy to cleanexit chare array for launch
CProxy_ATMSGRP ATMSGRP_Proxy;  // proxy to atmsGrp chare array for launch
CProxy_PAWrho PAWrho_Proxy;   // proxy to PAWrho for launch
#endif
//============================================================================


//============================================================================
/*  Main module for PAW code - reads input, pups it out and exits            */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
#ifndef CHARM_OFF
main::main(CkArgMsg *msg) { // msg contains argc and argv[]
#else
int main (int argc, char *argv[]){
#endif
//============================================================================
/* I) Parse the command line input */
  char exe_name[80];
  char piny_input[80];
  char charm_input[80];

  int num_pes=1; 
#ifndef CHARM_OFF
  int argc = msg->argc;
  num_pes = CkNumPes();
#endif
  if (argc < 3) {
    PRINTF("Usage: pawmain.x PAW.in charm.in\n");
    EXIT(1);
  }//endif
#ifndef CHARM_OFF
  strcpy(exe_name,   msg->argv[0]);
  strcpy(piny_input, msg->argv[1]);
  strcpy(charm_input,msg->argv[2]);
#else
  strcpy(exe_name,   argv[0]);
  strcpy(piny_input, argv[1]);
  strcpy(charm_input,argv[2]);
#endif

//============================================================================
/* II. Tell the user you are what you are running and on how many processors */
  PRINT_LINE_STAR;
  PRINTF("Executing PAW: BINARY - %s - on %d processors\n", exe_name,num_pes);
  PRINT_LINE_STAR; 
  PRINTF("\n");

//============================================================================
/* III. Setup phase */
  //------------------------------------------------------------------------------
  /* III.1 Tell the user you are entering the setup phase */
  PRINT_LINE_STAR;
  PRINTF("Starting PAW Scalar Setup Phase\n");
  PRINT_LINE_DASH; PRINTF("\n");

  //------------------------------------------------------------------------------
  /* III.2 Invoke PINY input class on processor 0*/
  Config config;
#ifndef CHARM_OFF
  PAWINFO * pawinfo = PAWINFO::get();
#else
  PAWINFO pawinfo_master; PAWINFO *pawinfo = &pawinfo_master;
#endif
  config.readConfig(piny_input, pawinfo);
  //------------------------------------------------------------------------------
  /* III.3 Invoke Charm++ input class */
  config.readConfigCharm(charm_input, pawinfo);

  //------------------------------------------------------------------------------
  /* III.4 Invoke the atoms input, only on proc 0 */
  FASTATOMS fastatoms;
  CELL *cell = &(pawinfo->cell);
  ATOM_MAPS *atom_maps = &(pawinfo->atom_maps);
  fastatoms.readatoms(cell, atom_maps, pawinfo->atmFnameIn, pawinfo->NatmChunk);
  pawinfo->Natm = fastatoms.natm;

  //------------------------------------------------------------------------------
  /* III.5 Invoke the fgrid input, only on proc 0 */
  FGRID_CONTAINER *fgrid_container = &(pawinfo->fgrid_container);
  PAW_FGRID_INIT paw_fgrid_init;
  paw_fgrid_init.fill_fgrid(fgrid_container, atom_maps);

  //------------------------------------------------------------------------------
  /* III.6 Tell the user you are done */
  PRINT_LINE_DASH;
  PRINTF("PAW Scalar setup phase completed\n");
  PRINT_LINE_STAR; PRINTF("\n");

//============================================================================
/* IV. Tell the user you are starting PAW */
  PRINT_LINE_STAR;
  PRINTF("Starting PAW parallel setup phases\n");
  PRINT_LINE_DASH;PRINTF("\n");

  //---------------------------------------------------------------------------
  /* Launch a chare array of length num_pes which will then invoke a reduction which end the program*/
#ifndef CHARM_OFF  
    CkArrayOptions optsClean;                       // Container class of info about chare array to be launched
    optsClean.setNumInitial(num_pes);               // 1D array of length num_pes - making a group by hand
    optsClean.setStaticInsertion(true);             // autoinsert array instead of looping over and inserting 1 by 1
    optsClean.setAnytimeMigration(false);	         // don't migrate
    cleanexit_Proxy = CProxy_cleanexit::ckNew(optsClean); // launch the array which will do a round Robin and then wait 
#endif    

  //---------------------------------------------------------------------------
  /* Launch a chare which will handle the atoms */
#ifndef CHARM_OFF  
    int natm = fastatoms.natm;
    ATMSGRPMSG *atmMsg = new (natm, natm, natm, natm, natm, natm, natm, natm,8*sizeof(int)) ATMSGRPMSG;
    copyFastAtomsToMsg(atmMsg, &fastatoms); // copy fastatoms into msg and send the msg to allthe chare array elements
    ATMSGRP_Proxy = CProxy_ATMSGRP::ckNew(atmMsg); // launch the chare array with the atoms in the msg,
                                                   //each guy gets its own copy of the atoms 
    fastatoms.destroy();
#endif
	
  //---------------------------------------------------------------------------
  /* Launch a chare array of length NatmChunk to handle to PAW density inside each atom*/
#ifndef CHARM_OFF
    int NatmChunk = pawinfo->NatmChunk;
    CkArrayOptions optsPAWrho;                       // Container class of info about chare array to be launched
    optsPAWrho.setNumInitial(NatmChunk);             // 1D array of length num_pes - making a group by hand
    optsPAWrho.setStaticInsertion(true);             // autoinsert array instead of looping over and inserting 1 by 1
    optsPAWrho.setAnytimeMigration(false);	         // don't migrate
    PAWrho_Proxy = CProxy_PAWrho::ckNew(optsPAWrho); // launch the array which will process the PAW density and exit
#endif    
//============================================================================

/* V. Scalar Exit */
#ifdef CHARM_OFF
    PAWscalarCleanExit();
    return 1;
#endif
  }// end Main
//============================================================================

//============================================================================
/*  Scalar exit : in scalar just come here directly: in parallel reduction sends you here */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PAWscalarCleanExit(){
   PRINTF("\n");
   PRINT_LINE_DASH;
   PRINTF("Completed PAW parallel setup phases\n");
   PRINT_LINE_STAR;
   EXIT(1);
}//end routine

//============================================================================
/*  Scalar exit : in scalar just come here directly: in parallel reduction sends you here */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void copyFastAtomsToMsg(ATMSGRPMSG *msg, FASTATOMS *fastatoms){
	int natm = fastatoms->natm;
	msg->natm = natm;
	for(int i=0; i<natm; i++) {
		msg->x[i]    = fastatoms->x[i];
		msg->y[i]    = fastatoms->y[i];
		msg->z[i]    = fastatoms->z[i];
		msg->q[i]    = fastatoms->q[i];
		msg->qt[i]   = fastatoms->qt[i];
		msg->alp[i]  = fastatoms->alp[i];
		msg->beta[i] = fastatoms->beta[i];
		msg->Rpc[i]  = fastatoms->Rpc[i];
	} // end for
}//end routine
	
 
//============================================================================
/* Strange charm++ input that is needed at the bottom of files */
#ifndef CHARM_OFF
#include "pawmain.def.h"  // main module def's : cleanexit is declared extern in pawmain.ci taking care of it
#endif
//============================================================================
