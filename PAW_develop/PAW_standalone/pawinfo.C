//#define CHARM_ON
#ifdef CHARM_ON
#define PUP_ON
#endif
int numPes;   // Global: number of procs
#include "standard_include.h"
#include "configure.h"
#include "class_PAW_info.h"
#ifdef CHARM_ON
class main : public Chare {
  public:
    main(CkMigrateMessage *m) { }
    main(CkArgMsg *);
    ~main();
};
#else
int main(int, char *[]);
#endif
//============================================================================
#ifdef CHARM_ON
extern PAWINFO  readonly_pawinfo;
#endif
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

/**
 * \brief The Main of CPAIMD, it calls all the init functions.
 */
#ifdef CHARM_ON
main::main(CkArgMsg *msg) {
#else
int main (int argc, char *argv[]){
#endif
  //============================================================================
  /**
     # Sequential startup within Main */
  /* Check arguments : Verbose output about startup procedures */
	char exe_name[80];
	char piny_input[80];
	char charm_input[80];
#ifdef CHARM_ON
	strcpy(exe_name, msg->argv[0]);
	strcpy(piny_input, msg->argv[1]);
	strcpy(charm_input, msg->argv[2]);
	int argc = msg->argc;
#else
	strcpy(exe_name, argv[0]);
	strcpy(piny_input, argv[1]);
	strcpy(charm_input, argv[2]);
#endif

  if (argc < 3) {
    PRINTF("Usage: cpaimd.x cpaimd_config pinysystem.input");
		EXIT(1);
  }//endif

  PRINTF("Executing OpenAtom: BINARY - %s\n", exe_name);
  PRINTF("\n");
  PRINT_LINE_STAR;
  PRINTF("Starting Cpaimd-Charm-Driver Setup Phase\n");
  PRINT_LINE_DASH;
	int num_pes = 1;
#ifdef CHARM_ON
	num_pes = CkNumPes();
#endif
  PRINTF("  Cpaimd-Charm-Driver running on %d processors. \n", num_pes);

  //============================================================================
  /* Invoke PINY input class */

  //============================================================================
  /* Invoke parallel driver input class */
  /** \addtogroup startup
      2)  Read the cpaimd_config file which determines parallel
      decomposition.  This is mostly stored in the config object, but a
      bunch of readonly globals also get instantiated with this
      information.
  */
  /**@{*/

  PRINTF("  Reading Driver input from %s\n",piny_input);

//  config.readConfig(input_name);

  PRINT_LINE_DASH;
  PRINTF("Cpaimd-Charm-Driver input completed\n");
  PRINT_LINE_STAR; PRINTF("\n");

  //============================================================================
  PRINT_LINE_DASH;
  PRINTF("Cpaimd-Charm-Driver setup phase completed\n");
  PRINTF("Total time to execute set up in the main chare\n");
  PRINT_LINE_STAR; PRINTF("\n");
  PRINT_LINE_STAR;
  PRINT_LINE_DASH;PRINTF("\n");
  /**@}*/
	
#ifndef CHARM_ON
	return 1;
#endif
  //============================================================================
}// end Main
//============================================================================

//============================================================================
#ifdef CHARM_ON
#include "pawinfo.def.h"
#endif
//============================================================================
