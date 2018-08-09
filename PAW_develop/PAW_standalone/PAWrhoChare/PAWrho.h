//==========================================================================
//                PAW info                          
//             {Variables needed for mem allocation:                        
//                 vxctyplen,nsplin_g,nsplin_g_tot,n_ang_max,               
//                 num_nl_lst       }                                       
//                                                                          
//==========================================================================

#ifndef _PAWRHO_
#define _PAWRHO_

class PAWrho : public CBase_PAWrho {
public:
  //======================================
  // integer member variables
  int me;
	int myAtmStart;
	int myAtmEnd;
	int numMyAtm;
	int nElements;
	int Natm;
  //==========================
  // Constructor
  PAWrho(void); 
  //==========================
  // Memeber funcitons
  void PAWrhoConstructed(CkReductionMsg *);
  void PAWrhoConstructed0();
};
#endif
