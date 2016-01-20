//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
/** \file debug_flags.h
 *  Useful debugging flags
 *
 */
//=============================================================================

//=============================================================================
// Major League Debug Controllers

// do reductions to time each phase of the time step
#define _CP_SUBSTEP_TIMING_

//-----------------------------------------------------------------------
//#define  _CP_DEBUG_SCALC_ONLY_  // REALLY scalc only NOTHING ELSE
#ifdef _CP_DEBUG_SCALC_ONLY_
#define _CP_DEBUG_SFNL_OFF_
#define _CP_DEBUG_VKS_OFF_
#define _CP_DEBUG_RHO_OFF_  // not really necesary as no path to rho
#endif

//-----------------------------------------------------------------------
//#define  _CP_DEBUG_VKS_ONLY_  //what we called scalc only : states make rho but don't send
#ifdef _CP_DEBUG_VKS_ONLY_
#define _CP_DEBUG_SFNL_OFF_
#define _CP_DEBUG_RHO_OFF_     // blocks send to rho from state chares
#endif

//-----------------------------------------------------------------------
//#define  _CP_DEBUG_NON_LOCAL_ONLY_  // SFNL + scalc, no density creation 
#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
#define _CP_DEBUG_VKS_OFF_
#define _CP_DEBUG_RHO_OFF_  // not really necessary as there is no path to rho
#endif

//-----------------------------------------------------------------------
//#define _CP_DEBUG_HARTEEXT_OFF_   // this leaves everything on but eext/hartree

#define _HFCALCULATOR_VERBOSE_

//-----------------------------------------------------------------------
/// Debug the PairCalculator and its companion whistles
//#define _PAIRCALC_DEBUG_
//#define DEBUG_CP_PAIRCALC_ALL
#ifdef DEBUG_CP_PAIRCALC_ALL
	#define _PAIRCALC_CREATE_DEBUG_
	#define _PAIRCALC_DEBUG_PLACE_
	#define DEBUG_CP_PAIRCALC_CREATION
	#define DEBUG_CP_PAIRCALC_COMM
	#define DEBUG_CP_PAIRCALC_PSIV
	#define DEBUG_CP_PAIRCALC_RDMA
	#define DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
	#define DEBUG_MESSAGEDATACOLLATOR_ALL
#endif


//=============================================================================
// src_charm_driver/paircalc/ckPairCalculator.C and a dozen other places
//#define _NAN_CHECK_
//#define PRINT_DGEMM_PARAMS

//=============================================================================
// src_charm_driver/main/groups.C
//#define _CP_DEBUG_PSI_OFF_

//=============================================================================
// src_charm_driver/fft_slab_ctrl
//#define _CP_DEBUG_NON_LOCAL_ONLY_
//#define _CP_DEBUG_FFTR_VKSR_

//=============================================================================
// src_charm_driver/utilities/util.C
//#define _CP_DEBUG_UTIL_VERBOSE_
//#define _CP_DEBUG_LINE_

//=============================================================================
// include/CPcharmParaInfo.h
//#define _CP_DEBUG_PARAINFO_VERBOSE_

//=============================================================================
// src_charm_driver/cp_density_ctrl/CP_Rho_RealSpacePlane.C : 
//#define _CP_DEBUG_RHOR_VERBOSE_
//#define _CP_DEBUG_RHOR_RHO_
//#define _CP_DEBUG_RHOR_VKSA_
//#define _CP_DEBUG_RHOR_VKSB_
//#define _CP_DEBUG_RHOR_VKSC_
//#define _CP_DEBUG_RHOR_VKSD_
//#define _CP_DEBUG_RHOR_VKSE_
//#define _CP_DEBUG_HARTEEXT_OFF_

//=============================================================================
// src_charm_driver/cp_density_ctrl/CP_Rho_GSpacePlane.C : 
//#define _CP_DEBUG_VKS_GSPACE_
//#define _CP_DEBUG_VKS_FORCES_
//#define _CP_DEBUG_RHOG_VERBOSE_
//#define _CP_DEBUG_RHOG_RHOG_
//#define _CP_DEBUG_RHOG_VKSA_
//#define _CP_DEBUG_HARTEEXT_OFF_

//=============================================================================
// src_charm_driver/cp_state_ctrl/CP_State_GSpacePlane.C    
//#define _CP_DEBUG_NON_LOCAL_ONLY_
//#define _CP_DEBUG_OLDFORCE_
//#define _CP_DEBUG_STATEG_VERBOSE_
//#define _CP_DEBUG_STATE_GPP_VERBOSE_
//#define _CP_DEBUG_LMAT_
//#define _CP_DEBUG_DYNAMICS_
//#define _CP_DEBUG_ORTHO_OFF_
//#define _CP_DEBUG_SFNL_OFF_
//#define _CP_DEBUG_PSI_OFF_
#define _CP_DEBUG_COEF_SCREEN_
//#define _CP_DEBUG_UPDATE_OFF_
//#define DEBUG_CP_GSPACE_CREATION
//#define DEBUG_CP_GSPACE_PSIV
//#define  _CP_GSPACE_PSI_FORCE_OUTPUT_BEFORE_LAMBDA_
/**
 * Enables dumping of Psi and Lambda from GSpace
 * Psi and Lambda are dumped both before being sent to PC ("Bf" suffix)
 * and immediately after receiving back from PC ("Af" suffix). Note that
 * the dump files are overwritten at each iteration.
*/
//#define PAIRCALC_TEST_DUMP

//If paircalc test will be run, lambda and psi dumps must be activated
#ifdef PAIRCALC_TEST_DUMP
#define _CP_GS_DUMP_LAMBDA_
#define _CP_GS_DUMP_PSI_
#endif
/** @note: Only one of these barriers is guaranteed to work at any given time. 
 * Turning on multiple barriers might beat the already convoluted RTH code in GSpaceDriver and crash!!
 */
//#define BARRIER_CP_GSPACE_PSI
//#define BARRIER_CP_GSPACE_PSIV
//#define BARRIER_CP_PARTICLEPLANE_NONLOCAL
//#define BARRIER_CP_GSPACE_IFFT

//=============================================================================
// src_charm_driver/cp_state_ctrl/CP_State_RealSpacePlane.C : 
//#define _CP_DEBUG_STATER_VKS_
//#define _CP_DEBUG_STATER_VERBOSE_
//#define _CP_DEBUG_RHO_OFF_

//=============================================================================
// src_charm_driver/structure_factor/StructureFactor.C : 
//#define _CP_DEBUG_SF_CALC_

//=============================================================================
// src_charm_driver/structure_factor/StructFactorCache.C : 
//#define _CP_DEBUG_SF_CACHE_

//=============================================================================
// src_charm_driver/orthog_ctrl/ortho.C
//#define  _CP_DEBUG_ORTHO_VERBOSE_
//#define VERBOSE_ORTHO

//#define _CP_DEBUG_SMAT_
//#define _CP_DEBUG_TMAT_
//#define _CP_ORTHO_DUMP_LMAT_
//#define _CP_ORTHO_DUMP_SMAT_
//#define _CP_ORTHO_DUMP_TMAT_
//#define _CP_ORTHO_DEBUG_COMPARE_LMAT_
//#define _CP_ORTHO_DEBUG_COMPARE_SMAT_
//#define _CP_ORTHO_DEBUG_COMPARE_TMAT_ 

//=============================================================================
// src_charm_driver/load_balance/MapTable.C
//#define _MAP_DEBUG_
//#define _MAP_VERBOSE_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_orthog/cp_lowdin.C
//#define _CP_DEBUG_NODIAG_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_integrate/cp_min_std.C
//#define _CP_DEBUG_NEWFORCE_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_ions/cp_rspace_ion.C
// src_piny_physics_v1.0/classical_physics/integrate/control_integrate.C
//#define _CP_DEBUG_ATM_FORC_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_local/cp_hart_ext.C
//#define _CP_DEBUG_VKS_HART_EEXT_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_nonlocal/cp_eke.C
//#define _CP_DEBUG_OLDFORCE_

//=============================================================================
// src_piny_physics_v1.0/piny_to_driver/Parainfoinit.C
//#define _CP_DEBUG_LINES_PER_PLANE_

// ---- All matrix dumps ----
//#define _CP_GS_DUMP_VKS_
//#define _CP_GS_DUMP_LAMBDA_
//#define _CP_GS_DUMP_PSI_
//#define _CP_GS_DEBUG_COMPARE_PSI_
//#define _PAIRCALC_DEBUG_PARANOID_FW_
//#define _CP_ORTHO_DUMP_SMAT_
//#define _CP_ORTHO_DUMP_TMAT_
//#define _CP_ORTHO_DEBUG_COMPARE_TMAT_
//#define _CP_ORTHO_DUMP_LMAT_
//#define _CP_ORTHO_DUMP_GMAT_

