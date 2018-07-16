//==========================================================================
//                PAW info                          
//             {Variables needed for mem allocation:                        
//                 vxctyplen,nsplin_g,nsplin_g_tot,n_ang_max,               
//                 num_nl_lst       }                                       
//                                                                          
//==========================================================================

#ifndef _PAWINFO_
#define _PAWINFO_

//==========================================================================
class PAWINFO{

  //-------------------------------------------------------------------------
  public:
    int paw_on;                 // Opt: Is PAW ``on''
		int nf;             				// total number of f points
		int nr;              				// number of r grid points
		int ntheta;          				// number of theta grid points
		int nphi;            				// number of phi grid points
		int nrfull;          				// number of r grid points*2
		int nang;            				// number of angular grid points = ntheta*nphi
		double beta_unitless;       // Coulomb screening parameter
    //-------------------------------------------------------------------------
    //con-destruct:
    PAWINFO(){
			paw_on 							= 0;
			nr									= 1;
			ntheta							= 1;
			nphi								= 1;
			nf									= nr*ntheta*nphi;
			nang								= ntheta*nphi;
			beta_unitless				= 1.0;
    };
    ~PAWINFO(){};

    //-------------------------------------------------------------------------
#ifdef PUP_ON
    //pupping
    void pup(PUP::er &p){
      //pupping ints
			p | paw_on;
			p | nf;
			p | nr;
			p | ntheta;
			p | nphi;
			p | nang;
      //pupping dbles
			p | beta_unitless;

      //pupping arrays
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())state_class_out ();
#endif     
    } // end pup routine
#endif

    //----------------------------------------------------------------------
    void state_class_out(){
			int my_pe = 1;
#ifdef CHARM_ON
			my_pe = CkMyPe();
#endif
      char fileName [255];
      sprintf (fileName, "%d_pawinfo.state", my_pe);
      FILE *fp;  fp = fopen(fileName,"w");

      // ints
      fprintf(fp,"paw_on %d\n",paw_on);
      fprintf(fp,"nf %d\n",nf);
      fprintf(fp,"nr %d\n",nr);
      fprintf(fp,"ntheta %d\n",ntheta);
      fprintf(fp,"nphi %d\n",nphi);
      fprintf(fp,"nang %d\n",nang);
			// doubles
      fprintf(fp,"beta_unitless %g\n",beta_unitless);
      fclose(fp);

    }// end routine 

    //-------------------------------------------------------------------------
}; //CPPSEUDO
//==========================================================================

#ifdef PUP_ON
PUPmarshall(PAWINFO);
#endif

#endif
//==========================================================================
