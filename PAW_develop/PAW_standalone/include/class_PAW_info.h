//==========================================================================
//                PAW info                          
//             {Variables needed for mem allocation:                        
//                 vxctyplen,nsplin_g,nsplin_g_tot,n_ang_max,               
//                 num_nl_lst       }                                       
//                                                                          
//==========================================================================

#ifndef _PAW_INFO_
#define _PAW_INFO_

#ifndef CHARM_OFF
#define PUP_ON
#endif

#include "cell.h"
#include "atom_maps.h"
#include "fgrid.h"

#ifndef CHARM_OFF
#define PUP_ON
class PAWINFO; 
extern PAWINFO readonly_pawinfo;
#endif

//==========================================================================
class PAWINFO{

  //-------------------------------------------------------------------------
  public:
    int paw_on;                        // Opt: Is PAW ``on''
    int nf;            		             // total number of f points
    int nr;            		             // number of r grid points
    int ntheta;        		             // number of theta grid points
    int nphi;          		             // number of phi grid points
    int nang;          		             // number of angular grid points = ntheta*nphi
    int model;	   	                   // model type, bessel = 1, gauss = 2 
    int Natm;		                       // number of atoms
    int NatmChunk;		                 // number of chunks of atoms in parallel
    double beta_unitless;              // Coulomb screening parameter
    CELL cell;			                   // cell
		FGRID_CONTAINER fgrid_container;   // fgrid container
    ATOM_MAPS atom_maps;	             // atom_maps
    NAME atmFnameIn;		               // atom input file name, not pupped, only used on proc 0 for main chare
    //-------------------------------------------------------------------------
    //con-destruct:
    PAWINFO(){
			paw_on        = 0;
			nr            = 1;
			ntheta        = 1;
			nphi          = 1;
			nf	          = nr*ntheta*nphi;
			nang          = ntheta*nphi;
			model 	      = 1;
			Natm          = 1;
      NatmChunk     = 1;
			beta_unitless = 1.0;
			strcpy(atmFnameIn, "atoms.in");
    };
    ~PAWINFO(){};

#ifndef CHARM_OFF
    static PAWINFO *get(){
	    return &readonly_pawinfo;  // return the pointer of the global instance
    }
#endif

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
      p | model;
      p | Natm;
      p | NatmChunk;
      //pupping dbles
      p | beta_unitless;
      //pupping cell
      cell.pup(p);
      //pupping atom_maps
      atom_maps.pup(p);
      //pupping fgrid_container
      fgrid_container.pup(p);
      //pupping arrays
#ifdef _PARALLEL_DEBUG_
      state_class_out ();
#endif     
    } // end pup routine
#endif

    //----------------------------------------------------------------------
    void state_class_out(){
      int my_pe = 1;
#ifndef CHARM_OFF
      my_pe = CkMyPe();
#endif
      char fileName [255];
      sprintf (fileName, "%d_pawinfo.state.out", my_pe);
      FILE *fp;  fp = fopen(fileName,"w");

      // ints
      fprintf(fp,"paw_on %d\n",paw_on);
      fprintf(fp,"nf %d\n",nf);
      fprintf(fp,"nr %d\n",nr);
      fprintf(fp,"ntheta %d\n",ntheta);
      fprintf(fp,"nphi %d\n",nphi);
      fprintf(fp,"nang %d\n",nang);
      fprintf(fp,"model %d (1 is sphericalBessel, 2 is Gaussian)\n",model);
      fprintf(fp,"Natm %d\n",Natm);
      fprintf(fp,"NatmChunk %d\n",NatmChunk);
      // doubles
      fprintf(fp,"beta_unitless %g\n",beta_unitless);
      fclose(fp);

    }// end routine 

    //-------------------------------------------------------------------------
}; //PAWINFO
//==========================================================================

#ifdef PUP_ON
PUPmarshall(PAWINFO);
#endif

#endif
//==========================================================================
