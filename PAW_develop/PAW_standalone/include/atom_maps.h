//==========================================================================
//                  Simulation cell data                                    

#ifndef _ATOM_MAPS_
#define _ATOM_MAPS_
#include "pup_utils.h"

class ATOM_MAPS{

  //----------------
  public:
    int natm_typ;          // number of atom types
    int natm;              // number of atoms
    int natm_atm_typ_max;  // max number of atoms of any type
    int lmax;
    int lmax_rho;
    int *lmax_by_typ;
    int *index_atm_typ;    // index of atom type of each atom
    int *natm_atm_typ;     // the number of atoms of each type
    int **list_atm_by_typ; // list of atoms sorted by atom type
    double *Rpc_by_typ;    // list of Rpc for all atom types
    double *alp_by_typ;
    double *beta_by_typ;
		NAME *atm_typ;         // names of the atom types

    //----------------
    //con-destruct:
    ATOM_MAPS(){
			natm_typ         = 1;
			natm             = 1;
			natm_atm_typ_max = 1;
      lmax             = 1;
      lmax_by_typ      = NULL;
			index_atm_typ    = NULL;
			natm_atm_typ     = NULL;
			list_atm_by_typ  = NULL;
			atm_typ          = NULL;
			Rpc_by_typ       = NULL;
			alp_by_typ       = NULL;
			beta_by_typ      = NULL;
    };
    ~ATOM_MAPS(){};
    //----------------
    //Allocate routine
     void allocate(){
       natm_atm_typ    = new int [natm_typ];
       index_atm_typ   = new int [natm];
       lmax_by_typ     = new int [natm_typ];
       list_atm_by_typ = new int *[natm_typ];
       for (int i=0; i<natm_typ; i++) { list_atm_by_typ[i] = new int [natm_atm_typ_max]; }
       atm_typ         = new NAME [natm_typ];
			 Rpc_by_typ      = new double [natm_typ];
			 alp_by_typ      = new double [natm_typ];
			 beta_by_typ     = new double [natm_typ];
     }//end routine

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | natm_typ;
      p | natm;
      p | natm_atm_typ_max;
      p | lmax;
     //pupping 1d arrays
      if (p.isUnpacking()) { allocate(); }
      PUParray(p,index_atm_typ,natm);
      PUParray(p,natm_atm_typ,natm_typ);
      PUParray(p,lmax_by_typ,natm_typ);
      PUParray(p,Rpc_by_typ,natm_typ);
      PUParray(p,beta_by_typ,natm_typ);
      PUParray(p,alp_by_typ,natm_typ);
     //pupping 2d arrays
      pup2d_int(p,list_atm_by_typ, natm_typ, natm_atm_typ_max);
     //pupping atom names
      pup1d_name(p, atm_typ, natm_typ);
#ifdef _PARALLEL_DEBUG_        
      state_class_out ();
#endif      
    } // end pup
#endif

    //----------------
    // Out funtion for debugging:
    void state_class_out(){
      int i;
      char fileName [255];
      sprintf (fileName, "%d_atom_maps.state.out", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      // ints
      fprintf(fp,"natm_typ %d\n",natm_typ);
      fprintf(fp,"natm %d\n",natm);
      fprintf(fp,"lmax %d\n",lmax);
      fprintf(fp,"natm_atm_typ_max %d\n",natm_atm_typ_max);
      // 1d  arrays
      for(i=0;i<natm;i++) {fprintf(fp,"index_atm_typ[%d]= %d\n",i,index_atm_typ[i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"natm_atm_typ[%d]= %d\n",i,natm_atm_typ[i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"lmax_by_typ[%d]= %d\n", i,lmax_by_typ[i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"Rpc_by_typ[%d]= %g\n",  i,Rpc_by_typ[i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"alp_by_typ[%d]= %g\n",  i,alp_by_typ[i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"beta_by_typ[%d]= %g\n", i,beta_by_typ[i]);}
      // 2d array
      for(i=0;i<natm_typ;i++) {
	for(int j=0;j<natm_atm_typ[i];j++) {
  	  fprintf(fp,"list_atm_by_typ[%d][%d]= %d\n",i,j,list_atm_by_typ[i][j]);
      }}
      // name of atom types
      for(i=0;i<natm_typ;i++) {fprintf(fp,"atm_typ[%d]= %s\n",i,atm_typ[i]);}
      fclose(fp);
    }// end routine
}; // ATOM_MAPS;

#ifdef PUP_ON
PUPmarshall(ATOM_MAPS);
#endif

#endif

//==========================================================================
