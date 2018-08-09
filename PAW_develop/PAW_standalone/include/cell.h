//==========================================================================
//                  Simulation cell data                                    

#ifndef _CELL_
#define _CELL_

class CELL{
  //----------------
  public:
    int iperd;            // periodicity
    int nimg;             // number of images
    double animg;         // (double) nimg
    double alpb;          // Ewald alpha
    double Rcut;          // Ewald cutoff real space
    double Ecut;          // planewave cutoff in Ryd
    double gcut;          // state g space
    double Gcut;          // density g cutoff, twice gcut
    double Rpc_max;       // max PAW radius of any atom
    double hmat[10];      // the simulation box
    double hmati[10];     // inverse simulation box
    double volume;        // simulation box volume
    double beta_unitless; // Coulomb screening parameter
    //----------------
    //con-destruct:
    CELL(){
	iperd          = 3;
	nimg           = 0;
	animg          = 0.0;
	alpb           = 1.0;
	Rcut           = 1.0;
	gcut           = 1.0;
	Gcut           = 1.0;
	Rpc_max        = 0.0;
	for (int i=0; i<10; i++) {hmat[i] = 0.0; hmati[i] = 0.0;}
	hmat[1] = 1.0; hmat[5] = 1.0; hmat[9] = 1.0;
	hmati[1] = 1.0; hmati[5] = 1.0; hmati[9] = 1.0;
	volume         = 1.0;
    };
    ~CELL(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | iperd;
      p | nimg;
      //pupping dbles
      p | animg;
      p | alpb;
      p | Rcut;
      p | gcut;
      p | Gcut;
      p | volume;
      p | Rpc_max;
      p | beta_unitless;
      //pupping dbl  arrays
      PUParray(p,hmat,10);
      PUParray(p,hmati,10);
#ifdef _PARALLEL_DEBUG_        
      state_class_out ();
#endif      
   } // end pup
#endif

    void state_class_out(){
      int i;
      char fileName [255];
      sprintf (fileName, "%d_cell.state.out", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      // ints
      fprintf(fp,"iperd %d\n",iperd);
      fprintf(fp,"nimg %d\n",nimg);
      // dbles
      fprintf(fp,"animg %g\n",animg);
      fprintf(fp,"alpb %g\n",alpb);
      fprintf(fp,"Rcut %g\n",Rcut);
      fprintf(fp,"gcut %g\n",gcut);
      fprintf(fp,"Gcut %g\n",Gcut);
      fprintf(fp,"volume %g\n",volume);
      fprintf(fp,"Rpc_max %g\n",Rpc_max);
			fprintf(fp,"beta_unitless %g\n",beta_unitless);
      // dbl  arrays
      for(i=1;i<=9;i++){fprintf(fp,"hmat[%d], %lg\n",i,hmat[i]);}
      for(i=1;i<=9;i++){fprintf(fp,"hmati[%d], %g\n",i,hmati[i]);}
      fclose(fp);
    }// end routine
}; // CELL;

#ifdef PUP_ON
PUPmarshall(CELL);
#endif

#endif
//==========================================================================
