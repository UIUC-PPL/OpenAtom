//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                  GWBSE simulation options                                   
//
//                class definition for GW_PARALLEL
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
#ifndef _GW_PARALLEL_
#define _GW_PARALLEL_


class GW_PARALLEL{

  //---------------
  public:
    unsigned K, L, M;      // Number of k points, occupied, and unoccupied psis
    unsigned n_qpt;          // Number of q points
    unsigned *Q;              // Discrete q points
    unsigned n_elems;         // Number of elements in psi
    unsigned pipeline_stages; // Number of stages in the M pipeline
    unsigned rows_per_chare;  // Rows per PMatrix chare
    unsigned cols_per_chare;  // Columns per PMatrix chare
    unsigned transpose_stages; // Number of stages in transpose
    unsigned proc_rows;        //processor rows for diagonalizer tiling
    unsigned proc_cols;        //processor cols for diagonalizer tiling

    int fft_nelems[3];        // Num: size of FFT

    //----------------
    //con-destruct:
    GW_PARALLEL(){
      K = L = M = 0;
      n_elems = 0;
      n_qpt = 0;
      pipeline_stages = 0;
      rows_per_chare = 0;
      cols_per_chare = 0;
      transpose_stages = 0;
      proc_rows = 1;
      proc_cols = 1;
      Q = new unsigned[1];
      for (int i; i<3; i++){ fft_nelems[i] = 0; }
    };
    ~GW_PARALLEL(){};

#ifdef PUP_ON
    //-------------
    //pupping
    void pup(PUP::er &p){
      p | K; p | L; p | M;
      p | n_elems;
      p | n_qpt;
      p | pipeline_stages;
      p | rows_per_chare;
      p | cols_per_chare;
      p | transpose_stages;
      p | proc_rows;
      p | proc_cols;

      PUParray(p, fft_nelems, 3);
      PUParray(p, Q, n_qpt);
      
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif
    } // end pup
#endif

    void state_class_out(){
      char fileName [255];
      sprintf ( fileName, "%d_gw_parallel.out", CKMYPE());
      FILE *fp; fp = fopen(fileName,"w");
      //dbles
      fprintf(fp,"K %i\n", K);
      fprintf(fp,"L %i\n", L);
      fprintf(fp,"M %i\n", M);
      fprintf(fp,"Q %i\n", Q);
      fprintf(fp,"pipeline_stages %i\n", pipeline_stages);
      fprintf(fp,"rows_per_chare %i\n", rows_per_chare);
      fprintf(fp,"cols_per_chare %i\n", cols_per_chare);
      fprintf(fp,"transpose_stages %i\n", transpose_stages);
      fprintf(fp,"proc_rows %i\n", proc_rows);
      fprintf(fp,"proc_cols %i\n", proc_cols);
      fprintf(fp,"fft size %d  %d  %d\n", fft_nelems[0], fft_nelems[1], fft_nelems[2]);
      fclose(fp);
    }// end routine
    
}; // GW_PARALLEL

#ifdef PUP_ON
PUPmarshall(GW_PARALLEL);
#endif

#endif
//==========================================================================
