//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                  GWBSE simulation options                                   
//
//                class definition for GW_EPSILON
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
#ifndef _GW_IO_
#define _GW_IO_

struct IOConfig {
  bool read, write, verify;
  std::string read_prefix, write_prefix, verify_prefix;
  void clear() {
    read = write = verify = false;
    read_prefix = write_prefix = verify_prefix = "";
  }
  void pup(PUP::er& p) {
    p | read;
    p | write;
    p | verify;
    p | read_prefix;
    p | write_prefix;
    p | verify_prefix;
  }
};

class GW_IO {
  public:
    IOConfig p_matrix, epsilon, epsilon_inv;

    GW_IO() {
      p_matrix.clear();
      epsilon.clear();
      epsilon_inv.clear();
    };

    ~GW_IO(){};

    void pup(PUP::er &p) {
      p | p_matrix;
      p | epsilon;
      p | epsilon_inv;
      
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking()) {
        state_class_out();
      }
#endif
    } // end pup

    void state_class_out() {
      char fileName [255];
      sprintf(fileName, "%d_gw_io.out", CkMyPe());
      FILE* fp = fopen(fileName, "w");

      fprintf(fp,"PMatrix read: %i \"%s\"\n",
          p_matrix.read, p_matrix.read_prefix.c_str());
      fprintf(fp,"PMatrix write: %i \"%s\"\n",
          p_matrix.write, p_matrix.write_prefix.c_str());
      fprintf(fp,"PMatrix verify: %i \"%s\"\n",
          p_matrix.verify, p_matrix.verify_prefix.c_str());

      fprintf(fp,"Epsilon read: %i \"%s\"\n",
          epsilon.read, epsilon.read_prefix.c_str());
      fprintf(fp,"Epsilon write: %i \"%s\"\n",
          epsilon.write, epsilon.write_prefix.c_str());
      fprintf(fp,"Epsilon verify: %i \"%s\"\n",
          epsilon.verify, epsilon.verify_prefix.c_str());

      fprintf(fp,"Epsilon Inverse read: %i \"%s\"\n",
          epsilon_inv.read, epsilon_inv.read_prefix.c_str());
      fprintf(fp,"Epsilon Inverse write: %i \"%s\"\n",
          epsilon_inv.write, epsilon_inv.write_prefix.c_str());
      fprintf(fp,"Epsilon Inverse verify: %i \"%s\"\n",
          epsilon_inv.verify, epsilon_inv.verify_prefix.c_str());

      fclose(fp);
    } // end state_class_out
    
}; // GW_IO

PUPmarshall(GW_IO);

#endif
//==========================================================================
