#ifndef __MAIN_H__
#define __MAIN_H__

class Main : public CBase_Main {

 public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);
  
  /// Entry Methods ///
  void done();
  void StartHi(int elems);
};

class EpsMap : public CkArrayMap {
    int _x, _y, pe, rem, quo;
  public:
    EpsMap(MatrixConfig cfg2D){
      _x = cfg2D.mat_rows/cfg2D.tile_rows;
      _y = cfg2D.mat_cols/cfg2D.tile_cols;
      int total = _x*_y;
      int node_count = CkNumNodes();
      if(node_count > _x*_y) node_count = _x*_y;
      rem = total%node_count;
      quo = total/node_count;
    }
    ~EpsMap(){}
    int procNum(int, const CkArrayIndex &idx) {
      int count = idx.data()[0]*_y +idx.data()[1];

      int chares_on_node = quo;
      if(count < quo+rem) chares_on_node += rem;

      int node_size = CkNumPes()/CkNumNodes();

      int chares_per_pe = ceil((double)chares_on_node/(double)node_size);
      if(chares_per_pe == 0) chares_per_pe = 1;

      if(count < quo+rem){ //node = 0
        int rank = count/chares_per_pe;
        return rank;
      }

      int sub = count - (quo+rem);
      int node = sub/quo +1;
      int count_on_node = sub%quo;
      pe = node*node_size;
      int rank = count_on_node/chares_per_pe;
      return pe+rank;
    }
};


#endif //__MAIN_H__
