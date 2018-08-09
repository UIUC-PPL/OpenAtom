#include "standard_include.h"
#include "pup_utils.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 2d char arrays                                                           */
/*==========================================================================*/
void pup1d_name(PUP::er &p,NAME *vec, int nlen){
  int i;
  char *scratch;

  if(nlen > 0){
    scratch = new char [MAXWORD];
    for(i=0;i<nlen;i++){
      if (p.isPacking()) {strcpy(scratch,vec[i]);}
      p(scratch,MAXWORD);
      if (p.isUnpacking()) {strcpy(vec[i],scratch);}
    }/*endfor*/
    delete [] scratch;
  }/* endif */

}/*end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* 2d integer arrays                                                        */
/*==========================================================================*/
void pup2d_int(PUP::er &p,int **vec, int nlen_1,int nlen_2){

  if(nlen_1 > 0 && nlen_2 > 0){
    int ir,ic;
    int index;
    int  psize    = nlen_1*nlen_2;
    int *pscratch = new int [psize];
    if (p.isPacking()) {
      // packing 2d array into 1d array to be pupped
      index = 0;
      for(ir=0; ir < nlen_1; ir++){
        for(ic=0; ic < nlen_2; ic++){
          pscratch[index] = vec[ir][ic];
          ++index;
        }}
    }// endif
    p(pscratch,psize);
    // unpacking 1d array into 2d array
    if (p.isUnpacking()) {
      index = 0;
      for(ir=0; ir < nlen_1; ir++){
        for(ic=0; ic < nlen_2; ic++){
          vec[ir][ic] = pscratch[index];
          ++index;
        }}
    }//endif
    delete [] pscratch;
  }// endif nlen 

}// end routine 
/*==========================================================================*/

