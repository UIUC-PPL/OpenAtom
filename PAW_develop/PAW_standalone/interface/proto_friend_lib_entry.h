#ifndef _FRIEND_LIB_
#define _FRIEND_LIB_

#include "ckcomplex.h"

/*---------------------------------------------------------------------*/
/*     cmalloc.c                                                       */

void *cmalloc(size_t,const char * );
void cfree(void *,const char *);
void *crealloc(void *,size_t,const char * );
int **cmall_int_mat(long,long,long,long,const char *);
double **cmall_mat(long,long,long,long,const char *);
double ***cmall_tens3(long,long,long,long,long,long,const char *);
double ****cmall_tens4(long,long,long,long,long,long,long,long,const char *);
int ****cmall_itens4(long,long,long,long,long,long,long,long,const char *);
int **creall_int_mat(int **,long,long,long,long ,long ,long ,long ,long ,const char *);
double **creall_mat(double **,long,long,long,long ,long ,long ,long ,long ,const char *);
void cfree_mat(double **,long , long , long , long );
void cfree_int_mat(int **,long , long , long , long );
void cfree_tens3(double ***,long ,long ,long ,long ,long ,long );
int ***cmall_itens3(long ,long ,long ,long ,long ,long , const char *);
void cfree_itens3(int ***,long ,long ,long ,long ,long ,long );

/*---------------------------------------------------------------------*/
/*     friend_lib.c                                                    */

void spline_fit(double *,double * ,double * ,double * ,double * ,int );
FILE *cfopen(const char [], const char *);

#endif
