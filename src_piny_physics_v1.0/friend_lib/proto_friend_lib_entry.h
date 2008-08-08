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

/*---------------------------------------------------------------------*/
/* piny_pup.C */

#ifdef PUP_ON
void pup1d_int(PUP::er &,int **, int );
void pup1d_dbl(PUP::er &,double  **, int );
void pup1d_cpl(PUP::er &,complex **, int );
void pup1d_char(PUP::er &,char **, int );
void pup2d_int_test(PUP::er &,int ***, int ,int);
void pup2d_int(PUP::er &,int ***, int ,int);
void pup2d_dbl(PUP::er &,double ***, int ,int,const char *);
void pup1d_name(PUP::er &,NAME **, int );
void pup3d_dbl(PUP::er &,double ****, int ,int ,int );
#endif

