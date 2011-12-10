/*==========================================================================*/

void set_cpmass(int ,int *,int *,int *,double *,double *,
		double *,double ,int *);

void calc_cutoff(int , double *,double *,int ,int *, int *,
                 double *, double ,int *,int );

/*==========================================================================*/
/*               Proto for check_kvec:                                      */ 

void check_kvec(int ,int [],int [],int [],
		int ,int [],int [],int []);

/*==========================================================================*/
/*               Proto for setkvec3d:                                      */ 

void setkvec3d(int ,double ,int *,double *,
		int *, int *, int *, 
                int *, int *, int , 
                double *, double *,double *);

void setkvec3d_res(int , double *, 
		  int *, int *, int *, int *, 
                  int , int );

/*==========================================================================*/
/*               Proto for setkvec3d_sm:                                    */ 

void setkvec3d_sm(int ,double ,int *,double *,int *, int *, int *, 
		  int *, int *, double *, double *,int );

void setkvec3d_sm_kpt(int ,double ,int *,double *,int *, int *, int *, 
		  int *, int *, double *, double *,int );

/*==========================================================================*/
/*               Proto for radixme:                                         */ 

void radixme(int , int , int ,int *, int *, int *,int);

/*==========================================================================*/
/*               Proto for countkvec3d:                                     */ 

void countkvec3d(int *,double ,int *,double *,double *,double *,double *);

/*==========================================================================*/
/*               Proto for countkvec3d_sm:                                  */ 

void countkvec3d_sm(int *, double , int *, double * , double *, double *);
void countkvec3d_sm_kpt(int *, double , int *, double * , double *, double *);

/*==========================================================================*/
/*               Proto for makemap_cp:                                      */ 
void makemap_cp(int *, double *, 
		int , int , int , int , int , int ,
		double , int *, int *, int *, 
		int *, int *, int *, int *, 
		int );

void makemap_cp_f_(int *, double *, 
		int *, int *, int *, int *, int *, int *,
		double *, int *, int *, int *, 
		int *, int *, int *, int *, 
		int *);

/*==========================================================================*/
/*               Proto for countmap_cp:                                    */ 
void countmap_cp(int *, double *, 
		 int , int , int , int , 
		 int , int , double );


/*==========================================================================*/
/*               Proto for set_pme_grid:                                    */ 

void set_pme_grid(double ,double ,double *,int *,
                  int *,int *,int *,int ,int );

void set_pme_wght(int ,int *,int *,int *,
                   int ,int ,int ,int,
                   int ,int ,int ,
                   double *,double *,
                   double *,int ,
                   double *,double *,double *);


void get_bspline_wght1d(int ,int ,double *,double *,
                   double *,double *,double *,
                   double *,double *,double *);

/*==========================================================================*/


void init_nonlocal_ees(int *,double, PSNONLOCAL *,int);
void init_eext_ees(int *,CPPSEUDO *,int);
void set_fftsizes(int, int *, int *,int);

/*==========================================================================*/
