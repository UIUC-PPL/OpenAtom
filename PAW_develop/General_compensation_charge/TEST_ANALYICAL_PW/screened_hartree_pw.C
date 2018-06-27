//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files
//==========================================================================
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>
#include <time.h>
using namespace std;

//==========================================================================
// General constants and Functions

#define M_PI_G 3.14159265358979323846264338327950288419716939937510
#define RT_PI_G (sqrt(M_PI_G))
#define MAX(A,B) (((A)>(B))?(A):(B))
#define MIN(A,B) (((A)<(B))?(A):(B))

//==========================================================================
// Structures

typedef struct ATOM_INFO{
  int natm;                      // int  : number of atms
  int natm_typ;                  // int  : number of atm types
  int lmax_psi_max;              // int  : max number of psi angular momentum channels of any atom type
  double beta_unitless;          // dble : unitless inverse erf screening length
  int *iatm_atm_typ;             // list [natm]    : atom type of every atom
  int *lmax_psi;                 // list [natm_typ]: max channel of this atom type
  double *rmax;                  // list [natm_typ]: rmax of this atom type
  double *alpha;                 // list [natm_typ]: gauss width of this atom type
  double *beta;                  // list [natm_typ]: inverse erf screening length of this atom type
}ATOM_INFO;
  
typedef struct FGRID_GEN{
  int nr;                        // int  : number of r integration points
  int ntheta;                    // int  : number of cos(theta) integration points
  int nphi;                      // int  : number of phi integration points
  int ntot;                      // int  : tot num integration points (nr+1)*ntheta*nphi
  int nang;                      // int  : angular integration points ntheta*nphi
  double *r,*wr,*j0;             // list[nr+1]   : r integration points and weights
  double *cost,*sint,*wCost;     // list[ntheta] : cos(theta) integration points and weights
  double *phi,*wPhi;             // list[nphi]   : phi integration points and weights
  double *xf,*yf,*zf,*wf;        // list[ntot]   :x,y,z,w on full grid -- r inner index
  double *wCostPhi;              // list[nang]   : weight of angular integration points
}FGRID_GEN;

typedef struct YLM_DATA{
  int lang;                      // int  : angular momentum component of this guy
  int ntheta;                    // int  : number of cos(theta) integration points
  int nphi;                      // int  : number of phi integration points
  int nang;                      // int  : angular integration points ntheta*nphi
  double **ylm;                  // int [2l+1][nang] : Y_lang,m(cost,phi)
}YLM_DATA;
  
typedef struct PW_ERF_CONS {
  int lmax;                      // int : max number of erf pw channels
  int ng;                        // int : number of GL integation points
  double *i_n;                   // list[lmax+1] : mod spher bessel fnc 1st kind small arg
  double *i_n_0;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_1;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_2;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_3;                 // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind small arg
  double *i_n_1r;                // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind large arg
  double *i_n_2r;                // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind large arg
  double *i_n_3r;                // list[lmax+1] : expansion coef of mod spher bessel fnc 1st kind large arg
  double *pref;                  // list[lmax+1] : partial wave pre facotr
  double *pw_erf_num;            // list[lmax+1] : partial wave for one given {r>,r<}
  double *pw_erf_math;           // list[lmax+1] : partial wave for one given {r>,r<}
  double *pw_erf_big;            // list[lmax+1] : partial wave for one given {r>,r<}
  double *gx,*wx;                // list[ng] : GL integation pts and nodes
}PW_ERF_CONS ;
  
typedef struct PW_ERF_L{
  int nr;                        // int : number of r points
  int lang;                      // int : angular mo of this expansion coef
  double pre_l;                  // dble : 4*pi/(2*lang+1)
  double **pw_erf;               // list [nr][nr] : partial wave coefs
}PW_ERF_L;
  
typedef struct PW_ERF_DATA{
  int lmax;                      // int  : number of partial waves
  int nr;                        // int  : number of r integration points
  double beta;                   // dble : inverse erf screening length
  double rmax;                   // dble : upper limit of r integration
  double *r,*wr;                 // list[nr+1] : r integration points
  PW_ERF_L *pw_erf_L;            // list[lmax+1] : partial wave expansion coefs
}PW_ERF_DATA;
  
typedef struct RHO_SCR{
  int nr;                        // int  : number of r integration points
  int ntheta;                    // int  : number of cos(theta) integration points
  int nphi;                      // int  : number of phi integration points
  int nang;                      // int  : angular integration points ntheta*nphi
  double *rho_lm;                // list[nr+1] : lm component in "r" of density  
  double *wCostPhi;              // list[nang] : weight of dcos(theta)dphi integration
  double **rho;                  // list[nr+1][nang] : density stored nicely for pw method
}RHO_SCR;
  
typedef struct RHO_ATM_DATA{
  int nr;                        // int  : number of r integration points
  int ntheta;                    // int  : number of cos(theta) integration points
  int nphi;                      // int  : number of phi integration points
  int ntot;                      // int  : tot num integration points (nr+1)*ntheta*nphi
  int iatm;                      // int  : whose atom density am I
  int lmax_psi;                  // int  : lmax_psi for this atom
  int lmax_rho;                  // int  : lmax_rho=2lmax_psi for this atom < lmax=lmax_pw
  int lmax_pw;                   // int  : partial wave lmax
  double rmax;                   // dble : maximum r
  double *rho;                   // list[ntot] : density around the atom
}RHO_ATM_DATA;

//==========================================================================
// Local funcitons
  int  main (int , char *[]);
  void init_atom_info(ATOM_INFO *,double );
  void init_fgrid_gen(FGRID_GEN *,int ,int ,int );
  void get_Ylm(YLM_DATA *,FGRID_GEN*,int);
  void init_pw_erf_cons(PW_ERF_CONS *,int); 
  void init_pw_erf_data(PW_ERF_CONS *,PW_ERF_DATA *,FGRID_GEN *,double ,double);
  void init_rho_scr(RHO_SCR *,FGRID_GEN *);
  void get_rho_atm(RHO_ATM_DATA *,FGRID_GEN *,YLM_DATA *,int ,double,int,int,double);
  void compute_self_screen_Hartree(RHO_ATM_DATA *,YLM_DATA *,PW_ERF_DATA *,RHO_SCR *,double *);
  void ctrl_set_erf_partial_wave(double , double ,  double , PW_ERF_CONS *, int *, int *);
  void pw_erf_analytic(double *, double , double , double , int , int *);
  void readtoendofline(FILE *fp);

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller to manage compensation charge energy
//==========================================================================
int main (int argc, char *argv[]){
//==========================================================================  
// Grid and lmax_pw inputs and constants

  printf("\n");
  if(argc<5){
    printf("=================================================================\n");
    printf(" Command line argument: code.exe nr ntheta lmax_pw beta_unitless \n");
    printf("=================================================================\n");
    exit(1);
  }//endif
  int nr               = atoi(argv[1]);
  int ntheta           = atoi(argv[2]);  // ntheta >= lmax_pw+1 to integrate lmax+pw waves in cos(theta)
  int nphi             = 2*ntheta;       // nphi >=2*lmax_pw+1 so 2*ntheta - 1 works and hence so would 2*ntheta
  int lmax_pw          = atoi(argv[3]);  // partial wave lmax
  double beta_unitless = atof(argv[4]);  // dimensionless inverse screening length of erf

//==================================================================================
// Initialize atom stuff

  ATOM_INFO atom_info;
  init_atom_info(&atom_info,beta_unitless);
  int natm          = atom_info.natm;
  int natm_typ      = atom_info.natm_typ;
  int lmax_psi_max  = atom_info.lmax_psi_max;
  int *iatm_atm_typ = atom_info.iatm_atm_typ;

//==================================================================================
// Need enough partial wave for square of wave function so 2x more than you think
  if(lmax_pw<2*lmax_psi_max){ 
    printf("=============================================================\n");
    printf(" lmax_pw of partial wave expansion must be >= 2*lmax_psi of states\n");
    printf("=============================================================\n");
    exit(1);
  }//endif
//==================================================================================
// Cos(theta): need to integrate Plmax_pw(cost) * Plmax_pw(cost)*= poly degree 2*lmax_pw
//    2*ntheta-1 = 2*lmax_pw  : ntheta = (2*lmax_pw+1)/2 : ntheta = lmax_pw + 1   smallest integer satisfies
  if(ntheta <= lmax_pw){
    printf("=============================================================\n");
    printf(" cos(theta) GL integration points must be > lmax_pw\n");
    printf("=============================================================\n");
    exit(1);
  }//endif

//==================================================================================
// Initialize genric grid  (remove scaling factor of r integration and its generic)
   FGRID_GEN fgrid_gen;
   init_fgrid_gen(&fgrid_gen,nr,ntheta,nphi);

//==================================================================================
// Initialize spherical harmonic stuff (real linear combos) does not depend on atm typ
   YLM_DATA *ylm_data = new YLM_DATA [(lmax_pw+1)];
   get_Ylm(ylm_data,&fgrid_gen,lmax_pw);

//==================================================================================
// Initialize partial wave stuff including scratch for density

  PW_ERF_CONS pw_erf_cons;
  init_pw_erf_cons(&pw_erf_cons,lmax_pw); // does not depend on atm typ

  PW_ERF_DATA *pw_erf_data = new PW_ERF_DATA [natm_typ];
  double *beta = atom_info.beta;
  double *rmax = atom_info.rmax;
  for(int iatm_typ=0;iatm_typ<natm_typ;iatm_typ++){
    pw_erf_data[iatm_typ].pw_erf_L = new PW_ERF_L [(lmax_pw+1)];
    init_pw_erf_data(&pw_erf_cons,&pw_erf_data[iatm_typ],&fgrid_gen,beta[iatm_typ],rmax[iatm_typ]);
  }//endfor
  double *gx = pw_erf_cons.gx;  double *wx = pw_erf_cons.wx;
  delete [] gx;  delete [] wx;


//==================================================================================
// Initialize electron density around each atom :

  RHO_SCR rho_scr; // same for everyone
  init_rho_scr(&rho_scr,&fgrid_gen);

  RHO_ATM_DATA *rho_atm_data = new RHO_ATM_DATA [natm];
  int *lmax_psi              = atom_info.lmax_psi;
  for(int iatm=0;iatm<natm;iatm++){
    int iatm_typ = iatm_atm_typ[iatm];
    get_rho_atm(&rho_atm_data[iatm],&fgrid_gen,ylm_data,lmax_psi[iatm_typ],rmax[iatm_typ],iatm,lmax_pw,beta[iatm_typ]);
  }//endfor

//==================================================================================
// Compute Hartree self by partial wave expansion
  double vself = 0.0;
  for(int iatm=0;iatm<natm;iatm++){
    int iatm_typ   = iatm_atm_typ[iatm];
    double vself_i = 0.0;
    compute_self_screen_Hartree(&rho_atm_data[iatm],ylm_data,&pw_erf_data[iatm_typ],&rho_scr,&vself_i);
    printf("Vself(iatm %d iatm_typ %d) %.10g\n",iatm,iatm_typ,vself_i);
    vself  += vself_i;
  }//endfor
  printf("\n");

//------------------------------------------------------------------------------
   return 1;
 }// end routine : main
//==================================================================================


//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void init_fgrid_gen(FGRID_GEN *fgrid_gen,int nr,int ntheta,int nphi){
//===============================================================================================
// Allocate and compute
  int nr1  = nr+1;
  int ntot = nr1*ntheta*nphi;
  int nang = ntheta*nphi;
  double *r        = new double [nr1];
  double *j0       = new double [nr1];
  double *wr       = new double [nr1];
  double *cost     = new double [ntheta];
  double *sint     = new double [ntheta];
  double *wCost    = new double [ntheta];
  double *phi      = new double [nphi];
  double *wPhi     = new double [nphi];
  double *xf       = new double [ntot];
  double *yf       = new double [ntot];
  double *zf       = new double [ntot];
  double *wf       = new double [ntot];
  double *wCostPhi = new double [nang];

//===============================================================================================
// r grid  : generic rmax = 1.0
  double dr = 1.0/((double)nr);
  for(int ir=0;ir<=nr;ir++){
    double dir = (double)ir;
    r[ir]      = dr*dir;
    wr[ir]     = dr*r[ir]*r[ir];
    double arg = r[ir]*M_PI_G;
    j0[ir]     = (ir!= 0 ? sin(arg)/arg : 1.0);
    j0[ir]    *= sqrt(2.0)*M_PI_G;
  }//endfor
  wr[0]    *= 0.5;
  wr[nr]   *= 0.5;

//===============================================================================================
// theta grid
  int ng;
  char fname[100];
  sprintf(fname,"gl_w_x_%d.dat",ntheta);
  FILE *fp = fopen(fname,"r");
  if(fp==NULL){
    printf("=============================================================\n");
    printf(" File %s not found\n",fname);
    printf("=============================================================\n");
    exit(1);
  }//endif
  int iii = fscanf(fp,"%d",&ng); 
  if(iii!=1){
    printf("=============================================================\n");
    printf(" 1st line of file %s invalid\n",fname);
    printf("=============================================================\n");
    exit(1);
  }//endif'
  readtoendofline(fp);
  if(ng!=ntheta){
    printf("=============================================================\n");
    printf(" File %s does not have %d points but %d\n",fname,ntheta,ng);
    printf("=============================================================\n");
    exit(1);
  }
  for(int ig=0;ig<ng;ig++){
    int jjj = fscanf(fp,"%lg %lg",&cost[ig],&wCost[ig]); 
    if(jjj!=2){
      printf("=============================================================\n");
      printf(" %d th line of file %s invalid\n",ig+2,fname);
      printf("=============================================================\n");
      exit(1);
    }//endif
    readtoendofline(fp);
  }//endfor
  fclose(fp);

  for(int icost=0;icost<ntheta;icost++){
    sint[icost] = sqrt(1.0-cost[icost]*cost[icost]);
  }//endfor

//===============================================================================================
// phi grid  
  double dphi = 2.0*M_PI_G/((double)nphi);
  for(int iphi=0;iphi<nphi;iphi++){
    double diphi = (double)iphi;
    phi[iphi]    = dphi*diphi;
    wPhi[iphi]   = dphi;
  }//endfor

//===============================================================================================
// angular grid  
  int iang = 0;
  for(int icost=0;icost<ntheta;icost++){
  for(int iphi=0;iphi<nphi;iphi++,iang++){
    wCostPhi[iang] = wPhi[iphi]*wCost[icost];
  }}//endfor

//===============================================================================================
// total grid
  int iff = 0;
  for(int icost=0;icost<ntheta;icost++){
  for(int iphi=0;iphi<nphi;iphi++){
  for(int ir=0;ir<=nr;ir++,iff++){
    xf[iff] = r[ir]*sint[icost]*cos(phi[iphi]);
    yf[iff] = r[ir]*sint[icost]*sin(phi[iphi]);
    zf[iff] = r[ir]*cost[icost];
    wf[iff] = wPhi[iphi]*wCost[icost]*wr[ir];
  }}}//endfor
  
//===============================================================================================
//  store in the structure

  fgrid_gen->nr       = nr;
  fgrid_gen->ntheta   = ntheta;
  fgrid_gen->nphi     = nphi;
  fgrid_gen->ntot     = ntot;
  fgrid_gen->nang     = nang;
  fgrid_gen->r        = r;
  fgrid_gen->j0       = j0;
  fgrid_gen->wr       = wr;
  fgrid_gen->cost     = cost;
  fgrid_gen->sint     = sint;
  fgrid_gen->wCost    = wCost;
  fgrid_gen->phi      = phi;
  fgrid_gen->wPhi     = wPhi;
  fgrid_gen->xf       = xf;
  fgrid_gen->yf       = yf;
  fgrid_gen->zf       = zf;
  fgrid_gen->wf       = wf;
  fgrid_gen->wCostPhi = wCostPhi;

//------------------------------------------------------------------------------
  }// end routine
//==================================================================================

//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void get_Ylm(YLM_DATA *ylm_data,FGRID_GEN *fgrid_gen,int lmax){
//===============================================================================================
// Local pointers etc.
  int lmax1        = lmax + 1;
  int nphi         = fgrid_gen->nphi;
  int ntheta       = fgrid_gen->ntheta;
  int nang         = fgrid_gen->nang;
  double *phi      = fgrid_gen->phi;
  double *cost     = fgrid_gen->cost;
  double *sint     = fgrid_gen->sint;
  double *wCostPhi = fgrid_gen->wCostPhi;

//===============================================================================================
// Alloc and fill the ylm structure

  for(int il=0;il<=lmax;il++){
    ylm_data[il].lang   = il;
    ylm_data[il].ntheta = ntheta;
    ylm_data[il].nphi   = nphi;
    ylm_data[il].nphi   = nang;
    int nm = 2*il+1;
    double **ylm   = new double *[nm];
    double *contig = new double [nm*nang];
    for(int m=0;m<nm;m++){ylm[m] = &contig[m*nang];}
    ylm_data[il].ylm = ylm;
  }//endfor

//===============================================================================================
// Allocate some local memory to store speparable parts of ylm

  double **chi    = new double *[lmax1];
  double **chip   = new double *[lmax1];
  double *ccontig  = new double [lmax1*nphi];
  double *ccontigp = new double [lmax1*nphi];
  for(int m=0;m<=lmax;m++){chi[m] = &ccontig[m*nphi]; chip[m] = &ccontigp[m*nphi];}

  double **plm    = new double *[lmax1];
  double *pcontig  = new double [lmax1*ntheta];
  for(int m=0;m<=lmax;m++){plm[m] = &pcontig[m*ntheta];}

//===============================================================================================
// Get the real linear combo chi_m(phi) chi_mp(phi) 0 1 1' 2 2' ... lmax lmax  2 lmax*phi matrices

  for(int iphi=0;iphi<nphi;iphi++){
    chi[0][iphi]  = 1.0;
    chip[0][iphi] = 1.0;
  }//enfor
  for(int m=1;m<=lmax;m++){
    double dm = (double)m;
    for(int iphi=0;iphi<nphi;iphi++){
      chi[m][iphi]  = cos(dm*phi[iphi])*sqrt(2.0);
      chip[m][iphi] = sin(dm*phi[iphi])*sqrt(2.0);
    }//enfor
  }//endfor

//===============================================================================================
// Get the normalized Plm (cost) and then the ylm
// all the P_l m=0 first
  for(int icost=0;icost<ntheta;icost++){
     plm[0][icost] = 1.0;
     plm[1][icost] = cost[icost];
  }//endfor
  for(int il=1;il<=lmax-1;il++){
     double dil  = (double)il;
     double cl   = (2.0*dil+1.0)/(dil+1.0);
     double clm1 = dil/(dil+1.0);
     for(int icost=0;icost<ntheta;icost++){
       plm[il+1][icost] = cl*cost[icost]*plm[il][icost] - clm1*plm[il-1][icost];
     }//endfr cost
  }//endfor : l

  for(int m=0;m<=lmax;m++){
     double dm = (double) m;
     int ind1  = 2*m;
     int ind2  = (m!=0 ? 2*m-1 : 0);
     for(int il=m;il<=lmax;il++){
        double dil = (double)il;
        double rat = 1.0; for(double dj = 1.0;dj<=dm;dj+=1.0){rat *= (dil+dj)*MAX((dil-dj+1.0),1.0);}
        double anorm = sqrt(0.5*(2.0*dil+1.0)/(2.0*M_PI_G*rat));
        int iang  = 0;
        for(int icost=0;icost<ntheta;icost++){
        for(int iphi=0;iphi<nphi;iphi++,iang++){
           ylm_data[il].ylm[ind1][iang] = chip[m][iphi]*plm[il][icost]*anorm;
           ylm_data[il].ylm[ind2][iang] =  chi[m][iphi]*plm[il][icost]*anorm;
        }}//endfor : iang
     }//endfor: il
     for(int il=lmax;il>=m+1;il--){ // next round of m we need one less
        double dil  = (double)il;
        double cm   = dil-dm;
        double cp   = dil+dm;
        for(int icost=0;icost<ntheta;icost++){
           plm[il][icost] = (cm*cost[icost]*plm[il][icost] - cp*plm[il-1][icost])/sint[icost];
        }//endfr cost
     }//endfor : l
  }//endofr m

//==================================================================================
// Test norm of diagonal elements
  for(int il=0;il<=lmax;il++){
    int nm = 2*il+1;
    for(int m=0;m<nm;m++){
        double test = 0.0;
        int iang  = 0;
        for(int icost=0;icost<ntheta;icost++){
        for(int iphi=0;iphi<nphi;iphi++,iang++){
           double tmp1 = ylm_data[il].ylm[m][iang];
           double tmp2 = wCostPhi[iang];
           test += (tmp1*tmp1*tmp2);
       }}//endfor
       if(fabs(test-1.0)>1.0e-10){printf("Norm |ylm(%d %d)|^2 = %g\n",il,m,test);}
    }//endfor: m
  }//endfor : l

  for(int il=0;il<=lmax;il++){
  for(int ilp=0;ilp<=lmax;ilp++){
     int nm  = 2*il+1; int nmp = 2*ilp+1;
     for(int m=0;m<nm;m++){
     for(int mp=0;mp<nmp;mp++){
        double test = 0.0;
        int iang  = 0;
        for(int icost=0;icost<ntheta;icost++){
        for(int iphi=0;iphi<nphi;iphi++,iang++){
           double tmp1 = ylm_data[il].ylm[m][iang];
           double tmp2 = ylm_data[ilp].ylm[mp][iang];
           double tmp3 = wCostPhi[iang];
           test += (tmp1*tmp2*tmp3);
        }}//endfor
        if(il!=ilp && m!=mp && fabs(test)>1.e-10){
           printf("Overlap ylm(%d %d) ylm(%d %d) = %g\n",il,m,ilp,mp,test);
        }//endif
     }}//endfor: m,mp
  }}//endfor : l,lp

//===============================================================================================
// Free local memory - only care about the doubles. 

  delete [] ccontig;
  delete [] ccontigp;
  delete [] pcontig;

//------------------------------------------------------------------------------
  }// end routine
//==================================================================================

//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void init_rho_scr(RHO_SCR *rho_scr,FGRID_GEN *fgrid_gen){
//===============================================================================================
  int nr             = fgrid_gen->nr;
  int nr1            = nr+1;
  int ntheta         = fgrid_gen->ntheta; 
  int nphi           = fgrid_gen->nphi;    
  int nang           = fgrid_gen->nang;    
  double *wCostPhi_f = fgrid_gen->wCostPhi;

  double *rho_lm     = new double [nr1];
  double *wCostPhi   = new double [nang];
  double **rho       = new double *[nr1];
  double *contig     = new double [nr1*nang];
  for(int ir=0;ir<=nr;ir++){rho[ir] = &contig[(ir*nang)];}

  for(int iang=0;iang<nang;iang++){wCostPhi[iang]=wCostPhi_f[iang];}

//===============================================================================================  
// Pack it away

  rho_scr->nr       = nr;
  rho_scr->ntheta   = ntheta;
  rho_scr->nphi     = nphi;
  rho_scr->nang     = nang;
  rho_scr->rho_lm   = rho_lm;
  rho_scr->wCostPhi = wCostPhi;
  rho_scr->rho      = rho;

//------------------------------------------------------------------------------
  }// end routine
//==================================================================================
 

//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void get_rho_atm(RHO_ATM_DATA *rho_atm_data,FGRID_GEN *fgrid_gen, YLM_DATA *ylm_data,
                 int lmax_psi,double rmax,int iatm, int lmax_pw, double beta){
//===============================================================================================
// Local pointers and some allocation
  int nr        = fgrid_gen->nr; 
  int ntheta    = fgrid_gen->ntheta;
  int nphi      = fgrid_gen->nphi;
  int ntot      = fgrid_gen->ntot;
  double *r     = fgrid_gen->r;
  double *cost  = fgrid_gen->cost;
  double *wr    = fgrid_gen->wr;
  double *xf    = fgrid_gen->xf;
  double *yf    = fgrid_gen->yf;
  double *zf    = fgrid_gen->zf;
  double *j0    = fgrid_gen->j0;
  double *wf    = fgrid_gen->wf;

  int lmax_rho  = 2*lmax_psi;
  double *rho   = new double [ntot];

//===============================================================================================
//  Fill rho

  double rmax3 = rmax*rmax*rmax;
  double rt3   = sqrt(3.0);
  double rt5   = sqrt(5.0);
  double anorm  = 0.0;
  double anorm1 = 0.0;
  double anorm2 = 0.0;
  double anorm3 = 0.0;
  double anorm4 = 0.0;
  int iff = 0;
  for(int icost=0;icost<ntheta;icost++){
  for(int iphi=0;iphi<nphi;iphi++){
  for(int ir=0;ir<=nr;ir++,iff++){
     double s_wave   = 1.0;                                        
     double px_wave  = (ir!=0 ? rt3*xf[iff]/r[ir] : 0.0);    
     double py_wave  = (ir!=0 ? rt3*yf[iff]/r[ir] : 0.0);
     double dz2_wave   = 0.5*(3.0*cost[icost]*cost[icost]-1.0)*rt5;
     rho[iff]        = j0[ir]*j0[ir]*(  s_wave*s_wave + px_wave*px_wave + py_wave*py_wave
                                      + dz2_wave*dz2_wave)/(16.0*rmax3*M_PI_G);
     //     rho[iff]        = j0[ir]*j0[ir]*(s_wave*s_wave)/(4.0*rmax3*M_PI_G);
     //     rho[iff]        = j0[ir]*j0[ir]*(dz2_wave*dz2_wave)/(4.0*rmax3*M_PI_G);
     //     rho[iff]        = j0[ir]*j0[ir]*(  s_wave*s_wave + px_wave*px_wave + py_wave*py_wave)/(12.0*rmax3*M_PI_G);
     anorm          += wf[iff]*rho[iff]*rmax3;
     anorm1         += wf[iff]*j0[ir]*j0[ir]*s_wave *s_wave/(4.0*M_PI_G);
     anorm2         += wf[iff]*j0[ir]*j0[ir]*px_wave*px_wave/(4.0*M_PI_G);
     anorm3         += wf[iff]*j0[ir]*j0[ir]*py_wave*py_wave/(4.0*M_PI_G);
     anorm4         += wf[iff]*j0[ir]*j0[ir]*dz2_wave*dz2_wave/(4.0*M_PI_G);
  }}}//endfor
  //  printf("anorm %g %g %g %g %g\n",anorm,anorm1,anorm2,anorm3,anorm4);

#ifdef _COMPUTE_EASY_SLOW_
  double vself_slow = 0.0;
  for(int iff=0;iff<ntot;iff++){
  printf("%d ",iff);
  if((iff % 10)==0){printf("\n");}
  for(int iffp=0;iffp<ntot;iffp++){
    double dx = xf[iff]-xf[iffp];
    double dy = yf[iff]-yf[iffp];
    double dz = zf[iff]-zf[iffp];
    double r  = sqrt(dx*dx+dy*dy+dz*dz)*rmax;
    double vinter = 0.0;
    if(r>1e-9){
      vinter = erf(beta*r)/r;
    }else{
      vinter = 2.0*beta/sqrt(M_PI_G);
    }//endif
    vself_slow  += (wf[iff]*rho[iff]*rmax3)*(wf[iffp]*rho[iffp]*rmax3)*vinter;
  }}
  vself_slow *= 0.5;
  printf("Vself the easy but slow way %.10g\n",vself_slow);
#endif

  double vself_swave = 0.0;
  double beta2 = beta*beta;
  for (int ir=0; ir<=nr; ir++) {
  for (int jr=0; jr<=nr; jr++) {
     double rgt = MAX(r[ir],r[jr])*rmax;
     double rlt = MIN(r[ir],r[jr])*rmax;
     double rd = rgt - rlt;
     double rs = rgt + rlt;
     double complicated = 0.0;
     if (ir != 0 && jr != 0) {
        double part1 = (exp(-beta2*rs*rs) - exp(-beta2*rd*rd))/(2.0*beta*sqrt(M_PI_G)*rgt*rlt);
        double part2 = (rd*erfc(beta*rd) - rs*erfc(beta*rs))/(2.0*rgt*rlt);
        double part3 = 1.0/rgt;
        complicated = part1 + part2 + part3;
     } else {
        if (ir !=0 || jr != 0) {
           complicated = erf(beta*rgt)/rgt;
        } else {
           complicated = 2.0*beta/sqrt(M_PI_G);
       }// end if
     } // end if
     vself_swave += (wr[ir]*j0[ir]*j0[ir])*(wr[jr]*j0[jr]*j0[jr])*complicated;
  }}//endfor
  vself_swave *= 0.5;
  printf("Vself for s-wave only %.10g\n",vself_swave);

//===============================================================================================
// Pack into the stucture

  rho_atm_data->nr       = nr;
  rho_atm_data->ntheta   = ntheta;
  rho_atm_data->nphi     = nphi;
  rho_atm_data->ntot     = ntot;
  rho_atm_data->iatm     = iatm;
  rho_atm_data->lmax_psi = lmax_psi;
  rho_atm_data->lmax_rho = lmax_rho;
  rho_atm_data->lmax_pw  = lmax_pw;
  rho_atm_data->rmax     = rmax;
  rho_atm_data->rho      = rho;

//-----------------------------------------------------------------------------------------------
  }// end routine
//===============================================================================================

//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void init_pw_erf_data(PW_ERF_CONS *pw_erf_cons, PW_ERF_DATA *pw_erf_data,FGRID_GEN *fgrid_gen,
                      double beta,double rmax){
//===============================================================================================
// Local pointers 

  int lmax            = pw_erf_cons->lmax;
  int nr              = fgrid_gen->nr;
  double *rgen        = fgrid_gen->r;
  double *wrgen       = fgrid_gen->wr;
  double *pw_erf_num  = pw_erf_cons->pw_erf_num;
  PW_ERF_L *pw_erf_L  = pw_erf_data->pw_erf_L; // allocated before coming here.

//===============================================================================================
// Allocate memory and setup some of the simle variables

  double *r  = new double [nr+1];
  double *wr = new double [nr+1];
  for(int ir=0;ir<=nr;ir++){
    r[ir]  =  rgen[ir]*rmax;
    wr[ir] = wrgen[ir]*rmax*rmax*rmax;
  }//endfor

  for(int il=0;il<=lmax;il++){
    double dil = (double)il;
    pw_erf_L[il].nr    = nr;
    pw_erf_L[il].lang  = il;
    pw_erf_L[il].pre_l = (4.0*M_PI_G)/(2.0*dil+1.0);    //  (4*pi)/(2*l+1)
  }//endfor

  int nr1 = nr+1;
  for(int il=0;il<=lmax;il++){
    double **pw_erf = new double *[nr1];
    double *contig  = new double [nr1*nr1];
    for(int ir=0;ir<=nr;ir++){pw_erf[ir] = &contig[(ir*nr1)];}
    pw_erf_L[il].pw_erf = pw_erf;
  }//endfor

//===============================================================================================
// Get the partial wave coefs

  double total     =  0.0;
  double math      =  0.0;
  double limit     =  0.0;
  for(int ir=0;ir<=nr;ir++){
  for(int irp=ir;irp<=nr;irp++){
    double rlt = r[ir];
    double rgt = r[irp];
    int math_used  = 0;
    int limit_used = 0;
    if(ir!=0){
      ctrl_set_erf_partial_wave(rlt,rgt,beta,pw_erf_cons,&math_used,&limit_used);
      for(int il=0;il<=lmax;il++){
        pw_erf_L[il].pw_erf[ir][irp] = pw_erf_num[il];     
        pw_erf_L[il].pw_erf[irp][ir] = pw_erf_num[il];     
      }//endfor
    }else{
      limit_used = 0;  math_used = lmax+1; // we know these exactly from math
      for(int il=1;il<=lmax;il++){
        pw_erf_L[il].pw_erf[0][irp] = 0.0;     
        pw_erf_L[il].pw_erf[irp][0] = 0.0;     
      }//endif
      if(irp!=0){
        pw_erf_L[0].pw_erf[0][irp] = erf(beta*rgt)/rgt;
        pw_erf_L[0].pw_erf[irp][0] = erf(beta*rgt)/rgt;
      }else{
        pw_erf_L[0].pw_erf[0][0]  = 2.0*beta/sqrt(M_PI_G);
      }//endif
    }//endif
    limit += (double) limit_used;
    math  += (double) math_used;
    total += (double) (lmax+1);
  }}//endfor : ir,irp
  math  /= total;
  limit /= total;

  printf("Analtyical solution applied to %g percent of partial wave coefficients\n",math*100.0);
  printf("Numerical  solution applied to %g percent of partial wave coefficients\n",(1.0-limit-math)*100.0);
  printf("Large beta limit    applied to %g percent of partial wave coefficients\n\n",limit*100.0);

//==================================================================================
// Pack away

  pw_erf_data->lmax     = lmax;
  pw_erf_data->nr       = nr;
  pw_erf_data->beta     = beta;
  pw_erf_data->rmax     = rmax;
  pw_erf_data->r        = r;
  pw_erf_data->wr       = wr;
  pw_erf_data->pw_erf_L = pw_erf_L;

  //------------------------------------------------------------------------------
 }// end routine
//==================================================================================


//==============================================================================================
// Gauss-Legendre integration to get the pw at r< r>
//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void ctrl_set_erf_partial_wave(double rlt, double rgt, double beta, PW_ERF_CONS *pw_erf_cons, int *math_used,
                               int *limit_used){
//===============================================================================================
// Local pointers
  int lmax           = pw_erf_cons->lmax;
  int ng             = pw_erf_cons->ng;  
  double *i_n        = pw_erf_cons->i_n;
  double *i_n_0      = pw_erf_cons->i_n_0;
  double *i_n_1      = pw_erf_cons->i_n_1; 
  double *i_n_2      = pw_erf_cons->i_n_2;
  double *i_n_3      = pw_erf_cons->i_n_3; 
  double *i_n_1r     = pw_erf_cons->i_n_1r; 
  double *i_n_2r     = pw_erf_cons->i_n_2r;
  double *i_n_3r     = pw_erf_cons->i_n_3r; 
  double *pref       = pw_erf_cons->pref;
  double *gx         = pw_erf_cons->gx;
  double *wx         = pw_erf_cons->wx;
  double *pw_erf_num = pw_erf_cons->pw_erf_num;
  double *pw_erf_math= pw_erf_cons->pw_erf_math;
  double *pw_erf_big = pw_erf_cons->pw_erf_big;
  
//===============================================================================================
// Compute partial waves in the large beta limit
  for(int il=1;il<=lmax;il++){pw_erf_big[il] =1.0;}
  if(rlt>0 && rgt >0){
   double rat = rlt/rgt;
   double arg = rat;
   pw_erf_big[0] = 1.0/rgt;
   for(int il=1;il<=lmax;il++){
     pw_erf_big[il] = rat/rgt;
     rat *= arg;
   }//endfor
  }//endif

//===============================================================================================
// Compute partial waves analytically
  int lmax_math = -1;
  pw_erf_analytic(pw_erf_math,beta,rlt,rgt,lmax,&lmax_math);
  int lmax_use = MIN(lmax_math,lmax);

//===============================================================================================
// Compute the erf(beta*r)/r partial wave using GL integration
  for(int il=0;il<=lmax;il++){pw_erf_num[il] = 0.0; }
    for(int ig=0;ig<ng;ig++){
       double x       = gx[ig]*beta;
       double x2      = x*x;
       double a2      = (rlt*rlt+rgt*rgt);
       double eee     = exp(-a2*x2);
       double arg     = 2.0*x2*rlt*rgt;
       double arg2    = arg*arg;
       double arg4    = arg2*arg2;
       double arg6    = arg2*arg4;
       double sinch_a = (arg > 0.0 ? sinh(arg)/arg : 1.0);
       double cosh_a  = cosh(arg);
       if(arg<30.0){
        //  use general recurion relation at intermediate argument
         if(arg>0.01) {i_n[0] = sinch_a;}
         if(arg2>0.01){i_n[1] = (cosh_a - sinch_a)/arg;}
         double argl_p1 = arg2;
         for(int il=1;il<=lmax-1;il++){
           double al = (double) il;
           if(argl_p1*arg > 0.01){i_n[il+1] = i_n[il-1] - (2.0*al+1.0)*i_n[il]/arg;}
           argl_p1 *= arg;
         }//endfor
        //  Use expansion at small argument
         double argl = 1.0;
         for(int il=0;il<=lmax;il++){
            if(argl*arg<=0.01){
              i_n[il] = (1.0 + arg2*i_n_1[il] + arg4*i_n_2[il] + arg6*i_n_3[il])*argl*i_n_0[il];
            }//endif
            argl *= arg;
         }//endfor
       }else{
         //  Use asymptotic expansion at large argument to compute renormalized bessel: i_l*2*z*exp(-arg)
         eee = 0.5*exp(-a2*x2+arg)/arg;      // absorb renormalization into a renormalized eee
         double argi = 1.0/arg; double argi2 = argi*argi; double argi3 = argi2*argi;
         for(int il=0;il<=lmax;il++){i_n[il] = (1.0 + argi*i_n_1r[il] + argi2*i_n_2r[il] + argi3*i_n_3r[il]);}
       }//endif
       for(int il=0;il<=lmax;il++){pw_erf_num[il] += pref[il]*eee*i_n[il]*wx[ig]*beta;}
     }//endfor: GL points

//=========================================================================
// Check for nans and report

  int nan_cnt_num = 0;
  for(int il=0;il<=lmax;il++){
     if(isnan(pw_erf_num[il]) ||!isfinite(pw_erf_num[il])){nan_cnt_num++;}
  }//endfo

  int nan_cnt_math = 0;
  for(int il=0;il<=lmax_use;il++){
    if(isnan(pw_erf_math[il])|| !isfinite(pw_erf_math[il])){nan_cnt_math++;}
  }//endfo

  if(nan_cnt_num>0) {printf("Numercial partial wave nan detected for (%g %g)\n",rlt,rgt);}
  if(nan_cnt_math>0){printf("Analytic  partial wave nan detected for (%g %g)\n",rlt,rgt);}

#ifdef _COLLECT_NUM_MATH_PW_ERR_
  FILE *fp = fopen("math_num_pw_erf_diff.out","a");
  for(int il=0;il<=lmax_use;il++){
    if(isfinite(pw_erf_math[il]) && isfinite(pw_erf_num[il])){
      if(fabs(pw_erf_math[il]-pw_erf_num[il])>1.0e-11){
        double junk1 = beta*beta*(rgt*rgt+rlt*rlt);
        double junk2 = 2.0*beta*beta*rgt*rlt;
        fprintf(fp,"math diff num %g : pw_%d (%g %g) %.10g %.10g %.10g %.10g %g\n",
                fabs(pw_erf_math[il]-pw_erf_num[il]),il,rlt,rgt,pw_erf_math[il],pw_erf_num[il],junk1,junk2,beta*rlt);
      }//endif
    }//endif
  }//endif
  fclose(fp);
#endif

//=========================================================================
// Use the analtyic results if they are finite
  int iii = 0;
  for(int il=0;il<=lmax_use;il++){
    if(isfinite(pw_erf_math[il])){pw_erf_num[il]=pw_erf_math[il];iii++;}
  }//endif
  math_used[0]=iii;

//===============================================================================================
// Replace any crazy values with the asymptotic result for lack of a better option
  iii = 0;
  for(int il=0;il<=lmax;il++){
    if(!isfinite(pw_erf_num[0])){pw_erf_num[il] = pw_erf_big[il];iii++;}
  }//endfor
  limit_used[0]=iii;

//-----------------------------------------------------------------------------------------------
  }//end routine
//===============================================================================================


//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void  init_pw_erf_cons(PW_ERF_CONS *pw_erf_cons,int lmax){
//===============================================================================================
// constants for computation
  int lmax1          = lmax+1;
  double *i_n        = new double [lmax1];
  double *i_n_0      = new double [lmax1];
  double *i_n_1      = new double [lmax1];
  double *i_n_2      = new double [lmax1];
  double *i_n_3      = new double [lmax1];
  double *i_n_1r     = new double [lmax1];
  double *i_n_2r     = new double [lmax1];
  double *i_n_3r     = new double [lmax1];
  double *pref       = new double [lmax1];
  double *pw_erf_num = new double [lmax1];
  double *pw_erf_math= new double [lmax1];
  double *pw_erf_big = new double [lmax1];
//=======================================================================
// prefactor for pw expansion
  for(int il=0;il<=lmax;il++){
    double al = (double) il;
    pref[il]  = 2.0*(2.0*al+1.0)/sqrt(M_PI_G);
  }// end for

//=======================================================================
// small argument expansion coefs (i.e. in z)
  for(int il=0;il<=lmax;il++){
    double al = (double) il;
    pref[il]  = 2.0*(2.0*al+1.0)/sqrt(M_PI_G);
    i_n_1[il] = 0.5/(2.0*al+3.0);
    i_n_2[il] = 0.125/((2.0*al+3.0)*(2.0*al+5.0));
    i_n_3[il] = 1.0/(48.0*(2.0*al+3.0)*(2.0*al+5.0)*(2.0*al+7.0));
  }//endfor
  i_n_0[0] = 1.0;
  for(int il=1;il<=lmax;il++){
    double al = (double) il;
    i_n_0[il] = i_n_0[il-1]/(2*al+1.0);
  }//endfor

//=======================================================================
// large argument expansion coefs (i.e in 1/z)
  for(int il=0;il<=lmax;il++){
    double nu   = 0.5 + ((double) il); // need bessel order not spherical bessel order so add 1/2
    double fnu2 = 4.0*nu*nu;
    i_n_1r[il]  = -(fnu2-1.0)/8.0; 
    i_n_2r[il]  = -i_n_1r[il]*(fnu2-9.0)/(8.0*2.0);
    i_n_3r[il]  = -i_n_2r[il]*(fnu2-25.0)/(8.0*3.0);
  }//endfor

//===============================================================================================
// grid for GL computation of partial waves

  int ng;
  FILE *fp = fopen("gl_w_x_200.dat","r");
  if(fp==NULL){
    printf("=============================================================\n");
    printf(" File gl_w_x_200.dat not found\n");
    printf("=============================================================\n");
    exit(1);
  }//endif
  int iii = fscanf(fp,"%d",&ng); 
  if(iii!=1){
    printf("=============================================================\n");
    printf(" 1st line of file gl_w_x_200.dat invalid\n");
    printf("=============================================================\n");
    exit(1);
  }//endif
  readtoendofline(fp);
  double *gx = new double [ng];
  double *wx = new double [ng];
  for(int ig=0;ig<ng;ig++){
    int jjj = fscanf(fp,"%lg %lg",&gx[ig],&wx[ig]);
    if(jjj!=2){
      printf("=============================================================\n");
      printf(" %d th line of file gl_w_x_200.dat invalid\n",ig+2);
      printf("=============================================================\n");
      exit(1);
    }//endif
    readtoendofline(fp);
    gx[ig] = (gx[ig]+1.0)*0.5;
    wx[ig] *= 0.5;
  }//endfor
  fclose(fp);

//===============================================================================================
//  Store in structure
  pw_erf_cons->lmax       = lmax;
  pw_erf_cons->i_n        = i_n;
  pw_erf_cons->i_n_0      = i_n_0;
  pw_erf_cons->i_n_1      = i_n_1;
  pw_erf_cons->i_n_2      = i_n_2;
  pw_erf_cons->i_n_3      = i_n_3;
  pw_erf_cons->i_n_1r     = i_n_1r;
  pw_erf_cons->i_n_2r     = i_n_2r;
  pw_erf_cons->i_n_3r     = i_n_3r;
  pw_erf_cons->pref       = pref;
  pw_erf_cons->pw_erf_num = pw_erf_num;
  pw_erf_cons->pw_erf_math= pw_erf_math;
  pw_erf_cons->pw_erf_big = pw_erf_big;

  pw_erf_cons->ng         = ng;
  pw_erf_cons->gx         = gx;
  pw_erf_cons->wx         = wx;

//-----------------------------------------------------------------------------------------------
 }//end routine
//===============================================================================================


//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void compute_self_screen_Hartree(RHO_ATM_DATA *rho_atm_data, YLM_DATA *ylm_data, PW_ERF_DATA *pw_erf_data,
                                 RHO_SCR *rho_scr, double *vself_out){
//================================================================================================
// local pointers
  int nr           = rho_atm_data->nr;
  int ntheta       = rho_atm_data->ntheta;
  int nphi         = rho_atm_data->nphi;
  //  int iatm         = rho_atm_data->iatm;
  double *rho_in   = rho_atm_data->rho;

  int nr_scr       = rho_scr->nr;
  int nang         = rho_scr->nang;
  double *rho_lm   = rho_scr->rho_lm;
  double **rho     = rho_scr->rho;
  double *wCostPhi = rho_scr->wCostPhi; // indendent of rmax, beta etc. so is universal

  int lmax         = pw_erf_data->lmax; // must use lmax of pw expansion
  double *wr       = pw_erf_data->wr;   // depends on atom type
  //  double *r        = pw_erf_data->r;    // depends on atom type

//================================================================================================
// Some checking

  if(nang!=ntheta*nphi || nr != nr_scr){
      printf("=============================================================\n");
      printf(" Mismatched integration points in PAW hartree self (%d!=%d) || (%d!=%d)\n",
                 nang,ntheta*nphi,nr,nr_scr);
      printf("=============================================================\n");
      exit(1);
  }//endif

//================================================================================================
// Nicely store rho_in : used (lmax+1)^2 times below so no cost to arrange nicely for use below.

  int iii = 0;
  for(int itp=0;itp<nang;itp++){
    for(int ir=0;ir<=nr;ir++,iii++){rho[ir][itp] = rho_in[iii];} // assumes ir is inner index of rho_in
  }//endfor

//================================================================================================
// Screened Hatree by partial waves
  double vself = 0.0;
  for(int il=0;il<=lmax;il++){
    double **ylm = ylm_data[il].ylm;
    int nm = 2*il+1;
    for(int m=0;m<nm;m++){
     // Compute lm components of rho inside atom by theta phi integration of rho*ylm = rho_lm(r) (real)
      for(int ir=0;ir<=nr;ir++){rho_lm[ir] = 0.0;}
      for(int ir=0;ir<=nr;ir++){
        for(int itp=0;itp<nang;itp++){rho_lm[ir] += ylm[m][itp]*rho[ir][itp]*wCostPhi[itp];}
      }//endfor
     // Compute Hartree self by inegration over r,r' and summing over lm : rho_lm(r)^2 pw_erf(r,r')
      double **pw_erf  = pw_erf_data->pw_erf_L[il].pw_erf;
      double  pre_l    = pw_erf_data->pw_erf_L[il].pre_l; // (4*pi)/(2*l+1)
      double pre_l_rt  = sqrt(pre_l);
      for(int ir=0;ir<=nr;ir++){rho_lm[ir]*=(pre_l_rt*wr[ir]);} // absorb weights and constant prefactor
      double vself_now = 0.0;
      for(int ir=0;ir<=nr;ir++){ 
      for(int irp=0;irp<=nr;irp++){
        vself_now += (rho_lm[ir]*pw_erf[ir][irp]*rho_lm[irp]);
      }}//endfor: ir, irp
      vself += (vself_now*0.5);
      //      if(vself_now>1.e-10){printf("atm(%d) Ylm(%d %d) %g\n",iatm,il,m,vself_now*0.5);}
    }//endfor : m
  }//enfor: il
  vself_out[0] = vself; // return vself

//-----------------------------------------------------------------------------------------------
  }//end routine
//===============================================================================================


//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void  init_atom_info(ATOM_INFO *atom_info, double beta_unitless){
//===============================================================================================
// Set up 1 atm and 1 atm typ by hand
  int natm_typ      = 1;
  int natm          = 1;
  int *iatm_atm_typ = new int [natm];
  int *lmax_psi     = new int [natm_typ];
  double *rmax      = new double [natm_typ];
  double *alpha     = new double [natm_typ];
  double *beta      = new double [natm_typ];

  for(int iatm=0;iatm<natm;iatm++){iatm_atm_typ[iatm] = 0;}

  int lmax_psi_max = 0;
  for(int iatm_typ=0;iatm_typ<natm_typ;iatm_typ++){
    lmax_psi[iatm_typ] = 2;
    lmax_psi_max       = MAX(lmax_psi_max,lmax_psi[iatm_typ]);
    rmax[iatm_typ]     = 4.0;
    alpha[iatm_typ]    = 2.8/rmax[iatm_typ];
    beta[iatm_typ]     = alpha[iatm_typ]*beta_unitless;
  }//endfor

//===============================================================================================
// pack away the data
  atom_info->natm          = natm;
  atom_info->natm_typ      = natm_typ;
  atom_info->beta_unitless = beta_unitless;
  atom_info->iatm_atm_typ  = iatm_atm_typ;
  atom_info->lmax_psi      = lmax_psi;
  atom_info->lmax_psi_max  = lmax_psi_max;
  atom_info->rmax          = rmax;
  atom_info->alpha         = alpha;
  atom_info->beta          = beta;

//-----------------------------------------------------------------------------------------------
  }//end routine
//===============================================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// readtoendofline: Function to read to end of line in read_coord files     
//==========================================================================
void readtoendofline(FILE *fp){
  int eol,ch;
  eol = (int )'\n';
  ch = eol+1;
  while(ch!=eol){ch=fgetc(fp);}
  if(ch==EOF){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("error: unexpected end of file reached          \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }//endif
}// end routine 
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void pw_erf_analytic(double *pw_erf, double beta, double rlt, double rgt, int lmax, int *lmax_out){
//==========================================================================
   double rgt2    = rgt*rgt;           double rlt2     = rlt*rlt;
   double rgt4    = rgt2*rgt2;         double rlt4     = rlt2*rlt2;
   double brgt    = beta*rgt;          double brlt     = beta*rlt;
   double brgt2   = brgt*brgt;         double brlt2    = brlt*brlt;
   double brgt4   = brgt2*brgt2;       double brlt4    = brlt2*brlt2;
   double brgt6   = brgt2*brgt4;       double brlt6    = brlt2*brlt4;
   double brgt8   = brgt4*brgt4;       double brlt8    = brlt4*brlt4;

   double argP    = brgt+brlt;         double argM     = brgt-brlt;
   double argH    = 2.0*brgt*brlt;     double argE     = brgt2+brlt2;

   double erfP    = erf(argP);         double erfM     = erf(argM);
   double derfP   = (erfP+erfM)*0.5;   double derfM    = (erfP-erfM)*0.5;

   double sinchH  = sinh(argH)/argH;   double coshH   = cosh(argH);
   double norm    = 2.0*beta/RT_PI_G;  double gaussH  = norm*exp(-argE);

   double poly    = 1.0 + 2.0*brlt2*(1.0+brgt2)/3.0;
   double derfM_rlt_S = gaussH;        double derfM_rlt_F = gaussH*poly;

   double combo   = brgt4 + brlt4 + brgt2*brlt2;
   double eee  = norm*exp((-argE+argH));
//==========================================================================
// Treat derfM at small beta*r< with 6th order Taylor expansion (good to O(brlt^8))
   double dderfM_rlt_S;
   double dderfM_rlt_F;
   double brlt_cut = 0.085; 
   if(brlt>brlt_cut){
     dderfM_rlt_S = (derfM-rlt*derfM_rlt_S)/rlt;
     dderfM_rlt_F = (derfM-rlt*derfM_rlt_F)/rlt;
   }else{
     // remove gaussian diff for precision improvement 
     double poly_m1      = 2.0*brlt2*(1.0+brgt2)/3.0;
     double derfM_rlt_S_G= 0.0;        double derfM_rlt_F_G = gaussH*poly_m1;
     double eeegt        = norm*exp(-brgt2);
     double diffG        = eeegt*(1.0-exp(-brlt2));
     double derfM_rlt_sm = eeegt*(+ ( 2.0*brgt2 -   1.0)*(brlt2/3.0) 
                                  + ( 4.0*brgt4 -  12.0*brgt2 + 3.0)*(brlt4/30.0)
                                  + ( 8.0*brgt6 -  60.0*brgt4 +  90.0*brgt2 - 15.0)*(brlt6/630.0)  
                                  + (16.0*brgt8 - 224.0*brgt6 + 840.0*brgt4 - 840.0*brgt2 + 105.0)*(brlt8/22680.0)
                                  );
     dderfM_rlt_S = derfM_rlt_sm - derfM_rlt_S_G + diffG;
     dderfM_rlt_F = derfM_rlt_sm - derfM_rlt_F_G + diffG;
   }//endif
//==========================================================================
// The first 4 partial waves from analytic work
  //-----------------------------------------------------------------------
  //Intermediate argH :
   if(argH<20.0 && argH>0.2){ 
     pw_erf[0] = (derfP)/rgt + (dderfM_rlt_S) - (gaussH)*(sinchH-1.0);

     if(lmax>=1){
       pw_erf[1] = (derfP)*(rlt/rgt2) + (dderfM_rlt_S)*(rgt/rlt)
                 - (gaussH)*( (2.0*argE-1.0)*sinchH + coshH - 2.0*brgt2)/argH;
     }//endif

     if(lmax>=2){
       pw_erf[2] = (derfP)*(rlt2/(rgt2*rgt)) + (dderfM_rlt_F)*(rgt2/(rlt2))
                 - (gaussH)*((4.0*combo -2.0*argE + 3.0)*sinchH 
                                     +(2.0*argE - 3.0)*coshH - 4.0*brgt4*poly)/(argH*argH);
     }//endif

     if(lmax>=3){
       pw_erf[3] = (derfP)*(rlt2*rlt/rgt4)  + (dderfM_rlt_F)*(rgt2*rgt/(rlt2*rlt))
                 - (gaussH)*( (-15.0 + 6.0*argE - 4.0*(brgt4 + brlt4 + 6.0*brgt2*brlt2)
                                     + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6) )*sinchH 
                             +(+15.0 - 6.0*argE + 4.0*(brgt4 + brlt4 + 1.0*brlt2*brgt2))*coshH
                             -8.0*poly*brgt6
                        )/(argH*argH*argH);
     }//endif
   }//endif : intermediate argH
  //-----------------------------------------------------------------------
  //Large argH :
   if(argH>=20.0){
     double eeeH   = exp(-argH);
     gaussH = norm*exp(-argE+argH);
     sinchH = 0.5*(1.0-eeeH*eeeH)/argH;
     coshH  = 0.5*(1.0+eeeH*eeeH);
     pw_erf[0] = (derfP)/rgt + (dderfM_rlt_S) - (gaussH)*(sinchH-1.0*eeeH);

     if(lmax>=1){
       pw_erf[1] = (derfP)*(rlt/rgt2) + (dderfM_rlt_S)*(rgt/rlt)
                 - (gaussH)*( (2.0*argE-1.0)*sinchH + coshH - 2.0*brgt2*eeeH)/argH;
     }//endif

     if(lmax>=2){
       pw_erf[2] = (derfP)*(rlt2/(rgt2*rgt)) + (dderfM_rlt_F)*(rgt2/(rlt2))
                 - (gaussH)*((4.0*combo -2.0*argE + 3.0)*sinchH 
                                       +(2.0*argE - 3.0)*coshH - 4.0*brgt4*poly*eeeH)/(argH*argH);
     }//endif

     if(lmax>=3){
       pw_erf[3] = (derfP)*(rlt2*rlt/rgt4)  + (dderfM_rlt_F)*(rgt2*rgt/(rlt2*rlt))
                 - (gaussH)*( (-15.0 + 6.0*argE - 4.0*(brgt4 + brlt4 + 6.0*brgt2*brlt2)
                                     + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6) )*sinchH 
                             +(+15.0 - 6.0*argE + 4.0*(brgt4 + brlt4 + 1.0*brlt2*brgt2))*coshH
                             -8.0*poly*brgt6*eeeH
                        )/(argH*argH*argH);
    }//endif
  }//endif : large argH
 //-----------------------------------------------------------------------
 // Small argH :
   if(argH<=0.2){
     double argH2   = argH*argH;
     double rats    = argH2/6.0;
     double ratc    = argH2/2.0;
     double sinchH1 = 0.0;
     double coshH1  = 0.0;
     for(int i = 3; i<=13;i+=2){
       sinchH1 += rats;
       coshH1  += ratc;
       double di = (double) i;
       rats *= (argH2/((di+2.0)*(di+1.0)));
       ratc *= (argH2/(di*(di+1.0)));
     }//endfor

     pw_erf[0] = (derfP)/rgt + (dderfM_rlt_S) - (gaussH)*(sinchH1);

     if(lmax>=1){
       pw_erf[1] = (derfP)*(rlt/rgt2) + (dderfM_rlt_S)*(rgt/rlt)
                 - (gaussH)*( (2.0*argE-1.0)*sinchH1 + coshH1 + 2.0*argE - 2.0*brgt2)/argH;
     }//endif

     if(lmax>=2){
       pw_erf[2] = (derfP)*(rlt2/(rgt2*rgt)) + (dderfM_rlt_F)*(rgt2/(rlt2))
                 - (gaussH)*((4.0*combo -2.0*argE + 3.0)*sinchH1  + 4.0*combo
                                       +(2.0*argE - 3.0)*coshH1 - 4.0*brgt4*poly)/(argH*argH);
     }//endif

     if(lmax>=3){
       pw_erf[3] = (derfP)*(rlt2*rlt/rgt4)  + (dderfM_rlt_F)*(rgt2*rgt/(rlt2*rlt))
                 - (gaussH)*( (-15.0 + 6.0*argE - 4.0*(brgt4 + brlt4 + 6.0*brgt2*brlt2)
                                     + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6) )*sinchH1
                                     - 20.0*brgt2*brlt2 + 8.0*(brgt6 + brlt4*brgt2 + brlt2*brgt4 + brlt6)
                             +(+15.0 - 6.0*argE + 4.0*(brgt4 + brlt4 + 1.0*brlt2*brgt2))*coshH1
                             -8.0*poly*brgt6
                            )/(argH*argH*argH);
     }//endif
   }//endif : small argH
//==========================================================================
// nan warning

  for(int il=0;il<=3;il++){if(isnan(pw_erf[il])){printf("analtyical_pw_%d(%g %g) nan\n",il,rlt,rgt);}}

//==========================================================================
// lmax = 3 analytic solutions

  lmax_out[0] = 3;
 
//==========================================================================
 }// end routine
//==========================================================================
