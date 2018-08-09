//==========================================================================
//					PAW f-grid

#ifndef _FGRID_CONTAINER_
#define _FGRID_CONTAINER_

//==========================================================================
// YLM_DATA structure: 
typedef struct YLM_DATA{
  int lang;                      // int  : angular momentum component of this guy
  int ntheta;                    // int  : number of cos(theta) integration points
  int nphi;                      // int  : number of phi integration points
  int nang;                      // int  : angular integration points ntheta*nphi
  int nm;                        // 2*lang+1
  double *contig;                // contiguous memory - size[lang*nm]
  double **ylm;                  // int [nm][nang] : Y_lang,m(cost,phi)
}YLM_DATA;

//==========================================================================
// PSI_COMP structure: it is the normalized radial wavefunction for each angular momentum 
typedef struct PSI_COMP_DATA{
  int nr;
  int lang;           // the current angular momentum
  double *r_psi0;     // size[nr]
}PSI_COMP_DATA;

//==========================================================================
// PW_ERF structure: 
typedef struct PW_ERF_DATA{
  int nr;
  int lang;         // the current angular momentum
  double alp;       // inverse convergence length for this partial wave
  double *contig;   // contiguous memory - size[nr*nr]
  double **pw_erf;  // [nr] by [nr]
}PW_ERF_DATA;

//==========================================================================
// FGRID structure: spherical polar with DVR
typedef struct FGRID{
  int nf;                  // total number of f points
  int nr;                  // number of r grid points
  int my_lmax;             // the maximum number of angular momentum channel for this atom type
  int my_lmax_rho;         // the maximum number of angular momentum channel for the density around this atom type
  int ityp;                // current atom type
  int lmax_pw;
  double alp;              // Gaussian parameter
  double alpb;             // Ewald inverse convergence length
  double beta;             // screen parameter
  double Rpc;              // Rpc
  double *wf;              // weight including the Jacobian - size[nf]
  double *xf;              // x coordinate - size[nf]
  double *yf;              // y - size[nf]
  double *zf;              // z - size[nf]
  double *xr;              // this is the r grid only - size[nr]
  double *wr;              // this is the r grid only - size[nr]
  double *r2_rho_tot;      // the total compensation charge density for this type - size[nf]
  PW_ERF_DATA *pw_erfA;    // Ewald alpha - size[lmax_pw+1] may converge at my_lmax_rho, to be tested
  PW_ERF_DATA *pw_erfB;    // Coulomb cutoff beta - size[lmax_pw+1]
  PW_ERF_DATA *pw_erfC;    // pure Coulomb - size[lmax_pw+1]
  PSI_COMP_DATA *psi_comp; // compensation charge wavefunction - size[my_lmax+1]
}FGRID;

class FGRID_CONTAINER{
  //----------------
public:
	int nf;              // total number of f points
	int nr;              // number of r grid points
	int ntheta;          // number of theta grid points
	int nphi;            // number of phi grid points
	int nang;            // number of angular grid points = ntheta*nphi
  int natm_typ;        // number of atom types
  int lmax;            // the maximum number of angular momentum channel for any atom type
  int lmax_rho;        // the maximum number of angular momentum channel for the density - just 2*lmax
  int lmax_pw;         // the maximum number of partial waves in the screened Coulomb expansion
  int *lmax_by_typ;    // the lmax of each atom type - size[natm_typ]
  double alpb;         // Ewald inverse convergence length
  double *alp_by_typ;  // alpha by atom type
  double *beta_by_typ; // beta by atom type
  double *Rpc_by_typ;  // Rpc by atom type
	double *wang;        // weight for the angular integral - size[nang]
	double *xcostheta;   // xcostheta - used to generate Ylmf - size[ntheta]
	double *xphi;        // xphi - used to generate Ylmf - size[nphi]
	double *wr_bare;     // - size[nr]
  double *xr_bare;     // everybody has the same r, but they are scaled by Rpc - size[nr]
  FGRID *fgrid;        // all the atom_type specific stuff - size[natm_typ]
  YLM_DATA *ylm_data;  // the spherical harmonics data with - size [lmax_pw+1] we may only need lmax_rho

  //----------------
  //con-destruct:
  FGRID_CONTAINER(){
    nf          = 0;
    nr          = 0;
    ntheta      = 0;
    nphi        = 0;
    nang        = 0;
    natm_typ    = 0;
    lmax        = 0;
    lmax_rho    = 0;
    lmax_by_typ = NULL;
    alpb        = 0.0;
    alp_by_typ  = NULL;
    beta_by_typ = NULL;
    Rpc_by_typ  = NULL;
    wang        = NULL;
    xcostheta   = NULL;
    xphi        = NULL;
    wr_bare     = NULL;
    xr_bare     = NULL;
    ylm_data    = NULL;
    fgrid       = NULL;
  };
  ~FGRID_CONTAINER(){};

  void allocate(){
    // allocating general arrays
    wang        = new double [nang];
    xcostheta   = new double [ntheta];
    xphi        = new double [nphi];
    wr_bare     = new double [nr];
    xr_bare     = new double [nr];
    fgrid       = new FGRID [natm_typ];
    ylm_data    = new YLM_DATA [lmax_pw+1];
    
    // allocating inside fgrid
    for(int i=0;i<natm_typ;i++) {
      int my_lmax = lmax_by_typ[i];
      fgrid[i].lmax_pw = lmax_pw;
      fgrid[i].nf = nf;
      fgrid[i].nr = nr;
      fgrid[i].my_lmax = my_lmax;
      fgrid[i].my_lmax_rho = 2*my_lmax;
      fgrid[i].ityp = i;
      fgrid[i].alp  = alp_by_typ[i];
      fgrid[i].beta = beta_by_typ[i];
      fgrid[i].Rpc  = Rpc_by_typ[i];
      fgrid[i].alpb = alpb;
      fgrid[i].wf = new double [nf];
      fgrid[i].xf = new double [nf];
      fgrid[i].yf = new double [nf];
      fgrid[i].zf = new double [nf];
      fgrid[i].xr = new double [nr];
      fgrid[i].wr = new double [nr];
      fgrid[i].r2_rho_tot= new double [nf];
      // allocating PW_ERF_DATA
      fgrid[i].pw_erfA = new PW_ERF_DATA [lmax_pw+1];
      fgrid[i].pw_erfB = new PW_ERF_DATA [lmax_pw+1];
      fgrid[i].pw_erfC = new PW_ERF_DATA [lmax_pw+1];
      for (int j=0;j<=lmax_pw;j++) {
        double * contigA = new double [nr*nr];
        double * contigB = new double [nr*nr];
        double * contigC = new double [nr*nr];
        double ** pwA = new double * [nr];
        double ** pwB = new double * [nr];
        double ** pwC = new double * [nr];
        for(int m=0;m<nr;m++){pwA[m] = &contigA[m*nr];}
        for(int m=0;m<nr;m++){pwB[m] = &contigB[m*nr];}
        for(int m=0;m<nr;m++){pwC[m] = &contigC[m*nr];}
        fgrid[i].pw_erfA[j].lang = j;
        fgrid[i].pw_erfB[j].lang = j;
        fgrid[i].pw_erfC[j].lang = j;
        fgrid[i].pw_erfA[j].nr = nr;
        fgrid[i].pw_erfB[j].nr = nr;
        fgrid[i].pw_erfC[j].nr = nr;
        fgrid[i].pw_erfA[j].alp = alp_by_typ[i];
        fgrid[i].pw_erfB[j].alp = beta_by_typ[i];
        fgrid[i].pw_erfC[j].alp = 1000.0; // a big number for Coulomb, never gonna be used
        fgrid[i].pw_erfA[j].contig = contigA;
        fgrid[i].pw_erfB[j].contig = contigB;
        fgrid[i].pw_erfC[j].contig = contigC;
        fgrid[i].pw_erfA[j].pw_erf = pwA;
        fgrid[i].pw_erfB[j].pw_erf = pwB;
        fgrid[i].pw_erfC[j].pw_erf = pwC;
      } // end for j 
      // allocating PSI_COMP_DATA
      fgrid[i].psi_comp = new PSI_COMP_DATA [my_lmax+1];
      for (int j=0;j<=my_lmax;j++) {
        fgrid[i].psi_comp[j].lang   = j;
        fgrid[i].psi_comp[j].nr     = nr;
        fgrid[i].psi_comp[j].r_psi0 = new double [nr];
      } // end for j
    } // end for i

    // allocating inside ylm_data
    for(int i=0; i<=lmax_pw; i++) {
      int nm = 2*i+1;
      ylm_data[i].lang   = i;
      ylm_data[i].nm     = nm;
      ylm_data[i].ntheta = ntheta;
      ylm_data[i].nphi   = nphi;
      ylm_data[i].nang   = nang;
      double **ylm   = new double *[nm];
      double *contig = new double [nm*nang];
      for(int m=0;m<nm;m++) {ylm[m] = &contig[m*nang];}
      ylm_data[i].contig = contig;
      ylm_data[i].ylm    = ylm;
    } // end for i
  } // end allocate()

#ifdef PUP_ON
  //----------------
  //pupping
  void pup(PUP::er &p){
    //pupping ints
    p | nf;
    p | nr;
    p | ntheta;
    p | nphi;
    p | nang;
    p | natm_typ;
    p | lmax;
    p | lmax_rho;
    p | lmax_pw;
    //pupping dbls
    p | alpb;
    //pupping by_typ arrays
    if (p.isUnpacking()) {lmax_by_typ = new int    [natm_typ];} PUParray(p, lmax_by_typ, natm_typ);
    if (p.isUnpacking()) {alp_by_typ  = new double [natm_typ];} PUParray(p, alp_by_typ,  natm_typ);
    if (p.isUnpacking()) {beta_by_typ = new double [natm_typ];} PUParray(p, beta_by_typ, natm_typ);
    if (p.isUnpacking()) {Rpc_by_typ  = new double [natm_typ];} PUParray(p, Rpc_by_typ,  natm_typ);
    // allocating all the memory required
    if (p.isUnpacking()) {allocate();}
    //pupping dbl arrays
    PUParray(p,wang     , nang  );
    PUParray(p,xcostheta, ntheta);
    PUParray(p,xphi     , nphi  );
    PUParray(p,wr_bare  , nr    );
    PUParray(p,xr_bare  , nr    );
		//pupping the fgrid
    for (int i=0;i<natm_typ;i++) {
      // the memory of these pointers is already allocated, so we just pup
      double *wf = fgrid[i].wf; PUParray(p,wf,nf);
      double *xf = fgrid[i].xf; PUParray(p,xf,nf);
      double *yf = fgrid[i].yf; PUParray(p,yf,nf);
      double *zf = fgrid[i].zf; PUParray(p,zf,nf);
      double *xr = fgrid[i].xr; PUParray(p,xr,nr);
      double *wr = fgrid[i].wr; PUParray(p,wr,nr);
      double *r2_rho_tot = fgrid[i].r2_rho_tot; PUParray(p,r2_rho_tot,nf);
      // pupping the PW_ERFs
      int my_lmax_rho = fgrid[i].my_lmax_rho;
      for (int j=0;j<=lmax_pw;j++) {
        double * contigA = fgrid[i].pw_erfA[j].contig; PUParray(p, contigA, nr*nr); 
        double * contigB = fgrid[i].pw_erfB[j].contig; PUParray(p, contigB, nr*nr); 
        double * contigC = fgrid[i].pw_erfC[j].contig; PUParray(p, contigC, nr*nr); 
      } // end for j
      // pupping the PSI_COMP
      int my_lmax = fgrid[i].my_lmax;
      for (int j=0;j<=my_lmax;j++) {
        double * r_psi0 = fgrid[i].psi_comp[j].r_psi0; PUParray(p, r_psi0, nr); 
      } // end for j
    } // end for i
    //pupping the ylm_data
    for(int i=0; i<=lmax_pw; i++) {
      int nm = 2*i+1;
      double * contig = ylm_data[i].contig; PUParray(p, contig, nm*nang);
    } // end for i
#ifdef _PARALLEL_DEBUG_        
    state_class_out ();
#endif      
  } // end pup
#endif

  void state_class_out(){
    int i;
    char fileName [255];
    sprintf (fileName, "%d_fgrid_container.state.out", CkMyPe());
    FILE *fp;
    fp = fopen(fileName,"w");
      // ints
      fprintf(fp,"nf       %d\n",nf      );
      fprintf(fp,"nr       %d\n",nr      );
      fprintf(fp,"ntheta   %d\n",ntheta  );
      fprintf(fp,"nphi     %d\n",nphi    );
      fprintf(fp,"nang     %d\n",nang    );
      fprintf(fp,"natm_typ %d\n",natm_typ);
      fprintf(fp,"lmax     %d\n",lmax    );
      fprintf(fp,"lmax_rho %d\n",lmax_rho);
      fprintf(fp,"lmax_pw  %d\n",lmax_pw );
      // int arrays
      for(i=0;i<natm_typ;i++) {fprintf(fp,"lmax_ty_typ[%d],  %d\n",i,lmax_by_typ[i]);}
      // dbles
      fprintf(fp,"alpb %g\n",alpb);
      // dbl arrays
      for(i=0;i<natm_typ;i++) {fprintf(fp,"alp_ty_typ [%d], %lg\n",i,alp_by_typ [i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"beta_ty_typ[%d], %lg\n",i,beta_by_typ[i]);}
      for(i=0;i<natm_typ;i++) {fprintf(fp,"Rpc_ty_typ [%d], %lg\n",i,Rpc_by_typ [i]);}
      for(i=0;i<nang    ;i++) {fprintf(fp,"wang       [%d], %lg\n",i,wang       [i]);}
      for(i=0;i<ntheta  ;i++) {fprintf(fp,"xcostheta  [%d], %lg\n",i,xcostheta  [i]);}
      for(i=0;i<nphi    ;i++) {fprintf(fp,"xphi       [%d], %lg\n",i,xphi       [i]);}
      for(i=0;i<nr      ;i++) {fprintf(fp,"wr_bare    [%d], %lg\n",i,wr_bare    [i]);}
      for(i=0;i<nr      ;i++) {fprintf(fp,"xr_bare    [%d], %lg\n",i,xr_bare    [i]);}
      // fgrid
      for(int j=0;j<natm_typ;j++){
        fprintf(fp,"  fgrid[%d].nf:          %d\n",j,fgrid[j].nf);
        fprintf(fp,"  fgrid[%d].nr:          %d\n",j,fgrid[j].nr);
        fprintf(fp,"  fgrid[%d].my_lmax:     %d\n",j,fgrid[j].my_lmax);
        fprintf(fp,"  fgrid[%d].lmax_pw:     %d\n",j,fgrid[j].lmax_pw);
        fprintf(fp,"  fgrid[%d].my_lmax_rho: %d\n",j,fgrid[j].my_lmax_rho);
        fprintf(fp,"  fgrid[%d].ityp:        %d\n",j,fgrid[j].ityp);
        fprintf(fp,"  fgrid[%d].alp:         %g\n",j,fgrid[j].alp);
        fprintf(fp,"  fgrid[%d].alpb:        %g\n",j,fgrid[j].alpb);
        fprintf(fp,"  fgrid[%d].beta:        %g\n",j,fgrid[j].beta);
        fprintf(fp,"  fgrid[%d].Rpc:         %g\n",j,fgrid[j].Rpc);
        for(i=0;i<nf;i++) {fprintf(fp,"  fgrid[%d].wf[%d]: %lg\n",j,i,fgrid[j].wf[i]);}
        for(i=0;i<nf;i++) {fprintf(fp,"  fgrid[%d].xf[%d]: %lg\n",j,i,fgrid[j].xf[i]);}
        for(i=0;i<nf;i++) {fprintf(fp,"  fgrid[%d].yf[%d]: %lg\n",j,i,fgrid[j].yf[i]);}
        for(i=0;i<nf;i++) {fprintf(fp,"  fgrid[%d].zf[%d]: %lg\n",j,i,fgrid[j].zf[i]);}
        for(i=0;i<nr;i++) {fprintf(fp,"  fgrid[%d].xr[%d]: %lg\n",j,i,fgrid[j].xr[i]);}
        for(i=0;i<nr;i++) {fprintf(fp,"  fgrid[%d].wr[%d]: %lg\n",j,i,fgrid[j].wr[i]);}
        // pw_erfA
        int my_lmax_rho = fgrid[j].my_lmax_rho;
        for(i=0;i<=lmax_pw;i++) {
          fprintf(fp,"  fgrid[%d].pw_erfA[%d]:     %d\n", j,i,fgrid[j].pw_erfA[i].nr);
          fprintf(fp,"  fgrid[%d].pw_erfA[%d]:     %d\n", j,i,fgrid[j].pw_erfA[i].lang);
          fprintf(fp,"  fgrid[%d].pw_erfA[%d]:     %lg\n",j,i,fgrid[j].pw_erfA[i].alp);
          for(int l=0;l<nr*nr;l++) {fprintf(fp,"  fgrid[%d].pw_erfA[%d].contig[%d]:  %lg\n",j,i,l,fgrid[j].pw_erfA[i].contig[l]);}
        } // end for i
        // pw_erfB
        for(i=0;i<=lmax_pw;i++) {
          fprintf(fp,"  fgrid[%d].pw_erfB[%d]:     %d\n", j,i,fgrid[j].pw_erfB[i].nr);
          fprintf(fp,"  fgrid[%d].pw_erfB[%d]:     %d\n", j,i,fgrid[j].pw_erfB[i].lang);
          fprintf(fp,"  fgrid[%d].pw_erfB[%d]:     %lg\n",j,i,fgrid[j].pw_erfB[i].alp);
          for(int l=0;l<nr*nr;l++) {fprintf(fp,"  fgrid[%d].pw_erfB[%d].contig[%d]:  %lg\n",j,i,l,fgrid[j].pw_erfB[i].contig[l]);}
        } // end for i
        // pw_erfC
        for(i=0;i<=lmax_pw;i++) {
          fprintf(fp,"  fgrid[%d].pw_erfC[%d]:     %d\n", j,i,fgrid[j].pw_erfC[i].nr);
          fprintf(fp,"  fgrid[%d].pw_erfC[%d]:     %d\n", j,i,fgrid[j].pw_erfC[i].lang);
          fprintf(fp,"  fgrid[%d].pw_erfC[%d]:     %lg\n",j,i,fgrid[j].pw_erfC[i].alp);
          for(int l=0;l<nr*nr;l++) {fprintf(fp,"  fgrid[%d].pw_erfC[%d].contig[%d]:  %lg\n",j,i,l,fgrid[j].pw_erfC[i].contig[l]);}
        } // end for i
        // psi_comp
        for(i=0;i<=my_lmax;i++) {
          fprintf(fp,"  fgrid[%d].psi_comp[%d]:    %d\n", j,i,fgrid[j].psi_comp[i].nr);
          fprintf(fp,"  fgrid[%d].psi_comp[%d]:    %d\n", j,i,fgrid[j].psi_comp[i].lang);
          for(int l=0;l<nr;l++) {fprintf(fp,"  fgrid[%d].psi_comp[%d].r_psi0[%d]:  %lg\n",j,i,l,fgrid[j].psi_comp[i].r_psi0[l]);}
        } // end for i
      } // end for j (atm_typ)
      // ylm_data
      for(i=0; i<=lmax_pw; i++) {
        int nm = ylm_data[i].nm;
        fprintf(fp,"  ylm_data[%d].lang:     %d\n",i,ylm_data[i].lang);
        fprintf(fp,"  ylm_data[%d].ntheta:   %d\n",i,ylm_data[i].ntheta);
        fprintf(fp,"  ylm_data[%d].nphi:     %d\n",i,ylm_data[i].nphi);
        fprintf(fp,"  ylm_data[%d].nang:     %d\n",i,ylm_data[i].nang);
        fprintf(fp,"  ylm_data[%d].nm:       %d\n",i,ylm_data[i].nm);
        for(int j=0;j<nm*nang;i++) {fprintf(fp,"  ylm_data[%d].contig[%d]:     %lg\n",i,j,ylm_data[i].contig[j]);}
      } // end for i
    fclose(fp);
  }// end routine
}; // FGRID_CONTAINER;

#ifdef PUP_ON
PUPmarshall(FGRID_CONTAINER);
#endif

#endif
//==========================================================================
