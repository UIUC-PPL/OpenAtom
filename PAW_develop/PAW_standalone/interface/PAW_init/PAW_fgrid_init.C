//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//            interface/PAW_init/PAW_fgrid_init.C
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "standard_include.h"
#include "atom_maps.h"
#include "fgrid.h"
#include "PAW_fgrid_init.h"

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void PAW_FGRID_INIT::fill_fgrid(FGRID_CONTAINER * fgrid_container, ATOM_MAPS * atom_maps, int model)
  //===================================================================================
{//begin routine
  //===================================================================================
  // Local pointers and consistency checks
  natm_typ = atom_maps->natm_typ      ;
  lmax     = atom_maps->lmax          ;
  lmax_rho = atom_maps->lmax_rho      ;
  lmax_pw  = fgrid_container->lmax_pw ;
  nr       = fgrid_container->nr      ;
  ntheta   = fgrid_container->ntheta  ;
  nphi     = fgrid_container->nphi    ;
  nf       = fgrid_container->nf      ;
  nang     = fgrid_container->nang    ;
  alpb     = fgrid_container->alpb    ;
  fgrid_container->lmax     = lmax    ;
  fgrid_container->natm_typ = natm_typ;

  //----------------------------------------------------------------------------------
  // Need enough partial wave for square of wave function so 2x more than you think
  if(lmax_pw < lmax_rho){
    printf("=============================================================\n");
    printf(" lmax_pw of partial wave expansion must be >= 2*lmax_psi of states = lmax_rho\n");
    printf("=============================================================\n");
    exit(1);
  }//endif
  //----------------------------------------------------------------------------------
  // Cos(theta): need to integrate Plmax_pw(cost) * Plmax_pw(cost)*= poly degree 2*lmax_pw
  //    2*ntheta-1 = 2*lmax_pw  : ntheta = (2*lmax_pw+1)/2 : ntheta = lmax_pw + 1   smallest integer satisfies
  if(ntheta <= lmax_pw){
    printf("=============================================================\n");
    printf(" cos(theta) GL integration points must be > lmax_pw\n");
    printf("=============================================================\n");
    exit(1);
  }//endif

  //----------------------------------------------------------------------------------
  //reserving memory for the by type arrays
  lmax_by_typ = new int    [natm_typ];
  alp_by_typ  = new double [natm_typ];
  beta_by_typ = new double [natm_typ];
  Rpc_by_typ  = new double [natm_typ];

  for (int i=0;i<natm_typ;i++) {
    lmax_by_typ[i] = atom_maps->lmax_by_typ[i];
    alp_by_typ[i]  = atom_maps->alp_by_typ[i];
    beta_by_typ[i] = atom_maps->beta_by_typ[i];
    Rpc_by_typ[i]  = atom_maps->Rpc_by_typ[i];
  } // end for i

  fgrid_container->lmax_by_typ = lmax_by_typ;
  fgrid_container->alp_by_typ  = alp_by_typ;
  fgrid_container->beta_by_typ = beta_by_typ;
  fgrid_container->Rpc_by_typ  = Rpc_by_typ;

  //===================================================================================
  // allocate all the required memory and get some more local pointers
  fgrid_container->allocate();

  //===================================================================================
  // set up the generic grid 
  double *r_psi0    = new double [nr];
  double *wcostheta = new double [ntheta];
  double *wphi      = new double [nphi];
  init_fgrid_gen(fgrid_container, wcostheta, wphi, r_psi0, model);
  double *xr_bare   = fgrid_container->xr_bare;
  double *wr_bare   = fgrid_container->wr_bare;
  double *xphi      = fgrid_container->xphi; 
  double *xcostheta = fgrid_container->xcostheta; 
  double *wang      = fgrid_container->wang;

  //=====================================================================
  // use the quadrature to create the fgrid struct for each atom type
  for (int i=0; i<natm_typ; i++) {
    FGRID *fgrid = &(fgrid_container->fgrid[i]);
    PSI_COMP_DATA *psi_comp = fgrid_container->psi_comp;
    double *wf = fgrid->wf;
    double *xf = fgrid->xf;
    double *yf = fgrid->yf;
    double *yf = fgrid->yf;
    double *xr = fgrid->xr;
    double *wr = fgrid->wr;
 /*
 56   double *r2_rho_tot;      // the total compensation charge density for this type - size[nf]
 57   PW_ERF_DATA *pw_erfA;    // Ewald alpha - size[lmax_pw+1] may converge at my_lmax_rho, to be tested
 58   PW_ERF_DATA *pw_erfB;    // Coulomb cutoff beta - size[lmax_pw+1]
 59   PW_ERF_DATA *pw_erfC;    // pure Coulomb - size[lmax_pw+1]
 60   PSI_COMP_DATA *psi_comp; // compensation charge wavefunction - size[my_lmax+1]
*/
    double scale = 1.0;
    if (model == 1) {scale = fgrid[i].Rpc;}
    if (model == 2) {scale = fgrid[i].alp;}
    double *xr = fgrid[i].xr;
    double *wr = fgrid[i].wr;
    for (int ir=0; ir<nr; ir++) {
      xr[ir] = xr_bare[ir]*scale;
      wr[ir] = wr_bare[ir]*scale;
    } // end for ir

    // over simplified model
    int my_lmax = fgrid[i].my_lmax;
    for (int j=0; j<=my_lmax; j++) {
      for (int k=0; k<nr; k++) { fgrid[i].psi_comp[j].r_psi0[k] = r_psi0[k]/sqrt(scale); }
    } // end for j

    // total grid
    double *wf = fgrid_container->fgrid
    int iff = 0;
    for(int icost=0;icost<ntheta;icost++){
    for(int iphi=0;iphi<nphi;iphi++){
    for(int ir=0;ir<=nr;ir++,iff++){
      double sint = sqrt(1.0-xcostheta[icost]*xcostheta[icost]);
      xf[iff] = r[ir]*sint*cos(xphi[iphi]);
      yf[iff] = r[ir]*sint*sin(xphi[iphi]);
      zf[iff] = r[ir]*xcostheta[icost];
      wf[iff] = wphi[iphi]*wcostheta[icost]*wr[ir];
    }}}//endfor

    double gamma2_inv = 0.5/(alp_tmp*alp_tmp) + 0.25/(beta_tmp*beta_tmp);
    double gamma = 1.0/sqrt(gamma2_inv);
#ifdef _INTEGRATION_CHECK_
    PRINTF("2*Hartreeself_screen: %.10g %.10g beta_tmp = %g\n", gamma/sqrt(M_PI_QI), result4, beta_tmp);
#endif
    for (int itheta=0; itheta<ntheta; itheta++) { xcostheta[itheta] = xcostheta[itheta];}
    for (int iphi=0; iphi<nphi; iphi++) { xphi[iphi] = xphi[iphi];}

    int f = 0; int iii = 0;
    double result5 = 0.0, result6 = 0.0;
    for (int itheta=0; itheta<ntheta; itheta++) {
       for (int iphi=0; iphi<nphi; iphi++) {
          wang[iii] = wcostheta[itheta]*wphi[iphi];
          ylm[iii]  = 1.0/sqrt(4.0*M_PI_QI);
          for (int ir=0; ir<nr; ir++) {
                double xsintheta = sqrt(1.0-(xcostheta[itheta])*(xcostheta[itheta]));
                wf[f] = wr[ir]*wcostheta[itheta]*wphi[iphi];
                xf[f] = xr[ir]*xsintheta*cos(xphi[iphi]);
                yf[f] = xr[ir]*xsintheta*sin(xphi[iphi]);
                zf[f] = xr[ir]*xcostheta[itheta];
                rf[f] = xr[ir];
                rho[f] = r2_psi0_2[ir]/(4.0*M_PI_QI);
                result6 += rho[f]*wr_bare[ir]*wang[iii];
                result5 += wf[f];
                f++;
          } // end for iphi
          iii++;
        } // end for itheta
    } // end for ir
#ifdef _INTEGRATION_CHECK_
    PRINTF("result5: %.10g\n", result5/(4.0*M_PI_QI));
    PRINTF("result6: %.10g\n", result6);
#endif

  } // end for i: atom type

  for (int i=0; i<natm_typ; i++) {
    int lmax_tmp = 0; int nr = fgrid[i].nr; int nang = fgrid[i].nang;
    double *ylm = fgrid[i].ylm; double *wang = fgrid[i].wang;
    double *wr_bare = fgrid[i].wr_bare;
    double **pw_coul = fgrid[i].pw_coul;
    double **pw_erfA = fgrid[i].pw_erfA;
    double **pw_erfB = fgrid[i].pw_erfB;
    double *rho_in = fgrid[i].rho; double **rho_scr = fgrid[i].rho_scr;
    double *rho_lm = fgrid[i].rho_lm; double vself_coul = 0.0, vselfA = 0.0, vselfB = 0.0;
    screen_self_Hartree(lmax_tmp, nr, nang, ylm, wang, wr_bare, pw_coul, rho_in, rho_scr, rho_lm, &vself_coul);
    screen_self_Hartree(lmax_tmp, nr, nang, ylm, wang, wr_bare, pw_erfA, rho_in, rho_scr, rho_lm, &vselfA);
    screen_self_Hartree(lmax_tmp, nr, nang, ylm, wang, wr_bare, pw_erfB, rho_in, rho_scr, rho_lm, &vselfB);
#ifdef _INTEGRATION_CHECK_
    PRINTF("vself_coul, vself_screenA, vself_screenB: %.14g %.14g %.14g\n", 2.0*vself_coul, 2.0*vselfA, 2.0*vselfB);
#endif
  } // end for i: atom type
} // end routine
  
  //===================================================================================
}//end routine
//========================================================================


//===============================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===============================================================================================
void init_fgrid_gen(FGRID_CONTAINER *fgrid_container, double *wcostheta, double *wphi, double *r_psi0, int model){
//===============================================================================================
  
//===============================================================================================
// Local pointers
  double *xr_bare   = fgrid_container->xr_bare;
  double *wr_bare   = fgrid_container->wr_bare;
  double *xphi      = fgrid_container->xphi; 
  double *xcostheta = fgrid_container->xcostheta; 
  double *wang      = fgrid_container->wang;

//===================================================================================
// set up the radial quadrature, which depends on model 
  //-----------------------------------------------------------------------------------
  // Gaussian model with half space Hermite quadrature
  if (model == 2) {
    int ng;
    char fname[100];
    sprintf(fname,"hh_w_x_%d.dat",nr);
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
      }//endif
      readtoendofline(fp);
      if(ng!=nr){
        printf("=============================================================\n");
        printf(" File %s does not have %d points but %d\n",fname,ng,nr);
        printf("=============================================================\n");
        exit(1);
      }
      for(int ig=0;ig<ng;ig++){
        int jjj = fscanf(fp,"%lg %lg",&xr_bare[ig],&wr_bare[ig]);
        if(jjj!=2){
          printf("=============================================================\n");
          printf(" %d th line of file %s invalid\n",ig+2,fname);
          printf("=============================================================\n");
          exit(1);
        }//endif
        readtoendofline(fp);
      }//endfor
    fclose(fp);
    double norm  = 4.0/sqrt(M_PI_QI);
    double check = 0.0;
    for (int i=0; i < nr; i++) {
      double x2    = xr_bare[i]*xr_bare[i];
      r_psi0[i]    = xr_bare[i]*sqrt(norm);
      r2_psi0_2[i] = r_psi0[i]*r_psi0[i];  // the weight function of the quadrature has the Gaussian
      wr_bare[i]   = wr_bare[i];
      check       += wr_bare[i]*r2_psi0_2[i];
    } // end for
#ifdef _INTEGRATION_CHECK_
    PRINTF("check = %g\n", check);
#endif
  } // end if

  //-----------------------------------------------------------------------------------
  // spherical Bessel model with equally spaced quadrature
  if (model == 1) {
    double Rpc_bare   = 1.0;
    double delta_bare = Rpc_bare/((double) (nr - 1));
    double check      = 0.0;
    double norm       = 2.0/Rpc_bare;
    for (int i=0; i < nr; i++) {
      double r     = ((double) i)*delta_bare;
      xr_bare[i]   = r; // depends on Rpc - only guy that needs to change if Rpc != 1.0
      double ak1   = M_PI_QI/Rpc_bare;
      r_psi0[i]    = sin(ak1*r)*sqrt(norm); // ak1*r is unitless
      r2_psi0_2[i] = r_psi0[i]*r_psi0[i];
      wr_bare[i]   = delta_bare;
      check       += wr_bare[i]*r2_psi0_2[i];
    } // end for
    wr_bare[0] *= 0.5; wr_bare[nr-1] *= 0.5;
    wr_bare[0] *= 0.5; wr_bare[nr-1] *= 0.5;
#ifdef _INTEGRATION_CHECK_
    PRINTF("check = %g\n", check);
#endif
  } // end if

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
    }//endif
    readtoendofline(fp);
    if(ng!=ntheta){
      printf("=============================================================\n");
      printf(" File %s does not have %d points but %d\n",fname,ng,ntheta);
      printf("=============================================================\n");
      exit(1);
    }
    for(int ig=0;ig<ng;ig++){
      int jjj = fscanf(fp,"%lg %lg",&xcostheta[ig],&wcostheta[ig]);
      if(jjj!=2){
        printf("=============================================================\n");
        printf(" %d th line of file %s invalid\n",ig+2,fname);
        printf("=============================================================\n");
        exit(1);
      }//endif
      readtoendofline(fp);
    }//endfor
  fclose(fp);

//===============================================================================================
// phi grid  
  double dphi = 2.0*M_PI_G/((double)nphi);
  for(int iphi=0;iphi<nphi;iphi++){
    double diphi = (double)iphi;
    xphi[iphi]   = dphi*diphi;
    wphi[iphi]   = dphi;
  }//endfor

//===============================================================================================
// angular grid  
  int iang = 0;
  for(int icost=0;icost<ntheta;icost++){
  for(int iphi=0;iphi<nphi;iphi++,iang++){
    wang[iang] = wphi[iphi]*wcostheta[icost];
  }}//endfor

//------------------------------------------------------------------------------
  }// end routine
//==================================================================================

