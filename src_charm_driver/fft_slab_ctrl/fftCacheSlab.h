//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file fftCacheSlab.h
 * Add type declarations for simulationConstants class (readonly vars) 
 * and once class for each type of object we currently have
 *(GStateXYslab, RealStateXZslab, RealRhoSlab, GRhoSlab, etc.) 
 */
//==============================================================================

#ifndef _fftcacheslab_h_
#define _fftcacheslab_h_

#include "../../include/RunDescriptor.h"

#ifdef _CP_USE_BLAS_GATHER_
include "../../src_mathlib/mathlib.h"
#endif

#define LEN_NHC_CP 4;
#define NUM_NHC_CP 20;


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class GStateSlab {
//==============================================================================
// data 
public:

   void pup(PUP::er &);
   int cp_min_opt;
   int numNonZeroPlanes;
   int numRuns, numLines;      // numLines=numRun/2=number of lines in collection
   int numPoints;              // number of non-zero pts in the collection
   int numFull;                // expanded data : nfftz*numLines
   int numFullNL;              // expanded data : nfftz*numLines
   int ees_nonlocal;
   size2d planeSize;           // size of the ffT in Y*Z dimension
   bool fftReqd;               // flags to indicate whether this data set of pencils
   int mysizeX;                // size of the ffT in X dimension 
   int S_grainSize;            // PC decomposition
   int xdim, ydim, zdim;       // FFT sizes
   int ngridaNL,ngridbNL,ngridcNL;       // FFT sizes for non-local
   int iplane_ind;             // collection index
   int istate_ind;             // state index of this collection
   int ihave_kx0;              // plane zero is in the collection (stored consecutively)
   int ihave_g000;
   int ind_g000;
   int kx0_strt;               // starting pt
   int kx0_end;                // ending pt
   int nkx0,nkx0_uni,nkx0_red; // split kx=0 into unique and redundant parts
   int nkx0_zero;              // ncoef_true=numPoints-nkx0_red

   double eke_ret;            // kinetic energy
   double fictEke_ret;        // fictitious kinetic energy
   double ekeNhc_ret;         // NHC energies
   double potNHC_ret;
   double degfree;            // Degrees of freedom (ncoef_true+num_nhc-1)
   double degfreeNHC;         // Degrees of freedom (num_nhc-1)*len_nhc
   double gammaNHC;           // Degrees of freedom degfree/(degfree+1.0)

   complex *packedPlaneData;   // Non-zero data pts [numPoints]
   complex *packedPlaneDataTemp; 
   complex *packedPlaneDataScr; 
   complex *packedForceData; 
   complex *packedVelData; 
   complex *packedRedPsi;
   int     len_nhc_cp;
   int     num_nhc_cp;
   double  kTCP;
   double  tauNHCCP;
   double  **xNHC;
   double  **xNHCP;
   double  *mNHC;
   double **vNHC;
   double **vNHC_scr;
   double **xNHC_scr;
   double **xNHCP_scr;
   double **fNHC;
   double *v0NHC,*a2NHC,*a4NHC;

//==============================================================================
// Constuctor, Destructor and utilities

   GStateSlab() {packedPlaneData=NULL; packedPlaneDataTemp=NULL; 
                 packedForceData=NULL; packedPlaneDataScr=NULL; 
                 packedVelData=NULL;}
   ~GStateSlab();

   void copyNHC(){
     for(int i=0;i<num_nhc_cp;i++){
     for(int j=0;j<len_nhc_cp;j++){
       vNHC_scr[i][j]  = vNHC[i][j];
       xNHC_scr[i][j]  = xNHC[i][j];
       xNHCP_scr[i][j] = xNHCP[i][j];
     }}//endfor
   }//end routine

   void initNHC(){
     xNHC      = new double *[20];
     xNHCP     = new double *[20];
     vNHC      = new double *[20];
     xNHC_scr  = new double *[20];
     xNHCP_scr = new double *[20];
     vNHC_scr  = new double *[20];
     fNHC      = new double *[20];
     for(int i=0;i<20;i++){
       xNHC[i]     = new double[4];
       xNHC_scr[i] = new double[4];
       xNHCP[i]     = new double[4];
       xNHCP_scr[i] = new double[4];
       vNHC[i]     = new double[4];
       vNHC_scr[i] = new double[4];
       fNHC[i]     = new double[4];
     }//endfor
     mNHC  = new double[4];
     v0NHC = new double[20];
     a2NHC = new double[20];
     a4NHC = new double[20];
     for(int i=0;i<20;i++){
     for(int j=0;j<4;j++){
       xNHC[i][j]     = 0.0;
       xNHC_scr[i][j] = 0.0;
       xNHCP[i][j]     = 0.0;
       xNHCP_scr[i][j] = 0.0;
       vNHC[i][j]     = 0.0;
       vNHC_scr[i][j] = 0.0;
       fNHC[i][j]     = 0.0;
     }}//endfor
   }//end routine

   void destroyNHC(){
     for(int i=0;i<num_nhc_cp;i++){
       delete []xNHC[i];
       delete []xNHC_scr[i];
       delete []xNHCP[i];
       delete []xNHCP_scr[i];
       delete []vNHC[i];
       delete []vNHC_scr[i];
       delete []fNHC[i];
     }//endfor
     delete []xNHC;
     delete []xNHC_scr;
     delete []xNHCP;
     delete []xNHCP_scr;
     delete []vNHC;
     delete []vNHC_scr;
     delete []fNHC;
     delete []mNHC;
     delete []v0NHC;
     delete []a2NHC;
     delete []a4NHC;
   }//end routine

   void setKRange(int , int *, int *, int *);
   void addForces(complex *,const int *);

//-----------------------------------------------------------------------------
  };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RealStateSlab {
public:
	size2d planeSize;         // size of the state in Z*X 
	complex *planeArr;
        double *planeArrR;   //planeArr cast to double
	int thisState;    
	int thisPlane;      
	int numPlanesToExpect;
        int nsize;
        int rsize;
        int size;
        double e_gga; 
	RealStateSlab() {}
	~RealStateSlab(); 
	void zeroOutPlanes();
	void allocate();
	void destroy();
	void pup(PUP::er &);

//-----------------------------------------------------------------------------
  };
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RhoRealSlab {
public:
   int size;    //plane size
   int trueSize;
   int sizeX, sizeY, sizeZ; //fft size
   int rhoRsubplanes;
   int csizeInt;
   int rsizeInt;
   double exc_ret, muxc_ret, exc_gga_ret;  //energy

   double *Vks;      // we have to keep him around
   double *density;  // we have to keep him around     
   double *rhoIRX,*rhoIRY,*rhoIRZ; //needed to receive stuff as it comes in
   double *VksHart;   //needed to receive stuff as it comes in
   double *rhoIRXint; 
   double *rhoIRYint; 
   double *rhoIRZint; 
   double *VksHartint; 

   // complex pointers to the same memory as the corresponding double array
   complex *VksC;     
   complex *densityC; 
   complex *rhoIRXC,*rhoIRYC,*rhoIRZC; 
   complex *VksHartC; 
   complex *rhoIRXCint; 
   complex *rhoIRYCint; 
   complex *rhoIRZCint; 
   complex *VksHartCint; 
 
   RhoRealSlab() {
       sizeX=sizeY=sizeZ=size=trueSize=0;
       exc_ret=muxc_ret=exc_gga_ret=0.0;
       Vks       = NULL;
       VksC      = NULL;
       density   = NULL;
       densityC  = NULL;
       rhoIRX    = NULL;
       rhoIRXC   = NULL;
       rhoIRY    = NULL;
       rhoIRYC   = NULL;
       rhoIRZ    = NULL;
       rhoIRZC   = NULL;
       Vks       = NULL;
       VksHartC  = NULL;
   }
   ~RhoRealSlab();

    void uPackScaleGrow(double *,double *,double );   // dest bigger  src : cp+scale
    void uPackScaleShrink(double *,double *,double ); // dest smaller src : cp+scale
    void uPackShrink(double *,double *);              // dest smaller src : cp only
    void uPackScale(double *, double *,double );      // dest size=   src : scale
    void scale(double *,double );                     // dest = src       : scale

    void pup(PUP::er &);

//-----------------------------------------------------------------------------
 };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class RhoGSlab {
 public:
	/* Stuff that is currently in CP_Rho_GSpacePlane's private variables */
	int sizeX, sizeY, sizeZ;
	int runsToBeSent;
	int numRuns;
	int numLines;
	int numFull;
	int numPoints;
	int nPacked;
	RunDescriptor *runs;

	complex *Rho;       // Is big enough to be expanded!
	complex *divRhoX;   // Is big enough to be expanded!
	complex *divRhoY;   // Is big enough to be expanded!
	complex *divRhoZ;   // Is big enough to be expanded!
	complex *packedRho;
        complex *packedVks;
        complex *Vks;       // Is big enough to be expanded!
	                    
 
        int *k_x, *k_y, *k_z;

	/* return values from rhoGSubroutine in subroutine.C */
	double ehart_ret, eext_ret, ewd_ret;
	
	/* used in the subroutine doRhoGSubroutine */
	int size;
	int xdim, ydim, zdim; 
	
	RhoGSlab() { //initialization paranoia
	    sizeX= sizeY= sizeZ=runsToBeSent= numRuns= numLines=numFull=0;
	    numPoints= nPacked=size=xdim=ydim= zdim=0; 
	    k_x=NULL;
	    k_y=NULL;
	    k_z=NULL;
	    Rho=NULL;       
	    divRhoX=NULL;   
	    divRhoY=NULL;   
	    divRhoZ=NULL;   
	    packedRho=NULL;
	    packedVks=NULL;
	    Vks=NULL;       
	    runs=NULL;
	}
	~RhoGSlab();
	void setKVectors(int *n);
        void divRhoGdot(double *,double ,complex *);
        void createWhiteByrd(double *, double );
	void pup(PUP::er &p);

//-----------------------------------------------------------------------------
  };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
//    Holder classes for the plans : Allows many fft libaries to be used
//==============================================================================

typedef struct fftplanholder {
  int option;             // 0= fftw, 1=essl
  int nfft;               // eesl stuff
  int ostride;            // eesl stuff
  int odist;              // eesl stuff
  int isign;              // eesl stuff
  int nwork1,    nwork2;  // eesl stuff
  double scale;           // eesl stuff
  double *work1, *work2;  // eesl stuff
  fftw_plan fftwPlan;     // fftw stuff
} FFTplanHolder;

typedef struct rfftplanholder {
  int option;             // 0= fftw, 1=essl
  int nfft;               // eesl stuff
  int ostride;            // eesl stuff
  int odist;              // eesl stuff
  int isign;              // eesl stuff
  int nwork1,    nwork2;  // eesl stuff
  double scale;           // eesl stuff
  double *work1, *work2;  // eesl stuff
  rfftwnd_plan rfftwPlan; // fftw stuff
} RFFTplanHolder;

//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class FFTcache: public Group {
 public:
    //-----------------------------------------------------------
    // FFT-Sizes and Ees method Options
     size2d planeSize;  
     int ngridaEext;    
     int ngridbEext;
     int ngridcEext;
     int ngridaNL;
     int ngridbNL;
     int ngridcNL;
     int ees_eext_on;
     int ees_NL_on;

    //-----------------------------------------------------------
    // Generic plane temporaries used everywhere possible to avoid memcpys
     complex *tmpData; 
     double  *tmpDataR;

    //-----------------------------------------------------------
    // Da Plans : Use the holder class to allow other fft libs
    //            All plans are double pack plans

     RFFTplanHolder fwdXPlanRho, bwdXPlanRho;     // state and density 
     FFTplanHolder  fwdYPlanRho, bwdYPlanRho;     
     FFTplanHolder  fwdYPlanRhoS, bwdYPlanRhoS;   // special subplane plan
     FFTplanHolder  fwdZPlanRho, bwdZPlanRho;

     RFFTplanHolder fwdXPlanNL, bwdXPlanNL;       // ees NL 
     FFTplanHolder  fwdYPlanNL, bwdYPlanNL;  
     FFTplanHolder  fwdZPlanNL, bwdZPlanNL;    

     RFFTplanHolder fwdXPlanEext, bwdXPlanEext;   // ees Eext
     FFTplanHolder  fwdYPlanEext, bwdYPlanEext;
     FFTplanHolder  fwdYPlanEextS, bwdYPlanEextS; // special subplane plan
     FFTplanHolder  fwdZPlanEext, bwdZPlanEext; 

    //-----------------------------------------------------------
    // The constructor 
     FFTcache(size2d planeSIZE, int , int , int , int , int , int , int , int );

    //-----------------------------------------------------------
    // Generic G-space expanders and contractors
     void expandGSpace(complex* data, complex *packedData, 
                       RunDescriptor *runs, int numRuns, int numFull,
 		       int numPoints, int nfftz);
     void packGSpace(complex* data, complex *packedPlaneData, 
                     RunDescriptor *runs, int numRuns, int numFull,
	             int numPoints, int nfftz);

    //-----------------------------------------------------------
    // Density FFTs
     void doHartFFTGtoR_Gchare(complex *,complex *,int , int ,int , int , 
				RunDescriptor *, int );
     void doRhoFFTRtoG_Gchare(complex *,complex *,int ,int ,int ,int ,RunDescriptor *, 
                              int ,int );
     void doRhoFFTGtoR_Gchare(complex *,complex *,int ,int ,int ,int ,RunDescriptor *, 
                              int ,int );
     void doRhoFFTRtoG_Rchare(complex *,double *,int , int ,int );
     void doRhoFFTGtoR_Rchare(complex *,double *,int , int ,int );

     void doRhoFFTRxToGx_Rchare(complex *,double *,int , int ,int );
     void doRhoFFTRyToGy_Rchare(complex *,double *,int , int ,int );
     void doRhoFFTGxToRx_Rchare(complex *,double *,int , int ,int );
     void doRhoFFTGyToRy_Rchare(complex *,double *,int , int ,int );

    //-----------------------------------------------------------
    // State FFTs
     void doStpFFTRtoG_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int);
     void doStpFFTGtoR_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int);
     void doStpFFTGtoR_Rchare(complex *,double *,int , int ,int );
     void doStpFFTRtoG_Rchare(complex *,double *,int , int ,int );

   //-----------------------------------------------------------
   // non-local fft
     void doNlFFTRtoG_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int);
     void doNlFFTGtoR_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int);
     void doNlFFTRtoG_Rchare(complex *,double *,int ,int ,int );
     void doNlFFTGtoR_Rchare(complex *,double *,int ,int ,int );

   //-----------------------------------------------------------
   // eext fft
     void doEextFFTRtoG_Gchare(complex *,int, int ,int ,int, RunDescriptor *,int);
     void doEextFFTGtoR_Gchare(complex *,complex *,int, int ,int ,int, RunDescriptor *,int);

     void doEextFFTRtoG_Rchare(complex *,double *,int ,int ,int );
     void doEextFFTGtoR_Rchare(complex *,double *,int ,int ,int );

     void doEextFFTRxToGx_Rchare(complex *,double *,int ,int ,int );
     void doEextFFTRyToGy_Rchare(complex *,double *,int ,int ,int );
     void doEextFFTGxToRx_Rchare(complex *,double *,int ,int ,int );
     void doEextFFTGyToRy_Rchare(complex *,double *,int ,int ,int );

//-----------------------------------------------------------------------------
  };
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
class GSlabInfo {
 public:
    CkIndex2D index;
    int numPoints;
    int numPlanes; 
    GSlabInfo() {}
    GSlabInfo(const GSlabInfo& info) { 
	index.x = info.index.x; index.y = info.index.y;
	numPoints = info.numPoints;
	numPlanes = info.numPlanes;
    }
    GSlabInfo(CkIndex2D idx, int _numPoints, int _numPlanes){
	index.x = idx.x; index.y = idx.x;
	numPoints = _numPoints;
	numPlanes = _numPlanes;
    }
};
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// slab initialization helpers

void initGStateSlab(GStateSlab *gs, int sizeX, size2d size, int gSpaceUnits, 
                    int realSpaceUnits, int s_grain,int iplane_ind,int istate_ind);
void initRealStateSlab(RealStateSlab *rs, size2d planeSize, int gSpaceUnits, 
                       int realSpaceUnits, int stateIndex, int thisPlane);
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim,
                     int xdimA, int ydimA, int myIndexX,int myIndexY,
                     int rhoRsubplanes);

//==============================================================================
// Eric's really cool BG/L progress callers

void fft_split(FFTplanHolder *fftplanholder, int howmany, 
               fftw_complex *in,  int istride, int idist, 
	       fftw_complex *out, int ostride, int odist, int split);

void rfftwnd_complex_to_real_split(RFFTplanHolder *rfftplanholder, int howmany, 
               fftw_complex *in, int istride, int idist, 
               fftw_real *out,   int ostride, int odist, int split);

void rfftwnd_real_to_complex_split(RFFTplanHolder *rfftplanholder, int howmany, 
    	       fftw_real *in,     int istride, int idist, 
               fftw_complex *out, int ostride, int odist, int split);

void initFFTholder  ( FFTplanHolder *,int *,int *,int *,double *,int *,int *,int *,int *);
void initRCFFTholder(RFFTplanHolder *,int *,int *,int *,double *,int *,int *,int *,int *);
void initCRFFTholder(RFFTplanHolder *,int *,int *,int *,double *,int *,int *,int *,int *);

//==============================================================================
#endif
