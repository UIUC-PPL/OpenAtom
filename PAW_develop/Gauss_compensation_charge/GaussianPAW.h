//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  The function prototypes for GaussianPAW.C
//
//==========================================================================
// PINY constants, inline functions and typedefs
 
typedef struct Eterm{
//	string name;
	NAME   name;
	double E;
	double EGrid;
	void pres () {
	//	PRINTF("%-17s  |   % 09.7f  |  % 09.7f  |  %9.7f  |  %4.2f\n",name.c_str(),E,EGrid,fabs(E-EGrid),fabs(E-EGrid)/MAX(fabs(E),1e-20)*100);
		PRINTF("%-17s  |   % 09.7f  |  % 09.7f  |  %9.7f  |  %4.2f\n",name,E,EGrid,fabs(E-EGrid),fabs(E-EGrid)/MAX(fabs(E),1e-20)*100);
	}
}Eterm;

typedef struct ESTRUCT{
	Eterm ENN;
	Eterm EeNself;
	Eterm EeN;
	Eterm EHarself;
	Eterm EHar;
	Eterm ENNshort;
	Eterm EeNshortself;
	Eterm EeNshort;
	Eterm EHarshortself;
	Eterm EHarshort;
	Eterm Elong;
	Eterm ENNselflong;
	Eterm Etot0D;
	Eterm Etot3D;
	Eterm EHarselfscr;
}ESTRUCT;
 
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Function prototypes
void gaussran(int , long *, double *);
double altRandom(long *); 
void gethinv(double *, double *, double *, int );
void computePAWreal(ATOM_MAPS *, ATOM_POS *, CELL *, ESTRUCT *);
void computePAWGrid(int, ATOM_MAPS *, ATOM_POS *, CELL *, ESTRUCT *, FGRID *);
void computePAWlong(ATOM_MAPS *, ATOM_POS *, CELL *, ESTRUCT *, FGRID *); 
double dist(double, double, double);	
void gen_Ylmf (int, int, double *, int, double *, int, int, complex *);
double erfc_a_r_over_r(double, double);
//==========================================================================
