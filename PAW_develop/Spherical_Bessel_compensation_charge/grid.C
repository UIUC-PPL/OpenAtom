//=========================================================================t
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//  This is the spherical polar grid for PAW
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "standard_include.h"
#include "ckcomplex.h"
#include "fgrid.h"

//==========================================================================
// Function prototypes

#include "grid.h"

//==========================================================================

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
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("ERROR: Unexpected end of file reached          \n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif
}// end routine 
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// SphericalHarmonicOnGrid: Function to calculate the Spherical Harmonics 
//==========================================================================
void SphericalHarmonicOnGrid(int l, int m, int thetaorder, double *xcostheta, int phiorder, double *xphi, complex **Ylm){
	double xsintheta[thetaorder];
	for (int i=0; i < thetaorder; i++) xsintheta[i] = sqrt((1.0-xcostheta[i]*xcostheta[i]));
	if (abs(m) > l) {
		PRINTF("Invalid Spherical Harmonics! abs(m) > l !!\n");
		EXIT(1);
	}
	else if (l < 0) {
		PRINTF("Invalid Spherical Harmonics! l MUST be > 0 !!\n");
		EXIT(1);
	}
	switch (l) {
		case 0: {
			for (int i=0; i<thetaorder; i++) {
				for (int j=0; j<phiorder; j++) {
					Ylm[i][j] = complex (1.0/2.0/sqrt(M_PI_QI),0.0);
				}
			}
			break;
		} 
		case 1: {
			if (m == -1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/2.0*sqrt(3.0/(2.0*M_PI_QI))*xsintheta[i]*CkExpIm(-xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 0) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/2.0*sqrt(3.0/(M_PI_QI))*xcostheta[i];
					} // end for i
				} //end for j
			}// end if
			if (m == 1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = -1.0/2.0*sqrt(3.0/(2.0*M_PI_QI))*xsintheta[i]*CkExpIm(xphi[j]);
					} // end for i
				} //end for j
			}// end if
			break;
		} 
		case 2: {
			if (m == -2) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/4.0*sqrt(15.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*CkExpIm(-2.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == -1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/2.0*sqrt(15.0/(2.0*M_PI_QI))*xsintheta[i]*xcostheta[i]*CkExpIm(-1.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 0) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/4.0*sqrt(5.0/(1.0*M_PI_QI))*(3.0*xcostheta[i]*xcostheta[i]-1.0);
					} // end for i
				} //end for j
			}// end if
			if (m == 1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = -1.0/2.0*sqrt(15.0/(2.0*M_PI_QI))*xsintheta[i]*xcostheta[i]*CkExpIm(1.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 2) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/4.0*sqrt(15.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*CkExpIm(2.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			break;
		} 
		case 3: {
			if (m == -3) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/8.0*sqrt(35.0/(1.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xsintheta[i]*CkExpIm(-3.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == -2) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/4.0*sqrt(105.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xcostheta[i]*CkExpIm(-2.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == -1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/8.0*sqrt(21.0/(1.0*M_PI_QI))*xsintheta[i]*(5*xcostheta[i]*xcostheta[i]-1.0)*CkExpIm(-1.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 0) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/4.0*sqrt(7.0/(1.0*M_PI_QI))*xcostheta[i]*(5*xcostheta[i]*xcostheta[i]-3.0);
					} // end for i
				} //end for j
			}// end if
			if (m == 1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = -1.0/8.0*sqrt(21.0/(1.0*M_PI_QI))*xsintheta[i]*(5*xcostheta[i]*xcostheta[i]-1.0)*CkExpIm(1.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 2) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 1.0/4.0*sqrt(105.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xcostheta[i]*CkExpIm(2.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 3) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = -1.0/8.0*sqrt(35.0/(1.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xsintheta[i]*CkExpIm(3.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			break;
		} 
		case 4: {
			if (m == -4) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/16.0*sqrt(35.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xsintheta[i]*xsintheta[i]*CkExpIm(-4.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == -3) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/8.0*sqrt(35.0/(1.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xsintheta[i]*xcostheta[i]*CkExpIm(-3.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == -2) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/8.0*sqrt(5.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*(7*xcostheta[i]*xcostheta[i]-1.0)*CkExpIm(-2.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == -1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/8.0*sqrt(5.0/(1.0*M_PI_QI))*xsintheta[i]*(7*xcostheta[i]*xcostheta[i]*xcostheta[i]-1.0*xcostheta[i])*CkExpIm(-1.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 0) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/16.0*sqrt(1.0/(1.0*M_PI_QI))*(35*pow(xcostheta[i],4) - 30*pow(xcostheta[i],2) + 3);
					} // end for i
				} //end for j
			}// end if
			if (m == 1) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = -3.0/8.0*sqrt(5.0/(1.0*M_PI_QI))*xsintheta[i]*(7*xcostheta[i]*xcostheta[i]*xcostheta[i]-1.0*xcostheta[i])*CkExpIm(1.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 2) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/8.0*sqrt(5.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*(7*xcostheta[i]*xcostheta[i]-1.0)*CkExpIm(2.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 3) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = -3.0/8.0*sqrt(35.0/(1.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xsintheta[i]*xcostheta[i]*CkExpIm(3.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			if (m == 4) {
				for (int i=0; i<thetaorder; i++) {
					for (int j=0; j<phiorder; j++) {
						Ylm[i][j] = 3.0/16.0*sqrt(35.0/(2.0*M_PI_QI))*xsintheta[i]*xsintheta[i]*xsintheta[i]*xsintheta[i]*CkExpIm(4.0*xphi[j]);
					} // end for i
				} //end for j
			}// end if
			break;
		} 
		default: {	
			PRINTF("Currently we only support Spherical Harmonics l <= 4!\n");
			EXIT(1);
		}
	} //end switch
}// end routine 
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// calcHermite: Calculate the value of the Gauss-Hermite quadrature
//==========================================================================
void calcHermiteHalfSpace(int neworder, int norder, double *w, double *x, double *neww, double *newx) {
	if (norder == 1) {
		if (neworder != 1) {
    		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			PRINTF("ERROR in calcHermiteHalfSpace! neworder != 1\n");
		    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			EXIT(1);
		}
		neww[0] = w[0]/2.0;
		newx[0] = x[0];
	} // end if
	if (norder % 2 == 0) {
		if (neworder != norder/2) {
		    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			PRINTF("ERROR in calcHermiteHalfSpace! neworder != %d\n", norder/2);
    		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			EXIT(1);
		}
  		for (int i=0; i<neworder; i++) {
			neww[i] = w[norder/2+i];
			newx[i] = x[norder/2+i];
  		} //endfor
	}else {
		if (neworder != norder/2 + 1) {
    		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			PRINTF("ERROR in calcHermiteHalfSpace! neworder != %d\n", norder/2 + 1);
		    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
			EXIT(1);
		}
		neww[0] = w[neworder-1]/2.0;
		newx[0] = x[neworder-1];
  		for (int i=1; i<neworder; i++) {
			neww[i] = w[neworder+i-1];
			newx[i] = x[neworder+i-1];
  		} //endfor
	} //end else
} //end routine

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// genphigrid: generate the phi grid
//==========================================================================
void genphigrid(int phiorder, double *w, double *x) {
	if (phiorder%2 != 0) {
   		PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		PRINTF("ERROR in genphigrid! phiorder MUST be even!\n");
	    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		EXIT(1);
	}//end if
	double aorder = (double) phiorder;
	for (int i=0; i<phiorder; i++) {
		double ai = (double) i;
		w[i] = (2.0*M_PI_QI)/aorder;
		x[i] = (2.0*M_PI_QI)/aorder*ai;
	} //endfor
} //end routine

