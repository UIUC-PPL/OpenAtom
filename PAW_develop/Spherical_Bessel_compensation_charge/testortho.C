//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Standard include files

#include "standard_include.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Main Program : Controller to manage compensation charge energy
//==========================================================================
int main (int argc, char *argv[]){
	int k = atoi(argv[1]);
	int kp = atoi(argv[2]);
	int n = atoi(argv[3]);
	double Rpc = 0.45;
	double ak = ((double) k)*M_PI_QI/Rpc;
	double ak1 = M_PI_QI/Rpc;
	double akp = ((double) kp)*M_PI_QI/Rpc;
	double delta = Rpc/((double) n);
	double * xr = new double [n+1];
	double * wr = new double [n+1];
	double * Bessel2 = new double [n+1];
	double result = 0.0;
	
	for (int i=0; i <= n; i++) {
		double r = ((double) i)*delta;
		xr[i] = r;
		Bessel2[i] = sin(ak1*r)*sin(ak1*r);
		wr[i] = delta*Bessel2[i]*2.0/Rpc;
		double argk = r*ak;
		double argkp = r*akp;
		result += delta*sin(argk)*sin(argkp)*2.0/Rpc;
	} // end for
	wr[0] *= 0.5;
	wr[n] *= 0.5;

	PRINTF("result1: %.10g\n", result);
	
	double alp_tmp = 1.8/Rpc;
	double beta_unitless = atof(argv[4]);
	double beta_tmp = alp_tmp*beta_unitless;

    double result4 = 0.0;
	double result3 = 0.0;
    double beta2 = beta_tmp*beta_tmp;
    for (int ir=0; ir<=n; ir++) {
        for (int jr=0; jr<=n; jr++) {
            double rgt = MAX(xr[ir],xr[jr]);
            double rlt = MIN(xr[ir],xr[jr]);
            double rd = rgt - rlt;
            double rs = rgt + rlt;
			double complicated;
			double simple;
			if (ir != 0 && jr != 0) {
          	 	double part1 = (exp(-beta2*rs*rs) - exp(-beta2*rd*rd))/(2.0*beta_tmp*sqrt(M_PI_QI)*rgt*rlt);
            	double part2 = (rd*erfc(beta_tmp*rd) - rs*erfc(beta_tmp*rs))/(2.0*rgt*rlt);
            	double part3 = 1.0/rgt;
	            complicated = part1 + part2 + part3;
				simple = 1.0/rgt;
			} else {
				if (ir !=0 || jr != 0) {
					complicated = erf(beta_tmp*rgt)/rgt;
					simple = 1.0/rgt;
				} else {
					complicated = 2.0*beta_tmp/sqrt(M_PI_QI);
					simple = 0.0;
				} // end if
			} // end if
			result3 += wr[ir]*wr[jr]*simple;
            result4 += wr[ir]*wr[jr]*complicated;
        } // end for jr
    } // end for ir
    PRINTF("result3: %.14g\n", result3);
    PRINTF("result4: %.14g\n", result4);

	return 1;	
}//end routine
//==========================================================================
