rm model_PAW.x grid.o quad_rule.o GaussianPAW.o gen_fgrid.o gen_Ylmf.o gen_Gauss_quad.o gen_Gauss_quad_driver.o
g++ -Wall -g -c -O2 grid.C 
g++ -Wall -g -c -O2 gen_Gauss_quad.C -lquadmath -lblas -llapack -lm
g++ -Wall -g -c -O2 gen_Gauss_quad_driver.C -lquadmath
g++ -Wall -g -c -O2 gen_Ylmf.C 
g++ -Wall -g -c -O2 quad_rule.C 
g++ -Wall -g -c -O2 GaussianPAW.C
g++ -std=c++11 -Wall -g -O2 -o model_PAW.x model_PAW.C gen_Gauss_quad.o gen_Gauss_quad_driver.o GaussianPAW.o gen_Ylmf.o grid.o quad_rule.o -lblas -llapack -lquadmath -lm
rm -rf grid.o quad_rule.o GaussianPAW.o gen_fgrid.o gen_Ylmf.o gen_Gauss_quad.o gen_Gauss_quad_driver.o *.dSYM
