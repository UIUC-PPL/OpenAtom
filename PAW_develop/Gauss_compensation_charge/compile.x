# for mac, its compile.x FORT UNDER_OFF
# for others, compile.x FORT UNDER
rm model_PAW.x grid.o quad_rule.o GaussianPAW.o gen_fgrid.o gen_Ylmf.o gen_Gauss_quad.o gen_Gauss_quad_driver.o
g++ -Wall -g -c -O2 -D$1 grid.C 
g++ -Wall -g -c -O2 -D$1 gen_Gauss_quad.C -lquadmath -lblas -llapack -lm
g++ -Wall -g -c -O2 -D$1 gen_Gauss_quad_driver.C -lquadmath
g++ -Wall -g -c -O2 -D$1 gen_Ylmf.C 
g++ -Wall -g -c -O2 -D$1 quad_rule.C 
g++ -Wall -g -c -O2 -D$1 GaussianPAW.C
g++ -std=c++11 -Wall -g -O2 -D$1 model_PAW.C -o model_PAW.x gen_Gauss_quad.o gen_Gauss_quad_driver.o GaussianPAW.o gen_Ylmf.o grid.o quad_rule.o -lblas -llapack -lquadmath -lm
rm -rf grid.o quad_rule.o GaussianPAW.o gen_fgrid.o gen_Ylmf.o gen_Gauss_quad.o gen_Gauss_quad_driver.o *.dSYM
