#=============================================
#         Instructions
# Mac or others sans fort_under: compile.x FORT UNDER_OFF
# Compilers with fort under: compile.x FORT UNDER
#=============================================
#         Clean up
#rm -rf model_PAW.x *.o *.a
rm -rf model_PAW.x *.o
#=============================================
#         Quad precision ggev
#cd ../ggev_quad_pkg/compile/
#./compile.x
#rm -rf *.o
#cd -
#/bin/cp ./../ggev_quad_pkg/compile/libqggev.a .
#/bin/rm ./../ggev_quad_pkg/compile/libqggev.a 

#=============================================
g++-5 -Wall -g -c -O2 grid.C 
g++-5 -Wall -g -c -O2 gen_Gauss_quad.C 
g++-5 -Wall -g -c -O2 gen_Gauss_quad_driver.C
g++-5 -Wall -g -c -O2 gen_Ylmf.C 
g++-5 -Wall -g -c -O2 quad_rule.C 
g++-5 -Wall -g -c -O2 GaussianPAW.C
g++-5 -Wall -g -O2 -D$1 -o model_PAW.x model_PAW.C \
	 gen_Gauss_quad.o gen_Gauss_quad_driver.o GaussianPAW.o \
	 gen_Ylmf.o grid.o quad_rule.o \
	 libqggev.a -lblas -llapack -lm -lgfortran -lquadmath

#=============================================
rm -rf *.o *.dSYM
#=============================================
