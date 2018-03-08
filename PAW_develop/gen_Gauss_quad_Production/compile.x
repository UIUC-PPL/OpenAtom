#=============================================
#         Instructions
# Mac or others sans fort_under: compile.x FORT UNDER_OFF
# Compilers with fort under: compile.x FORT UNDER
#=============================================
#         Clean up
rm -rf driver.x *.o *.a
#=============================================
#         Quad precision ggev
cd ../ggev_quad_pkg/compile/
./compile.x
rm -rf *.o
cd -
/bin/cp ./../ggev_quad_pkg/compile/libqggev.a .
#=============================================
#        Quad precision gen gauss quad 
g++ -c -Wall -g -O2 -D$1 gen_Gauss_quad.C
g++ -c -Wall -g -O2 gen_Gauss_quad_driver.C
#=============================================
#        Link
g++ -Wall -g -O2 -o main.x main.C \
        gen_Gauss_quad.o gen_Gauss_quad_driver.o \
        libqggev.a -lblas -llapack -lm -lgfortran -lquadmath
#=============================================
/bin/rm -rf *dSYM
/bin/rm -rf *.o *.a
#=============================================

