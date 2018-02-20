rm -rf driver.x *.o
g++ -c -std=c++11 -Wall -g -O2  gen_Gauss_quad.C gen_Gauss_quad_driver.C
g++ -std=c++11 -Wall -g -O2 driver.C -o driver.x gen_Gauss_quad.o gen_Gauss_quad_driver.o -lblas -llapack -lquadmath -lm
rm -rf *dSYM 
