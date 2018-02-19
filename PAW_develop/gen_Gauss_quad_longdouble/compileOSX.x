rm -rf driver.x *.o
g++ -c -std=c++11 -Wall -g -O2  gen_Gauss_quad.C
g++ -std=c++11 -Wall -g -O2 driver.C -o driver.x gen_Gauss_quad.o -lblas -llapack -lm
rm -rf *dSYM 
