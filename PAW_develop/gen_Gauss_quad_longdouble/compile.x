# for mac, its compile.x FORT UNDER_OFF
# for others, compile.x FORT UNDER
rm -rf driver.x *.o
g++ -c -std=c++11 -Wall -g -O2 -D$1 gen_Gauss_quad.C
g++ -std=c++11 -Wall -g -O2 -D$1 driver.C -o driver.x gen_Gauss_quad.o -lblas -llapack -lm
rm -rf *dSYM 
