rm -rf driver.x
g++ -std=c++11 -Wall -g -O2 driver.C -o driver.x -lblas -llapack -lm
rm -rf *dSYM 
