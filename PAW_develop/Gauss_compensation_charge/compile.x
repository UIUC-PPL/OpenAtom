rm model_PAW.x grid.o quad_rule.o GaussianPAW.o gen_fgrid.o gen_Ylmf.o
g++ -Wall -g -c -O2 grid.C 
g++ -Wall -g -c -O2 gen_fgrid.C 
g++ -Wall -g -c -O2 gen_Ylmf.C 
g++ -Wall -g -c -O2 quad_rule.C 
g++ -Wall -g -c -O2 GaussianPAW.C
g++ -Wall -g -O2 -o model_PAW.x model_PAW.C GaussianPAW.o gen_fgrid.o gen_Ylmf.o grid.o quad_rule.o -lm
rm grid.o quad_rule.o GaussianPAW.o gen_fgrid.o gen_Ylmf.o
