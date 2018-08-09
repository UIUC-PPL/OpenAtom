/bin/rm -rf ../runable/pawmain_scalar.x *.o
g++ -O2 -c -g -Wall -DCHARM_OFF ../interface/configure.C -I../interface -I../include
g++ -O2 -c -g -Wall -DCHARM_OFF ../interface/interface_hand.C -I../interface -I../include
g++ -O2 -c -g -Wall -DCHARM_OFF ../main/pawmain.C -I../interface -I../include -I../main
g++ -O2    -g -Wall -DCHARM_OFF -o ../runable/pawmain_scalar.x pawmain.o configure.o interface_hand.o -lm
/bin/rm -rf *.o
