/bin/rm -rf ../runable/pawmain_charm.x *.o

INCLUDE='-I../include -I../main -I../interface -I../PAWrhoChare -I../math -I../PAWrhoPhysics -I../AtmsGrp -I./'
charmc='/home/gmartyna/charm/bin/charmc'
$charmc -O2 -Wall -c -g -D_PARALLEL_DEBUG_ $1 $INCLUDE

/bin/rm -rf *.o
