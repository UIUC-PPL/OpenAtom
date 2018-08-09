/bin/rm -rf ../runable/pawmain_charm.x *.o

INCLUDE='-I../include -I../main -I../interface -I../PAWrhoChare -I../math -I../PAWrhoPhysics -I../AtmsGrp -I./'
charmc='/Users/Qi-MAC/Work/openatom/charm/netlrts-darwin-x86_64-new/bin/charmc'
$charmc -c $1 $INCLUDE

/bin/rm -rf *.o
