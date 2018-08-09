/bin/rm -rf ../runable/pawmain_charm.x *.decl.h *.def.h *.o
INCLUDE='-I../include -I../main -I../interface -I../PAWrhoChare -I../math -I../PAWrhoPhysics -I../AtmsGrp -I./'
charmc='/Users/Qi-MAC/Work/openatom/charm/netlrts-darwin-x86_64-new/bin/charmc'
sanitize='-fsanitize=address -O0 -fno-omit-frame-pointer'

$charmc $INCLUDE -c -D_PARALLEL_DEBUG_ ../main/pawmain.ci   
$charmc $INCLUDE -c -D_PARALLEL_DEBUG_ ../main/cleanexit.ci 
$charmc $INCLUDE -c -D_PARALLEL_DEBUG_ ../PAWrhoChare/PAWrho.ci $INCLUDE
$charmc $INCLUDE -c -D_PARALLEL_DEBUG_ ../AtmsGrp/ATMSGRP.ci $INCLUDE

$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../interface/interface_hand.C 
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../interface/configure.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../math/pup_utils.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../math/mathlib.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../PAWrhoChare/PAWrho.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../main/cleanexit.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../main/pawmain.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../PAWrhoPhysics/fastAtoms.C
$charmc $INCLUDE $sanitize -c -g -D_PARALLEL_DEBUG_ ../AtmsGrp/ATMSGRP.C

$charmc $sanitize -O0 -Wall -g -D_PARALLEL_DEBUG_ -o ../runable/pawmain_charm.x interface_hand.o configure.o pup_utils.o mathlib.o PAWrho.o cleanexit.o pawmain.o fastAtoms.o ATMSGRP.o -lm
/bin/rm -rf *.decl.h *.def.h *.o

/bin/cp ./charmrun ../runable/.
