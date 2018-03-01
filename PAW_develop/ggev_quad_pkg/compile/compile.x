/bin/rm *.o
/bin/rm *.a
#=========================================================
#compile all the quad precision routines and make library
#---------------------------------------------------------
# qblas routines
gfortran -c -Wall  -O2 -g ../qblas/iqamax.f
gfortran -c -Wall  -O2 -g ../qblas/qaxpy.f
gfortran -c -Wall  -O2 -g ../qblas/qcopy.f
gfortran -c -Wall  -O2 -g ../qblas/qdot.f
gfortran -c -Wall  -O2 -g ../qblas/qgemv.f
gfortran -c -Wall  -O2 -g ../qblas/qgemm.f
gfortran -c -Wall  -O2 -g ../qblas/qger.f
gfortran -c -Wall  -O2 -g ../qblas/qnrm2.f
gfortran -c -Wall  -O2 -g ../qblas/qrot.f
gfortran -c -Wall  -O2 -g ../qblas/qscal.f
gfortran -c -Wall  -O2 -g ../qblas/qswap.f
gfortran -c -Wall  -O2 -g ../qblas/qtrmm.f
gfortran -c -Wall  -O2 -g ../qblas/qtrmv.f
#---------------------------------------------------------
# quad (laplack) routines
gfortran -c -Wall  -O2 -g ../quad/qgeqr2.f
gfortran -c -Wall  -O2 -g ../quad/qgeqrf.f
gfortran -c -Wall  -O2 -g ../quad/qggbak.f
gfortran -c -Wall  -O2 -g ../quad/qggbal.f
gfortran -c -Wall  -O2 -g ../quad/qggev.f
gfortran -c -Wall  -O2 -g ../quad/qgghrd.f
gfortran -c -Wall  -O2 -g ../quad/qhgeqz.f
gfortran -c -Wall  -O2 -g ../quad/qisnan.f
gfortran -c -Wall  -O2 -g ../quad/qlabad.f
gfortran -c -Wall  -O2 -g ../quad/qlacpy.f
gfortran -c -Wall  -O2 -g ../quad/qladiv.f
gfortran -c -Wall  -O2 -g ../quad/qlag2.f
gfortran -c -Wall  -O2 -g ../quad/qlaisnan.f
gfortran -c -Wall  -O2 -g ../quad/qlaln2.f
gfortran -c -Wall  -O2 -g ../quad/qlange.f
gfortran -c -Wall  -O2 -g ../quad/qlanhs.f
gfortran -c -Wall  -O2 -g ../quad/qlapy2.f
gfortran -c -Wall  -O2 -g ../quad/qlapy3.f
gfortran -c -Wall  -O2 -g ../quad/qlarfb.f
gfortran -c -Wall  -O2 -g ../quad/qlarf.f
gfortran -c -Wall  -O2 -g ../quad/qlarfg.f
gfortran -c -Wall  -O2 -g ../quad/qlarft.f
gfortran -c -Wall  -O2 -g ../quad/qlartg.f
gfortran -c -Wall  -O2 -g ../quad/qlascl.f
gfortran -c -Wall  -O2 -g ../quad/qlaset.f
gfortran -c -Wall  -O2 -g ../quad/qlassq.f
gfortran -c -Wall  -O2 -g ../quad/qlasv2.f
gfortran -c -Wall  -O2 -g ../quad/qorg2r.f
gfortran -c -Wall  -O2 -g ../quad/qorgqr.f
gfortran -c -Wall  -O2 -g ../quad/qorm2r.f
gfortran -c -Wall  -O2 -g ../quad/qormqr.f
gfortran -c -Wall  -O2 -g ../quad/qtgevc.f
#---------------------------------------------------------
# quad (laplack) util routines
gfortran -c -Wall  -O2 -g ../qutil/ilaqlc.f
gfortran -c -Wall  -O2 -g ../qutil/ilaqlr.f
gfortran -c -Wall  -O2 -g ../qutil/qlmach.f
#---------------------------------------------------------
# make the library
ar rvs libqggev.a \
iqamax.o \
qaxpy.o \
qcopy.o \
qdot.o \
qgemm.o \
qgemv.o \
qger.o \
qnrm2.o \
qrot.o \
qscal.o \
qswap.o \
qtrmm.o \
qtrmv.o \
qgeqr2.o \
qgeqrf.o \
qggbak.o \
qggbal.o \
qggev.o \
qgghrd.o \
qhgeqz.o \
qisnan.o \
qlabad.o \
qlacpy.o \
qladiv.o \
qlag2.o \
qlaisnan.o \
qlaln2.o \
qlange.o \
qlanhs.o \
qlapy2.o \
qlapy3.o \
qlarfb.o \
qlarf.o \
qlarfg.o \
qlarft.o \
qlartg.o \
qlascl.o \
qlaset.o \
qlassq.o \
qlasv2.o \
qorg2r.o \
qorgqr.o \
qorm2r.o \
qormqr.o \
qtgevc.o \
ilaqlc.o \
ilaqlr.o \
qlmach.o 
#=========================================================
