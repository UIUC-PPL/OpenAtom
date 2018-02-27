#=========================================================
#compile all the quad precision routines and make library
#---------------------------------------------------------
# qblas routines
f77 -c -O2 -g qblas/qaxpy.f
f77 -c -O2 -g qblas/qcopy.f
f77 -c -O2 -g qblas/qdot.f
f77 -c -O2 -g qblas/qgemv.f
f77 -c -O2 -g qblas/qgemm.f
f77 -c -O2 -g qblas/qger.f
f77 -c -O2 -g qblas/qnrm2.f
f77 -c -O2 -g qblas/qrot.f
f77 -c -O2 -g qblas/qscal.f
f77 -c -O2 -g qblas/qswap.f
f77 -c -O2 -g qblas/qtrmm.f
f77 -c -O2 -g qblas/qtrmv.f
#---------------------------------------------------------
# quad (laplack) routines f77 -c -O2 -g qgeqr2.f
f77 -c -O2 -g quad/qgeqrf.f
f77 -c -O2 -g quad/qggbak.f
f77 -c -O2 -g quad/qggbal.f
f77 -c -O2 -g quad/qggev.f
f77 -c -O2 -g quad/qgghrd.f
f77 -c -O2 -g quad/qhgeqz.f
f77 -c -O2 -g quad/qisnan.f
f77 -c -O2 -g quad/qlabad.f
f77 -c -O2 -g quad/qlacpy.f
f77 -c -O2 -g quad/qladiv.f
f77 -c -O2 -g quad/qlag2.f
f77 -c -O2 -g quad/qlaisnan.f
f77 -c -O2 -g quad/qlaln2.f
f77 -c -O2 -g quad/qlange.f
f77 -c -O2 -g quad/qlanhs.f
f77 -c -O2 -g quad/qlapy2.f
f77 -c -O2 -g quad/qlapy3.f
f77 -c -O2 -g quad/qlarfb.f
f77 -c -O2 -g quad/qlarf.f
f77 -c -O2 -g quad/qlarfg.f
f77 -c -O2 -g quad/qlarft.f
f77 -c -O2 -g quad/qlartg.f
f77 -c -O2 -g quad/qlascl.f
f77 -c -O2 -g quad/qlaset.f
f77 -c -O2 -g quad/qlassq.f
f77 -c -O2 -g quad/qlasv2.f
f77 -c -O2 -g quad/qorg2r.f
f77 -c -O2 -g quad/qorgqr.f
f77 -c -O2 -g quad/qorm2r.f
f77 -c -O2 -g quad/qormqr.f
f77 -c -O2 -g quad/qtgevc.f
#---------------------------------------------------------
# quad (laplack) util routines
f77 -c -O2 -g qutil/ilaqlc.f
f77 -c -O2 -g quitl/ilaqlr.f
f77 -c -O2 -g quitl/qlmach.f
#---------------------------------------------------------
# make the library
ar rvs libqggev.a \
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
