#=========================================================
#compile all the quad precision routines and make library
#---------------------------------------------------------
# qblas routines
f77 -c -O2 -g qaxpy.f
f77 -c -O2 -g qcopy.f
f77 -c -O2 -g qdot.f
f77 -c -O2 -g qgemv.f
f77 -c -O2 -g qgemm.f
f77 -c -O2 -g qger.f
f77 -c -O2 -g qnrm2.f
f77 -c -O2 -g qrot.f
f77 -c -O2 -g qscal.f
f77 -c -O2 -g qswap.f
f77 -c -O2 -g qtrmm.f
f77 -c -O2 -g qtrmv.f
#---------------------------------------------------------
# quad (laplack) routines
f77 -c -O2 -g qgeqr2.f
f77 -c -O2 -g qgeqrf.f
f77 -c -O2 -g qggbak.f
f77 -c -O2 -g qggbal.f
f77 -c -O2 -g qggev.f
f77 -c -O2 -g qgghrd.f
f77 -c -O2 -g qhgeqz.f
f77 -c -O2 -g qisnan.f
f77 -c -O2 -g qlabad.f
f77 -c -O2 -g qlacpy.f
f77 -c -O2 -g qladiv.f
f77 -c -O2 -g qlag2.f
f77 -c -O2 -g qlaisnan.f
f77 -c -O2 -g qlaln2.f
f77 -c -O2 -g qlange.f
f77 -c -O2 -g qlanhs.f
f77 -c -O2 -g qlapy2.f
f77 -c -O2 -g qlapy3.f
f77 -c -O2 -g qlarfb.f
f77 -c -O2 -g qlarf.f
f77 -c -O2 -g qlarfg.f
f77 -c -O2 -g qlarft.f
f77 -c -O2 -g qlartg.f
f77 -c -O2 -g qlascl.f
f77 -c -O2 -g qlaset.f
f77 -c -O2 -g qlassq.f
f77 -c -O2 -g qlasv2.f
f77 -c -O2 -g qorg2r.f
f77 -c -O2 -g qorgqr.f
f77 -c -O2 -g qorm2r.f
f77 -c -O2 -g qormqr.f
f77 -c -O2 -g qtgevc.f
#---------------------------------------------------------
# quad (laplack) util routines
f77 -c -O2 -g ilaqlc.f
f77 -c -O2 -g ilaqlr.f
f77 -c -O2 -g qlmach.f
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
