// constants for approx to erfc
#ifndef M_PI_QI
#define M_PI_QI 3.14159265358979323846264338327950288419716939937510
#endif
#define PERFC  (0.3614)
#define CERFC1 (0.2041422096422003)
#define CERFC2 (0.1997535956961481)
#define CERFC3 (0.2213176596405576)
#define CERFC4 (0.03360430734640255)
#define CERFC5 (0.4732592578721755) 
#define CERFC6 (-0.509078520069735)
#define CERFC7 (0.6772631491947646)
#define CERFC8 (-0.369912979092217)
#define CERFC9 (0.06965131976970335)
#define DCERFC1 (1.0*CERFC1)
#define DCERFC2 (2.0*CERFC2)
#define DCERFC3 (3.0*CERFC3)
#define DCERFC4 (4.0*CERFC4)
#define DCERFC5 (5.0*CERFC5)
#define DCERFC6 (6.0*CERFC6)
#define DCERFC7 (7.0*CERFC7)
#define DCERFC8 (8.0*CERFC8)
#define DCERFC9 (9.0*CERFC9)
#define PRE_ERFC (2.0/sqrt(M_PI_QI))