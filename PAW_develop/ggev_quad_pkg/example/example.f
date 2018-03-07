      program main
      double precision junk,junk2
      real*16 junk3, junk4
      parameter(junk2=3.d0)      
      parameter(junk4=3.q0)      
      junk  = 1.0d0
      junk3 = 1.0q0
      write(*,*)junk,sizeof(junk)
      write(*,*)junk2,sizeof(junk2)
      write(*,*)junk3,sizeof(junk3)
      write(*,*)junk4,sizeof(junk4)
      stop
      end
