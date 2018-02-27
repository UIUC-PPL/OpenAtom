      program main
c      real*8 junk
      double precision junk
      real*16 junk2
      junk  = 1.0d0
      junk2 = 1.0q0
      write(*,*)junk,junk2, sizeof(junk), sizeof(junk2)
      stop
      end
