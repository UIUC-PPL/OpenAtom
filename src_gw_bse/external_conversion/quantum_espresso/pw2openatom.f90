!
! Sep. 2014   minjung.kim@yale.edu
!
! Last modified: Nov. 2016 by MK

! note:
! qexml.f90 had been moved to Modules directory from PP/src directory
! and some subroutines had also been changed.
! qexml_read_bands does not exist anymore
! qexml_read_bands_pw should be used instead.
! this affects read_eig_occ subroutine.



Program converter

   use qexml_module

   implicit none

   character(len=5) :: codename = 'PW2OA'
   integer, parameter :: iunit = 10
   integer, parameter :: dp = kind(1.0d0)

   type cell_info
      real(dp) :: alat
      real(dp) :: a1(3), a2(3), a3(3)
      real(dp) :: b1(3), b2(3), b3(3)
      character(len=256) :: aunit, bunit
   end type
   
   type(cell_info) :: cell
  
   integer :: ngk_tot, nb
   integer, allocatable :: ipwk(:), igk_all(:)

   integer :: i
   integer :: ierr
   
   character(256) :: work_dir, prefix, dirname, filename, fsysname
   logical :: shift_flag ! whether wavefunctions are shifted or not
   integer :: nk ! it should be decided at first
   integer :: nspin

   ! total number of g vectors
   integer :: ngm
   ! master gvector list
   integer, allocatable :: master_gv(:,:)
   ! eigenvalues and its occupancies   
   real(dp), allocatable :: eig(:,:,:), occ(:,:,:)
   ! wavefunction
   complex(dp), allocatable :: wfn(:)
   ! k points and their weight
   real(dp), allocatable :: xk(:,:), wk(:)
   ! fftgrid for density (dense fft grid)
   integer :: fftsize(3)
   ! wavefunction cutoff
   real(dp) :: ecutwfc
   ! array to collect all g list in all k points
   !integer, allocatable :: gkidx_allkpt(:)
   ! kmax_cp
   integer :: kmax_cp(3)
   ! number of g vectors
   integer :: ngdoublepack, ngkpt, ngoa
   ! doublepack has only half sphere
   logical :: doublepack
   ! new g list that fits to OpenAtom
   integer, allocatable :: glist_doublepack(:,:), glist_kpt(:,:)
   ! mapping between OA g list and QE g list
   integer, allocatable :: idxmap(:)
   ! number of planewaves at each k points (same numbers!)
   !integer :: ncoeff 
   
   integer :: ib, ispin, ik, npwk
   integer, allocatable :: idxgk(:)

   ! GPP calculation related variables
   integer :: nr(3)
   real(dp), allocatable :: rho(:,:,:)
   logical :: gpp_flag
   
   ! Vxc variables
   logical :: Vxc_flag
   character(256) :: fname_vxc

   ! optional variables
   integer :: nbandschunk ! number of bands to be read at each loop
                         ! this number sometimes matter since this converter is serial!
                         ! this variable will be used in read_write_wfn subroutine

   integer :: stdin = 5

   NAMELIST /INPUT/ prefix, work_dir, fsysname, doublepack, shift_flag, gpp_flag, Vxc_flag, nbandschunk

!---------------------------------------------------------------------
! starts here
!---------------------------------------------------------------------

   ! default setting
   work_dir = './'
   doublepack = .false.
   shift_flag = .false.
   gpp_flag = .false.
   Vxc_flag = .false.
   fname_vxc = 'Vxcr.dat'
   nbandschunk = 0

   read( stdin, INPUT, iostat=ierr )

   ! initialize QEXML library
   
   dirname = trim(work_dir) // '/' // trim(prefix) // '.save/'
   call qexml_init( iunit, Dir=dirname )

   filename = trim(dirname)//"data-file.xml"
   
   call qexml_openfile( filename, "read", IERR=ierr )


   ! read cell information
   call get_cell_info( cell )

   ! read number of k points = nk
   call qexml_read_bands_info( NUM_K_POINTS=nk, IERR=ierr)

   ! read BZ: k points and their weight factor
   allocate( xk(3,nk), wk(nk) )
   call qexml_read_bz( XK=xk, WK=wk, IERR=ierr )

   ! get master g vector
   call qexml_read_planewaves( NGM=ngm, IERR=ierr )
   allocate( master_gv(3,ngm) )
   call qexml_read_planewaves( IGV=master_gv, IERR=ierr )

   ! number of gvectors in each k
   allocate( ipwk( nk ) )
   call read_pw_index( nk, ipwk )

   ! Read number of bands and number of spins
   call qexml_read_bands_info( NBND=nb, NSPIN=nspin, IERR=ierr)

   ! read eigenvalues and occupancies
   allocate( eig(nb, nk, nspin), occ(nb, nk, nspin) )
   call read_eig_occ( nspin, nk, nb, eig, occ )

   ! read FFTsize (dense FFT grid for density)
   call qexml_read_planewaves( NR1=fftsize(1), NR2=fftsize(2), NR3=fftsize(3), IERR=ierr )

   ! read wavefunction cutoff
   call qexml_read_planewaves( ECUTWFC=ecutwfc, IERR=ierr )

   !---------------------------------------------------------------------------
   !   here it starts setting the same number of k points at each k points
   !---------------------------------------------------------------------------
   ! Adjust ecutwfc so that all g vectors at each k points are inside of the new ecutwfc
   ! of course, 0.5*|g+k|^2 < ecutwfc(original) but g itself can be outside of ecutwfc(original)
   call adjust_ecutwfc( cell, ecutwfc, nk, xk )

! --- lets not set gkidx_allkpt
!  ! set gkidx_allkpt and initialize
!  allocate( gkidx_allkpt( ngm ) )
!  gkidx_allkpt(:) = 0 ! 0 - do not include, 1 - include
!  
!  do ik = 1, nk
!     ! read number of plane-waves at this k point
!     call qexml_read_gk( ik, NPWK=npwk, IERR=ierr)
!     allocate( idxgk( npwk ) )
!     ! read plane-wave index to idxgk. idxgk connects to master_gv
!     call qexml_read_gk( ik, index=idxgk, IERR=ierr )
!     ! now, update gkidx_allkpt
!     call set_gkidx_allkpt( npwk, idxgk, ngm, master_gv, gkidx_allkpt )

!     deallocate( idxgk )
!  enddo
!  ! finally, we got the number of planewaves at each k point (each k point has the same number of planewaves)
!  ncoeff = sum(gkidx_allkpt)
!
!  ! finished
   
   !print*, "(new) number of planewaves at each k points are:", ncoeff

   ! Here, we set the glist for OpenAtom
   ! There are two similar subroutines, one for doublepack, one of k points calculations
   ! doublepack means we save the wavefunctions only for the half sphere
   ! ngdoublepack = number of g vectors for doublepack
   ! ngkpt = number of g vectors for kpoints calcualtions
   ! ngkpt >= ncoeff 
  
   ! set kmax_cp
   ! initialization
   kmax_cp = 0 
   do ispin = 1, nspin
      do ik = 1, nk
         ! read number of plane-waves at this k point
         call qexml_read_gk( ik, NPWK=npwk, IERR=ierr )
         allocate( idxgk( npwk ) )
         ! read planewave index to idxgk. idxgk connects to master_gv
         call qexml_read_gk( ik, index=idxgk, IERR=ierr )
         ! now, we want to convert wavefunctions to openatom output format
         ! all k points should include the same number of g vectors
         call get_kmax( ngm, master_gv, npwk, idxgk, kmax_cp )
         deallocate( idxgk )
      enddo
   enddo
   ! to be safe, let's add 1 in each direction
   !kmax_cp(:) = kmax_cp(:) + 1

   print*, 'now we are modifying fftsize for density by calling radixme routine'
   print*, 'original fft size: ', fftsize
   call radixme(kmax_cp,fftsize)
   print*, 'new fft size after radixme: ', fftsize

   ! let's set the list of G vectors here
   ! warning: gamma point will contain only half sphere if doublepack flag is on
   if (doublepack) then
      call countkvec3d_sm( kmax_cp, ngdoublepack, ecutwfc, cell )
      allocate( glist_doublepack(3,ngdoublepack) )
      call setkvec3d_sm_simple( kmax_cp, ngdoublepack, ecutwfc, cell, glist_doublepack )
   else
      call countkvec3d_sm_kpt( kmax_cp, ngkpt, ecutwfc, cell )
      allocate( glist_kpt(3,ngkpt) )
      call setkvec3d_sm_kpt( kmax_cp, ngkpt, ecutwfc, cell, glist_kpt)
   endif

   ! now we need to map the two gvector lists to shuffle things around
   ! mapping between OA and QE g vector list


   ! read and write wfc into state.out file
   do ispin = 1, nspin
      do ik = 1, nk

         print*, 'reading and writing... spin index:',ispin,'   k point index:',ik
         ! get number of plane-waves
         call qexml_read_gk( ik, NPWK=npwk, IERR=ierr)
         allocate( idxgk( npwk ) ) 
         ! read planewave index to idxgk, idxgk connects to master_gv
         call qexml_read_gk( ik, index=idxgk, IERR=ierr)

         ! now we need to map the two gvector lists to shuffle things around
         ! mapping between OA and QE g vector list
         if (doublepack) then
            allocate( idxmap(ngdoublepack) )
            call mapping_glist( doublepack, ngdoublepack, glist_doublepack, ngm, master_gv, npwk, idxgk, idxmap )
         else
            allocate( idxmap(ngkpt) )
            call mapping_glist( doublepack, ngkpt, glist_kpt, ngm, master_gv, npwk, idxgk, idxmap )
         endif
         ! if k is at gamma point and doublepack flag is on
         ! it saves only half sphere
         if ( ik==1 .and. doublepack .eqv. .true. ) then
            call read_write_wfn( nb, ik, ispin, npwk, idxgk, master_gv, ngm,&
                 fftsize, shift_flag, doublepack, ngdoublepack, glist_doublepack, idxmap, nbandschunk)
         ! if not gamma point or not doublepack
         else
            call read_write_wfn( nb, ik, ispin, npwk, idxgk, master_gv, ngm,&
                 fftsize, shift_flag, doublepack, ngkpt, glist_kpt, idxmap, nbandschunk)
         endif
         deallocate( idxmap, idxgk )
      enddo! ik loop
   enddo! is loop
   ! wfn order: is-ik-ib-ig

   call qexml_closefile( "read", IERR=ierr )

   if (doublepack) ngoa = ngdoublepack
   if ( .not. doublepack) ngoa = ngkpt


   if ( shift_flag .eqv. .false.) then
      ! Write system information
      call write_system_info( cell, ngoa, xk, wk, &
                       nspin, nk, nb, fftsize, fsysname )
   else
      continue
   endif

   ! write eigenvalues and occupation numbers
   do ispin = 1, nspin
      do ik = 1, nk
         call write_eig_occ(ispin, ik, nb, eig(:,ik,ispin), occ(:,ik,ispin), shift_flag)
      enddo
   enddo
   

! If RHO flag is on, create rho.dat file for GPP calculations
   if ( gpp_flag ) then
      call qexml_read_rho( NR1=nr(1), NR2=nr(2), NR3=nr(3), IERR=ierr )
      allocate( rho(nr(1),nr(2),nr(3)) )
      call write_rho( nr, rho )
   endif

! If Vxc_flag is on, create vxc.dat file for sigma calculations
!  if ( Vxc_flag ) then
!     call write_vxc( fname_vxc, codename  )
!  endif


   deallocate(xk, wk, eig, occ, master_gv )





!****************************************************************!
!
!                        SUBROUTINES
!
!****************************************************************!

contains


!************ subroutines ********************************************

!subroutine write_vxc( fname_vxc, codename )
!  ! adopted from pw2bgw.f90 
!  use environment,          only: environment_start, environment_end
!  use fft_base,             only: dfftp
!  use fft_interfaces,       only: fwfft
!  use ener,                 only: etxc, vtxc
!  use gvect,                only: ngm, ngm_g, ig_l2g, nl, mill
!  use lsda_mod,             only: nspin
!  use scf,                  only: rho, rho_core, rhog_core
!  use wavefunctions_module, only: psic

!  character(len=100), intent(in) :: fname_vxc
!  character(len=5), intent(in) :: codename
!  integer, parameter :: iu = 13

!  ! local variables
!  integer :: id, ig, is, ir
!  integer :: nr, ns, nd
!  integer, allocatable :: gvec(:,:)
!  
!  real(dp), allocatable :: vxcr(:,:)
!  complex(dp), allocatable :: vxcg(:,:)
!  
!  
!  ! assign variables
!  ns = nspin       ! number of spin
!  nr = dfftp%nnr   ! number of r grid
!  nd = 3           ! dimension


!  open(unit=iu,file=trim(fname_vxc),status='replace',form='formatted')

!  call environment_start( codename )
!  
!  call read_file()
!  
!  allocate( vxcr( nr, ns ) )
!  vxcr(:,:) = 0.d0
!  rho_core(:) = 0.d0
!  rhog_core(:) = 0.d0
!  call v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxcr )
!  
!  ! ngm is the number of g vectors in this process (with gamma tricks, only G>=0)
!  ! ngm_g is the total number of g vectors
!  ! if serial, ngm = ngm_g
!  allocate( gvec( 3, ngm ) )
!  gvec = 0.d0
!  do ig = 1, ngm
!     do id = 1, nd
!        ! ig_l2g converts the local G-vector index into the global index
!        ! mill is a miller index
!        gvec( id, ig_l2g(ig) ) = mill( id, ig )
!     enddo
!  enddo
!  
!  ! now we need to get vxcg
!  ! nl: fft index to G index
!  allocate( vxcg( ngm_g, ns ) )
!  do is = 1, ns
!     do ir = 1, nr
!        psic(ir) = complex( vxcr(ir,is), 0.d0 )
!     enddo
!     ! forward Fourier transform
!     call fwfft( 'Dense', psic, dfftp )
!     do ig = 1, ngm
!        vxcg( ig_l2g(ig), is ) = psic( nl(ig) )
!     enddo
!  enddo  
!  
!  do is = 1, ns
!     do ig = 1, ngm_g
!        write(iu,*) vxcg( ig, is ), ( gvec(i,ig), i=1,3 )
!     enddo
!  enddo
!  
!  deallocate( gvec, vxcr, vxcg )
!  
!  call environment_end( codename ) 
!  
!end subroutine




   

subroutine write_eig_occ( ispin, ik, nstate, eig, occ, shift_flag )

   integer, intent(in) :: ispin, ik, nstate
   real(dp), dimension(nstate), intent(in) :: eig, occ
   logical, intent(in) :: shift_flag
   character(len=100) :: fdir, fname
   
   integer :: iunit = 30
   
   ! directory name
   if (shift_flag .eqv. .false.) then
      if ( ik .le. 10 ) then
         write( fdir, '( "./STATES/Spin.", I1, "_Kpt.", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fdir, '( "./STATES/Spin.", I1, "_Kpt.", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   elseif (shift_flag .eqv. .true. ) then
      if ( ik .le. 10 ) then
         write( fdir, '( "./STATES/Spin.", I1, "_Kpt.0", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fdir, '( "./STATES/Spin.", I1, "_Kpt.0", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   endif   
   
   
   fname = 'eigenvalues.in'
   
   
   ! open file
   
   open(iunit, file=trim(fdir)//trim(fname), form='formatted', status='unknown')
   
   do i = 1, nstate
      write(iunit,*) eig(i)  !, occ(i) ! openatom won't read occupation 
   enddo
   close(iunit)

end subroutine





  
  


!------------------------------------------------------------------------------
! name of the output file: state\\ib\\.out.gz
subroutine read_write_wfn( nb, ik, ispin, npwk, idxgk, master_gv, ngm, &
                           fftsize, shift_flag, doublepack, ngoa, glist_oa, idxmap, nbandschunk)

   integer, intent(in) :: nb, ik, ispin, npwk, fftsize(3), ngm
   integer, intent(in) :: master_gv(3,ngm)
   integer, intent(in) :: idxgk(npwk)
   logical, intent(in) :: shift_flag
   logical, intent(in) :: doublepack
   integer, intent(in) :: ngoa ! number of g vectors in openatom g list
   integer, intent(in) :: glist_oa(3,ngoa)
   integer, intent(in) :: idxmap(ngoa) ! mapping array to write wavefunctions in OA g vector order
   integer, intent(in) :: nbandschunk ! this sets nblockband

   complex(dp), allocatable :: wfn(:,:)
   integer :: iunit = 30
   character (len=100) :: fplace, fname
   integer :: open_status

   integer :: i, j
   ! total number of gvector to be written
   integer :: ngktot
   logical :: notfound
   real(dp) :: ZERO
   integer :: nblockband  ! number of bands to be read at each qexml_read_wfc function call
   integer :: iloop, nloop
   integer :: ib, ibstart, ibend, ibcounter, wfcounter, jstart
   integer :: counter
 
   ! variables to decide nblockband if nbandschunk is not given
   complex(dp) :: complexnumber
   integer :: sizeofcomplex
   integer :: memoryinGB = 10 ! let assume that the single node has 10GB memory.


   ZERO = dble(0)

   ! set up the directory where state files are stored
   if ( shift_flag .eqv. .false. ) then
      if ( ik .le. 10 ) then
         write( fplace, '( "./STATES/Spin.", I1, "_Kpt.", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fplace, '( "./STATES/Spin.", I1, "_Kpt.", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   elseif( shift_flag .eqv. .true. ) then
      if ( ik .le. 10 ) then
         write( fplace, '( "./STATES/Spin.", I1, "_Kpt.0", I1, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1
      elseif ( ik .le. 100 ) then
         write( fplace, '( "./STATES/Spin.", I1, "_Kpt.0", I2, "_Bead.0_Temper.0/" )' ) ispin-1, ik-1   
      else
         print*, 'Cannot locate directory because you have more than 100 k points'
         stop
      endif
   endif
   
   ! I assume that the login node (or computing node) has at least 10GB of memory (usually they offer much bigger memory)

   if ( nbandschunk .ne. 0 ) then
      nblockband = nbandschunk
   else
      sizeofcomplex = sizeof(complexnumber)
      ! the single state file takes npwk*sizeofcomplex bytes of memeory
      nblockband = memoryinGB*10**9/(npwk*sizeofcomplex) 
   endif 

   ! if nblockband is less than total number of bands
   if ( nblockband .ge. nb ) then 
      nblockband = nb
      nloop = 0
   endif

   print*, 'number of bands to be read at each loop:', nblockband

   do iloop = 1, nloop+1
      if (iloop .lt. nloop+1) then
         ibstart = 1 + (iloop-1)*nblockband
         ibend = iloop*nblockband
      elseif ( iloop .eq. nloop+1 ) then
         ibstart = 1 + nloop*nblockband
         ibend = nb
      else
         print*, 'No band index found in routine read_write_wfc, program exits'
         stop
      endif

      print*, 'reading wavefunctions from',ibstart,'to',ibend
      
      ! allocate memory to save savefunctions
      ! total number of bands for iloop is (ibend-ibstart+1)
      allocate( wfn( npwk, ibend-ibstart+1 ) )

      if ( ispin .eq. 1 ) then
         call qexml_read_wfc ( IBNDS=ibstart, IBNDE=ibend, IK=ik, WF=wfn, IERR=ierr)
      else
         call qexml_read_wfc ( IBNDS=ibstart, IBNDE=ibend, IK=ik, ISPIN=ispin, WF=wfn, IERR=ierr)
      endif

      ! writing
      do ib = ibstart, ibend
         ! set-up the file name
         if ( ib .lt. 10) write( fname, '( "state", I1, ".out" )' )ib
         if ( ib .ge. 10 .and. ib .lt. 100) &
            write( fname, '( "state", I2, ".out" )' ) ib
         if ( ib .ge. 100 .and. ib .lt. 1000) &
            write( fname, '( "state", I3, ".out" )' ) ib
         if ( ib .ge. 1000 .and. ib .lt. 10000) &
            write( fname, '( "state", I4, ".out" )' ) ib
         if ( ib .ge. 10000 .and. ib .lt. 100000) &
            write( fname, '( "state", I2, ".out" )' ) ib
         ! open states file name
         open(iunit, file=trim(fplace)//trim(fname), form='formatted', status='replace', iostat=open_status)
         if (open_status .ne. 0 ) then
            print*, 'Could not open file ', trim(fplace)//trim(fname), '. Please check your STATES directory'
         endif
         !---------------------------------------------
         ! writing starts here
         write(iunit,"(I10,3I4)") ngoa, fftsize(1:3)

         ! we need to set some counter 
         ! wfn(1:npwk,1:nblockband)
         ibcounter = ib - ibstart ! this will set the starting point

         ! initialization
         jstart = 1 ! initialize
         counter = 1

         print*, 'writing band #', ib

         do i = 1, ngoa
            if ( idxmap(i) .ne. -1 ) then
               write(iunit,99) real( wfn(idxmap(i),ibcounter+1) ), imag( wfn(idxmap(i),ibcounter+1) ), &
                    glist_oa(1:3,i)
            elseif ( idxmap(i) .eq. -1 ) then
               write(iunit,99) ZERO, ZERO, glist_oa(1:3,i)
            endif
         enddo
         



! let's not use below---------- BEGIN DO NOT USE
!        do i = 1, ngm
!            notfound = .true.
!           if ( gkidx_allkpt(i) == 0 ) then
!           ! do nothing
!           elseif ( gkidx_allkpt(i) == 1) then
!              ! let's find it
!              do j = jstart, npwk
!                 if( idxgk(j) == i ) then
!                    jstart = j + 1 ! updating searching index
!                    notfound = .false.
!                    write(iunit,99) real( wfn(j,ibcounter+1) ), imag( wfn(j,ibcounter+1) ),&
!                      master_gv(1:3,idxgk(j))
!                    ! --- modification starts here
!                    ! save wavefunctions at this band index
!                    thisBandwfn(counter) = wfn(j,ibcounter+1)
!                    qeglist(1:3,counter) = master_gv(1:3,idxgk(j))
!                    counter = counter + 1
!                    ! finish do j loop
!                    exit
!                 endif
!              enddo
!              if( notfound ) then
!                 write(iunit,99) ZERO, ZERO, master_gv(1:3,i)
!              endif
!           else
!              ! something wrong
!              print*, '@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@'
!              print*, 'gkidx_allkpt has illegal value for the point', i
!              print*, '@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@'
!              stop
!           endif
!        enddo ! end i
! END DO NOT USE

         close(iunit)

      enddo ! end ib

      deallocate(wfn)

   enddo ! end iloop
   
   ! final call for qexml_read_wfc
   
!  FORMAT WARNING: If your G vector index is smaller than -999 (e.g. -1000), you need to change the format
99 FORMAT(2E13.5,3I5)
! old format for comparison
!99 FORMAT(2E25.15,3I10)
 
end subroutine


!------------------------------------------------------------------------------
subroutine write_system_info( cell, ncoeff, xk, wk, &
                        nspin, nk, nb, fftsize,  outname )

   integer, parameter :: iunit = 20

   type(cell_info) :: cell
   integer, intent(in) :: ncoeff
   integer, intent(in) :: nspin, nk, nb 
   integer, intent(in) :: fftsize(3)

   real(dp), intent(in) :: xk(3,nk), wk(nk)
   integer :: i, j, ik, ib, ispin
   character(len=256), intent(in) :: outname

   open(iunit,file=trim(outname),form='formatted', status='replace')
 
   ! cell lattice
   write(iunit,'(f12.6)') cell%alat
   
   ! lattice vectors
   write(iunit,'(3f12.6)') ( cell%a1(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%a2(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%a3(i), i=1,3 )


   ! reciprocal lattice vectors
   write(iunit,'(3f12.6)') ( cell%b1(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%b2(i), i=1,3 )
   write(iunit,'(3f12.6)') ( cell%b3(i), i=1,3 )
   


   write(iunit,*) nspin, nk, nb
   
   ! number of planewaves at each k point 
   ! update: the number of planewaves at each k point is now the same
   do ik = 1, nk
      write(iunit,*) ncoeff
   enddo

   ! k points and weights
   do ik = 1, nk
      write(iunit,'(4f12.6)') ( xk(i,ik), i=1,3 ), wk(ik)
   enddo
   
   ! dense fftsize
   write(iunit,*) fftsize(1), fftsize(2), fftsize(3)
   
   

   close(iunit)

end subroutine



!---------------------------------------------------------------------
! get information of simulation cell
subroutine get_cell_info( c )
   
   type(cell_info), intent(inout) :: c
   integer :: ierr
   
   call qexml_read_cell( ALAT=c%alat, &
                        A1=c%a1(:), A2=c%a2(:), A3=c%a3(:), &
                        B1=c%b1(:), B2=c%b2(:), B3=c%b3(:), &
                        A_UNITS=c%aunit, B_UNITS=c%bunit, IERR=ierr)
   
end subroutine




!-------------------------------------------------------------------
subroutine read_eig_occ( nspin, nk, nb, eigen, occp )

   integer, intent(in) :: nb, nk, nspin
   real(dp), intent(inout) :: eigen(nb,nk,nspin), occp(nb,nk,nspin)

   ! work variables
   real(dp), allocatable :: eigtmp(:,:), occtmp(:,:)
   integer :: nkstot
   integer :: iks
   integer :: ierr
   logical :: lsda, lkpoint_dir
   character(len=256) :: filename
   
   nkstot = nspin * nk
   
   if ( nspin==1 ) then
      lsda = .false.
   elseif( nspin==2 ) then
      lsda = .true.
   endif

   lkpoint_dir = .true.
   ! FIXME: I am not sure what this exactly is.
   ! if this is false, filename must be given.
   ! It seems like if this is true, qexml searches "DATAFILE" from data-file.xml
   ! If it is false, qexml searches "DATA_EIG" from the file whose name was given by "filename"

   filename = 'data-file.xml'
   ! FIXME: so, if I set lkpoint_dir = .true., filename does not matter
   ! but FILENAME is not an optional argument, so we need this.
   
   allocate(eigtmp(nb,nkstot),occtmp(nb,nkstot))
   
   call qexml_read_bands_pw( NUM_K_POINTS=nk, NBND=nb, NKSTOT=nkstot, LSDA=lsda, LKPOINT_DIR=lkpoint_dir,&
       FILENAME=filename, ET=eigtmp, WG=occtmp, IERR=ierr)

   iks = 0
   do ispin = 1, nspin
      do ik = 1, nk
         iks = iks+1
         eigen(:,ik,ispin) = eigtmp(:,iks)
         occp(:,ik,ispin) = occtmp(:,iks)
      enddo
   enddo

end subroutine

!---------------------------------------------------------------------
! get number of plane waves at each k points
Subroutine read_pw_index( nk, ipwk )

   integer, intent(in) :: nk
   integer :: npwk, ik, ierr

   integer, intent(inout) :: ipwk(nk)

   do ik = 1, nk
      call qexml_read_gk( ik, NPWK=npwk, IERR=ierr )
      ipwk( ik ) = npwk
   enddo

end subroutine



!------------------------------------------------------------------------------
!  we need to save rho for GPP calculations
subroutine write_rho( nr, rho )
  
   integer, intent(in) :: nr(3)
   real(dp), intent(inout) :: rho(nr(1),nr(2),nr(3))
   integer,parameter :: iu = 31
   integer :: i,j,k
  
   call qexml_read_rho( RHO=rho, IERR=ierr )

   ! for FORTRAN CODE
  !open(iu,file='rho.dat',form='unformatted',status='unknown')

  !write(iu) nr(1:3)   
  !do k = 1, nr(3)
  !   write(iu) (rho(1:nr(1),j,k),j=1,nr(2))
  !enddo
  !
  !close(iu)


   ! for C++ CODE
   print*, "Rho is printed to a file. This rho is for C++ serial code"
   open(iu,file='rho.dat',form='formatted',status='unknown')
   write(iu,*) nr(1:3)
   do i = 1, nr(1)
      do j = 1, nr(2)
         write(iu,*) rho(i,j,1:nr(3))
      enddo
   enddo
   
 end subroutine write_rho





!-----------------------------------------------------------------------------
! find new ecutwfc. It will be slightly bigger than original ecutwfc
subroutine adjust_ecutwfc( cell, ecutwfc, nk, xk )

  ! cell information (i.e., lattice vectors and reciprocal lattice vectors)
  type(cell_info) :: cell
  ! wavefunction cutoff
  real(dp), intent(inout) :: ecutwfc
  ! number of k points
  integer, intent(in) :: nk
  ! k points in 2pi/alat inverse bohr unit
  real(dp), dimension(3,nk) ,intent(in) :: xk
  
  ! size of the k vector
  real(dp) :: trykmax, kmax
  ! define PI
  real(dp), parameter :: PI = 3.14159265359
  integer :: ik
  

  ! 1. find the largest k vector in inverse bohr unit
  ! initializie
  kmax = 0
  do ik = 1, nk
     trykmax = xk(1,ik)**2 + xk(2,ik)**2 + xk(3,ik)**2
     kmax = sqrt( kmax )
     if ( trykmax .ge. kmax ) then
        kmax = trykmax
     endif
  enddo

  print*, 'original ecutwfc: (Hartree)', ecutwfc
  ! new cutwfc
  ecutwfc = 0.5 * (sqrt(2*ecutwfc) + kmax*2*PI/cell%alat ) **2
  print*, 'new ecutwfc: (Hartree)', ecutwfc

end subroutine adjust_ecutwfc





!-----------------------------------------------------------------------------
! find new ecutwfc. It will be slightly bigger than original ecutwfc
  subroutine set_gkidx_allkpt( npwk, idxgk, ngm, master_gv, gkidx_allkpt )
     ! number of plane waves at this k point
     integer, intent(in) :: npwk
     ! connect to master_gv list
     integer, intent(in) :: idxgk(npwk)
     ! number of plane waves in the master g vector list
     integer, intent(in) :: ngm
     ! master g vector list
     integer, intent(in) :: master_gv(3,ngm)
     ! gkidx_allkpt. it keeps updating
     integer, intent(inout) :: gkidx_allkpt(ngm)

     integer :: i, this_gidx

     ! initially, gkidx_allkpt(:) is ZERO
     ! loop over npwk
     do i = 1, npwk
        this_gidx = idxgk(i)
        if( gkidx_allkpt( this_gidx ) .ne. 1 ) then
           gkidx_allkpt( this_gidx ) = 1
        endif
     enddo

     print*, 'Total g vector list up to this point:', sum( gkidx_allkpt ) 

  end subroutine set_gkidx_allkpt



!-----------------------------------------------------------------------------
! find universial g list
! kmax_cp is updated throughout spin/kpt loop
subroutine get_kmax( ngm, master_gv, npwk, idxgk, kmax_cp)
  ! number of g in the master gvec list
  integer, intent(in) :: ngm
  ! master gvector list
  integer, intent(in) :: master_gv(3,ngm)
  ! number of planewave at this k point
  integer, intent(in) :: npwk
  ! mapping to the master_gv
  integer, intent(in) :: idxgk(npwk)
  ! kmax_cp is used to count and set new g list
  integer, intent(inout) :: kmax_cp(3)

  integer :: ipw, i, tmpkmax(3)

  do ipw = 1, npwk
     do i = 1, 3
        tmpkmax(i) = master_gv(i,idxgk(ipw))
        ! if tmpkmax(i) is negative, make it positive
        if (tmpkmax(i) .lt. 0) then
           tmpkmax(i) = -1 * tmpkmax(i)
        endif
        ! compare if tmpkmax(i) is bigger than kmax_cp(i)
        if ( tmpkmax(i) .ge. kmax_cp(i) ) then
           ! if bigger, then we set new kmax_cp
           kmax_cp(i) = tmpkmax(i)
        endif
     enddo
  enddo
  
end subroutine get_kmax


!-----------------------------------------------------------------------------
! Below subroutines come from OpenAtom program
! It is actually NOT called inside of the main program
! But I leave it here in case we may modify this converter later
! Okay, now I'm modifying the converter, 
! so below routines are now called in the main program (Nov. 2017)
!-----------------------------------------------------------------------------
subroutine radixme( kmax_cp, fftsize )

   integer, intent(inout) :: kmax_cp(3)
   integer, intent(inout) :: fftsize(3) 
   ! work variables
   integer :: i1, i2, i3
   integer :: i, k1, k2, k3
   integer :: kk1, kk2, kk3
   integer :: nrad
   integer :: krad(0:180)
   integer :: nrad_in = 180
   integer :: n(3)

   call set_fftsizes(nrad_in, nrad, krad )

   kk1 = 2*(2*kmax_cp(1) + 1);
   kk2 = 2*(2*kmax_cp(2) + 1);
   kk3 = 2*(2*kmax_cp(3) + 1);

   i1 = 0
   do i=nrad,1,-1
      if ( krad(i) .ge. kk1 ) then
         k1 = krad(i)
         i1 = i
      endif
   enddo
   if ( i1 == 0 ) then 
      print*, 'Bad Radix' 
      stop
   endif
   
   i2 = 0
   do i=nrad,1,-1
      if ( krad(i) .ge. kk2 ) then
         k2 = krad(i)
         i2 = i
      endif
   enddo
   if ( i2 == 0 ) then
      print*, 'Bad Radix'
      stop
   endif

   i3 = 0
   do i=nrad,1,-1
      if ( krad(i) .ge. kk3 ) then
         k3 = krad(i)
         i3 = i
      endif
   enddo
   if ( i3 == 0 ) then
      print*, 'Bad Radix'
      stop
   endif

   print*, 'Radix data : k1, kmax1, kk1', k1 ,kmax_cp(1), kk1
   print*, 'Radix data : k2, kmax2, kk2', k2, kmax_cp(2) ,kk2
   print*, 'Radix data : k3, kmax3, kk3', k3, kmax_cp(3), kk3
  

   n(1) = k1 
   n(2) = k2
   n(3) = k3
   
   fftsize(1:3) = n(1:3)


end subroutine radixme


subroutine set_fftsizes( nrad_in, nrad, krad )

   integer, intent(in) :: nrad_in
   integer, intent(inout) :: nrad
   integer, intent(inout) :: krad(0:180)
   !integer, intent(in) :: fft_opt
  
   integer :: nrad_tmp = 179
   
   nrad = nrad_tmp

   if( nrad_in < nrad_tmp) then
     print*,'@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@'
     print*,'Internal Error in hardcoded radix size array.'
     print*,'@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@' 
   endif

   krad(1)   = 4;
   krad(2)   = 6;
   krad(3)   = 8;
   krad(4)   = 10;
   krad(5)   = 12;
   krad(6)   = 14;
   krad(7)   = 16;
   krad(8)   = 18;
   krad(9)   = 20;
   krad(10)  = 22;
   krad(11)  = 24;
   krad(12)  = 28;
   krad(13)  = 30;
   krad(14)  = 32;
   krad(15)  = 36;
   krad(16)  = 40;
   krad(17)  = 42;
   krad(18)  = 44;
   krad(19)  = 48;
   krad(20)  = 56;
   krad(21)  = 60;
   krad(22)  = 64;
   krad(23)  = 66;
   krad(24)  = 70;
   krad(25)  = 72;
   krad(26)  = 80;
   krad(27)  = 84;
   krad(28)  = 88;
   krad(29)  = 90;
   krad(30)  = 96;
   krad(31)  = 112;
   ! if(fft_opt==0){krad[31]  = 100;}
   krad(32)  = 112;
   krad(33)  = 120;
   krad(34)  = 128;
   krad(35)  = 128;
   krad(36)  = 132;
   krad(37)  = 140;
   krad(38)  = 144;
   krad(39)  = 154;
   krad(40)  = 160;
   krad(41)  = 168;
   krad(42)  = 176;
   krad(43)  = 180;
   krad(44)  = 192;
   krad(45)  = 198;
   krad(46)  = 210;
   krad(47)  = 220;
   krad(48)  = 224;
   krad(49)  = 240;
   krad(50)  = 252;
   krad(51)  = 256;
   krad(52)  = 264;
   krad(53)  = 280;
   krad(54)  = 288;
   krad(55)  = 308;
   krad(56)  = 320;
   krad(57)  = 330;
   krad(58)  = 336;
   krad(59)  = 352;
   krad(60)  = 360;
   krad(61)  = 384;
   krad(62)  = 396;
   krad(63)  = 420;
   krad(64)  = 440;
   krad(65)  = 448;
   krad(66)  = 462;
   krad(67)  = 480;
   krad(68)  = 504;
   krad(69)  = 512;
   krad(70)  = 528;
   krad(71)  = 560;
   krad(72)  = 576;
   krad(73)  = 616;
   krad(74)  = 630;
   krad(75)  = 640;
   krad(76)  = 660;
   krad(77)  = 672;
   krad(78)  = 704;
   krad(79)  = 720;
   krad(80)  = 768;
   krad(81)  = 770;
   krad(82)  = 792;
   krad(83)  = 840;
   krad(84)  = 880;
   krad(85)  = 896;
   krad(86)  = 924;
   krad(87)  = 960;
   krad(88)  = 990;
   krad(89)  = 1008;
   krad(90)  = 1024;
   krad(91)  = 1056;
   krad(92)  = 1120;
   krad(93)  = 1152;
   krad(94)  = 1232;
   krad(95)  = 1260;
   krad(96)  = 1280;
   krad(97)  = 1320;
   krad(98)  = 1344;
   krad(99)  = 1386;
   krad(100) = 1408;
   krad(101) = 1440;
   krad(102) = 1536;
   krad(103) = 1540;
   krad(104) = 1584;
   krad(105) = 1680;
   krad(106) = 1760;
   krad(107) = 1792;
   krad(108) = 1848;
   krad(109) = 1920;
   krad(110) = 1980;
   krad(111) = 2016;
   krad(112) = 2048;
   krad(113) = 2112;
   krad(114) = 2240;
   krad(115) = 2304;
   krad(116) = 2310;
   krad(117) = 2464;
   krad(118) = 2520;
   krad(119) = 2560;
   krad(120) = 2640;
   krad(121) = 2688;
   krad(122) = 2772;
   krad(123) = 2816;
   krad(124) = 2880;
   krad(125) = 3072;
   krad(126) = 3080;
   krad(127) = 3168;
   krad(128) = 3360;
   krad(129) = 3520;
   krad(130) = 3584;
   krad(131) = 3696;
   krad(132) = 3840;
   krad(133) = 3960;
   krad(134) = 4032;
   krad(135) = 4096;
   krad(136) = 4224;
   krad(137) = 4480;
   krad(138) = 4608;
   krad(139) = 4620;
   krad(140) = 4928;
   krad(141) = 5040;
   krad(142) = 5120;
   krad(143) = 5280;
   krad(144) = 5376;
   krad(145) = 5544;
   krad(146) = 5632;
   krad(147) = 5760;
   krad(148) = 6144;
   krad(149) = 6160;
   krad(150) = 6336;
   krad(151) = 6720;
   krad(152) = 6930;
   krad(153) = 7040;
   krad(154) = 7168;
   krad(155) = 7392;
   krad(156) = 7680;
   krad(157) = 7920;
   krad(158) = 8064;
   krad(159) = 8192;
   krad(160) = 8448;
   krad(161) = 8960;
   krad(162) = 9216;
   krad(163) = 9240;
   krad(164) = 9856;
   krad(165) = 10080;
   krad(166) = 10240;
   krad(167) = 10560;
   krad(168) = 10752;
   krad(169) = 11088;
   krad(170) = 11264;
   krad(171) = 11520;
   krad(172) = 12288;
   krad(173) = 12320;
   krad(174) = 12672;
   krad(175) = 13440;
   krad(176) = 13860;
   krad(177) = 14080;
   krad(178) = 14336;
   krad(179) = 14784;   
 
end subroutine set_fftsizes



!------------------------------------------------------------------------------
!  count number of g vectors for doublepack (only half sphere is saved)
 subroutine countkvec3d_sm( kmax_cp, ng , ecut, cell )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(inout) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell

   integer :: i1, i2, i3
   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: xk, yk, zk, tryme
   real(dp) :: factor
   integer :: i, kamax, kbmax, kcmax, kbmin, kcmin, ka, kb, kc, icount
   real(dp) :: aka, akb, akc, g, gmin, gmax

   ! parameter set-up
   factor = 2*PI/cell%alat
   ! initialization
   gmin = 10000000
   gmax = 0
   icount = 0

! FIXME(?) Below doesn't work for non-cubic cell (like diamond structure), so I'll skip this for the moment
!  ! starts from b1 direction
!  i1 = kmax_cp(1)

!  ! expands to only b1 direction
!  do i = 1, i1
!     xk = i * cell%b1(1) * factor
!     yk = i * cell%b1(2) * factor
!     zk = i * cell%b1(3) * factor
!     tryme = (xk*xk + yk*yk + zk*zk)*0.5
!     print*, i, tryme
!     if (tryme .gt. ecut) exit
!  enddo

!  kamax = i-1  ! kamax =  maximum b1 
!  i1 = kamax

!  ! ka goes through 0 to the maximum b1  half sphere includes only ka >= 0
!  ! b1, b2, b3 are real reciprocal lattice vectors
!  do ka = 0, i1
!     aka = dble( ka )
!     kbmin = -kmax_cp(2)
!     if (ka==0) kbmin = 0
!     
!     ! 1. find the real minimum for b2 (kbmin) at this "ka(=aka)"
!     do i = kbmin, 0 ! kbmin = -kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk )*0.5
!        if (tryme .le. ecut ) exit
!     enddo
!     kbmin = i ! this is the minimum for b2

!     ! 2. find the real maximum for b2 (kbmax) at this "ka(=aka)"
!     i2 = kmax_cp(2)
!     do i = 1, i2  ! 1 to kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!        if (tryme .gt. ecut ) exit
!     enddo
!     kbmax = i - 1

!     ! then loop over kbmin =< kb =< kbmax
!     i2 = kbmax
!     do kb = kbmin, i2
!        akb = dble( i )
!        kcmin = -kmax_cp(3)
!        if ( ka==0 .and. kb==0 ) kcmin = 1
!        ! 3. Find kcmin at this "ka(=aka)" and "kb(=akb)"
!        do i = kcmin, 0
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .le. ecut) exit
!        enddo
!        kcmin = i

!        ! 4. Fine kcmax at this "ka(=aka)" and "kb(=akb)"
!        i3 = kmax_cp(3)
!        do i = 1, i3
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .gt. ecut) exit
!        enddo
!        kcmax = i - 1

!        i3 = kcmax
!        ! then loop over kcmin =< kc =< kcmax
!        do kc = kcmin, i3
!           akc = dble( kc )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           g = sqrt( xk*xk + yk*yk + zk*zk )
!           if( gmin .ge. g ) gmin = g
!           if( gmax .le. g ) gmax = g
!           icount = icount + 1
!        enddo! kc loop

!     enddo! kb loop
!  enddo! ka loop
!  
!  ! total number of g vectors except g = 0
!  ng = icount
!  ! add one to include g = 0
!  ng = ng + 1

   ! count again to check using kmax_cp
   icount = 0
   do ka = 0, kmax_cp(1)
      aka = dble( ka )
      kbmin = -kmax_cp(2)
      if( ka==0 ) kbmin = 0
      do kb = kbmin, kmax_cp(2)
         akb = dble( kb )
         kcmin = -kmax_cp(3)
         if( ka==0 .and. kb==0 ) kcmin = 1
         do kc = kcmin, kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) icount = icount + 1
         enddo! kc
      enddo! kb
   enddo! ka

   ng = icount
   ng = ng + 1 ! to include g=0

   ! print results
   print*, 'gvector counts: half sphere without g=0 has', icount , 'vectors'
   print*, 'Ecutwfc:', ecut, 'kmax_cp:', kmax_cp(1), kmax_cp(2), kmax_cp(3)
   !print*, 'kamax, kbmax, kbmin, kcmax, kcmin:', kamax, kbmax, kbmin, kcmax, kcmin
   !print*, 'gmax:', gmax, 'gmin:', gmin
            
 end subroutine countkvec3d_sm


 ! set g vectors for when doublepack = .true.
 subroutine setkvec3d_sm_simple( kmax_cp, ng, ecut, cell, glist )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(in) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell
   integer, intent(inout) :: glist(3,ng)

   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: xk, yk, zk, tryme
   real(dp) :: factor
   integer :: i, kamax, kbmax, kcmax, kbmin, kcmin, ka, kb, kc, icount
   real(dp) :: aka, akb, akc
   
   ! parameter set-up
   factor = 2*PI/cell%alat
   icount = 0
   do ka = 0, kmax_cp(1)
      aka = dble( ka )
      kbmin = -kmax_cp(2)
      if( ka==0 ) kbmin = 0
      do kb = kbmin, kmax_cp(2)
         akb = dble( kb )
         kcmin = -kmax_cp(3)
         if( ka==0 .and. kb==0 ) kcmin = 1
         do kc = kcmin, kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) then
               icount = icount + 1
               glist(1,icount) = ka
               glist(2,icount) = kb
               glist(3,icount) = kc
            endif
            
         enddo! kc
      enddo! kb
   enddo! ka
   icount = icount + 1
   if ( ng .ne. icount ) then
      print*, '@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@'
      print*, ' number of g vectors do not match'
      print*, ' ng', ng , 'vs', icount
      print*, '@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@'
      stop
   endif
   ! g=0 is added
   glist(:,icount) = 0
      
 end subroutine setkvec3d_sm_simple
 


 ! count g vectors for the full sphere (doublepack = .false.)
 subroutine countkvec3d_sm_kpt( kmax_cp, ng, ecut, cell )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(inout) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell

   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: factor
   integer :: icount, ka, kb, kc
   real(dp) :: xk, yk, zk, aka, akb, akc, g, gmin, gmax, tryme
   
   factor = 2 * PI / cell%alat
   ! initialization
   gmin = 10000000
   gmax = 0

   icount = 0
   do ka = -kmax_cp(1), kmax_cp(1)
      aka = dble( ka )
      do kb = -kmax_cp(2), kmax_cp(2)
         akb = dble( kb )
         do kc = -kmax_cp(3), kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            g = sqrt( xk*xk + yk*yk + zk*zk )
            if ( g .ge. gmax ) gmax = g
            if ( g .le. gmin ) gmin = g
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) icount = icount + 1
         enddo
      enddo
   enddo

   ! number of k points
   ng = icount
   
   ! printint results
   print*, 'gvector counts: full sphere has', icount ,'vectors'
   print*, 'ecutwfc:', ecut, 'kmax_cp:', kmax_cp(1), kmax_cp(2), kmax_cp(3)
   !print*, 'gmax:', gmax, 'gmin:', gmin

 end subroutine countkvec3d_sm_kpt



 ! set g vectors for the full sphere ( doublepack = .false. )
 subroutine setkvec3d_sm_kpt( kmax_cp, ng, ecut, cell, glist )
   integer, intent(in) :: kmax_cp(3)
   integer, intent(in) :: ng
   real(dp), intent(in) :: ecut
   type(cell_info) :: cell
   integer, intent(inout) :: glist(3,ng)

   real(dp), parameter :: PI = 3.14159265359
   real(dp) :: factor
   integer :: icount, ka, kb, kc
   real(dp) :: xk, yk, zk, aka, akb, akc, tryme
   
   factor = 2 * PI / cell%alat
   icount = 0

   do ka = -kmax_cp(1), kmax_cp(1)
      aka = dble( ka )
      do kb = -kmax_cp(2), kmax_cp(2)
         akb = dble( kb )
         do kc = -kmax_cp(3), kmax_cp(3)
            akc = dble( kc )
            xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
            yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
            zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
            tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
            if( tryme .le. ecut ) then
               icount = icount + 1
               glist(1,icount) = ka
               glist(2,icount) = kb
               glist(3,icount) = kc
            endif
         enddo
      enddo
   enddo

   ! printint errocr
   if ( icount .ne. ng ) then
      print*, '@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@'
      print*, 'incorrect number of small kvectors in full sphere'
      print*, icount, 'vs', ng
      print*, '@@@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@@@'
   endif

 end subroutine setkvec3d_sm_kpt

! we don't use this

!subroutine setkvec3d_sm( kmax_cp, ng , ecut, cell, glist )
!  integer, intent(in) :: kmax_cp(3)
!  integer, intent(in) :: ng
!  real(dp), intent(in) :: ecut
!  type(cell_info) :: cell
!  integer, intent(inout) :: glist(3,ng)

!  integer :: i1, i2, i3
!  real(dp), parameter :: PI = 3.14159265359
!  real(dp) :: xk, yk, zk, tryme
!  real(dp) :: factor
!  integer :: i, kamax, kbmax, kcmax, kbmin, kcmin, ka, kb, kc, icount
!  real(dp) :: aka, akb, akc, g, gmin, gmax

!  ! parameter set-up
!  factor = 2*PI/cell%alat
!  ! initialization
!  gmin = 10000000
!  gmax = 0
!  icount = 0

!  ! starts from b1 direction
!  i1 = kmax_cp(1)

!  ! expands to only b1 direction
!  do i = 1, i1
!     xk = i * cell%b1(1) * factor
!     yk = i * cell%b1(2) * factor
!     zk = i * cell%b1(3) * factor
!     tryme = (xk*xk + yk*yk + zk*zk)*0.5
!     if (tryme .gt. ecut) exit
!  enddo

!  kamax = i-1  ! kamax =  maximum b1 - 1
!  i1 = kamax

!  ! ka goes through 0 to the maximum b1  half sphere includes only ka >= 0
!  ! b1, b2, b3 are real reciprocal lattice vectors
!  do ka = 0, i1
!     aka = dble( ka )
!     kbmin = -kmax_cp(2)
!     if (ka==0) kbmin = 0
!     
!     ! 1. find the real minimum for b2 (kbmin) at this "ka(=aka)"
!     do i = kbmin, 0 ! kbmin = -kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk )*0.5
!        if (tryme .le. ecut ) exit
!     enddo
!     kbmin = i ! this is the minimum for b2

!     ! 2. find the real maximum for b2 (kbmax) at this "ka(=aka)"
!     i2 = kmax_cp(2)
!     do i = 1, i2  ! 1 to kmax_cp[2]
!        akb = dble( i )
!        xk = ( aka * cell%b1(1) + akb * cell%b2(1) ) * factor
!        yk = ( aka * cell%b1(2) + akb * cell%b2(2) ) * factor
!        zk = ( aka * cell%b1(3) + akb * cell%b2(3) ) * factor
!        tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!        if (tryme .gt. ecut ) exit
!     enddo
!     kbmax = i - 1

!     ! then loop over kbmin =< kb =< kbmax
!     i2 = kbmax
!     do kb = kbmin, i2
!        akb = dble( i )
!        kcmin = -kmax_cp(3)
!        if ( ka==0 .and. kb==0 ) kcmin = 1
!        ! 3. Find kcmin at this "ka(=aka)" and "kb(=akb)"
!        do i = kcmin, 0
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .le. ecut) exit
!        enddo
!        kcmin = i

!        ! 4. Fine kcmax at this "ka(=aka)" and "kb(=akb)"
!        i3 = kmax_cp(3)
!        do i = 1, i3
!           akc = dble( i )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           tryme = ( xk*xk + yk*yk + zk*zk ) * 0.5
!           if (tryme .gt. ecut) exit
!        enddo
!        kcmax = i - 1

!        i3 = kcmax
!        ! then loop over kcmin =< kc =< kcmax
!        do kc = kcmin, i3
!           akc = dble( kc )
!           xk = ( aka * cell%b1(1) + akb * cell%b2(1) + akc * cell%b3(1) ) * factor
!           yk = ( aka * cell%b1(2) + akb * cell%b2(2) + akc * cell%b3(2) ) * factor
!           zk = ( aka * cell%b1(3) + akb * cell%b2(3) + akc * cell%b3(3) ) * factor
!           g = sqrt( xk*xk + yk*yk + zk*zk )
!           if( gmin .ge. g ) gmin = g
!           if( gmax .le. g ) gmax = g
!           icount = icount + 1
!           glist(1,icount) = ka
!           glist(2,icount) = kb
!           glist(3,icount) = kc
!        enddo! kc loop

!     enddo! kb loop
!  enddo! ka loop

!  if ( ng-1 .ne. icount ) then
!     print*, "@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@"
!     print*, "Mismatch number of small kvectors:"
!     print*, icount ,'vs', ng
!     print*, "@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@"
!  endif
!  
!  ! total number of g vectors
!  icount = icount + 1
!  glist(:,icount) = 0

!end subroutine setkvec3d_sm

 ! this subroutine maps OA glist to QE glist
 ! this is called at each k point
 
subroutine mapping_glist( doublepack, ngoa, glist_oa, ngm, master_gv, npwk, gidx_qe, idxmap )
   logical, intent(in) :: doublepack
   integer, intent(in) :: ngoa ! number of g vectors for openatom
   integer, intent(in) :: glist_oa(3,ngoa) ! g list of openatom
   integer, intent(in) :: ngm ! number of g vectors in the master_gv list
   integer, intent(in) :: master_gv(3,ngm)
   integer, intent(in) :: npwk
   integer, intent(in) :: gidx_qe(ngm) ! index array that indicates which gvector we include and whatnot
   integer, intent(inout) :: idxmap(ngoa) ! indexing map from oa to qe

   integer :: i, j, k
   integer, allocatable :: glist_qe(:,:)


   allocate (glist_qe(3,npwk))
   do i=1,npwk
      glist_qe(:,i) = master_gv(:,gidx_qe(i))
   enddo


   ! initialize
   idxmap(:) = -1
   do i=1, ngoa
      do j=1, npwk
         if ( glist_oa(1,i) == glist_qe(1,j) .and. glist_oa(2,i) == glist_qe(2,j) .and. &
              glist_oa(3,i) == glist_qe(3,j) ) then
            idxmap(i) = j
            ! finish search
            exit
         endif
      enddo
   enddo

   ! mapping done

   deallocate( glist_qe )
   
end subroutine mapping_glist
 

end program
