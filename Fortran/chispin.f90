!================================================================
! BerkeleyGW utility: chispin.f90
!
! Wrapper for the interacting spin susceptibility chi calculation
! 
! An example input file is given in chispin.inp_example
! Change to chispin.inp for use.
!
! Authors: TBZ, JJL
! Last updated: Dec 2012
!================================================================

program chispin
  !--------------------------------------------------------------
  !     VARIABLES DECLARATIONS
  !--------------------------------------------------------------
  ! .. types ..
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  ! .. scalars ..
  integer, dimension(3) :: S
  integer :: lenS, error, Nq, ii, jj, outunit, debugunit 
  integer :: Nfreq, iq, ind1, ind2, nmtx, ifreq, ncol
  real(DP) :: Vcell,dV, latt_const
  real(DP) :: a, b ! local variables used to read chi0
  ! .. strings ..
  character*256 :: filecd, filegv, filechi0, filenmtx, filechi, filedebug
  character*256, parameter :: inputfile = "chispin.inp"
  ! .. integer arrays ..
  integer, allocatable, dimension(:) :: nmtxdat,indx
  ! .. real arrays ..
  real(DP), dimension(3,3) :: R,inverseR
  real(DP), allocatable, dimension(:) :: G2,Ixc
  real(DP), allocatable, dimension(:,:) :: G,gvecs,gvecsdat
  ! .. complex arrays ..
  complex(DPC), allocatable, dimension(:) :: row_chi0, chi0
  complex(DPC), allocatable, dimension(:,:) :: epsmat, chi, chi0w, chi0f, Ident
  !--------------------------------------------------------------
  !     READING THE INPUT FILE & INITIALIZING VARIABLES
  !--------------------------------------------------------------
  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(inputfile)
  
  open(55,file=trim(inputfile),form='formatted',status='old')
  read(55,'(a)') filenmtx
  read(55,'(a)') filecd
  read(55,'(a)') filegv
  read(55,'(a)') filechi0
  read(55,'(a)') filechi
  read(55,*) (S(ii),ii=1,3)
  read(55,*) latt_const
  read(55,*) ((R(ii,jj),jj=1,3),ii=1,3)
  read(55,*) Nq
  read(55,*) Nfreq
  close(55)
  
  write(6,'(2a)')       "     nmtx    file  = ", trim(filenmtx)
  write(6,'(2a)')       "     cd      file  = ", trim(filecd)
  write(6,'(2a)')       "     gv      file  = ", trim(filegv)
  write(6,'(2a)')       "     chi0    file  = ", trim(filechi0)
  write(6,'(2a)')       "     chi     file  = ", trim(filechi)
  write(6,'(a, 3i6)') "     Fourier grid  = ", S(1:3)
  write(6,'(a, f10.5)')  "     latt_const    = ", latt_const
  write(6,'(a, 3f10.5)')   "     R(1,:)        = ", R(1,1:3)
  write(6,'(a, 3f10.5)')   "     R(2,:)        = ", R(2,1:3)
  write(6,'(a, 3f10.5)')   "     R(3,:)        = ", R(3,1:3)
  write(6,'(a, i6)')    "     Nq            = ", Nq
  write(6,'(a, i6)')    "     Nfreq         = ", Nfreq
  write(6,*)
  ! get inverse lattice tensor
  !-------------- Calling invert -----------
  call invert_real(R, inverseR, size(R, dim=1))
  ! total size of FFTbox and lattice tensor in atomic units
  lenS = S(1)*S(2)*S(3)
  R = R * latt_const
  ! unit cell volume, volume element
  Vcell = R(1,1)*R(2,2)*R(3,3)
  dV    = Vcell/lenS
  !--------------------------------------------------------------
  !     I/O initialization
  !--------------------------------------------------------------
  ! output file for chi
  outunit = 99
  open(unit=outunit,file=filechi,form='formatted',status='replace')
  ! output file for debug output
  !debugunit = 66
  !open(unit=debugunit,file=filedebug,form='formatted',status='replace')
  ! allocate sizes of chi0 matrices (number of Gvecs)
  allocate( nmtxdat(Nq), stat = error)
  open(10,file=filenmtx,form='formatted')
  do ii=1,Nq
     read(10,*) nmtxdat(ii)
  enddo
  close(10)
  ! allocate arrays of reciprocal lattice vectors (G) and their length squared (G2)
  allocate( G(lenS,3) , stat=error)
  allocate( G2(lenS)  , stat=error)
  allocate( chi0(Nfreq*sum(nmtxdat**2)), stat=error)
  allocate( gvecsdat(sum(nmtxdat**2),3) , stat=error)
  allocate( Ixc(lenS) , stat = error)
  !----------------- call setup ----------------
  call setup(inverseR,S,G,G2)
  !--------------------------------------------------------------
  !     READ THE DATAFILES
  !--------------------------------------------------------------
  ! print *, "size of chi0", Nfreq*sum(nmtxdat**2), size(chi0)
  open(11,file=filechi0,form='formatted')
  do ii = 1, Nfreq*sum(nmtxdat**2)
     read(11,*) a, b
     chi0(ii) = cmplx(a, b)/2.0
   enddo
  close(11)

  open(12,file=filegv,form='formatted')
  do ii=1,sum(nmtxdat)
     read(12,*) gvecsdat(ii,1),gvecsdat(ii,2),gvecsdat(ii,3)
  enddo
  close(12)
  ! ---------------- call getIxc ---------------
  call getIxc(lenS,Ixc, filecd)

  !**************************************************************
  !     MAIN BODY OF THE PROGRAM
  !**************************************************************
  !---------------------------------------------
  ! LOOP THROUGH THE q-vectors
  !---------------------------------------------
  do iq = 1,Nq
     print *, "qpoint #    ", iq, "out of", Nq
     ! get the number of matrix elements for a current q-point iq
     nmtx = nmtxdat(iq)
     ! Check if the ind1 is zero and assign it and ind2
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1))
     endif
     ind2 = sum(nmtxdat(1:iq))
  !-----------------------------------------------
  ! Allocate arrays
  !-----------------------------------------------
     allocate(gvecs(nmtx,3), stat=error)
     allocate(indx(nmtx), stat=error)
     allocate( chi0w(Nfreq,nmtx**2), stat=error)
     allocate( chi0f(nmtx,nmtx), stat=error)
     allocate(epsmat (nmtx,nmtx), stat=error)
     allocate(chi (nmtx,nmtx), stat=error)
     allocate(Ident(nmtx,nmtx), stat=error)
     Ident = 0d0
     do ii = 1, nmtx
           Ident(ii,ii) = cmplx(1d0,0d0)
     enddo
     ! allocate row_chi0
     allocate(row_chi0 (lenS), stat=error) 

     ! assigning values to gvecs
     gvecs = gvecsdat( ind1+1:ind2, 1:3)
     
     !--------- call fftbox----------------------
     call get_fftbox(gvecs,S,nmtx,indx)

     ! check if the ind1 is zero and assign it and ind2 again for chi0
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1)**2)
     endif
     ind2 = sum(nmtxdat(1:iq)**2)

     ! get the part of chi0 that corresponds to the current iq
     ! and reshape it into a matrix
     chi0w = reshape(chi0( ind1*Nfreq + 1 : ind2*Nfreq), (/ Nfreq, nmtx**2 /) )

     !-------------------------------------------
     ! LOOP THROUGH THE FREQUENCIES
     !-------------------------------------------
     do ifreq = 1,Nfreq
        print *, "frequency # ", ifreq, "out of", Nfreq
        ! reshaping the chi0w array
        chi0f = reshape (chi0w(ifreq,:), (/ nmtx,nmtx/))

        epsmat = 0d0
        chi = 0d0
        !----------------------------------------
        ! LOOP THROUGH THE G-VECTORS
        !----------------------------------------
        do ncol = 1,nmtx
           ! .. output progress ..
           if ( mod(ncol,50).eq.0) print *, "column #    ", ncol, "out of", nmtx
           ! .. nullify row_chi0 and assign it
           row_chi0 = 0d0
           row_chi0(indx) = chi0f(:,ncol)
           !-------------- Calling cJ -----------
           call cJ( row_chi0, S)
           ! multiply by the interaction
           row_chi0 = row_chi0 * Ixc
           !-------------- Calling cI -----------
           call cI( row_chi0, S)
           ! init epsilon
           epsmat( : ,ncol) = row_chi0(indx)
        enddo ! loop over g-vectors

        ! get the epsilon matrix using identity matrix
        epsmat = Ident - epsmat
        !-------------- Calling invert -----------
        call invert_cmpx(epsmat, chi, nmtx)
        ! multiply by chi0
        chi = matmul(chi, chi0f)
        ! write output into a file
        write(unit=outunit, fmt="(A, I10, A, I10)", iostat=error, advance='YES') &
           "chi for iq # ", iq, " and ifreq # ", ifreq 
        do ii=1, nmtx
          do jj=1,nmtx
               write(unit=outunit, fmt="(2F25.15)", iostat=error, advance='YES') &
                  dble(chi(ii,jj)), aimag(chi(ii,jj))
            if ( error /= 0 ) stop "Write error in file unit ounit"
          enddo
        enddo   

     enddo ! loop over frequencies
     !-----------------------------------------------
     ! Deallocate arrays
     !-----------------------------------------------
     deallocate(row_chi0, stat=error)
     deallocate(epsmat, stat=error) 
     deallocate(chi, stat=error)
     deallocate(Ident, stat=error)
     deallocate(chi0f, stat=error)
     deallocate(chi0w, stat=error)
     deallocate(gvecs, stat=error)
     deallocate(indx, stat=error)

  enddo ! loop over q-vectors
  
  !-----------------------------------------------
  ! Deallocate preliminary arrays
  !-----------------------------------------------
  deallocate( nmtxdat, stat = error)
  deallocate( G, stat=error)
  deallocate( G2, stat=error)
  deallocate( chi0, stat=error)
  deallocate( gvecsdat, stat=error)
  deallocate( Ixc, stat = error)
    
  ! close the output file
  close(outunit)
  
  print *, "PROGRAM FINISHED"

end program chispin
