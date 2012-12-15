!================================================================
!
! BerkeleyGW utility:
!
! chispin
!
! wrapper for the interacting spin susceptibility chi calculation
!
! Authors: TBZ, JJL
!
! Last updated: Dec 2012
!================================================================

program chispin
  !--------------------------------------------------------------
  !     VARIABLES DECLARATIONS
  !--------------------------------------------------------------
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3) :: S
  integer :: lenS, error, Nq, ii, jj, outunit, Nfreq, iq, ind1, ind2, nmtx, ifreq, ncol
  integer, allocatable, dimension(:) :: nmtxdat,indx

  real(DP), dimension(3,3) :: R,inverseR
  real(DP), allocatable, dimension(:) :: G2,Ixc
  real(DP), allocatable, dimension(:,:) :: G,gvecs,gvecsdat, Identity
  complex(DPC), allocatable, dimension(:) :: row_chi0, chi0
  complex(DPC), allocatable, dimension(:,:) :: epsmat, chi, chi0w, chi0f
 
  real(DP) :: Vcell,dV, latt_const
  real(DP) :: a, b ! local variables used to read chi0
  !--------------------------------------------------------------
  !     VARIABLES INITIALIZATION
  !--------------------------------------------------------------
  ! Fourier grid size
  S(1) = 60
  S(2) = 60
  S(3) = 360
  lenS = S(1)*S(2)*S(3)

  ! lattice tensor (matrix of lattice vectors)
  latt_const = 6.97 ! atomic units
  R(1,1) =  1.0;  R(1,2) =  0.0;  R(1,3) =  0.0
  R(2,1) =  0.0;  R(2,2) =  1.0;  R(2,3) =  0.0
  R(3,1) =  0.0;  R(2,3) =  0.0;  R(3,3) = 5.86
  R = R * latt_const

  ! inverse lattice tensor
  inverseR(1,1) = 1./R(1,1); inverseR(1,2) = 0.0;       inverseR(1,3) = 0.0
  inverseR(2,1) = 0.0;       inverseR(2,2) = 1./R(2,2); inverseR(2,3) = 0.0
  inverseR(3,1) = 0.0;       inverseR(3,2) = 0.0;       inverseR(3,3) = 1./R(3,3)

  ! unit cell volume, volume element
  Vcell = R(1,1)*R(2,2)*R(3,3)
  dV    = Vcell/lenS
  
  ! number of qpoints
  Nq = 1

  ! number of frequencies
  Nfreq = 18
 !--------------------------------------------------------------
 !     I/O initialization
 !--------------------------------------------------------------
  ! Output file for chi
  outunit = 99
  open(unit=outunit,file='chi.dat',form='formatted',status='replace')

  ! allocate sizes of chi0 matrices (number of Gvecs)
  allocate( nmtxdat(Nq), stat = error)
  open(10,file="data/nmtx.dat",form='formatted')
  do ii=1,Nq
     read(10,*) nmtxdat(ii)
  enddo
  close(10)

  ! allocate arrays of reciprocal lattice vectors (G) and their length squared (G2)
  allocate( G(lenS,3) , stat=error)
  allocate( G2(lenS)  , stat=error)
  ! allocate space for chi0
  allocate( chi0(Nfreq*sum(nmtxdat**2)), stat=error)
  ! allocate space for gvecs
  allocate( gvecsdat(sum(nmtxdat**2),3) , stat=error)
  ! allocate Ixc 
  allocate( Ixc(lenS) , stat = error)
  
  !----------------- call setup ----------------
  call setup(inverseR,S,G,G2)
 !--------------------------------------------------------------
 !     READ THE DATAFILES
 !--------------------------------------------------------------
  print *, "size of chi0", Nfreq*sum(nmtxdat**2), size(chi0)
  
  open(11,file="data/chimat.dat",form='formatted')
  do ii = 1, Nfreq*sum(nmtxdat**2)
     read(11,*) a, b
     chi0(ii) = cmplx(a, b)
     !print *, "chi0", ii, chi0(ii)
  enddo
  close(11)

  open(12,file="data/gvecs.dat",form='formatted')
  do ii=1,sum(nmtxdat)
     read(12,*) gvecsdat(ii,1),gvecsdat(ii,2),gvecsdat(ii,3)
  enddo
  close(12)

  ! ---------------- call getIxc ---------------
  call getIxc(lenS,Ixc)

  !--------------------------------------------------------------
  !--------------------------------------------------------------
  !     MAIN BODY OF THE PROGRAM
  !--------------------------------------------------------------
  !--------------------------------------------------------------
  
  !---------------------------------------------
  ! LOOP THROUGH THE q-vectors
  !---------------------------------------------
  do iq = 1,Nq
     
     print *, "doing qpoint #", iq, "out of total", Nq

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
     print *, "allocating arrays"

     allocate(gvecs(nmtx,3), stat=error)
     allocate(indx(nmtx), stat=error)
     ! allocate chi0w and individual arrays for real and imag parts     
     allocate( chi0w(Nfreq,nmtx**2), stat=error)
     ! allocate chi0f individual arrays for real and imag parts     
     allocate( chi0f(nmtx,nmtx), stat=error)
     ! allocate epsmat, chi, Identity
     allocate(epsmat (nmtx,nmtx), stat=error)
     allocate(chi (nmtx,nmtx), stat=error)
     ! Create identity matrix
     allocate(Identity(nmtx,nmtx), stat=error)
     Identity = 0.0
     do ii = 1, nmtx; Identity(ii,ii) = 1.0; enddo
     ! allocate row_chi0
     allocate(row_chi0 (lenS), stat=error) 

     ! assigning values to gvecs
     gvecs = gvecsdat( ind1+1:ind2, 1:3)
     
     !--------- call fftbox----------------------
     print *, "calling fftbox"
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
!     print *, "reshaping chi0 ", "ind1 =", ind1, "Nfreq=", Nfreq, "ind2=", ind2
     chi0w = reshape (chi0( ind1*Nfreq + 1 : ind2*Nfreq), (/ Nfreq, nmtx**2 /) )

     !-------------------------------------------
     ! LOOP THROUGH THE FREQUENCIES
     !-------------------------------------------
     do ifreq = 1,Nfreq
        
        print *, "doing frequency #", ifreq, "out of total", Nfreq

        ! reshaping the chi0w array
        chi0f = reshape (chi0w(ifreq,:), (/ nmtx,nmtx/))
        
        epsmat = 0d0
        chi = 0d0
        
        !----------------------------------------
        ! LOOP THROUGH THE G-VECTORS
        !----------------------------------------
        do ncol = 1,nmtx
           
           if ( mod(ncol,50).eq.0) then
             print *, "doing column #", ncol, "out of total", nmtx
           endif
           
           row_chi0 = 0d0

           row_chi0(indx) = chi0f(:,ncol)
           
           !-------------- Calling cJ -----------
  !         print *, "calling cJ"
           call cJ( row_chi0, S)
  !         print *, "cJ done"
           
           ! multiply by the interaction
           row_chi0 = row_chi0 * Ixc
           
           !-------------- Calling cI -----------
 !          print *, "calling cI"
           call cI( row_chi0, S)
!           print *, "cI done"

           ! init epsilon
           epsmat( : ,ncol) = row_chi0(indx)

           if ( mod(ncol,50).eq.0) then
             print *, "DONE! column #", ncol, "out of total", nmtx    
           endif
           
        enddo ! loop over g-vectors

        ! get the epsilon matrix using identity matrix
        epsmat = Identity - epsmat
        
        ! invert epsilon
        !-------------- Calling invert -----------
        call invert(epsmat, chi, nmtx, error)

        print *, "epsilon matrix inverse DONE!"

        ! multiply by chi0
        chi = chi *  chi0f
        
        ! write output into file
        write(unit=outunit, fmt="(A, I10, A, I10)", iostat=error, advance='YES') "chi for iq # ", iq, " and ifreq # ", ifreq 
        do ii=1, nmtx
          do jj=1,nmtx
            write(unit=outunit, fmt="(2F25.15)", iostat=error, advance='YES') dble(chi(ii,jj)), aimag(chi(ii,jj))
            if ( error /= 0 ) stop "Write error in file unit ounit"
          enddo
        enddo   

        print *, "DONE! frequency #", ifreq, "out of total", Nfreq    

     enddo ! loop over frequencies

  !-----------------------------------------------
  ! Deallocate arrays
  !-----------------------------------------------

     ! deallocate row_chi0
     deallocate(row_chi0, stat=error)

     ! deallocate epsmat, chi
     deallocate(epsmat, stat=error) 
     deallocate(chi, stat=error)
     deallocate(Identity, stat=error)
     ! deallocate chi0f individual arrays for real and imag parts     
     deallocate(chi0f, stat=error)
     ! deallocate chi0w and individual arrays for real and imag parts     
     deallocate(chi0w, stat=error)
 
     ! deallocating gvecs and indx
     deallocate(gvecs, stat=error)
     deallocate(indx, stat=error)

     print *, "DONE! q-vector #", iq, "out of total", Nq    

  enddo ! loop over q-vectors
  
  !-----------------------------------------------
  ! Deallocate preliminary arrays
  !-----------------------------------------------

  ! deallocate sizes of chi0 matrices (number of Gvecs)
  deallocate( nmtxdat, stat = error)
  ! deallocate arrays of reciprocal lattice vectors (G) and their length squared (G2)
  deallocate( G, stat=error)
  deallocate( G2, stat=error)
  ! deallocate space for chi0
  deallocate( chi0, stat=error)
  ! deallocate space for gvecs
  deallocate( gvecsdat, stat=error)
  ! deallocatede Ixc 
  deallocate( Ixc, stat = error)
    
  ! close the output file
  close(outunit)

end program chispin
