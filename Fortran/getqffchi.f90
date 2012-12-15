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
!================================================================

program chispin


!--------------------------------------------------------------
!     VARIABLES DECLARATIONS
!--------------------------------------------------------------
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3) :: S
  integer :: lenS, error, Nq, ii, jj, Nfreq, iq, ind1, ind2, nmtx, ifreq, ncol
  integer, allocatable, dimension(:) :: nmtxdat,indx

  real(DP), dimension(3,3) :: R,inverseR, a, b
  real(DP), allocatable, dimension(:) :: G2,Ixc
  real(DP), allocatable, dimension(:,:) :: G,chi0,gvecs,gvecsdat, chi0w, chi0w_real, chi0w_cmpx, chi0f_real, chi0f_cmpx, Identity
  complex(DPC), allocatable, dimension(:) :: row_chi0
  complex(DPC), allocatable, dimension(:,:) :: epsmat, chi
 
  real(DP) :: Vcell,dV, latt_const


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
!     ALLOCATE ARRAYS AND READ THE DATAFILES
!--------------------------------------------------------------

  open(unit=99,file='chi.dat',form='formatted',status='replace')

  ! sizes of chi0 matrices (number of Gvecs)
  allocate( nmtxdat(Nq), stat = error)

  ! get reciprocal lattice vectors (G) and their length squared (G2)
  allocate( G(lenS,3) , stat=error)
  allocate( G2(lenS)  , stat=error)
  
  !----------------- call setup ----------------
  call setup(inverseR,S,G,G2)
  !---------------------------------------------
  
  open(10,file="nmtx.dat",form='formatted')
  do ii=1,Nq
     read(10,*) nmtxdat(ii)
  enddo
  close(10)

  ! allocate chi0
  allocate( chi0(Nfreq*sum(nmtxdat**2),2), stat=error)
  print *, "size of chi0", Nfreq*sum(nmtxdat**2)
  
  ! allocate space for gvecs
  allocate( gvecsdat(sum(nmtxdat**2),3) , stat=error)

  ! open chi0 file
  open(11,file="chimat.dat",form='formatted')
  do ii=1,Nfreq*sum(nmtxdat**2)
     read(11,*) chi0(ii,1),chi0(ii,2)
     !print *, "chi0", ii,chi0(ii,1), chi0(ii,2)
  enddo
  close(11)

  open(12,file="gvecs.dat",form='formatted')
  do ii=1,sum(nmtxdat)
     read(12,*) gvecsdat(ii,1),gvecsdat(ii,2),gvecsdat(ii,3)
  enddo
  close(12)

  ! allocate Ixc 
  allocate( Ixc(lenS) , stat = error)

  ! ---------------- call getIxc ---------------
  call getIxc(lenS,Ixc)
  ! --------------------------------------------
  

!--------------------------------------------------------------
!
!     MAIN BODY OF THE PROGRAM
!
!--------------------------------------------------------------

  !---------------------------------------------
  ! LOOP THROUGH THE q-vectors
  !---------------------------------------------
  do iq = 1,Nq
     
     print *, "doing qpoint #", iq, "out of total", Nq

     !do ii=1,Nfreq*sum(nmtxdat**2)
     !    read(11,*) a, b
     !    chi0(ii
     !enddo

     ! get the number of matrix elements for a current q-point iq
     nmtx = nmtxdat(iq)

     ! Check if the ind1 is zero and assign it and ind2
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1))
     endif
     
     ind2 = sum(nmtxdat(1:iq))
     
     ! allocating gvecs and indx
     print *, "allocating gvecs and indx"
     allocate(gvecs(nmtx,3), stat=error)
     allocate(indx(nmtx), stat=error)

     ! assigning values to gvecs
     gvecs = gvecsdat( ind1+1:ind2, 1:3)
     
     !--------- call fftbox------------
     print *, "calling fftbox"
     call get_fftbox(gvecs,S,nmtx,indx)
     !--------------------------------
     
     ! check if the ind1 is zero and assign it and ind2 again for chi0
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1)**2)
     endif
     ind2 = sum(nmtxdat(1:iq)**2)

     ! allocate chi0w and individual arrays for real and imag parts     
     allocate( chi0w(Nfreq*nmtx**2,2), stat=error)
     allocate( chi0w_real(Nfreq*nmtx**2,1), stat=error)
     allocate( chi0w_cmpx(Nfreq*nmtx**2,1), stat=error)

     ! get the part of chi0 that corresponds to the current iq
     ! and reshape it into a matrix
     ! TODO: understand how the complex algebra works and get rid of 2 arrays or real am cmpx
     print *, "reshaping chi0 ", "ind1 =", ind1, "Nfreq=", Nfreq, "ind2=", ind2
     chi0w = chi0( ind1*Nfreq + 1 : ind2*Nfreq, : )
     chi0w_real = reshape (chi0w(:,1), (/ Nfreq, nmtx**2 /))
     chi0w_cmpx = reshape (chi0w(:,2), (/ Nfreq, nmtx**2 /))

     !-------------------------------------------
     ! LOOP THROUGH THE FREQUENCIES
     !-------------------------------------------
     do ifreq = 1,Nfreq
        
        print *, "doing frequency #", ifreq, "out of total", Nfreq
        print *, "reshaping the arrays chi0w. nmtx= ", nmtx
        ! allocate chi0f individual arrays for real and imag parts     
        allocate( chi0f_real(nmtx,nmtx), stat=error)
        print *, "reshaping the arrays chi0w1"
        allocate( chi0f_cmpx(nmtx,nmtx), stat=error)

        print *, "reshaping the arrays chi0w2"

        ! reshaping the chi0w array
        chi0f_real = reshape (chi0w_real(ifreq,:), (/ nmtx,nmtx/))
        chi0f_cmpx = reshape (chi0w_cmpx(ifreq,:), (/ nmtx,nmtx/))
        
        ! allocate epsmat, chi and nullify also
        allocate(epsmat (nmtx,nmtx), stat=error)
        epsmat = 0d0
        allocate(chi (nmtx,nmtx), stat=error)
        chi = 0d0
        
        !-------------------------------------------
        ! LOOP THROUGH THE G-VECTORS
        !-------------------------------------------
        do ncol = 1,nmtx
           
           print *, "doing column #", ncol, "out of total", nmtx
           
           ! allocate row_chi0 and nullify it also
           allocate(row_chi0 (lenS), stat=error) 
           row_chi0 = 0d0

           row_chi0(indx) = cmplx ( chi0f_real(:,ncol) , chi0f_cmpx(:,ncol) )
           
           !-------------- Calling cJ ----------------------
           print *, "calling cJ"
           call cJ( row_chi0, S)
           print *, "cJ done"
           !------------------------------------------------
           
           ! multiply by the interaction
           row_chi0 = row_chi0 * Ixc
           
           !-------------- Calling cI ----------------------
           print *, "calling cI"
           call cI( row_chi0, S)
           print *, "cI done"
           !------------------------------------------------

           !
           epsmat( : ,ncol) = row_chi0(indx)
         
           ! deallocate row_chi0
           deallocate(row_chi0, stat=error)
           
        enddo ! loop over g-vectors
      
        print *, "DONE! column #", ncol, "out of total", nmtx    

        ! Create identity matrix
        allocate(Identity(nmtx,nmtx), stat=error)
        Identity = 0.0
        do ii = 1, nmtx; Identity(ii,ii) = 1.0; enddo

        print *, "Identity matrix DONE!"

        ! get the epsilon matrix using identity matrix
        epsmat = Identity - epsmat
        
        ! invert epsilon and multiply by chi0
        !----------------- Calling findinv ------------
        call invert(epsmat, chi, nmtx, error)
        !----------------------------------------------

        print *, "epsilon matrix inverse DONE!"

        chi = chi *  cmplx ( chi0f_real , chi0f_cmpx )
        
        do ii = 1,size(chi,1)
           write(99,*) "Chi0 for ifreq # %d", ifreq, "\n \n"
           write(99,100) (chi(ii,jj), jj=1,nmtx)
        enddo

        ! deallocate epsmat, chi
        deallocate(epsmat, stat=error) 
        deallocate(chi, stat=error) 

        ! deallocate chi0f individual arrays for real and imag parts     
        deallocate( chi0f_real, stat=error)
        deallocate( chi0f_cmpx, stat=error)

     enddo ! loop over frequencies

     ! deallocate chi0w and individual arrays for real and imag parts     
     deallocate( chi0w, stat=error)
     deallocate( chi0w_real, stat=error)
     deallocate( chi0w_cmpx, stat=error)
 
     ! deallocating gvecs and indx
     deallocate(gvecs, stat=error)
     deallocate(indx, stat=error)

  enddo ! loop over q-vectors
  
100 format(3f25.15)
  
close(99)
end program chispin
