!================================================================
!
! Utilities:
!
! chispin
!
! wrapper for the interacting spin susceptibility chi calculation
!
!================================================================

program chispin

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3) :: S
  integer :: lenS, error, Nq, ii, Nfreq, iq, ind1, ind2, nmtx, ifreq, ncol
  integer, allocatable, dimension(:) :: nmtxdat,indx

  real(DP), dimension(3,3) :: R,inverseR
  real(DP), allocatable, dimension(:) :: G2,Ixc
  real(DP), allocatable, dimension(:,:) :: G,chi0,gvecs,gvecsdat, chi0w, chi0f, epsmat
  real(DP) :: Vcell,dV, latt_const

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

  ! sizes of chi0 matrices (number of Gvecs)
  allocate( nmtxdat(Nq), stat = error)

  ! get reciprocal lattice vectors (G) and their length squared (G2)
  allocate( G(lenS,3) , stat=error)
  allocate( G2(lenS)  , stat=error)
  
  call setup(inverseR,S,G,G2)
  
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
   
  !--------------------------------------------
  ! LOOP THROUGH THE q-vectors
  !--------------------------------------------
  do iq = 1,Nq
     
     print *, "doing qpoint #", iq, "out of total", Nq
     
     ! get the number of matrix elements for a current q-point iq
     nmtx = nmtxdat(iq)

     ! Check if the ind1 is zero and assign it and ind2
     !---
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1))
     endif
     
     ind2 = sum(nmtxdat(1:iq))
     !---
     
     ! allocating gvecs and indx
     allocate(gvecs(nmtx,3), stat=error)
     allocate(indx(nmtx), stat=error)

     ! assigning values to gvecs
     gvecs = gvecsdat( ind1+1:ind2, 1:3)
     
     !--------- call fftbox------------
     call get_fftbox(gvecs,S,nmtx,indx)
     !--------------------------------
     
     ! check if the ind1 is zero and assign it and ind2 again for chi0
     !---
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1)**2)
     endif
     
     ind2 = sum(nmtxdat(1:iq)**2)
     !---
     
     ! get the part of chi0 that corresponds to the current iq
     ! and reshape in into a matrix
     chi0w = chi0( ind1*Nfreq + 1 : ind2*Nfreq, : )
     chi0w = reshape (chi0w, (/ Nfreq, nmtx**2 /))
     
     !-------------------------------------------
     ! LOOP THROUGHT THE FREQUENCIES
     !-------------------------------------------
     do ifreq = 1,Nfreq
        
        print *, "doing frequency #", ifreq, "out of total", Nfreq
        chi0f = reshape (chi0w(ifreq,:), (/ nmtx,nmtx/))
        
        ! allocate epsmat TODO: nullify it also
        allocate(epsmat (nmtx,nmtx), stat=error) 
        
        !-------------------------------------------
        ! LOOP THROUGHT THE FREQUENCIES
        !-------------------------------------------
        do ncol = 1,nmtx
           
           print *, "doing column #", ncol, "out of total", nmtx
         
           !do ii = 1,size(chi0w,1)
           !   print *, "chi0f", ii, chi0f(ii,1)
           !enddo
        
        enddo
           
        ! deallocate epsmat
        deallocate(epsmat, stat=error) 
     enddo

     ! deallocating gvecs and indx
     deallocate(gvecs, stat=error)
     deallocate(indx, stat=error)

  enddo
  
  

end program chispin
