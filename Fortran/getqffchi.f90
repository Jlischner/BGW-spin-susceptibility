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
  integer :: lenS,error,Nq,ii,Nfreq
  integer, allocatable, dimension(:) :: nmtxdat

  real(DP), dimension(3,3) :: R,inverseR
  real(DP), allocatable, dimension(:) :: G2,Ixc
  real(DP), allocatable, dimension(:,:) :: G,chi0,gvecs
  real(DP) :: Vcell,dV

  ! Fourier grid size
  S(1) = 25
  S(2) = 25
  S(3) = 25
  lenS = S(1)*S(2)*S(3)

  ! lattice tensor (matrix of lattice vectors)
  R(1,1) = 20.0;  R(1,2) =  0.0;  R(1,3) =  0.0
  R(2,1) =  0.0;  R(2,2) = 20.0;  R(2,3) =  0.0
  R(3,1) =  0.0;  R(2,3) =  0.0;  R(3,3) = 20.0

  ! inverse lattice tensor
  inverseR(1,1) = 1./R(1,1); inverseR(1,2) = 0.0;       inverseR(1,3) = 0.0
  inverseR(2,1) = 0.0;       inverseR(2,2) = 1./R(2,2); inverseR(2,3) = 0.0
  inverseR(3,1) = 0.0;       inverseR(3,2) = 0.0;       inverseR(3,3) = 1./R(3,3)

  ! unit cell volume, volume element
  Vcell = R(1,1)*R(2,2)*R(3,3)
  dV    = Vcell/lenS
  
  ! number of qpoints
  Nq = 4

  ! number of frequencies
  Nfreq = 62

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
  allocate( gvecs(sum(nmtxdat**2),3) , stat=error)

  open(11,file="fullchi.dat",form='formatted')
  do ii=1,Nfreq*sum(nmtxdat**2)
     read(11,*) chi0(ii,1),chi0(ii,2)
     !print *, "chi0", ii,chi0(ii,1), chi0(ii,2)
  enddo
  close(11)

  open(12,file="gvecs.dat",form='formatted')
  do ii=1,sum(nmtxdat)
     read(12,*) gvecs(ii,1),gvecs(ii,2),gvecs(ii,3)
  enddo
  close(12)
  
  do ii = 1,sum(nmtxdat)
     print *, "gvecs", gvecs(ii,1),gvecs(ii,2),gvecs(ii,3)
  enddo
  
  allocate( Ixc(lenS) , stat = error)
  call getIxc(lenS,Ixc)
  
  do ii = 1,lenS
     !print *, "Ixc", ii, Ixc(ii)
  enddo

end program chispin
