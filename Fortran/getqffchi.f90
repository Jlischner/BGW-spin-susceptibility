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
  integer :: lenS,error,Nq,ii,Nfreq,iq,ind1,ind2,nmtx
  integer, allocatable, dimension(:) :: nmtxdat,indx

  real(DP), dimension(3,3) :: R,inverseR
  real(DP), allocatable, dimension(:) :: G2,Ixc
  real(DP), allocatable, dimension(:,:) :: G,chi0,gvecs,gvecsdat
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
  
  allocate( Ixc(lenS) , stat = error)
  call getIxc(lenS,Ixc)
   
  do iq = 1,Nq
     
     print *, "doing qpoint #", iq, "out of total", Nq

     nmtx = nmtxdat(iq)

     ! Check if the ind1 is zero and assign it
     if( iq-1 .eq. 0) then
        ind1 = 0
     else
        ind1 = sum(nmtxdat(1:iq-1))
     endif
     ind2 = sum(nmtxdat(1:iq))
     
     print *, "ind",ind1, ind2
     
     allocate(gvecs(nmtx,3), stat=error)
     allocate(indx(nmtx), stat=error)
     gvecs = gvecsdat( ind1+1:ind2, 1:3)
     
     print *,"calling fftbox"
     call get_fftbox(gvecs,S,nmtx,indx)
  
     do ii = 1,size(gvecs,1)
        !print *, "indx", ii, indx(ii)
     enddo

  enddo
  
  

end program chispin
