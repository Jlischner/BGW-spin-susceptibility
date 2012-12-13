subroutine do_FFT(fftbox,Nfft,sign)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  include 'fftw_f77.i'

  integer, dimension(3),INTENT(IN) :: Nfft
  integer, INTENT(IN) :: sign
  complex(DPC), dimension(Nfft(1),Nfft(2),Nfft(3)),INTENT(INOUT) :: fftbox
  integer*8, save :: plus_plan, minus_plan
    
  call fftwnd_f77_create_plan(plus_plan,3,Nfft,FFTW_BACKWARD, &
       FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
  call fftwnd_f77_create_plan(minus_plan,3,Nfft,FFTW_FORWARD, &
       FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)

   if (sign == 1) then
    call fftwnd_f77_one(plus_plan,fftbox,0)
  else if (sign == -1) then
    call fftwnd_f77_one(minus_plan,fftbox,0)
  else
    print *,"sign is not 1 or -1 in do_FFT"
    stop
  endif

  return
end subroutine do_FFT

!**************************************************************
subroutine cI(fft_in,Nfft)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3),INTENT(IN) :: Nfft
  complex(DPC), dimension(Nfft(1)*Nfft(2)*Nfft(3)),INTENT(INOUT) :: fft_in
  complex(DPC), dimension(Nfft(1),Nfft(2),Nfft(3)) :: fftbox
  integer :: sign

  fftbox = RESHAPE( fft_in, (/ Nfft(1), Nfft(2), Nfft(3) /))
  sign = 1
  call do_FFT(fftbox,Nfft,sign)
  fft_in = RESHAPE( fftbox, (/ Nfft(1)*Nfft(2)*Nfft(3) /))
  
  return
end subroutine cI

!*******************************************************************
subroutine cJ(fft_in,Nfft)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3),INTENT(IN) :: Nfft
  complex(DPC), dimension(Nfft(1)*Nfft(2)*Nfft(3)),INTENT(INOUT) :: fft_in
  complex(DPC), dimension(Nfft(1),Nfft(2),Nfft(3)) :: fftbox
  integer :: sign,fftsize
  
  fftsize = Nfft(1)*Nfft(2)*Nfft(3)
  fftbox = RESHAPE( fft_in, (/ Nfft(1), Nfft(2), Nfft(3) /))
  sign = -1
  call do_FFT(fftbox,Nfft,sign)
  fft_in = RESHAPE( fftbox, (/ Nfft(1)*Nfft(2)*Nfft(3) /))
  fft_in = fft_in / REAL(fftsize)

  return
end subroutine cJ

!*********************************************************************
subroutine dmult( vec, mat_in, mat_out, nline, ncol)
  
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  
  INTEGER :: nline,ncol,i
  REAL(DP), dimension(nline), INTENT(IN) :: vec
  COMPLEX(DPC), dimension(nline,ncol), INTENT(IN) :: mat_in
  COMPLEX(DPC), dimension(nline,ncol), INTENT(OUT) :: mat_out
  
  DO i=1,ncol
     mat_out(:,i) = vec * mat_in(:,i)
  END DO

  return
end subroutine dmult

!***********************************************************
! setup
! 
! This subroutine does all the preparations for reading
! the G-vectors
!
! parameters: 
! - inverseR --- inverse lattice tensor of dim(3,3) 
! - S --- dim(3) vector of Fourier grid
! - G --- dim(:,3) matrix of g-vectors
! - G2 --- dim (:) vector of G^2 values 
!
!***********************************************************
subroutine setup(inverseR,S,G,G2)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  real(DP), dimension(3,3), intent(in) :: inverseR
  integer, dimension(3), intent(in) :: S
  integer :: i,lenS,error
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ms,m1,m2,m3,n1,n2,n3
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: M,N
  real(DP), dimension(S(1)*S(2)*S(3),3), intent(out) :: G
  real(DP), dimension(S(1)*S(2)*S(3)), intent(out) :: G2
  real(DP) :: pi

  lenS = S(1)*S(2)*S(3)
  pi = 2.*ACOS(0.0)
  !------------------allocate statements------------------------
  ALLOCATE( ms(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for ms"
     STOP
  END IF

  ALLOCATE( m3(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for m3"
     STOP
  END IF

  ALLOCATE( m2(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for m2"
     STOP
  END IF

  ALLOCATE( m1(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for m1"
     STOP
  END IF

  ALLOCATE( n1(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for n1"
     STOP
  END IF

  ALLOCATE( n2(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for n2"
     STOP
  END IF

  ALLOCATE( n3(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for n3"
     STOP
  END IF

  ALLOCATE( M(lenS,3) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for M"
     STOP
  END IF

  ALLOCATE( N(lenS,3) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for N"
     STOP
  END IF
  !-----------------done allocating-----------------

  DO i=1,lenS
     ms(i) = i-1
  END DO

  m3 = MOD(ms,S(3))
  m2 = MOD(FLOOR( REAL(ms)/S(3)),S(2))
  m1 = MOD(FLOOR( REAL(ms)/(S(3)*S(2))),S(1))

  M(:,1) = m1; M(:,2) = m2; M(:,3) = m3

  DO i=1,lenS

     IF( m1(i) > REAL(S(1))/2 ) THEN
        n1(i) = 1
     ELSE
        n1(i) = 0
     END IF

     IF( m2(i) > REAL(S(2))/2 ) THEN
        n2(i) = 1
     ELSE
        n2(i) = 0
     END IF

     IF( m3(i) > REAL(S(3))/2 ) THEN
        n3(i) = 1
     ELSE
        n3(i) = 0
     END IF

  END DO

  n1 = m1-n1*S(1)
  n2 = m2-n2*S(2)
  n3 = m3-n3*S(3)
  N(:,1) = n1; N(:,2) = n2; N(:,3) = n3

  G  = 2.*pi*MATMUL(N,inverseR)
  G2 = SUM(G**2,2)

  return
end subroutine setup

!***********************************************************
! getIxc
! 
! This subroutine forms a matrix that represents an LSDA
! exchange and correlation interaction
!
! parameters: 
! - lenS --- dimensions of the Fourier grid in G-space 
!   multiplied e.g. lenS = dimX*dimY*dimZ
! - Ixc --- resulting output matrix
!
!***********************************************************
subroutine getIxc(lenS,Ixc)
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  
  integer, intent(in):: lenS
  integer :: error,ii
  real(DP) :: A,a1,b1,b2,b3,b4,pi
  real(DP), allocatable, dimension(:) :: n,rs,rrs,ac,expp
  real(DP), dimension(lenS), intent(out) :: Ixc

  character(len=200), parameter :: infile = "cd.dat"

  allocate( n(lenS) , stat=error)
!  allocate( Ixc(lenS) , stat=error) - this was giving memory error
  allocate( rs(lenS) , stat=error)
  allocate( rrs(lenS) , stat=error)
  allocate( ac(lenS) , stat=error)
  allocate( expp(lenS) , stat=error)
 
  pi = 2.0 * acos(0.0)
  !# parameter for -ac from Perdew Wang paper
  !# last column of TABLE I (atomic units)
  A  = 0.016887
  a1 = 0.11125
  b1 = 10.357
  b2 = 3.6231
  b3 = 0.88026
  b4 = 0.49671

  !# load in density in real space from espresso
  open(15,file=infile,form='formatted')
  do ii=1,lenS
     read(15,*) n(ii)
  enddo
  close(15)

  
  rs = (3.0/4.0/pi/n)**(1.0/3.0)
  rrs= sqrt(rs)
  
  !# calculate spin stiffness
  ac   = log(1.0+1.0/(2.0*A*( b1*rrs + b2*rrs**2.0 + b3*rrs**3.0 + b4*rrs**4.0)))
  ac   = 2.0*A*(1.0+a1*rs) * ac
 
  !# calculate 2nd derivative of exchange energy
  expp  = -3.0/4.0/pi/rs * (9.0*pi/4.0)**(1.0/3.0)
  expp  = 4.0/9.0 * expp

  !# final result for I(n)
  Ixc = 2.0/n * (ac + expp)
  
  !# Rydberg: multiply by 2
  Ixc = 2.0 * Ixc

end subroutine getIxc

!***********************************************************
! get_fftbox
! 
! This subroutine forms a vector that represents the index
! that a particular G-vector has in the Fourier grid S
!
! parameters:
!  input
! - gkvectors --- dim(nmtx,3) array of g-vectors
! - S --- dim(3) vector of the Fourier grid
! - nmtx --- the length of gkvectors array (number of g-vectors)
!  output
! - indx --- resulting output vector of dim(nmtx)
!
!***********************************************************
subroutine get_fftbox(gkvectors,S,nmtx,indx)
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  
  integer, dimension(3),intent(in):: S
  integer, intent(in) :: nmtx
  real(DP), dimension( nmtx,3),intent(in) :: gkvectors
  integer, dimension(nmtx),intent(out) :: indx
  real(DP), allocatable,dimension(:,:) :: ivec
  real(DP) :: isign
  integer :: error,Ngk,ii,jj

  allocate( ivec(nmtx,3) , stat=error)
  ivec = gkvectors + 1
  
  ! TODO: rewrite this loop using "where" later
  do ii=1, nmtx
     do jj=1,3
  
        isign = 0.0d0
        if( ivec(ii,jj) .le. 0.0d0) isign = 1.0d0
        ivec(ii,jj) = ivec(ii,jj) + real(S(jj)) * isign
     enddo
  enddo
  
  do ii=1,nmtx
     indx(ii) = S(2)*S(1)*(int(ivec(ii,3))-1) + S(1)*(int(ivec(ii,2))-1) + int(ivec(ii,1))
  enddo

end subroutine get_fftbox
