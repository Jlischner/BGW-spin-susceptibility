! ----------------------------------------------------------------
!     Subroutines that calculate inverse matrix using Lapack's z(d)gesv
!     Dec 2012
!
!     Author: TBZ
! ----------------------------------------------------------------
subroutine INVERT_CMPX(input, output, N)
!  .. Types ..   
   implicit none
   integer, parameter :: DP = kind(1.0d0)
   integer, parameter :: DPC = kind((1.0d0,1.0d0))
!  .. Parameters ....
   INTEGER, INTENT(IN) :: N
   COMPLEX(DPC), DIMENSION(N,N), INTENT(IN) :: input
   COMPLEX(DPC), DIMENSION(N,N), INTENT(OUT) :: output
   COMPLEX(DPC), allocatable, DIMENSION(:,:) :: A, temp
   COMPLEX(DPC), allocatable, DIMENSION(:) :: ipvt
   INTEGER          NOUT
   PARAMETER        (NOUT=6)
!  .. Local Scalars ..
   INTEGER          INFO, ii
!     .. Executable Statements ..
   WRITE (NOUT,*) 'invert_cmpx() subroutine is called'
!  .. Allocate temporary arrays ..
   allocate(temp(N,N), stat=INFO)
   allocate(ipvt(N), stat=INFO)
   allocate(A(N,N), stat=INFO)
!  .. Init Input ,,
   A = input
!  .. Create Identity matrix ..
   temp(:,:) = 0d0
   do ii=1,N
     temp(ii,ii) = 1d0
   enddo

   CALL ZGESV(N,N,A,N,ipvt,temp,N,INFO)

   IF (INFO.EQ.0) THEN
       output = temp
   ELSE
       WRITE (NOUT,*) 'The matrix inversion subroutine (CMPX) failed. Error #', INFO 
   END IF
!  .. Deallocate temporary arrays ..
   deallocate(temp, stat=INFO)
   deallocate(ipvt, stat=INFO)
   deallocate(A, stat=INFO)

end subroutine INVERT_CMPX
! ------------------------------------------------------------------
subroutine INVERT_REAL(input, output, N)
!  .. Types ..   
   implicit none
   integer, parameter :: DP = kind(1.0d0)
   integer, parameter :: DPC = kind((1.0d0,1.0d0))
!  .. Parameters ....
   INTEGER, INTENT(IN) :: N
   REAL(DP), DIMENSION(N,N), INTENT(IN) :: input
   REAL(DP), DIMENSION(N,N), INTENT(OUT) :: output
   REAL(DP), allocatable, DIMENSION(:,:) :: A, temp
   REAL(DP), allocatable, DIMENSION(:) :: ipvt
   INTEGER          NOUT
   PARAMETER        (NOUT=6)
!  .. Local Scalars ..
   INTEGER          INFO, ii
!     .. Executable Statements ..
   WRITE (NOUT,*) 'invert_real() subroutine is called'
!  .. Allocate temporary arrays ..
   allocate(temp(N,N), stat=INFO)
   allocate(ipvt(N), stat=INFO)
   allocate(A(N,N), stat=INFO)
!  .. Init Input ,,
   A = input
!  .. Create Identity matrix ..
   temp(:,:) = 0d0
   do ii=1,N
     temp(ii,ii) = 1d0
   enddo

   CALL DGESV(N,N,A,N,ipvt,temp,N,INFO)

   IF (INFO.EQ.0) THEN
       output = temp
   ELSE
       WRITE (NOUT,*) 'The matrix inversion subroutine (REAL) failed. Error #', INFO 
   END IF
!  .. Deallocate temporary arrays ..
   deallocate(temp, stat=INFO)
   deallocate(ipvt, stat=INFO)
   deallocate(A, stat=INFO)

end subroutine INVERT_REAL