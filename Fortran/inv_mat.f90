!     
!     Program that Caclulates inverse matrix using Lapack procedures
!     NMAX - maximum size of the matrix below
subroutine INVERT(input, output, N)
   
   implicit none
   integer, parameter :: DP = kind(1.0d0)
   integer, parameter :: DPC = kind((1.0d0,1.0d0))
!  .. Parameters ....
   INTEGER, INTENT(IN) :: N
   COMPLEX(DPC), DIMENSION(N,N), INTENT(IN) :: input
   COMPLEX(DPC), DIMENSION(N,N), INTENT(OUT) :: output

   INTEGER          NIN, NOUT
   PARAMETER        (NOUT=6)
   INTEGER          NMAX, LDA, LWORK
   PARAMETER        (NMAX=1000,LDA=NMAX,LWORK=64*NMAX)
!  .. Local Scalars ..
   INTEGER          I, IFAIL, INFO, J
!  .. Local Arrays ..
   COMPLEX *16      A(LDA,NMAX), WORK(LWORK)
   INTEGER          IPIV(NMAX)
   CHARACTER        CLABS(1), RLABS(1)
!  .. Passing variables ..
   A = input
   output = 0d0
!  .. External Subroutines ..
!   EXTERNAL         ZGETRF, ZGETRI
!     .. Executable Statements ..
   WRITE (NOUT,*) 'invert() subroutine is called'

   
!--------------------------------------------------------------
! ORIGINAL FILEREAD STATEMENTS COMMENTED

!     Skip heading in data file
!   open(55,file="f07awfe.d",form='formatted',status='old')

!   NIN = 55
!   READ (NIN,*)
!   READ (NIN,*) N
!   IF (N.LE.NMAX) THEN
!     
!  Read A from data file
!
!      READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
!--------------------------------------------------------------

!  Factorize A
!
      CALL ZGETRF(N,N,A,LDA,IPIV,INFO  )
!
      WRITE (NOUT,*)
      IF (INFO.EQ.0) THEN
!
!           Compute inverse of A
!
            CALL ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
!
!           Print inverse
!
            IFAIL = 0
            output = A
!            do i=1,N 
!               print *, (A(i,j), j=1,N)
!            enddo
         ELSE
            WRITE (NOUT,*) 'The factor U is singular'
      END IF
!

end subroutine INVERT
