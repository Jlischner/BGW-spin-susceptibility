*     
*     Program that Caclulates inverse matrix using Lapack procedures
*     NMAX - maximum size of the matrix below
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NOUT=6)
      INTEGER          NMAX, LDA, LWORK
      PARAMETER        (NMAX=1000,LDA=NMAX,LWORK=64*NMAX)
*     .. Local Scalars ..
      INTEGER          I, IFAIL, INFO, J, N
*     .. Local Arrays ..
      COMPLEX *16      A(LDA,NMAX), WORK(LWORK)
      INTEGER          IPIV(NMAX)
      CHARACTER        CLABS(1), RLABS(1)
*     .. External Subroutines ..
      EXTERNAL         ZGETRF, ZGETRI
*     .. Executable Statements ..
      WRITE (NOUT,*) 'F07AWF Example Program Results'
*     Skip heading in data file
      open(55,file="f07awfe.d",form='formatted',status='old')

      NIN = 55
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*     
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
*        Factorize A
*
         CALL ZGETRF(N,N,A,LDA,IPIV,INFO)
*
!         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Compute inverse of A
*
            CALL ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
*
*           Print inverse
*
            IFAIL = 0
            do i=1,N 
               print *, (A(i,j), j=1,N)
            enddo
         ELSE
            WRITE (NOUT,*) 'The factor U is singular'
         END IF
      END IF
*
      END
