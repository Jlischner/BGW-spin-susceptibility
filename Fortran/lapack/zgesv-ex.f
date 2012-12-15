*     ZGESV Example Program Text
*     NAG Copyright 2005.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX
      PARAMETER        (NMAX=8)
      INTEGER          LDA
      PARAMETER        (LDA=NMAX)
*     .. Local Scalars ..
      INTEGER          I, IFAIL, INFO, J, N
*     .. Local Arrays ..
      COMPLEX *16      A(LDA,NMAX), B(NMAX)
      INTEGER          IPIV(NMAX)
      CHARACTER        CLABS(1), RLABS(1)
*     .. External Subroutines ..
      EXTERNAL         X04DBF, ZGESV
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZGESV Example Program Results'
      WRITE (NOUT,*)
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A and B from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
         READ (NIN,*) (B(I),I=1,N)
*
*        Solve the equations Ax = b for x
*
         CALL ZGESV(N,1,A,LDA,IPIV,B,N,INFO)
*
         IF (INFO.EQ.0) THEN
*
*           Print solution
*
            WRITE (NOUT,*) 'Solution'
            WRITE (NOUT,99999) (B(I),I=1,N)
*
*           Print details of factorization
*
            WRITE (NOUT,*)
            IFAIL = 0
            CALL X04DBF('General',' ',N,N,A,LDA,'Bracketed','F7.4',
     +                  'Details of factorization','Integer',RLABS,
     +                  'Integer',CLABS,80,0,IFAIL)
*
*           Print pivot indices
*
            WRITE (NOUT,*)
            WRITE (NOUT,*) 'Pivot indices'
            WRITE (NOUT,99998) (IPIV(I),I=1,N)
         ELSE
            WRITE (NOUT,99997) 'The (', INFO, ',', INFO, ')',
     +        ' element of the factor U is zero'
         END IF
      ELSE
         WRITE (NOUT,*) 'NMAX too small'
      END IF
      STOP
*
99999 FORMAT ((3X,4(' (',F7.4,',',F7.4,')',:)))
99998 FORMAT (1X,7I11)
99997 FORMAT (1X,A,I3,A,I3,A,A)
      END
