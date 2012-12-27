!===============================================================================
!
! Utilities:
!
! chi0toascii      Last Modified Dec 2012
!
! Serial code that converts chi0 matrix into a 2-column vector in ascii format 
! containing the real and imaginary parts of chi0. Input parameters
! are read from file chi0toascii.inp in the working directory.
!
! chi0toascii.inp:
!
! chi0mat      ! chimat file, full-frequency or static
! wfn0         ! wavefunction file on any grid
!
!===============================================================================

#include "f_defs.h"

program chi0toascii
  ! ----------------------------
  !  DECLARATIONS
  !-----------------------------
  use global_m
 
  implicit none
  integer :: i,j,k,iq,mq,ig,igp,indx,jj,np
  
  character*256, parameter :: fninp = "chi02ascii.inp"
  character*256 :: fnchi,fnwfn,fnrho
  
  integer :: qtot,ntranq
  ! I/O units
  integer :: outunit_chi0, outunit_gvecs, outunit_chi_, outunit_nmtx, inputunit_chi0
  
  integer :: freq_dep,nFreq,ii,nq,ng,nmtx,qgrid(3),qmax(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta,qvec(3)
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: dFreqGrid(:),ekin(:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
 
  complex(DPC), allocatable :: dFreqBrd(:)
  complex(DPC), allocatable :: chi0R(:)
  !-----------------------------
  ! Local variables declaration
  !-----------------------------
  outunit_chi0 = 8
  outunit_gvecs = 9
  outunit_chi_ = 10
  outunit_nmtx = 11
  inputunit_chi0=12
  !-----------------------------
  ! initial  I/O HANDLING
  !-----------------------------
  call open_file(unit=outunit_chi0,file='chimat.dat',form='formatted',status='replace')
  call open_file(unit=outunit_gvecs,file='gvecs.dat',form='formatted',status='replace')
  call open_file(unit=outunit_nmtx,file='nmtx.dat',form='formatted',status='replace')
  ! read input file
  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(fninp)
  
  call open_file(55,file=trim(fninp),form='formatted',status='old')
  
  read(55,'(a)') fnchi
  read(55,'(a)') fnrho
  
  call close_file(55)
  
  write(6,'(2a)')       "     eps  file  = ", trim(fnchi)
  write(6,'(2a)')       "     rho  file  = ", trim(fnrho)

!!-----------------------------
!!  BODY
!!-----------------------------
  !write(6,*)g(1)
  write(6,'(1x,"reading",1x,a,1x,"file",/)') trim(fnchi)

  !-----------------------------------
  ! READING THE CHI0 FILE
  !-----------------------------------
  call open_file(unit=inputunit_chi0,file=trim(fnchi),form='unformatted',status='old')

  ! Skipping one line
  read(inputunit_chi0)
  ! Reading Frequency dependence and number of frequencies
  read(inputunit_chi0) freq_dep,ii
 
  write(6,*) 'Frequency dependence ', freq_dep
 
  if (freq_dep.eq.2) then
    nFreq=ii
    write(6,*) 'Number of Fequencies ', nFreq
    SAFE_ALLOCATE(chi0R, (nFreq))
  endif
  
  ! Reading q-grid 
  read(inputunit_chi0) (qgrid(i),i=1,3)
  SAFE_ALLOCATE(dFreqGrid,(nFreq))
  SAFE_ALLOCATE(dFreqBrd,(nFreq))
 
  ! If frequancy dependent -> reading frequency grid and broadenings
  if (freq_dep.eq.2) then
    read(inputunit_chi0) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    if (nFreq.gt.1) dDeltaFreq=dFreqGrid(2)-dFreqGrid(1)
    dBrdning=IMAG(dFreqBrd(1))
  else
    do i=1,nFreq
      dFreqGrid(i)=dDeltaFreq*dble(i-1)
      dFreqBrd(i)=dBrdning
    enddo
    ! Skipping one line
    read(inputunit_chi0)
  endif
 
  ! Skipping one line
  read(inputunit_chi0)

  ! Skipping one line
  read(inputunit_chi0)

  ! Reading cutoff energy
  read(inputunit_chi0) ecuts
 
  write(6,*) 'Screened Coulomb cutoff', ecuts

  ! Skipping one line
  read(inputunit_chi0) !nrk

  ! Reading number of G-vectors, number of q-points
  read(inputunit_chi0) ng,qtot,qmax(1:3)
  write(6,*) 'Number of G-vectors, Number of k points', ng,qtot
  
  ! Going back one line 
  backspace(inputunit_chi0)
  
  SAFE_ALLOCATE(gvec,(3,ng))
  SAFE_ALLOCATE(ekin,(ng))
  SAFE_ALLOCATE(isrtx,(ng))

  ! Reading all the g-vectors
  read(inputunit_chi0) ng,qtot,qmax(1:3),((gvec(jj,ig),jj=1,3),ig=1,ng)
  write(6,*) 'G1', gvec(1,1),gvec(2,1),gvec(3,1)
  write(6,*) 'G2', gvec(1,2),gvec(2,2),gvec(3,2)

  ! Reading the number of q-vectors
  read(inputunit_chi0) nq
  write(6,*) 'Number of q vectors', nq

  do ii = 1,nq

    ! Reading ????
     read(inputunit_chi0) ntranq
     write(6,*) 'ntranq ', ntranq

     ! Reading number of matrix elements, 
     read(inputunit_chi0) nmtx,np,(isrtx(jj),ekin(jj),jj=1,ng)
     write(6,*) 'Size of chi0 ', nmtx
     write(outunit_nmtx,*) nmtx

     !-----------------------------------
     ! WRITING G-VECTORS
     !-----------------------------------
     do jj = 1,nmtx
        write(outunit_gvecs,*)gvec(1:3,isrtx(jj))
     enddo
     
     !-----------------------------------
     ! WRITING chi0
     !-----------------------------------
     do j=1,nmtx
        do i=1,nmtx            
           read(inputunit_chi0) (chi0R(k),k=1,nFreq)

           !-----------------------------------
           ! WRITING chi0(0,0)
           !-----------------------------------
           if ((j.eq.1).and.(i.eq.1)) then
               write(6,*) chi0R
           endif
              
           do indx = 1,nFreq
              write(outunit_chi0,100)chi0R(indx)
           enddo
        enddo

#ifdef CPLX
        do i=1,nmtx
           read(inputunit_chi0)
        enddo
#endif
     enddo
  enddo

  
  call close_file(inputunit_chi0)
  
  write(6,'(a,i6)')     "     omega num  = ", nFreq
  write(6,'(a,f7.3)')   "     omega step = ", dDeltaFreq
  write(6,'(a,f7.3,/)') "     omega brd  = ", dBrdning    

  call close_file(outunit_chi0)
  call close_file(outunit_gvecs)
  call close_file(outunit_nmtx)
!-----------------------------
 ! deallocate and finish
  
  SAFE_DEALLOCATE(chi0R)
  SAFE_DEALLOCATE(dFreqGrid)
  SAFE_DEALLOCATE(dFreqBrd)
  
100 format(3f25.15)

end program chi0toascii
