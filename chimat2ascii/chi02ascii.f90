!===============================================================================
!
! Utilities:
!
! (1) chi0toascii      Originally By gsm      Last Modified 11/17/2010 (gsm)
!
! Serial code that plots frequency dependence of Plasmon-Pole, Retarded,
! Advanced epsilon inverse for a given q, G, G` vectors. Input parameters
! are read from file chi0toascii.inp in the working directory.
!
! chi0toascii.inp:
!
! eps0mat      ! epsmat file, full-frequency or static
! wfn0         ! wavefunction file on any grid
! RHO          ! RHO file for Plasmon-Pole
! 0.0 0.0 0.0  ! q-vector in crystal coordinates
! 1 0 0        ! G-vector in crystal coordinates
! 1 0 0        ! G`-vector in crystal coordinates
! 201          ! nFreq, for static epsmat. For full-freq., overwritten from file
! 0.5          ! dDeltaFreq in eV, in case of static epsmat. See Above.
! 2.0          ! dBrdning in eV, in case of static epsmat. See Above.
! 0.1          ! Exciton binding energy to use to calculate the effective eps^-1 head for BSE.
! epsInvPP.dat ! Plasmon-Pole epsilon inverse, Re and Im parts
! epsInvR.dat  ! Retarded epsilon inverse, Re and Im parts
! epsInvA.dat  ! Advanced epsilon inverse, Re and Im parts
!
! Note that epsInvPP constructed with full-frequency epsmat will be
! different from epsInvPP constructed with static epsmat because of
! finite broadening used in full-frequency epsmat.
! JL
!===============================================================================

#include "f_defs.h"

program chi0toascii



!-----------------------------
!   DECLARATIONS
!-----------------------------

  use global_m
  implicit none

  integer :: i,j,k,iq,mq,ig,igp,itape,indx,jj,np
  real(DP) :: omega
  
  character*256, parameter :: fninp = "chi0toascii.inp"
  character*256, parameter :: filechi0 = "chi0_(0,0)_nfreq="
  character*256 :: fnchi,fnwfn,fnrho,fngpp,fnffr,fnffa
  real(DP) :: q(3)
  
  !integer :: g(3),gp(3),gmgp(3),nband,ktot,ntranq
  
  integer :: freq_dep,nFreq,ii,nq,ng,nmtx,kgrid(3),kmax(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta,qvec(3)
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: dFreqGrid(:),ekin(:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
  !SCALAR, allocatable :: eps(:)
  complex(DPC), allocatable :: dFreqBrd(:)
  complex(DPC), allocatable :: epsPP(:)
  complex(DPC), allocatable :: epsR(:)
  complex(DPC), allocatable :: epsA(:)
  
  real(DP) :: bdot(3,3),celvol,ebind
  integer :: nvecs,nspin
  
  integer nproc_para,num_gvec,gx,gy,gz
  complex(DPC) :: epsStatic,xcdum(2)
  real(DP) :: rho0
  SCALAR :: rhogmgp
  
  ! real(DP) :: wp2,qg(3),qgp(3),qgqg,qgqgp,lambda,phi
  ! SCALAR :: Omega2,wtilde2,epsggp,I_epsggp,eps_static,eps_dynamic
  ! complex(DPC) :: wtilde2_temp




!-----------------------------
!   I/O HANDLING
!-----------------------------

  call open_file(unit=9,file='chimat.dat',form='formatted',status='replace')
  call open_file(unit=8,file='gvecs.dat',form='formatted',status='replace')
! read input file

  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(fninp)
  
  call open_file(55,file=trim(fninp),form='formatted',status='old')
  
  read(55,'(a)') fnchi
  read(55,'(a)') fnwfn
  read(55,'(a)') fnrho
  read(55,*) (q(i),i=1,3)
  read(55,*) (g(i),i=1,3)
  read(55,*) (gp(i),i=1,3)
  read(55,*) nFreq
  read(55,*) dDeltaFreq
  read(55,*) dBrdning
  read(55,*) ebind
  read(55,'(a)') fngpp
  read(55,'(a)') fnffr
  read(55,'(a)') fnffa
  
  call close_file(55)
  
  write(6,'(2a)')       "     eps  file  = ", trim(fnchi)
  write(6,'(2a)')       "     wfn  file  = ", trim(fnwfn)
  write(6,'(2a)')       "     rho  file  = ", trim(fnrho)
  write(6,'(a, 3f7.3)') "     q  vector  = ", q(1:3)
  write(6,'(a, 3i4)')   "     G  vector  = ", g(1:3)
  write(6,'(a, 3i4)')   "     G' vector  = ", gp(1:3)
  write(6,'(a, i6)')    "     nFreq      = ", nFreq
  write(6,'(a, f7.3)')  "     dDeltaFreq = ", dDeltaFreq
  write(6,'(a, f7.3)')  "     dBrdning   = ", dBrdning
  write(6,'(a, f7.3)')  "     ebind      = ", ebind
  write(6,'(2a)')       "     GPP file   = ", trim(fngpp)
  write(6,'(2a)')       "     FFR file   = ", trim(fnffr)
  write(6,'(2a,/)')     "     FFA file   = ", trim(fnffa)
  
  gmgp(:)=g(:)-gp(:)
  


!-----------------------------
!   BODY
!-----------------------------


  write(6,*)g(1)
  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fnchi)

  !-----------------------------------
  ! READING THE CHI0 FILE
  !-----------------------------------

  itape=12
  call open_file(unit=itape,file=trim(fnchi),form='unformatted',status='old')

  read(itape)
  read(itape) freq_dep,ii
 
  write(6,*) 'Frequency dependence ', freq_dep
 
  if (freq_dep.eq.2) then
    nFreq=ii
    write(6,*) 'Number of Fequencies ', nFreq
    SAFE_ALLOCATE(epsR, (nFreq))
  endif
 
  read(itape) (kgrid(i),i=1,3)
  SAFE_ALLOCATE(dFreqGrid,(nFreq))
  SAFE_ALLOCATE(dFreqBrd,(nFreq))
 
  if (freq_dep.eq.2) then
    read(itape) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    if (nFreq.gt.1) dDeltaFreq=dFreqGrid(2)-dFreqGrid(1)
    dBrdning=IMAG(dFreqBrd(1))
  else
    do i=1,nFreq
      dFreqGrid(i)=dDeltaFreq*dble(i-1)
      dFreqBrd(i)=dBrdning
    enddo
    read(itape)
  endif
 
  read(itape)
  read(itape)
  read(itape) ecuts
 
  write(6,*) 'Screened Coulomb cutoff', ecuts

  read(itape) !nrk
  read(itape) ng,ktot,kmax(1:3)
  write(6,*) 'Number of G-vectors, Number of k points', ng,ktot
 
  backspace(itape)
  
  SAFE_ALLOCATE(gvec,(3,ng))
  SAFE_ALLOCATE(ekin,(ng))
  SAFE_ALLOCATE(isrtx,(ng))

  read(itape) ng,ktot,kmax(1:3),((gvec(jj,ig),jj=1,3),ig=1,ng)
  write(6,*) 'G1', gvec(1,1),gvec(2,1),gvec(3,1)
  write(6,*) 'G2', gvec(1,2),gvec(2,2),gvec(3,2)

  !-----------------------------------
  ! READING THE NUMBER OF Q-VECTORS
  !-----------------------------------
  read(itape) nq
  write(6,*) 'Number of q vectors', nq

  do ii = 1,nq

     read(itape) ntranq
     write(6,*) 'ntranq ', ntranq

     read(itape) nmtx,np,(isrtx(jj),ekin(jj),jj=1,ng)
     write(6,*) 'Size of epsilon ', nmtx

      !-----------------------------------
      ! WRITING G-VECTORS
      !-----------------------------------
     do jj = 1,nmtx
        write(8,*)gvec(1:3,isrtx(jj))
     enddo
     
      !-----------------------------------
      ! WRITING chi0
      !-----------------------------------
     do j=1,nmtx
        do i=1,nmtx            
           read(itape) (chi0R(k),k=1,nFreq)

      !-----------------------------------
      ! WRITING chi0(0,0)
      !-----------------------------------

           if (j.eq.1) && (i.eq.1) then

              !call open_file(unit=13,file=trim(filechi0),form='unformatted',status='replace') 
              write(6,'GOTCHA!')
           
           endif
              
           do indx = 1,nFreq
              write(9,100)chi0R(indx)
           enddo
        enddo
#ifdef CPLX
        do i=1,nmtx
           read(itape)
        enddo
#endif
     enddo
  enddo

  
  call close_file(itape)
  
  write(6,'(a,i6)')     "     omega num  = ", nFreq
  write(6,'(a,f7.3)')   "     omega step = ", dDeltaFreq
  write(6,'(a,f7.3,/)') "     omega brd  = ", dBrdning

!-----------------------------
! write full frequency files

  if (freq_dep.eq.2) then
    write(6,'(1x,"writing",1x,a,1x,"and",1x,a,1x,"files",/)') &
      trim(fnffr),trim(fnffa)
    

    call close_file(9)
    call close_file(8)
  endif

!-----------------------------
! deallocate and finish
  
  SAFE_DEALLOCATE(epsPP)
  SAFE_DEALLOCATE(epsR)
  SAFE_DEALLOCATE(epsA)
  SAFE_DEALLOCATE(dFreqGrid)
  SAFE_DEALLOCATE(dFreqBrd)
  
100 format(3f25.15)
  
end program chi0toascii
