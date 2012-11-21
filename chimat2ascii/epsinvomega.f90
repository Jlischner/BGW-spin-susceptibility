!===============================================================================
!
! Utilities:
!
! (1) epsinvomega      Originally By gsm      Last Modified 11/17/2010 (gsm)
!
! Serial code that plots frequency dependence of Plasmon-Pole, Retarded,
! Advanced epsilon inverse for a given q, G, G` vectors. Input parameters
! are read from file epsinvomega.inp in the working directory.
!
! epsinvomega.inp:
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

program epsinvomega

  use global_m
  implicit none

  integer :: i,j,k,iq,mq,ig,igp,itape,indx,jj,np
  real(DP) :: omega
  
  character*256, parameter :: fninp = "epsinvomega.inp"
  character*256 :: fneps,fnwfn,fnrho,fngpp,fnffr,fnffa
  real(DP) :: q(3)
  integer :: g(3),gp(3),gmgp(3),nband,ktot,ntranq
  
  integer :: freq_dep,nFreq,ii,nq,ng,nmtx,kgrid(3),kmax(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta,qvec(3)
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: dFreqGrid(:),ekin(:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
  SCALAR, allocatable :: eps(:)
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
  
  real(DP) :: wp2,qg(3),qgp(3),qgqg,qgqgp,lambda,phi
  SCALAR :: Omega2,wtilde2,epsggp,I_epsggp,eps_static,eps_dynamic
  complex(DPC) :: wtilde2_temp
!  real(DP) :: epsPPRe,epsPPIm
!  complex(DPC) :: wtilde,ampl

!-----------------------------
  call open_file(unit=9,file='fulleps.dat',form='formatted',status='replace')
  call open_file(unit=8,file='gvecs.dat',form='formatted',status='replace')
! read input file

  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(fninp)
  
  call open_file(55,file=trim(fninp),form='formatted',status='old')
  
  read(55,'(a)') fneps
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
  
  write(6,'(2a)')       "     eps  file  = ", trim(fneps)
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
! read eps file
  write(6,*)g(1)
  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fneps)

  itape=12
  call open_file(unit=itape,file=trim(fneps),form='unformatted',status='old')

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



  read(itape) nq
  write(6,*) 'Number of q vectors', nq

  do ii = 1,nq

     read(itape) ntranq
     write(6,*) 'ntranq ', ntranq

     read(itape) nmtx,np,(isrtx(jj),ekin(jj),jj=1,ng)
     write(6,*) 'Size of epsilon ', nmtx

     do jj = 1,nmtx
        write(8,*)gvec(1:3,isrtx(jj))
     enddo
     
     do j=1,nmtx
        do i=1,nmtx
           read(itape) (epsR(k),k=1,nFreq)
           do indx = 1,nFreq
              write(9,100)epsR(indx)
           enddo
        enddo
#ifdef CPLX
        do i=1,nmtx
           read(itape)
        enddo
#endif
     enddo
  enddo

!  read(itape) nq
!  read(itape) ng
!  
!  rewind(itape)
!  
!  SAFE_ALLOCATE(qpt, (3,nq))
!  SAFE_ALLOCATE(gvec, (3,ng))
!  
!  read(itape)
!  read(itape)
!  read(itape)
!  read(itape)
!  read(itape)
!  read(itape)
!  read(itape)
!  read(itape) nq,((qpt(j,i),j=1,3),i=1,nq)
!  read(itape) ng,((gvec(j,i),j=1,3),i=1,ng)
  
!  ig=0
!  igp=0
!  do i=1,ng
!    if (gvec(1,i).eq.g(1).and.gvec(2,i).eq.g(2).and. &
!      gvec(3,i).eq.g(3)) ig=i
!    if (gvec(1,i).eq.gp(1).and.gvec(2,i).eq.gp(2).and. &
!      gvec(3,i).eq.gp(3)) igp=i
!  enddo
!  if (ig.eq.0) then
!    call die("cannot find G vector in file " // trim(fneps))
!  endif
!  if (igp.eq.0) then
!    call die("cannot find G' vector in file " // trim(fneps))
!  endif
  
!  mq=-1
!  do iq=1,nq
!    if (abs(q(1)-qpt(1,iq)).lt.TOL_Zero.and. &
!      abs(q(2)-qpt(2,iq)).lt.TOL_Zero.and. &
!      abs(q(3)-qpt(3,iq)).lt.TOL_Zero) then
!      mq=iq-1
!      exit
!    endif
!  enddo
!  if (mq.eq.-1) then
!    call die("cannot find q vector in file " // trim(fneps))
!  endif
  
!  do iq=1,nq
!    read(itape) ng,nmtx
!    read(itape)
!    read(itape)
!    if (freq_dep.eq.0) then
!      do j=1,nmtx
!        read(itape)
!      enddo
!    endif
!    if (freq_dep.eq.2) then
!     do j=1,nmtx
!       do i=1,nmtx
!          read(itape) (epsR(k),k=1,nFreq)
!          do indx = 1,nFreq
!             write(9,100)epsR(indx)
!          enddo
!#ifdef CPLX
!          read(itape)
!#endif
!        enddo
!      enddo
!    endif
! enddo
  
!  SAFE_ALLOCATE(isrtx, (ng))
!  SAFE_ALLOCATE(isorti, (ng))
  
!  read(itape) ng,nmtx,(isrtx(i),isorti(i),i=1,ng)
!  read(itape)
!  read(itape) (qvec(i),i=1,3)
  
!  ig=isorti(ig)
!  igp=isorti(igp)
!  if (ig.eq.0) then
!    call die("cannot find G vector in file " // trim(fneps))
!  endif
!  if (igp.eq.0) then
!    call die("cannot find G' vector in file " // trim(fneps))
!  endif
  
!  SAFE_ALLOCATE(epsPP, (nFreq))
!  SAFE_ALLOCATE(epsA, (nFreq))
  
!  if (freq_dep.eq.0) then
!    SAFE_ALLOCATE(eps, (nmtx))
!    do j=1,nmtx
!      if (j.eq.igp) then
!        read(itape) (eps(i),i=1,nmtx)
!        epsStatic=eps(ig)
!      else
!        read(itape)
!      endif
!    enddo
!    SAFE_DEALLOCATE(eps)
!  endif
!  if (freq_dep.eq.2) then
!    do j=1,nmtx
!      do i=1,nmtx
!        if (j.eq.igp.and.i.eq.ig) then
!          read(itape) (epsR(k),k=1,nFreq)
!          do indx = 1,nFreq
!             write(9,100)epsR(indx)
!          enddo
!          epsStatic=epsR(1)
!        else
!          read(itape) (epsR(k),k=1,nFreq)
!          do indx = 1,nFreq
!             write(9,100)epsR(indx)
!          enddo
!        endif
!      enddo
!#ifdef CPLX
!      do i=1,nmtx
!        if (j.eq.igp.and.i.eq.ig) then
!          read(itape) (epsA(k),k=1,nFreq)
!        else
!          read(itape)
!        endif
!      enddo
!#endif
!    enddo
!#ifndef CPLX
!    do i=1,nFreq
!      epsA(i)=CONJG(epsR(i))
!    enddo
!#endif
!  endif
!  
!  SAFE_DEALLOCATE(isrtx)
!  SAFE_DEALLOCATE(isorti)
!  
!  SAFE_DEALLOCATE(qpt)
!  SAFE_DEALLOCATE(gvec)
  
  call close_file(itape)
  
  write(6,'(a,i6)')     "     omega num  = ", nFreq
  write(6,'(a,f7.3)')   "     omega step = ", dDeltaFreq
  write(6,'(a,f7.3,/)') "     omega brd  = ", dBrdning

!-----------------------------
! read wfn file
!
!  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fnwfn)
!  
!  itape=25
!  call open_file(unit=itape,file=trim(fnwfn),form='unformatted',status='old')
!
!  read(itape) ((bdot(i,j),i=1,3),j=1,3)
!  read(itape) celvol
!  
!  call close_file(itape)
!  
!  write(6,'(5x,"cel vol =",e20.12,/)') celvol
!
!-----------------------------
! read rho file

!  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fnrho)
!  
!  itape=95
!  call open_file(unit=itape,file=trim(fnrho),form='unformatted',status='old')
!  
!  rho0=0.0d0
!  rhogmgp=ZERO
  
!  read(itape)
!  read(itape) nvecs,nspin
!  read(itape) nproc_para
!  do i=1,nproc_para
!    read(itape) num_gvec
!    do j=1,num_gvec
!      read(itape) gx,gy,gz
!      read(itape) (xcdum(k),k=1,nspin)
!      if (gx.eq.0.and.gy.eq.0.and.gz.eq.0) then
!        do k=1,nspin
!          rho0=rho0+dble(xcdum(k))
!        enddo
!      endif
!      if (gx.eq.gmgp(1).and.gy.eq.gmgp(2).and.gz.eq.gmgp(3)) then
!        do k=1,nspin
!          rhogmgp=rhogmgp+ &
!#ifdef CPLX
!            xcdum(k)
!#else
!            dble(xcdum(k))
!#endif
!        enddo
!      endif
!    enddo
!  enddo
!    
!  if (abs(rho0).le.TOL_Zero) then
!    call die("cannot find rho(0) in file " // trim(fnrho))
!  endif
!  if (abs(rhogmgp).le.TOL_Zero) then
!    call die("cannot find rho(G-G') in file " // trim(fnrho))
!  endif
!  
!  call close_file(itape)
!  
!  write(6,'(5x,"rho(0) =",f10.3,/)') rho0

!-----------------------------
! construct generalized plasmon pole model

!  write(6,'(1x,"constructing GPP model",/)')
!  
!  wp2=ryd*ryd*16.0d0*PI_D*rho0/celvol
!  qg(:)=q(:)+dble(g(:))
!  qgp(:)=q(:)+dble(gp(:))
!  qgqg=dot_product(qg,matmul(bdot,qg))
!  qgqgp=dot_product(qg,matmul(bdot,qgp))
!  if (abs(qgqg) .lt. TOL_Zero) call die("GPP model diverges")
!  Omega2=wp2*qgqgp/qgqg*rhogmgp/rho0
!#ifdef CPLX
!  epsggp = epsStatic
!#else
!  epsggp = dble(epsStatic)
!#endif
!  if (all(g(1:3) .eq. gp(1:3))) then
!    delta = 1.0d0
!  else
!    delta = 0.0d0
!  endif
!  I_epsggp = delta - epsggp
!  if (abs(I_epsggp) .lt. TOL_Small) call die("GPP model diverges")
!#ifdef CPLX
! Complex GPP [PRB 40, 3162 (1989)]
!  wtilde2_temp = Omega2 / I_epsggp
!  lambda = abs(wtilde2_temp)
!  if (lambda .lt. TOL_Small) call die("GPP model diverges")
!  phi = atan2(IMAG(wtilde2_temp), dble(wtilde2_temp))
!  if (abs(cos(phi)) .lt. TOL_Small) call die("GPP model diverges")
!  wtilde2 = lambda / cos(phi)
!  Omega2 = Omega2 * CMPLX(1.0d0, -tan(phi))
!#else
! Real GPP [PRB 34, 5390 (1986)]
!  wtilde2 = Omega2 / I_epsggp
!  if (abs(wtilde2) .lt. TOL_Small) call die("GPP model diverges")
!#endif
!  wtilde=dble(sqrt(COMPLEXIFY(wtilde2)))
!  ampl=-0.5d0*PI_D*sqrt(COMPLEXIFY(Omega2/wtilde2))
!  do i=1,nFreq
!    omega=dFreqGrid(i)
!    epsPPRe=delta+dble(Omega2/(omega**2-wtilde2))
!    epsPPIm=dble(ampl)/sqrt(PI_D)/dBrdning* &
!     (exp(-(omega-wtilde)**2/dBrdning**2) &
!     -exp(-(omega+wtilde)**2/dBrdning**2))
!    epsPP(i)=CMPLX(epsPPRe,epsPPIm)
! Instead of using the above, we write real and imaginary in a compact 
! form by adding dbrdning in denominator. The imaginary part should be 
! a sharp (delta function) at the plasmon frequency.
!    epsPP(i)=delta+Omega2/(omega**2-wtilde2-CMPLX(0D0,dBrdning))
!  enddo
!
!  write(6,'(5x,"plasma frequency =",f10.3," eV",/)') sqrt(wp2)

!-----------------------------
! write generalized plasmon pole file

!  write(6,'(1x,"writing",1x,a,1x,"file",/)')trim(fngpp)
!    
!  call open_file(unit=7,file=trim(fngpp),form='formatted',status='replace')
!    
!  do i=1,nFreq
!    omega=dFreqGrid(i)
!    write(7,100)omega,epsPP(i)
!  enddo
!    
!  call close_file(7)

!-----------------------------
! write full frequency files

  if (freq_dep.eq.2) then
    write(6,'(1x,"writing",1x,a,1x,"and",1x,a,1x,"files",/)') &
      trim(fnffr),trim(fnffa)
    
!    call open_file(unit=8,file=trim(fnffr),form='formatted',status='replace')
!    call open_file(unit=9,file=trim(fnffa),form='formatted',status='replace')
    
!    eps_static=epsR(1)
!    eps_dynamic=eps_static
    
!    do i=1,nFreq
!      omega=dFreqGrid(i)
!      write(8,100)omega,epsR(i)
!      write(9,100)omega,epsA(i)
!      if (i .gt. 1) then
!        eps_dynamic=eps_dynamic-(2.0d0/PI_D)*(dFreqGrid(i)-dFreqGrid(i-1))* &
!          IMAG(epsR(i))*((1.0d0/omega)-1.0d0/(omega+ebind))
!      endif
 !   enddo
 !   write(6,*) 'Static Head:', eps_static, 'Dynamic Head:', eps_dynamic
 !   write(6,*)
      
 !   call close_file(8)
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
  
end program epsinvomega
