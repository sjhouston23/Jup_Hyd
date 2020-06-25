program InterpolateSIMXS
!*******************************************************************************
!* Created by Stephen J. Houston 06.24.20
!*******************************************************************************
!* This program reads in cross-sections for hydrogen from a .txt files and
!* interpolates them using both loglog linear and loglog spline interpolation.
!* All the cross-sections are from Schultz et al. (2019).
!* The initial cross-sections are stored in SIMXSall.dat, which is generated
!* with SIMXSCalculation.f08 in the SulfurXS directory.
!*******************************************************************************

implicit none!real*8(a-h,o-z)

!**************************** Variable Declaration *****************************

integer i,j,k,l,m,n
integer nProc,nEnergies,nInterpEnergies,nChS,ndum

integer ChS,Eng,E !Charge state and energy
integer Proc !Process
integer SI,DI,TI,SS,DS,SC,DC,TEX,PEX,ES! Processes

parameter(nProc=10) !Number of processes
parameter(nEnergies=15) !Number of inital energies
parameter(nInterpEnergies=25000) !Number of interpolated energies
parameter(nChS=3) !Number of charge states from 0-16
parameter(SI=1,DI=2,TI=3,SS=4,DS=5,SC=6,DC=7,TEX=8,PEX=9,ES=10) !Process Indices
! *****
! SI  = 1  - Single Ionization
! DI  = 2  - Double Ionization
! TI  = 3  - Transfer Ionization
! SS  = 4  - Single Stripping
! DS  = 5  - Double Stripping
! SC  = 6  - Single Capture
! DC  = 7  - Double Capture
! TEX = 8  - Target Excitation
! PEX = 9  - Projectile Excitation
! ES  = 10 - Elastic Scattering
! *****

real*8,dimension(nEnergies) :: Energy,xs_tmp2 !Each initial energy
real*8,dimension(nEnergies,2) :: xs_tmp
real*8 xs_tmpI,SpE,dum !Spline variables
real*8,dimension(nEnergies,nChS,2) :: SigTot
real*8,dimension(nInterpEnergies,nChS,2) :: SigTotInterp
real*8,dimension(nChS,nEnergies,nProc,2) :: xs !xs points
! 1 - CTMC XS
! 2 - Normailzed XS
real*8,dimension(nProc,nChS,nInterpEnergies,2) :: xsInterp !Interpolated xs

character(len=100) :: filename
character(len=50) :: tline
character(len=3),dimension(nProc) :: ProcNames !Target processes

! logical useNormalizedData

!****************************** Data Declaration *******************************
data Energy/1.0,2.0,5.0,10.0,25.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,&
            5000.0,10000.0,25000.0/
data ProcNames/'SI ','DI ','TI ','SS ','DS ','SC ','DC ','TEX','PEX','ES '/
!**************************** Initialize Variables *****************************
xs=0.0;xsInterp=0.0;SigTot=0.0;SigTotInterp=0.0;xs_tmp=0.0;xs_tmp2=0.0;
xs_tmpI=0.0
! useNormalizedData = .true.
!******************************** Main Program *********************************
!* Read in integral cross-sections
open(unit=100,file='./XS/OG_Integral_XS_CTMC.txt',status='old')
do Proc=1,nProc
  do i=1,4
    read(100,*) tline
    ! write(*,*) tline
  end do
  do Eng=1,nEnergies
    read(100,101) ndum,(xs(ChS,Eng,Proc,1),ChS=1,nChS)
    ! write(*,*) ndum,(xs(ChS,Eng,Proc),ChS=1,nChS),Proc,Eng
  end do
end do
close(100)
101 format(I6.1,3(4x,ES8.2E2))
open(unit=200,file='./XS/OG_Integral_XS_CTMC_Summed.dat')
!* Read in normalized integral cross-sections
open(unit=100,file='./XS/OG_Integral_XS_Normalized.txt',status='old')
do Proc=1,nProc
  do i=1,4
    read(100,*) tline
    ! write(*,*) tline
  end do
  do Eng=1,nEnergies
    read(100,102) ndum,(xs(ChS,Eng,Proc,2),dum,ChS=1,nChS)
    ! write(*,*) ndum,(xs(ChS,Eng,Proc),dum,ChS=1,nChS),Proc,Eng
  end do
end do
close(100)
102 format(I6.1,3(4x,ES8.2E2,4x,F4.3))
open(unit=201,file='./XS/OG_Integral_XS_Normalized_Summed.dat')
!* Calculate the sum of the integral cross-sections and write them out
write(200,*) 'Energy    H^-        H        H^+'
write(201,*) 'Energy    H^-        H        H^+'
do j = 1,2
  do i=1,nEnergies
    write(199+j,20400) int(Energy(i)),(sum(xs(ChS,i,:,j)),ChS=1,nChS)
    do ChS=1,nChS
      SigTot(i,ChS,j)=sum(xs(ChS,i,:,j))
    end do
  end do
  close(199+j)
end do
20400 format(I7,3(2x,ES8.2E2))

open(unit=206,file='./XS.dat')
do Proc=1,1!nProc !Loop through every process
  do ChS=1,1!nChS !Loop through every charge state (0-16)
    Eng=1 !Set Eng variable back to 1
    SpE=0.0
    !Calculate second derivative array
    ! do i=1,nEnergies
    !   xs_tmp(i) = xs(ChS,i,Proc,1)
    ! end do
    xs_tmp(:,:) = xs(ChS,:,Proc,:)
    call spline(log(Energy),log(xs_tmp(:,1)),nEnergies,xs_tmp2)
    do E=1,nInterpEnergies !Interpolation loop (1-25000 keV/u)
      SpE=SpE+1.0
      if(E.ge.Energy(Eng+1)) Eng=Eng+1 !Go to next Energy when appropriate
      if(E.eq.Energy(nEnergies)) Eng=nEnergies !Don't want Eng to go out of bounds

      if(E.lt.1.or.E.gt.10000)then
        xsInterp(Proc,ChS,E,1)=log(xs(ChS,Eng,Proc,1))+&
        (log(real(E))-log(Energy(Eng)))*&
        (log(xs(ChS,Eng+1,Proc,1))-log(xs(ChS,Eng,Proc,1)))/&
        (log(Energy(Eng+1))-log(Energy(Eng)))

        xsInterp(Proc,ChS,E,1)=exp(xsInterp(Proc,ChS,E,1))
      else
      ! if(E.ge.1.and.E.le.10000)then
      !   if(tProc.eq.1.and.pProc.gt.1.and.ChS.ge.8)goto 9000 !SI processes
      !   if(tProc.eq.1.and.pProc.eq.5)goto 9000 !SI+DPEX
      !   if(tProc.eq.2.and.pProc.gt.1.and.ChS.ge.9)goto 9000 !DI processes
      !   if(tProc.eq.2.and.pProc.eq.3.and.ChS.ge.8)goto 9000 !DI+DS
      !   if(tProc.eq.2.and.pProc.eq.5.and.ChS.ge.8)goto 9000 !DI+DPEX
      !   if(tProc.eq.3.and.pProc.eq.1.and.ChS.le.3)goto 9000 !TI
      !   if(tProc.eq.3.and.pProc.gt.1.and.ChS.ge.10)goto 9000 !TI processes
      !   if(tProc.eq.4.and.pProc.eq.1.and.ChS.le.7)goto 9000 !DCAI
      !   if(tProc.eq.4.and.pProc.gt.1)goto 9000 !DCAI processes
      !   if(tProc.eq.5.and.pProc.gt.1.and.ChS.ge.8)goto 9000 !SC processes
      !   if(tProc.eq.6.and.pProc.eq.1.and.ChS.le.6)goto 9000 !DC
      !   if(tProc.eq.6.and.pProc.eq.1.and.ChS.ge.9)goto 9000 !DC
      !   if(tProc.eq.6.and.pProc.gt.1.and.ChS.eq.3)goto 9000 !DC S^++
      !   if(tProc.eq.6.and.pProc.gt.1.and.ChS.ge.8)goto 9000 !DC S^++
      !   if(tProc.eq.7.and.pProc.gt.1.and.ChS.ge.6)goto 9000 !TEX processes
      !   if(E.eq.10)then
      !     do i=1,nEnergies
      !       xs_tmp(i)=xs(ChS,i,tProc,pProc) !Create a vector for spline
      !     end do
      !   end if

        call splineinterp(log(SpE),log(Energy),log(xs_tmp(:,1)),nEnergies,&
                          xs_tmp2,xs_tmpI)
        xsInterp(Proc,ChS,E,1)=exp(xs_tmpI)
      end if
      write(206,*) E,SpE,Energy(Eng),xs_tmp(Eng,1),nEnergies,xs_tmp2(Eng),xsInterp(Proc,ChS,E,1),xs_tmpI
      9000 continue
      ! if(xsInterp(tProc,pProc,ChS,E).ge.3E-14)then
      !   write(*,*) xsInterp(tProc,pProc,ChS,E),tProc,pProc,ChS,E
      !   stop
      ! end if
      ! write(*,*) E,xsInterp(tProc,pProc,ChS,E)
      ! if(tProc.eq.1.and.pProc.eq.1)then
      !   SigTotInterp(E,ChS)=log(SigTot(Eng,ChS))+&
      !   (log(real(E))-log(Energy(Eng)))*&
      !   (log(SigTot(Eng+1,ChS))-log(SigTot(Eng,ChS)))/&
      !   (log(Energy(Eng+1))-log(Energy(Eng)))
      !   SigTotInterp(E,ChS)=exp(SigTotInterp(E,ChS))
      !   ! write(*,*) E,SigTotInterp(E,ChS)!,Energy(Eng)
      ! end if
    end do !End interpolation loop (1-2000 keV/u)
  end do !End loop through every charge state (0-16)
end do !End loop through every target process
! ! xsInterp=exp(xsInterp) !Revert back from log
! ! do i=1,nInterpEnergies
! !   write(206,20400) i,(sum(xsInterp(:,:,ChS,i)),ChS=10,10)
! ! end do
!
! do tProc=1,nTargProc !Loop through every target process
!   do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
!     write(filename,"('./SIMXSInterp/',A,A,'.dat')")&
!     trim(TargProcCap(tProc)),trim(ProjProcCap(pProc))
!     open(unit=101,file=trim(filename))
!     write(101,1002) (i-1,i=1,17)
!     do E=1,nInterpEnergies !Energy interpolation loop (1-2000 keV/u)
!       do ChS=1,nChS !Loop through every charge state
!         if(isnan(xsInterp(tProc,pProc,ChS,E)))& !Get rid of any NaN values
!         xsInterp(tProc,pProc,ChS,E)=0.0
!       end do !End loop through every charge state
!       write(101,1001) real(E),(xsInterp(tProc,pProc,ChS,E),ChS=1,nChS)
!     end do !End energy interpolation loop (1-2000 keV/u)
!     close(101)
!     write(filename,"('./SIMXS/',A,A,'p.dat')")&
!     trim(TargProcCap(tProc)),trim(ProjProcCap(pProc))
!     open(unit=102,file=trim(filename))
!     write(102,1002) (i-1,i=1,17)
!     do Eng=1,nEnergies
!       write(102,1001) Energy(Eng),(xs(ChS,Eng,tProc,pProc),ChS=1,nChS)
!     end do
!     close(102)
!   end do !End loop through every projectile process plus 1
! end do !End loop through every target process
!
! open(unit=103,file='./SIMXSInterp/SIMXSInterpAll.txt')
! do tProc=1,nTargProc !Loop through every target process
!   do pProc=1,nProjProc+1 !Loop through every projectile process plus 1
!     write(103,1111) TargProcCap2(tProc),ProjProcCap(pProc)
!     write(103,*) 'Energy     S        S^+       S^++      S^3+      S^4+      &
!     &S^5+      S^6+      S^7+      S^8+      S^9+     S^10+     S^11+     S^12+&
!     &     S^13+     S^14+     S^15+     S^16+'
!     do E=1,nInterpEnergies !Energy interpolation loop (1-2000 keV/u)
!       do ChS=1,nChS !Loop through every charge state
!         if(isnan(xsInterp(tProc,pProc,ChS,E)))& !Get rid of any NaN values
!         xsInterp(tProc,pProc,ChS,E)=0.0
!       end do !End loop through every charge state
!       write(103,1001) real(E),&
!         (xsInterp(tProc,pProc,ChS,E),ChS=1,nChS)
!     end do !End energy interpolation loop (1-2000 keV/u)
!   end do !End loop through every projectile process plus 1
! end do !End loop through every target process
! close(103)
! ! open(unit=103,file='./SIMXSInterp/SIMXSInterpAll.dat')
! ! write(103,1100) exp(xsInterp) !Write out every cross-section to a single file
! ! close(103)
! do eng=1,nInterpEnergies
!   do ChS=1,nChS
!     SigTotInterp(Eng,ChS)=sum(xsInterp(:,:,ChS,Eng))
!   end do
! end do
! open(unit=205,file='./SIMXSInterp_TotalOG.dat')
! write(205,*) 'Energy     S        S^+       S^++      S^3+      S^4+      &
! &S^5+      S^6+      S^7+      S^8+      S^9+     S^10+     S^11+     S^12+&
! &     S^13+     S^14+     S^15+     S^16+'
! do i=1,nInterpEnergies
!   write(205,20400) i,(SigTotInterp(i,ChS),ChS=1,nChS)
! end do
! close(205)
!
! 1000 format(17(ES9.3E2,1x))
! 1001 format(F7.2,17(1x,ES9.3E2))
! 1002 format(' Energy',4x,9('S^',I0,'+',6x),8('S^',I0,'+',5x))
! 1100 format(20(ES9.3E2,1x))
! 1111 format(A4,A5,'----------')
! 1112 format(I7,17(1x,ES9.3E2))

end program
