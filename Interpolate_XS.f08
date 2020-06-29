program InterpolateSIMXS
!*******************************************************************************
!* Created by Stephen J. Houston 06.24.20
!*******************************************************************************
!* This program reads in cross-sections for hydrogen from a .txt files and
!* interpolates them using both loglog linear and loglog spline interpolation.
!* All the cross-sections are from Schultz et al. (2020).
!* The initial cross-sections are stored in OG_Integral_XS_CTMC.txt, which is
!* was hand made from the tables in the paper.
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
open(unit=200,file='./XS/Integral_XS_CTMC_Summed.dat')
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
open(unit=201,file='./XS/Integral_XS_Normalized_Summed.dat')
!* Calculate the sum of the integral cross-sections and write them out
write(200,*) 'Energy    H^-        H        H^+'
write(201,*) 'Energy    H^-        H        H^+'
do j = 1,2
  do i=1,nEnergies
    write(199+j,20400) int(Energy(i)),(sum(xs(ChS,i,:,j)),ChS=nChS,1,-1)
    do ChS=1,nChS
      SigTot(i,ChS,j)=sum(xs(ChS,i,:,j))
    end do
  end do
  close(199+j)
end do
20400 format(I7,3(2x,ES8.2E2))

! open(unit=206,file='./XS.dat')
do j=1,2 !CTMC and Normalized
  do Proc=1,nProc !Loop through every process
    do ChS=1,nChS !Loop through every charge state (+1 - -1)
      Eng=1 !Set Eng variable back to 1
      SpE=0.0
      !Calculate second derivative array
      ! do i=1,nEnergies
      !   xs_tmp(i) = xs(ChS,i,Proc,1)
      ! end do
      xs_tmp(:,:) = xs(ChS,:,Proc,:)
      call spline(log(Energy),log(xs_tmp(:,j)),nEnergies,xs_tmp2)
      do E=1,nInterpEnergies !Interpolation loop (1-25000 keV/u)
        SpE=SpE+1.0
        if(E.ge.Energy(Eng+1)) Eng=Eng+1 !Go to next Energy when appropriate
        if(E.eq.Energy(nEnergies)) Eng=nEnergies !Don't want Eng to go out of bounds

        if(E.lt.1.or.E.gt.10000)then
          xsInterp(Proc,ChS,E,j)=log(xs(ChS,Eng,Proc,j))+&
          (log(real(E))-log(Energy(Eng)))*&
          (log(xs(ChS,Eng+1,Proc,j))-log(xs(ChS,Eng,Proc,j)))/&
          (log(Energy(Eng+1))-log(Energy(Eng)))

          xsInterp(Proc,ChS,E,j)=exp(xsInterp(Proc,ChS,E,j))
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

          call splineinterp(log(SpE),log(Energy),log(xs_tmp(:,j)),nEnergies,&
                            xs_tmp2,xs_tmpI)
          xsInterp(Proc,ChS,E,j)=exp(xs_tmpI)
        end if
        if(isnan(xsInterp(Proc,ChS,E,j)))then
          if(E.gt.50)then
            xsInterp(Proc,ChS,E,j) = xsInterp(Proc,ChS,E-1,j)
          else
            xsInterp(Proc,ChS,E,j) = 0.0
          end if
        endif
        ! write(206,*) E,SpE,Energy(Eng),xs_tmp(Eng,1),nEnergies,xs_tmp2(Eng),xsInterp(Proc,ChS,E,1),xs_tmpI
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
end do

!* Write out the results
do j=1,2 ! Loop through CTMC and normalized
  if(j.eq.1)then
    open(unit=300,file='./XS/Integral_XS_CTMC_Interpolated.dat')
    write(300,3002) '! Interpolated Integral Cross-Sections for the CTMC Model'
  elseif(j.eq.2)then
    open(unit=300,file='./XS/Integral_XS_Normalized_Interpolated.dat')
    write(300,3002) '! Interpolated Integral Cross-Sections for the CTMC Model &
    &Normalized to ORNL Recommendations'
  end if
  do Proc=1,nProc !Loop through every process
    if(Proc.eq.SI)then
      write(300,3002) 'SI   ---------- H^{q+} + H_2; q = -1, 0, +1'
    elseif(Proc.eq.DI)then
      write(300,3002) '!'
      write(300,3002) 'DI   ---------- H^{q+} + H_2; q = -1, 0, +1'
    elseif(Proc.eq.TI)then
      write(300,3002) '!'
      write(300,3002) 'TI   ---------- H^{q+} + H_2; q = 0, +1'
    elseif(Proc.eq.SS)then
      write(300,3002) '!'
      write(300,3002) 'SS   ---------- H^{q+} + H_2; q = -1, 0'
    elseif(Proc.eq.DS)then
      write(300,3002) '!'
      write(300,3002) 'DS   ---------- H^{q+} + H_2; q = -1'
    elseif(Proc.eq.SC)then
      write(300,3002) '!'
      write(300,3002) 'SC   ---------- H^{q+} + H_2; q = 0, +1'
    elseif(Proc.eq.DC)then
      write(300,3002) '!'
      write(300,3002) 'DC   ---------- H^{q+} + H_2; q = +1'
    elseif(Proc.eq.TEX)then
      write(300,3002) '!'
      write(300,3002) 'TEX  ---------- H^{q+} + H_2; q = -1, 0, +1'
    elseif(Proc.eq.PEX)then
      write(300,3002) '!'
      write(300,3002) 'PEX   ---------- H^{q+} + H_2; q = 0'
    elseif(Proc.eq.ES)then
      write(300,3002) '!'
      write(300,3002) 'ElasticScattering - H^{q+} + H_2; q = -1, 0, +1'
    endif
    write(300,3002) 'Energy     Integral Cross Sections [cm^-2]'
    write(300,3002) ' [keV]       H^-          H          H^+'
    do E=1,nInterpEnergies !Loop through every energy
      write(300,3001) E,(xsInterp(Proc,ChS,E,j),ChS=nChS,1,-1)
    end do
  end do
  close(300)
end do
3001 format(I6,1x,3(3x,ES9.3e2))
3002 format(A)

do j=2,2 ! Loop through normalized factors
  if(j.eq.1)then
    ! open(unit=300,file='./XS/OG_Integral_XS_CTMC_Interpolated.dat')
    ! write(300,3002) '! Interpolated Integral Cross-Sections for the CTMC Model'
  elseif(j.eq.2)then
    open(unit=300,file='./XS/Normalized_Factors_Interpolated.dat')
    write(300,3002) '! Interpolated factors for the CTMC Model &
    &Normalized to ORNL Recommendations'
  end if
  do Proc=1,nProc !Loop through every process
    if(Proc.eq.SI)then
      write(300,3002) 'SI   ---------- H^{q+} + H_2; q = -1, 0, +1'
    elseif(Proc.eq.DI)then
      write(300,3002) '!'
      write(300,3002) 'DI   ---------- H^{q+} + H_2; q = -1, 0, +1'
    elseif(Proc.eq.TI)then
      write(300,3002) '!'
      write(300,3002) 'TI   ---------- H^{q+} + H_2; q = 0, +1'
    elseif(Proc.eq.SS)then
      write(300,3002) '!'
      write(300,3002) 'SS   ---------- H^{q+} + H_2; q = -1, 0'
    elseif(Proc.eq.DS)then
      write(300,3002) '!'
      write(300,3002) 'DS   ---------- H^{q+} + H_2; q = -1'
    elseif(Proc.eq.SC)then
      write(300,3002) '!'
      write(300,3002) 'SC   ---------- H^{q+} + H_2; q = 0, +1'
    elseif(Proc.eq.DC)then
      write(300,3002) '!'
      write(300,3002) 'DC   ---------- H^{q+} + H_2; q = +1'
    elseif(Proc.eq.TEX)then
      write(300,3002) '!'
      write(300,3002) 'TEX  ---------- H^{q+} + H_2; q = -1, 0, +1'
    elseif(Proc.eq.PEX)then
      write(300,3002) '!'
      write(300,3002) 'PEX   ---------- H^{q+} + H_2; q = 0'
    elseif(Proc.eq.ES)then
      write(300,3002) '!'
      write(300,3002) 'ElasticScattering - H^{q+} + H_2; q = -1, 0, +1'
    endif
    write(300,3002) 'Energy                 Factors'
    write(300,3002) ' [keV]       H^-          H          H^+'
    do E=1,nInterpEnergies !Loop through every energy
      write(300,3003) E,(xsInterp(Proc,ChS,E,2)/xsInterp(Proc,ChS,E,1),ChS=nChS,1,-1)
    end do
  end do
  close(300)
end do
! 3001 format(I6,3(4x,ES8.2e2))
! 3002 format(A)
3003 format(I6,3(4x,F8.4))

open(unit=400,file='./XS/Integral_XS_CTMC_Sum_Interpolated.dat')
open(unit=401,file='./XS/Integral_XS_Normalized_Sum_Interpolated.dat')
!* Calculate the sum of the integral cross-sections and write them out
write(400,3002) 'Energy       H^-          H          H^+'
write(401,3002) 'Energy       H^-          H          H^+'
do j = 1,2
  do i=1,nInterpEnergies
    write(399+j,3004) i,(sum(xsInterp(:,ChS,i,j)),ChS=nChS,1,-1)
    ! do ChS=1,nChS
    !   SigTot(i,ChS,j)=sum(xsInterp(ChS,i,:,j))
    ! end do
  end do
  close(399+j)
end do
3004 format(I6,2x,3(1x,ES11.5e2))

end program
