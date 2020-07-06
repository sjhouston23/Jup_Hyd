program Combine_Output_Recursive
!*******************************************************************************
!* Created by Stephen J. Houston 11.29.18
!*******************************************************************************
!* This program combines electron 2-stream data from hydrogen precipitation
!*******************************************************************************

use, intrinsic :: ISO_FORTRAN_ENV
use formatting
implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer energy,atmosLen,ChS,err,start,run

! PROCESSES (Proc)
integer Proc,nProc !Number of processes
integer SI,DI,TI,SS,DS,SC,DC,TEX,PEX,ES! Processes
parameter(nProc=10)
parameter(SI=1,DI=2,TI=3,SS=4,DS=5,SC=6,DC=7,TEX=8,PEX=9,ES=10)

parameter(nChS=13,atmosLen=1544)
parameter(nOutputFiles=2,MaxnTrials=1000,MaxnLines=100000)
parameter(nEnergiesNorm=15,nEnergiesJuno=36)
parameter(nE2strBins=260) !Number of 2 stream bins

integer trial(MaxnTrials),nLines(nOutputFiles) !Number of trials/lines in a file

real*8 norm
real*8,allocatable,dimension(:) :: IonEnergy !Initial ion energies once decided
real*8,dimension(nEnergiesJuno) :: IonEnergyJuno !Initial ion energies for Juno
real*8,dimension(nEnergiesNorm) :: IonEnergyNorm !Initial ion energies normally
real*8,dimension(atmosLen,nE2strBins) :: electFwd,electBwd
real*8,dimension(atmosLen,nE2strBins) :: electFwdComb,electBwdComb

character(len=100) filename,files(nOutputFiles) !Output file names
character(len=9) date
!****************************** Data Declaration *******************************
!* Initial ion enegy input:
data IonEnergyNorm/1.0,2.0,5.0,10.0,25.0,50.0,75.0,100.0,200.0,500.0,1000.0,&
     2000.0,5000.0,10000.0,25000.0/
!* Initial ion enegy input from interpoalted JEDI bins:
data IonEnergyJuno/5.312,6.062,6.893,7.759,8.714,9.766,11.140,12.271,13.518,&
     14.892,16.660,18.638,20.851,23.326,24.817,26.403,28.090,29.885,31.892,&
     34.035,36.321,38.761,43.293,48.355,54.009,60.324,69.950,81.112,94.054,&
     109.062,131.160,157.734,189.692,228.125,270.78,312.5/ !Interpolated energies
data files/'2Str_Elect_Fwd','2Str_Elect_Bwd'/
!********************************* Initialize **********************************
energy=0;trial=0;nEnergies=0
!*******************************************************************************
EnergySwitch=1 !1 for normal energy bins, 2 for Juno energy bins
if(EnergySwitch.eq.1)then !Normal energy bins
  nEnergies=nEnergiesNorm
  allocate(IonEnergy(nEnergies))
  IonEnergy=IonEnergyNorm
elseif(EnergySwitch.eq.2)then !JEDI interpolated energy bins
  nEnergies=nEnergiesJuno
  allocate(IonEnergy(nEnergies))
  IonEnergy=IonEnergyJuno
end if
!******************************** Main Program *********************************
date='01Jul2020'
do run=nEnergies,nEnergies
  nTrials=15;norm=0.0
  energy=nint(IonEnergy(run))
  do i=1,nTrials
    trial(i) = 99+i
  end do
  write(*,*) 'Reading in and combining files...'
  write(*,*) 'Number of files: ',nTrials,'At an energy of: ',energy,'keV/u.'
!********************************* Initialize **********************************
  nerr=0;start=1;electFwdComb=0.0;electBwdComb=0.0
!********** Open output data files for each set of initial energies ************
  1002 continue
  do n=start,nTrials
!********************************* Initialize **********************************
    nLines=0;electFwd=0.0;electBwd=0.0
!******************************** Read in Files ********************************
    do i=1,nOutputFiles !Open all of the files
      write(filename,'("./Output/",I0,"/",A,"/",A,"-",I0,".dat")') &
            energy,date,trim(files(i)),trial(n) !File output name
      filename=trim(filename)
      open(unit=100+i,file=filename,status='old',iostat=err) !Open the files
      if(err.gt.0)then !If there's an error opening the file
        write(*,*) 'File:',n,'Trial:',trial(n),'ERROR! File:',filename !Error note
        start=n+1 !If there's an error, want to go to next trial
        nerr=nerr+1 !Count the number of errors accrued
        goto 1002 !Go to the next trial
      end if !End error opening file if statemen
    end do !End do loop for all files of a specific trial
    write(*,*) 'File:',n,'Trial:',trial(n)
  !*** 2-Stream electrons
    do j=1,nE2strBins
      read(101,F2Str) (electFwd(i,j),i=atmosLen,1,-1)
      read(102,F2Str) (electBwd(i,j),i=atmosLen,1,-1)
    end do
    electFwdComb=electFwdComb+electFwd
    electBwdComb=electBwdComb+electBwd
    do i=1,nOutputFiles !Close all of the files
      close(100+i)
    end do
  end do !End number of trials loop
!********************************** Write Out **********************************
  write(*,*) 'Writing output files...'
  do i=1,nOutputFiles !Open the final combined files
    write(filename,'("./Output/",I0,"/",A,"/Totals/",A,"-Total.dat")') &
          energy,date,trim(files(i))
    filename=trim(filename)
    open(unit=200+i,file=filename,status='unknown')
  end do
  norm=real(nTrials)-real(nerr) !Normalization condition to per ion per cm
  write(*,*) norm,nTrials,nerr
  !*** 2-Stream electrons
  do j=1,nE2strBins
    write(201,F2Str) (electFwdComb(i,j)/norm,i=atmosLen,1,-1)
    write(202,F2Str) (electBwdComb(i,j)/norm,i=atmosLen,1,-1)
  end do
  do i=1,nOutputFiles !Close all of the files
    close(200+i)
  end do
end do !End of energies do loop


end program
