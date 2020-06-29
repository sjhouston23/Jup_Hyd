subroutine CollisionSim(E,xs,xs_Total,ChS,excite,elect,disso,PID)
!*******************************************************************************
!* Created by Stephen J. Houston 06.26.20
!*******************************************************************************
!* This subroutine reads in each individual cross-section and the total cross-
!* section for all the processes for each energy and charge state for 10
!* different collision types (see processes below). Then uses a random
!* number from 0-1 to and compare it to the probability for each
!* collision process (individual cross-section divided by total cross-section).
!* This probability will be calcultated from the cross-sections for the
!* different collision processes. The probability that is larger than the random
!* number will be the selected outcome process of the collision.
!*******************************************************************************
!*
!* 	Input:
!*		E --> Energy of the ion
!*		 Type: Integer
!*		 Units: keV/u
!*
!*		xs --> Cross-sections
!*		 Type: Real*8 Matrix
!*		 Units: cm^-2
!*
!*		xs_Total --> Sum of cross-sections vs. energy and charge state
!*		 Type: Real*8 Matrix
!*		 Units: cm^-2
!*
!*		ChS --> Charge state
!*		 Type: Integer
!*		 Units: None
!*
!*   Returns:
!*    ChS --> Replaced with new charge state
!*		 Type: Integer
!*		 Units: None
!*
!*    excite --> Number of ion excitations
!*     Type: Integer
!*     Units: None
!*
!*    elect --> Number of electrons ejected
!*     Type: Integer
!*     Units: None
!*
!*    disso --> Number indicating whether dissociation is possible (0-2)
!*     Type: Integer
!*     Units: None
!*
!*    PID --> Collision Type (1-10)
!*     Type: 2-Component Array, Integer
!*     Units: None
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer,intent(in) :: E!,PID(2)
integer,intent(inout) :: ChS
integer,intent(out) :: excite,elect,disso,PID

integer SI,DI,TI,DCAI,SC,DC,TEX !Target processes
integer SS,DS,SPEX,DPEX !Projectile processes
integer nProc,Proc

parameter(nProc=10) !Number of processes
parameter(nInterpEnergies=25000) !Number of interpolated energies
parameter(nChS=3) !Number of charge states from -1, 0, +1
parameter(SI=1,DI=2,TI=3,SS=4,DS=5,SC=6,DC=7,TEX=8,PEX=9,ES=10)
parameter(k1=0,k2=0,lux=3) !lux set to 3 for optimal randomness and timeliness

real*8,dimension(nChS,nInterpEnergies),intent(in) :: xs_Total
real*8,dimension(nProc,nChS,nInterpEnergies),intent(in) :: xs

real*8 sumProb
real*8,dimension(nProc) :: Prob,Ptmp

real ranVecB(1) !Number from RNG
character(len=3),dimension(nProc) :: ProcNames !Target processes

data ProcNames/'SI ','DI ','TI ','SS ','DS ','SC ','DC ','TEX','PEX','ES '/
!**************************** Initialize Variables *****************************
sumProb=0.0;Prob=0.0;excite=0;elect=0;diss=0;PID=0;Ptmp=0.0
!******************************** Main Program *********************************
!******************** Collision-Type Probability Calculation *******************
!* Calculate the transition probabilities by taking the collision process XS and
!* dividing it by the total cross-section.
!*******************************************************************************
! goto 2000
do Proc=1,nProc !Loop through every process
  Prob(Proc)=xs(Proc,ChS,E)/xs_Total(ChS,E)
end do !End loop through every process
if(sum(Prob).ge.1.00001.or.sum(Prob).le.0.9999)then !Warning if probability is bad
  write(*,*) "CollisionSim.f08: WARNING: Normalized collision probability is &
             &not close enough to 1. The value is: ", sum(Prob)
  write(*,*) 'Process  Proc Prob    Proc XS  Total XS   Prob SUM'
  do Proc=1,nProc !Loop through every process
    Ptmp(Proc)=xs(Proc,ChS,E)/xs_Total(ChS,E)
    write(*,100) ProcNames(Proc),Ptmp(Proc),xs(Proc,ChS,E),xs_Total(ChS,E),sum(Ptmp)
  end do !End loop through every process
  ! stop
  write(*,*) "tempQ: ", ChS, "Energy: ", E
end if
100 format(2x,A4,3x,F10.8,2x,ES9.3E2,1x,ES9.2E2,1x,F10.8)
!*******************************************************************************
!************************* Collision-Type Determination ************************
!*******************************************************************************
!* Now that the collision probability is determined, a simple monte carlo
!* technique is used to determine the collision type. After collision-type is
!* determined the final charge state and collision-type are recorded as
!* tempQ and process, respectively.
!*
!* Algorithm:
!*  Step 1) Generate a random number on the interval [0,1]
!*  Step 2) Parce a sub-interval to cover the probability of a single
!*          collision type; i.e. for Single Ionization from 0->P(SI) [P(SI)<1]
!*  Step 3) Ask if the randomly generated number is in the interval [0,P(SI)],
!*          if so write out the final charge state and process
!*  Step 4) if not, parce the next sub-interval of [0,1] from [P(SI),P(DI)]
!*  Step 5) Repeat step 3.
!*  Step 6) Repreat step 4&5 through all 10 collision types until collision is
!*          determined.
!*
!* PID    Target Processes      Abbreviation  Charge Change  Hydrogen excitation
!*   1    Single Ionization         (SI)            0               0
!*   2    Double Ionization         (DI)            0               0
!*   3    Transfer Ionization       (TI)           -1               0
!*   4    Single Stripping          (SS)           +1               0
!*   5    Double Stripping          (DS)           +2               0
!*   6    Single Capture            (SC)           -1               1
!*   7    Double Capture            (DC)           -2               2
!*   8    Target Excitation         (TEX)           0               0
!*   9    Projectile Excitation     (PEX)           0               1
!*  10    Elastic Scattering        (ES)            0               0
!*
!* Some processes will dissociate H2 always, sometimes, or never. This is
!* followed with the "disso" variable.
!* disso = 0, never dissociates
!* disso = 1, 10% chance of dissociation (determined in the main program)
!* disso = 2, always dissociates
!*
!* PID is used as a processes identification
!*******************************************************************************
call ranlux(ranVecB,1) !Only need 1 random number every time there's a collision
if(ranVecB(1).gt.0.99999)ranVecB(1)=ranVecB(1)-0.00001
do Proc=1,nProc !Loop through every process
  sumProb=sumProb+Prob(Proc) !Keep adding the probability
  if(ranVecB(1).le.sumProb)goto 1000 !If it falls within range, get out
end do !End loop through every process
!If we get into this portion, it means that ranVecB was .gt. sumProb
write(*,*)'CollisionSim.f08: ERROR: Random number greater than normalized &
&probability:', ranVecB(1), sumProb
STOP 'CollisionSim.f08: Stopping program...'
1000 continue
PID=Proc
! 2000 continue
if(Proc.eq.SI)then !Single Ionization
  ChS=ChS !Charge state stays the same
  excite=0 !Ion isn't excited
  elect=1 !One electron ejected
  disso=1 !10% chance of dissociation
elseif(Proc.eq.DI)then !Double Ionization
  ChS=ChS
  excite=0
  elect=2
  disso=2
elseif(Proc.eq.TI)then !Transfer Ionization
  ChS=ChS-1
  excite=0
  elect=1
  disso=2
elseif(Proc.eq.SS)then !Single Stripping
  ChS=ChS+1
  excite=excite
  elect=1
  disso=1
elseif(Proc.eq.DS)then !Double Stripping
  ChS=ChS+2
  excite=excite
  elect=2
  disso=1
elseif(Proc.eq.SC)then !Single Capture
  ChS=ChS-1
  excite=1
  elect=0
  disso=1
elseif(Proc.eq.DC)then !Double Capture
  ChS=ChS-2
  excite=2
  elect=0
  disso=2
elseif(Proc.eq.TEX)then !Target Excitation
  ChS=ChS
  excite=0
  elect=0
  disso=0
elseif(Proc.eq.PEX)then !Projectile Excitation
  ChS=ChS
  excite=1
  elect=0
  disso=1
elseif(Proc.eq.ES)then !Elastic Scattering
  ChS=ChS
  excite=0
  elect=0
  disso=0
end if !End of processes
end subroutine
