subroutine EnergyLoss(E,ChS,eEnergy,PID,eEnergySS,eAngleSS,eEnergyDS,eAngleDS,dE)
!*******************************************************************************
!* Created by Stephen J. Houston 06.26.20
!*******************************************************************************
!* This routine calculates the energy loss of a precpitating ion.
!* The energy loss obtained by energy loss models based on models
!* calculated by Schultz et. al. 2020 in: Data for secondary electron production
!* from ion precipitation at Jupiter III: Target and projectile processes in
!* H^+, H, and H^- + H_2 collisions.
!* Table A, pg. 10, Tables 1 and 2, pg. 31 and 32.
!* Also from Schultz et. al. 2016 Ionization of molecular hydrogen and stripping
!* of oxygen atoms and ions... Table 4, pg. 37.
!*******************************************************************************
!*
!* 	Input:
!*		E --> Energy of ions
!*		Type: Real
!*		Units: keV/u
!*
!*		ChS --> Charge state
!*		Type: Integer
!*		Units: None
!*
!*    eEnergy --> Energy of secondary electron(s)
!*    Type: Real*8
!*    Units: eV
!*
!*    PID --> Process ID (1-10) see CollisionSim.f08
!*    Type: Integer
!*    Units: None
!*
!*    eEnergySS --> Energy of single stripping electrons
!*    Type: Real*8
!*    Units: eV
!*
!*    eAngleSS --> Angle of single stripping electrons
!*    Type: Real*8
!*    Units: Degrees
!*
!*    eEnergyDS --> Energy of double stripping electrons
!*    Type: 2-Component array,, real*8
!*    Units: eV
!*
!*    eAngleDS --> angle of double stripping electrons
!*    Type: 2-Component array, real*8
!*    Units: Degrees
!*
!* Returns:
!*    dE --> Delta Energy
!*    Type: Real
!*    Units: eV
!*
!*******************************************************************************
!*
!* Note: For single stripping and double stripping, the electron energies need
!* to be boosted into the projectile frame rather than the target frame,
!* produced by the SDXS. This calculation can be found in the appendix of
!* Schultz et. al. 2018.
!*
!*******************************************************************************

implicit real*8(a-h,o-z)

!**************************** Variable Declaration *****************************
integer,intent(in) :: ChS,PID
real*8,intent(in) :: E,eEnergy
real*8,intent(in) :: eEnergySS,eAngleSS,eEnergyDS(2),eAngleDS(2)
real*8,intent(out) :: dE !Change in energy

integer SI,DI,TI,SS,DS,SC,DC,TEX,PEX,ES! Processes
real*8 IP1,IP2 !First two ionization potentials for H2
real*8 pi,c,hMass !Speed of light and hydrogen mass (9.389e5 keV/c^2)
real*8 eMass !Electron mass
real*8 auCon !eV to a.u. conversion (1 a.u. = 27.2116 eV)
! real*8 alpha !Fine-structure constant (1/137)
real*8 ThetaP !Theta in the projectile frame

parameter(nChS=3) !Number of charge states from 0-16
parameter(nEnergies=15) !Number of inital energies
parameter(SI=1,DI=2,TI=3,SS=4,DS=5,SC=6,DC=7,TEX=8,PEX=9,ES=10)
parameter(c=137.036) !Speed of light in a.u.
parameter(hMass=9.389e5,auCon=27.2116) !hMass conversion to keV/c^2
parameter(eMass=0.511e6) !Electron mass in eV/c^2
parameter(pi=4.0*atan(1.0d0))
parameter(IP1=15.4254,IP2=16.4287) !Ionization potentials for hydrogen
!IP Data from J. Liu et al, Determination of the ionization and dissociation
!energy of the hydrogen molecule, J. Chem. Phys. 130, 174306 (2009).

integer Eng
integer iBin !IonEnergy bin

real*8 Vproj !Projectile velocity (precipitating ion)
real*8 Vz !Electron velocity in projectile frame
real*8 Vsquared !The square of the ejected electron's vel. in projectile frame
real*8 electEnergySS,electEnergyDS1,electEnergyDS2,electEnergyAU1,electEnergyAU2
real*8,dimension(nEnergies) :: IonEnergy !Initial ion energy
real*8,dimension(nChS-1) :: IP !Ionization potentials for hydrogen (q=-1, 0)
real*8,dimension(nEnergies,nChS) :: dETEX,dEPEX,dEES !Process
!****************************** Data Declaration *******************************
!* Initial ion energy input:
data IonEnergy/1.0,2.0,5.0,10.0,25.0,50.0,75.0,100.0,200.0,500.0,1000.0,2000.0,&
     5000.0,10000.0,25000.0/
!* Ionization potentials of atomic hydrogen (eV):
data IP/0.7543,13.62/
!* New energy loss/gain is in Schultz et al. 2020 (eV)
 !* Energy loss/gain for Target Excitation (Schultz 2020 Table 12)
 data dETEX/&
8.7,8.0,8.0,8.0,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,&    !H^-
12.1,10.0,8.22,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,& !H
12.8,11.6,8.37,7.9,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7,7.7/  !H^+

 !* Energy loss/gain for Projectile Excitation (Schultz 2020 Table 13)
 data dEPEX/&
0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,&!H^-
12.2,11.0,10.5,10.2,10.0,9.94,9.93,9.91,9.88,9.81,9.74,9.67,9.58,9.54,9.47,&!H
0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/ !H^+

!* Energy loss for Elastic Scattering (Schultz 2020 Table A)
data dEES/&
1.0,1.1,1.2,1.3,1.5,1.7,2.5,3.5,8.0,15.0,22.0,30.0,42.0,52.0,65.0,& !H^-
1.0,1.1,1.2,1.3,1.5,1.7,2.5,3.5,8.0,15.0,22.0,30.0,42.0,52.0,65.0,& !H
1.0,1.1,1.2,1.3,1.5,1.7,2.5,3.5,8.0,15.0,22.0,30.0,42.0,52.0,65.0/  !H^+
!******************************** Main Program *********************************
!* Initialize:
iBin=0;k=0;f=0.0
do Eng=nEnergies,1,-1 !Loop through bins to get the correct ion energy bin
  if(E.ge.real(IonEnergy(Eng-1)+(IonEnergy(Eng)-IonEnergy(Eng-1))/2.0))then
    iBin = Eng !Get the ion energy bin number
    k = 1 !Used to get an interpolation function
    goto 1000
  elseif(E.ge.real(IonEnergy(Eng-1)))then
    iBin = Eng-1 !Get the ion energy bin number
    k = 2 !Used to get an interpolation function
    goto 1000
  endif
end do
1000 continue
if (iBin.eq.0) write(206,*) 'EnergyLoss.f08: Error! iBin=0'
!* Want to use f to somewhat interpolate the cross-section for ion energies that
!* lie between energy bins.
if (k.eq.1) f=(E-IonEnergy(iBin-1))/(IonEnergy(iBin)-IonEnergy(iBin-1))
if (k.eq.2) f=(E-IonEnergy(iBin))/(IonEnergy(iBin+1)-IonEnergy(iBin))

!* Initialize:
dE=0.0;electEnergyAU1=0.0;Vproj=0.0;Vz=0.0;Vsquared=0.0;electEnergySS=0.0
ThetaP=0.0;electEnergyDS1=0.0;electEnergyDS2=0.0;electEnergyAU2=0.0

Vproj=sqrt(2.0*E/hMass) !Units of c

!* Find the dE for the given process based on Schultz et al. 2018 Table 1
!* Go through the "target" processes
if(PID.eq.SI)then !Single Ionization
  dE=IP1+eEnergy
elseif(PID.eq.DI)then !Double Ionization
  dE=IP1+IP2+eEnergy !Both electron energies have been added together
elseif(PID.eq.TI)then !Transfer Ionization
  dE=eEnergy+IP(ChS)+0.5*eMass*Vproj**2
  ! if(f.ge.0.5)dE=IP1+eEnergy+(f*dETI(iBin,ChS)+(1-f)*dETI(iBin-1,ChS))
  ! if(f.lt.0.5)dE=IP1+eEnergy+(f*dETI(iBin+1,ChS)+(1-f)*dETI(iBin,ChS))
elseif(PID.eq.SS)then !Single Stripping
  electEnergyAU1=eEnergySS/auCon !Convert to a.u.
  Vproj=Vproj*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU1)*cos(eAngleSS*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU1)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electEnergySS=(0.5)*Vsquared*auCon !Convert back to eV
  ThetaP=acos(Vz/sqrt(Vsquared))*180/pi
  dE=IP(ChS)+electEnergySS
elseif(PID.eq.DS)then !Double Stripping
  electEnergyAU1=eEnergyDS(1)/auCon !Convert to a.u.
  Vproj=Vproj*c !In a.u.
  Vz=sqrt(2.0*electEnergyAU1)*cos(eAngleDS(1)*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU1)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electEnergyDS1=(0.5)*Vsquared*auCon !Convert back to eV
  electEnergyAU2=eEnergyDS(2)/auCon !Convert to a.u.
  Vz=sqrt(2.0*electEnergyAU2)*cos(eAngleDS(2)*pi/180.0)-Vproj !In a.u.
  Vsquared=(2*electEnergyAU2)-(2*Vz*Vproj)-Vproj**2 !Electron velocity squared
  electenergyDS2=(0.5)*Vsquared*auCon !Convert back to eV
  dE=IP(ChS)+IP(ChS+1)+electenergyDS1+electenergyDS2
elseif(PID.eq.SC)then !Single Capture
  dE=-IP1+IP(ChS)+0.5*eMass*Vproj**2
  ! if(f.ge.0.5)dE=(f*dESC(iBin,ChS)+(1-f)*dESC(iBin-1,ChS))
  ! if(f.lt.0.5)dE=(f*dESC(iBin+1,ChS)+(1-f)*dESC(iBin,ChS))
elseif(PID.eq.DC)then !Double Capture
  dE=-IP1-IP2+IP(ChS-2)+IP(ChS-1)+0.5*eMass*Vproj**2
  ! if(f.ge.0.5)dE=(f*dEDC(iBin,ChS)+(1-f)*dEDC(iBin-1,ChS))
  ! if(f.lt.0.5)dE=(f*dEDC(iBin+1,ChS)+(1-f)*dEDC(iBin,ChS))
elseif(PID.eq.TEX)then !Target Excitation
  if(f.ge.0.5)dE=(f*dETEX(iBin,ChS)+(1-f)*dETEX(iBin-1,ChS))
  if(f.lt.0.5)dE=(f*dETEX(iBin+1,ChS)+(1-f)*dETEX(iBin,ChS))
elseif(PID.eq.PEX)then !Single Projectile Excitation
  if(f.ge.0.5)dE=(f*dEPEX(iBin,ChS)+(1-f)*dEPEX(iBin-1,ChS))
  if(f.lt.0.5)dE=(f*dEPEX(iBin+1,ChS)+(1-f)*dEPEX(iBin,ChS))
elseif(PID.eq.ES)then !Double Projectile Excitation
  if(f.ge.0.5)dE=(f*dEES(iBin,ChS)+(1-f)*dEES(iBin-1,ChS))
  if(f.lt.0.5)dE=(f*dEES(iBin+1,ChS)+(1-f)*dEES(iBin,ChS))
end if


end subroutine
