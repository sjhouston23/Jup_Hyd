module formatting
implicit none

!****************************** Header Variables *******************************
character(len=*),parameter :: &
H01="(' Charge state equilibrium fractions.',/,' Energy [keV/u]',4x,'H-',9x,&
  'H',8x,'H+')",&
H02="(' Atmospheric H^+ production rate [cm^-3 s^-1] from hydrogen ion',&
  ' precipitation.',/,&
  ' Input of 1 ion/cm^2/s.',/,&
  ' Alt [km] Production Rate')",&
H03="(' Atmospheric H_2^+ production rate [cm^-3 s^-1] from hydrogen ion',&
  ' precipitation.',/,&
  ' Input of 1 ion/cm^2/s.',/,&
  ' Alt [km] Production Rate')",&
H04="(' Atmospheric H_2 excitation rate [cm^-3 s^-1] from hydrogen ion',&
  ' precipitation.',/,&
  ' Input of 1 ion/cm^2/s.',/,&
  ' H_2 is only excited from a single process, TEX.',/,&
  ' Alt [km] Production Rate')",&
H05="(' Altitude integrated photon production [photons cm^-2 s^-1] as a',&
  ' function of charge state.')",&
H06="(' ΔAlt [km]',7x,'H-',9x,'H',8x,'H+')",&
H07="(' Photon production [photons cm^-3 s^-1] as a function of altitude',&
  ' and charge state.')",&
H08="(' Alt [km]',8x,'H-',9x,'H',8x,'H+')",&
H09="(' Initial input of 1 ion cm^-2 s^-1')",&
H10="(' E',14x,'SP        Sig       dE     dN        SP1           Ions',/,&
  ' Various data to calculate stopping power for input of 1 ion/cm^2/s.',/,&
  ' Conversion from [eV*cm^2] to [(MeV/mg)*cm^2] is 3.345e-15.',/,&
  ' Ion Energy      dE/dN      Sigma      dE       dN     Sigma*dE     ',&
  ' Counts',/,' [keV/u]       [eV*cm^2]   [cm^2]     [eV]    [cm^-2]',&
  ' [eV*cm^2]     per bin')",&
!******************************* Notes Variables *******************************
N01="(' Note: These charge states are the resultant charge state from a',&
  ' collision. E.g. the production under H is the number',/,6x,&
  ' of photons produced from H+ gaining an electron that cascades',&
  ' resulting in H. H+ is,',/,6x,&
  ' meaningless here - it is only output as a check that the model is',&
  ' outputting correctly and should always,',/,6x,&
  ' be 0.00E+00.')",&
N02="(' Note: The photon production below is from direct excitation of a',&
  ' charge state. If an electron gets excited and then',/,6x,&
  ' relaxes back down to a more energy favorable state, it will emit a',&
  ' photon. H+ is meaningless here - it is only',/,6x,&
  ' output as a check that the model is outputting correctly and should',&
  ' always be 0.00E+00.')"
!************************** Data Formatting Variables **************************
!* Notes:
!*   1x to create a space before the data
!*   F7.2 is used for atmosphere ==> "3000.00" or " -88.00"
!*   F8.2 is used for energy ==> "25000.00" or "  150.00"
!*   ES9.2 is used for large/small data values ==> "-1.34E+04" or " 2.31E-11"
!*   Electron two-stream formatting - 912 (F2Str)
!*******************************************************************************
character(len=*),parameter :: &
F1="(1x,A10,43x,F7.2)",&
F2="(1x,A26,I4,':',I2,':',F5.2)",&
F3="(1x,60('*'))",&
F4="(135('-'))",&
F01="(1x,F8.2,5x,3(2x,ES8.2))",&
F02="(1x,F7.2,5x,ES8.2)",&
F03="(10(8x,A3,1x),12x,'Sum')",&
F04="(10(I11,1x),3x,I12)",&
F04a="(9x,5(I12,1x))",&
F05="(10(F11.3,1x))",&
F06="(1x,F8.2,3(2x,ES8.2))",&
F07="(1x,F8.2,6x,2(ES9.2,1x),F8.2,2(1x,ES9.2),1x,I11)",&
F2Str="(1x,1P10E11.3)"  !Electron 2-stream formatting (912)
!****************************** Processes Header *******************************
character(len=3) Coll(10)
! character(len=4) TargColl2(7),ProjColl2(5)
data Coll/' SI',' DI',' TI',' SS',' DS',' SC',' DC','TEX','PEX',' ES'/

end module
