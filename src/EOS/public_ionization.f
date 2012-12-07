      program teos
      include 'implno.dek'
      include 'vector_eos.dek'

c..
c..tests the timmes eos with full ionization 
c..
c..ionmax  = number of isotopes in the network
c..xmass   = mass fractions
c..ymass   = molar fractions
c..aion    = number of nucleons
c..zion    = number of protons

      integer          i,j,ionmax
      parameter        (ionmax=1)
      character*5      ionam(ionmax)
      double precision xmass(ionmax),ymass(ionmax),
     1                 aion(ionmax),zion(ionmax),temp,den,abar,zbar


c..set the mass fractions, z's and a's of the composition
c..hydrogen
      xmass(1) = 1.00d0
      aion(1)  = 28.0d0
      zion(1)  = 14.0d0
      ionam(1) = 'si28  '

c..helium
c     xmass(2) = 0.23d0
c     aion(2)  = 4.0d0
c     zion(2)  = 2.0d0
c     ionam(2) = 'he4 '

c..carbon 12
c     xmass(3) = 0.02d0
c     aion(3)  = 12.0d0
c     zion(3)  = 6.0d0
c     ionam(3) = 'c12 '


c..set the input vector
      temp_row(1) = 1.5d4
      den_row(1)  = 1.0d-13

c..here the pipeline is only 1 element long
      jlo_eos = 1
      jhi_eos = 1


c..set up for a call to the 

      niso = ionmax
      do j=jlo_eos,jhi_eos
       do i=1,ionmax
        xmass_row(i,j) = xmass(i)
        aion_row(i,j)  = aion(i)
        zion_row(i,j)  = zion(i)
       enddo
      enddo

      call eosion



c..write out the results
      call pretty_eos_out('eosion:',ionam)

      stop 'normal termination'
      end   








      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      include 'implno.dek'
      include 'vector_eos.dek'



c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax

c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar


c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),abar,zbar,zbarxx,ytot1

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end








      subroutine pretty_eos_out(whose,ionam)
      include 'implno.dek'
      include 'vector_eos.dek'


c..writes a pretty output for the eos tester


c..declare
      integer     i,j,k,kend
      character*5 ionam(1)
      character*7 whose


c..local variables
      character*6 roman(30)
      data roman /'I     ', 'II    ', 'III   ', 'IV    ', 'V     ', 
     1            'VI    ', 'VII   ', 'VIII  ', 'IX    ', 'X     ', 
     2            'XI    ', 'XII   ', 'XIII  ', 'XIV   ', 'XV    ', 
     3            'XVI   ', 'XVII  ', 'XVIII ', 'XIX   ', 'XX    ',
     4            'XXI   ', 'XXII  ', 'XXIII ', 'XXIV  ', 'XXV   ',
     5            'XXVI  ', 'XXVII ', 'XXVIII', 'XXIX  ', 'XXX   '/ 





c..popular formats
01    format(1x,t2,a,t11,'total',t24,'ion',t34,'electron',
     1       t46,'positron',t58,'radiation',t70,'coulomb',
     2       t82,'ionization')
02    format(1x,t2,a,1p7e12.4)
03    format(1x,t2,a6,1pe12.4,t22,a6,1pe12.4,
     1         t42,a6,1pe12.4,t62,a6,1pe12.4)
06    format(1x,t2,a,a,1pe12.4,
     1          t30,a,a,1pe12.4,
     2          t58,a,a,1pe12.4)




      do j=jlo_eos,jhi_eos


c..the input 
      write(6,03) 'temp =',temp_row(1),'den  =',den_row(1),
     1            'abar =',abar_row(1),'zbar =',zbar_row(1)
      write(6,*) ' ' 


c..and the output
c..first the totals from each of the components
      write(6,01)  whose
      write(6,02) 'pres =',
     1            ptot_row(j),pion_row(j),pele_row(j),
     2            ppos_row(j),prad_row(j),pcou_row(j),
     3            pip_row(j)
      write(6,02) 'ener =',
     1            etot_row(j),eion_row(j),eele_row(j),
     2            epos_row(j),erad_row(j),ecou_row(j),
     3            eip_row(j)
      write(6,02) 'entr =',
     1            stot_row(j),sion_row(j),sele_row(j),
     2            spos_row(j),srad_row(j),scou_row(j),
     3            sip_row(j)


c..derivatives of the totals with respect to the input variables
      write(6,*)  ' '
      write(6,03) 'dp/dd=',dpd_row(j),'dp/dt=',dpt_row(j),
     1            'dp/da=',dpa_row(j),'dp/dz=',dpz_row(j)
      write(6,03) 'de/dd=',ded_row(j),'de/dt=',det_row(j),
     1            'de/da=',dea_row(j),'de/dz=',dez_row(j)
      write(6,03) 'ds/dd=',dsd_row(j),'ds/dt=',dst_row(j),
     1            'ds/da=',dsa_row(j),'ds/dz=',dsz_row(j)


c..derivatives of the electron-positron compoenets with
c..respect to the input variables
      write(6,*) ' ' 
      write(6,03) 'dpepd=',dpepd_row(j),'dpept=',dpept_row(j),
     1            'dpepa=',dpepa_row(j),'dpepz=',dpepz_row(j)
      write(6,03) 'deepd=',deepd_row(j),'deept=',deept_row(j),
     1            'deepa=',deepa_row(j),'deepz=',deepz_row(j)
      write(6,03) 'dsepd=',dsepd_row(j),'dsept=',dsept_row(j),
     1            'dsepa=',dsepa_row(j),'dsepz=',dsepz_row(j)


c..the thermodynamic consistency relations, these should all be
c..at the floating poiint limit of zero
      write(6,*) ' ' 
      write(6,03) 'maxw1=',dse_row(j),'maxw2=',dpe_row(j),
     1            'maxw3=',dsp_row(j)


c..number density of electrons, poistrons, matter electrons, and ions
      write(6,03) 'xne  =',xne_row(j),'xnp  =',xnp_row(j),
     1            'xnem =',xnem_row(j),'xni  =',xni_row(j)


c..derivatibves of the electron number density with 
c..respect to the input variables
      write(6,03) 'dxned=',dxned_row(j),'dxnet=',dxnet_row(j),
     1            'dxnea=',dxnea_row(j),'dxnez=',dxnez_row(j)


c..electron chemical potential, positron chemical potential
c..and derivatives of electron chemical potential with respect
c..to the input variables
      write(6,03) 'eta  =',etaele_row(j),'etap =',etapos_row(j)
      write(6,03) 'detad=',detad_row(j),'detat=',detat_row(j),
     1            'detaa=',detaa_row(j),'detaz=',detaz_row(j)


c..specific heats, and ratio of electostatic to thermal energy
      write(6,03) 'cp   =',cp_row(j),'cv   =',cv_row(j),
     1            'plasg=',plasg_row(j)

c..the 3 gammas and the sound speed
      write(6,03) 'gam1 =',gam1_row(j),'gam2 =',gam2_row(j),
     1            'gam3 =',gam3_row(j),'csond=',cs_row(j)



c..write out the ionization fractions
      write(6,*) ' '
      do i=1,niso
       if (xmass_row(i,j) .gt. 1.0e-16) then
        kend = int(zion_row(i,j)) + 1
        write(6,06) (ionam(i),roman(k),frac_row(k,i,j), k=1,kend)
        write(6,*)
       end if
      enddo
      write(6,*) ' '




      enddo
      return
      end










c..here is the timmes eos with ionization: eosion

c..routine eosion computes a stellar eos with ionziation
c..routine etages provides a good guess to the chemical potential
c..routine geteta is used by a root finder to find the chemical potential
c..routine getxne computes the electron/positron thermodynamics
c..routine zbrac_eosion brackets the chemical potential root
c..routine zbrent_eosion finds the chemical potential root
c..routine bigsaha computes the fraction of electrons in an ionization stage
c..routine ionpot returns the ionization potentials
c..routine stat_weight returns the statistical weights
c..routine coulomb_ion implments coulomb corrections





      subroutine eosion
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar (average weight) and zbar (average charge), 
c..this routine returns all the other thermodynamic quantities. 

c..of interest is the pressure [erg/cm**3], specific thermal energy [erg/gr], 
c..the entropy [erg/g/K], with their derivatives with respect to temperature, 
c..density, abar, and zbar.

c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.

c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. the full fermi-dirac integrals and their derivatives 
c..with respect to eta and beta are computed to machine precision, and
c..all other derivatives are analytic.

c..references: cox & giuli (c&g) chapter 24, 
c..            timmes & arnett, apj supp. 125, 277, 1999
c..            timmes & swesty, apj supp. 126, 501, 2000 





c..a dictionary of terms used:
c..this routine has now been pipelined. 
c..all the input and output variables are in the file vector_eos.dek.
c..the vector name is the scaler name appended with an "_row",
c..for example, temp_row(i), den_row(i), and so on. 



c..input:
c..temp     = temperature
c..den      = density
c..niso     = number of isotopes
c..xmass    = array of mass fractions
c..aion     = array of nuclear weights
c..zion     = array of nuclear charges


c..output:

c..pres     = total pressure
c..dpresdd  = derivative of total pressure with respect to density
c..dpresdt  = derivative of total pressure with respect to temperature

c..ener     = total internal energy
c..denerdd  = derivative of total energy with respect to density
c..denerdt  = derivative of total energy with respect to temperature

c..entr     = total entropy
c..dentrdd  = derivative of total entropy with respect to density
c..dentrdt  = derivative of total entropy with respect to temperature


c..prad     = radiation pressure
c..dpraddd  = derivative of the radiation pressure with density
c..dpraddt  = derivative of the radiation pressure with temperature

c..erad     = radiation energy
c..deraddd  = derivative of the radiation energy with density
c..deraddt  = derivative of the radiation energy with temperature

c..srad     = radiation entropy
c..dsraddd  = derivative of the radiation entropy with density
c..dsraddt  = derivative of the radiation entropy with temperature

c..radmult  = radiation multiplier (useful for turning radiation off/on)




c..xni      = number density of ions
c..dxnidd   = derivative of the ion number density with density
c..dxnidt   = derivative of the ion number density with temperature

c..pion     = ion pressure
c..dpiondd  = derivative of the ion pressure with density
c..dpiondt  = derivative of the ion pressure with temperature

c..eion     = ion energy
c..deiondd  = derivative of the ion energy with density
c..deiondt  = derivative of the ion energy with temperature

c..sion     = ion entropy
c..dsiondd  = derivative of the ion entropy with density
c..dsiondt  = derivative of the ion entropy with temperature

c..ionmult  = ion multiplier (useful for turning ions off/on)


c..etaele   = electron chemical potential
c..detadd   = derivative of the electron chem potential with density
c..detadt   = derivative of the electron chem potential with temperature

c..etapos   = positron degeneracy parameter


c..xne       = number density of electrons
c..dxnedd    = derivative of the electron number density with density
c..dxnedt    = derivative of the electron number density with temperature

c..xnefer    = fermi integral electron number density
c..dxneferdd = derivative of the fermi electron number density with density
c..dxneferdt = derivative of the fermi electron number density with temperature

c..xnpfer    = fermi integral positron number density
c..dxnpferdd = derivative of the fermi positron number density with density
c..dxnpferdt = derivative of the fermi positron number density with temperature

c..pele      = electron pressure
c..dpeledd   = derivative of the electron pressure with density
c..dpeledt   = derivative of the electron pressure with temperature

c..eele     = electron energy
c..deeledd   = derivative of the electron energy with density
c..deeledt   = derivative of the electron energy with temperature

c..sele     = electron entropy
c..dseledd   = derivative of the electron entropy with density
c..dseledt   = derivative of the electron entropy with temperature


c..ppos     = positron pressure
c..dpposdd   = derivative of the positron pressure with density
c..dpposdt   = derivative of the positron pressure with temperature

c..epos     = electron energy
c..deposdd   = derivative of the positron energy with density
c..deposdt   = derivative of the positron energy with temperature

c..spos     = electron entropy
c..dsposdd   = derivative of the positron entropy with density
c..dsposdt   = derivative of the positron entropy with temperature

c..pep      = electron + positron pressure
c..dpepdd   = derivative of the electron+positron pressure with density
c..dpepdt   = derivative of the electron+positron pressure with temperature

c..eep      = electron + ositron energy
c..deepdd   = derivative of the electron+positron energy with density
c..deepdt   = derivative of the electron+positron energy with temperature

c..sep      = electron + positron entropy
c..dsepdd   = derivative of the electron+positron entropy with density
c..dsepdt   = derivative of the electron+positron entropy with temperature

c..elemult  = electron multiplier (useful for turning e-e+ off/on)


c..eip      = ionization potential energy
c..deipdd   = derivative of ionization energy with density
c..deipdt   = derivative of ionization energy with temperature

c..pip      = ionization potential pressure
c..dpipdd   = derivative of ionization pressure with density
c..dpipdt   = derivative of ionization pressure with temperature


c..sip      = ionization potential entropy
c..dsipdd   = derivative of ionization entropy with density
c..dsipdt   = derivative of ionization entropy with temperature

c..potmult  = ionization energy multiplier (useful for turning ionization off)


c..pcoul    = coulomb pressure correction
c..coulmult = coulomb component multiplier
c..dpcouldd = derivative of the coulomb pressure with density
c..dpcouldt = derivative of the coulomb pressure with temperature

c..ecoul    = coulomb energy correction
c..decouldd = derivative of the coulomb energy with density
c..decouldt = derivative of the coulomb energy with temperature

c..scoul    = coulomb entropy correction
c..dscouldd = derivative of the coulomb entropy with density
c..dscouldt = derivative of the coulomb entropy with temperature


c..kt       = kerg * temperature
c..beta     = dimensionless ratio of kerg*temp/me*c^2

c..chit     = temperature exponent in the pressure equation of state
c..chid     = density exponent in the pressure equation of state
c..cv       = specific heat at constant volume
c..cp       = specific heat at constant pressure
c..gam1     = first adiabatic exponent
c..gam2     = second adiabatic exponent
c..gam3     = third adiabatic exponent
c..nabad    = adiabatic gradient
c..sound    = relativistically correct adiabatic sound speed
c..plasg    = ratio of electrostatic to thermal energy


c..dse      = thermodynamic consistency check de/dt = t*ds/dt
c..dpe      = thermodynamic consistency check p = d**2 de/dd + t*dpdt
c..dsp      = thermodynamic consistency check dp/dt = - d**2 ds/dd





c..declare the input
      double precision zbar,abar

      double precision temp,den,
     1                 xmass(irowmax),ymass(irowmax), 
     2                 aion(irowmax),
     3                 zion(irowmax)

      common /bkdoor/  temp,den,xmass,ymass,aion,zion


      double precision fracs(jstagemax,irowmax)


      double precision pres,ener,entr,
     1                 dpresdd,dpresdt,
     2                 denerdd,denerdt,
     3                 dentrdd,dentrdt


c..radiation
      integer          radmult
      double precision prad,erad,srad,
     1                 dpraddd,dpraddt,
     2                 deraddd,deraddt,
     3                 dsraddd,dsraddt


c..ions
      integer          ionmult
      double precision pion,eion,sion,xni,etaion,
     1                 dpiondd,dpiondt,
     2                 deiondd,deiondt,
     3                 dsiondd,dsiondt,
     4                 dxnidd,dxnidt,
     5                 detaiondd,detaiondt


c..electron-positrons
      integer          elemult
      double precision etaele,detadd,detadt,
     1                 etapos,zeff,
     2                 xne,dxnedd,dxnedt,
     3                 xnefer,dxneferdd,dxneferdt,
     4                 xnpfer,dxnpferdd,dxnpferdt,
     5                 pele,dpeledd,dpeledt,
     6                 ppos,dpposdd,dpposdt,
     7                 pep,dpepdd,dpepdt,
     8                 eele,deeledd,deeledt,
     9                 epos,deposdd,deposdt,
     1                 eep,deepdd,deepdt,
     2                 sele,dseledd,dseledt,
     3                 spos,dsposdd,dsposdt,
     4                 sep,dsepdd,dsepdt



c..ionization contributions
      integer          potmult
      double precision eip,deipdd,deipdt,
     1                 pip,dpipdd,dpipdt,
     2                 sip,dsipdd,dsipdt


c..coulomb corrections
      integer          coulmult
      double precision plasg,
     1                 pcoul,dpcouldd,dpcouldt,
     2                 ecoul,decouldd,decouldt,
     3                 scoul,dscouldd,dscouldt



c..various physical quantities based on derivatives 
      double precision chit,chid,cv,cp,gam1,gam2,gam3,nabad,sound


c..for the maxwell relations
      double precision dse,dpe,dsp


c..miscelaneous local variables 
      integer          i,j,jj,k,kend,niter
      double precision kt,ktinv,x,y,z,xx,yy,zz,ages,agesav,agesnew,
     1                 ratio,ytot1,f,df,deninv,tempinv,zbarxx

      external         geteta
      logical          success
      double precision zbrent_eosion,geteta,etalo,etahi


c..various derived constants
      double precision third,sifac,eostol,fpmin
      parameter        (third  = 1.0d0/3.0d0,
     1                  sifac  = 8.6322745944370191d-45,
     2                  eostol = 1.0d-13,
     3                  fpmin  = 1.0d-14)

c..note: sifac = h**3/(2.0d0*pi*amu)**1.5d0



c..popular format statements for debugging
01    format(1x,5(a,1pe24.16))
02    format(1x,5(a,1pe16.8))
03    format(1x,1p5e16.8)
04    format(1x,1p6e12.3)



c..set the on/off switches
      radmult  = 1
      ionmult  = 1
      elemult  = 1
      coulmult = 1
      potmult  = 1




c..start pipeline loop
      do j=jlo_eos,jhi_eos


c..load the input
       temp  = temp_row(j)
       den   = den_row(j)
       do i=1,niso
        xmass(i) = xmass_row(i,j)
        aion(i)  = aion_row(i,j)
        zion(i)  = zion_row(i,j)
       enddo

c..and compute the number fractions ymass, abar and zbar 
       zbarxx  = 0.0d0
       ytot1   = 0.0d0
       do i=1,niso
        ymass(i) = xmass(i)/aion(i)
        ytot1    = ytot1 + ymass(i)
        zbarxx   = zbarxx + zion(i) * ymass(i)
       enddo
       abar   = 1.0d0/ytot1
       zbar   = zbarxx * abar
       abar_row(j) = abar
       zbar_row(j) = zbar



c..initialize
       prad     = 0.0d0
       dpraddd  = 0.0d0
       dpraddt  = 0.0d0

       erad     = 0.0d0
       deraddd  = 0.0d0
       deraddt  = 0.0d0

       srad     = 0.0d0
       dsraddd  = 0.0d0
       dsraddt  = 0.0d0

       xni      = 0.0d0
       dxnidd   = 0.0d0
       dxnidt   = 0.0d0

       pion     = 0.0d0
       dpiondd  = 0.0d0
       dpiondt  = 0.0d0

       eion     = 0.0d0
       deiondd  = 0.0d0
       deiondt  = 0.0d0

       sion     = 0.0d0
       dsiondd  = 0.0d0
       dsiondt  = 0.0d0

       xne      = 0.0d0
       dxnedd   = 0.0d0
       dxnedt   = 0.0d0

       etaele   = 0.0d0
       detadd   = 0.0d0
       detadt   = 0.0d0
       etapos   = 0.0d0

       xnefer    = 0.0d0
       dxneferdd = 0.0d0
       dxneferdt = 0.0d0

       xnpfer    = 0.0d0
       dxnpferdd = 0.0d0
       dxnpferdt = 0.0d0

       pele     = 0.0d0
       dpeledd  = 0.0d0
       dpeledt  = 0.0d0

       eele     = 0.0d0
       deeledd  = 0.0d0
       deeledt  = 0.0d0

       sele     = 0.0d0
       dseledd  = 0.0d0
       dseledt  = 0.0d0

       ppos     = 0.0d0
       dpposdd  = 0.0d0
       dpposdt  = 0.0d0

       epos     = 0.0d0
       deposdd  = 0.0d0
       deposdt  = 0.0d0

       spos     = 0.0d0
       dsposdd  = 0.0d0
       dsposdt  = 0.0d0

       pep      = 0.0d0
       dpepdd   = 0.0d0
       dpepdt   = 0.0d0

       eep      = 0.0d0
       deepdd   = 0.0d0
       deepdt   = 0.0d0

       sep      = 0.0d0 
       dsepdd   = 0.0d0
       dsepdt   = 0.0d0

       do i=1,niso
        kend = int(zion(i)) + 1 
        do k=1,kend
         fracs(k,i) = 0.0d0
        enddo
       enddo


       eip      = 0.0d0
       deipdd   = 0.0d0
       deipdt   = 0.0d0


       pip      = 0.0d0
       dpipdd   = 0.0d0
       dpipdt   = 0.0d0

       sip      = 0.0d0
       dsipdd   = 0.0d0
       dsipdt   = 0.0d0

       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0

       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0

       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0


       kt      = kerg * temp
       ktinv   = 1.0d0/kt
       deninv  = 1.0d0/den
       tempinv = 1.0d0/temp




c..radiation section:
       if (radmult .ne. 0) then

c..pressure in erg/cm**3
        prad    = asol * third * temp * temp * temp * temp
        dpraddd = 0.0d0 
        dpraddt = 4.0d0 * prad/temp

c..energy in erg/gr
        erad    = 3.0d0 * prad * deninv
        deraddd = -erad * deninv
        deraddt = 3.0d0 * dpraddt * deninv

c..entropy in erg/g/kelvin
        srad    = (prad*deninv + erad) * tempinv
        dsraddd = (dpraddd*deninv - prad*deninv**2 + deraddd) * tempinv 
        dsraddt = (dpraddt*deninv + deraddt - srad)  * tempinv
       end if




c..ion section:

c..number density in 1/cm**3,
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnidt  = 0.0d0

       if (ionmult .ne. 0) then

c..pressure in erg/cm**3
        pion    = xni * kt 
        dpiondd = dxnidd * kt 
        dpiondt = xni * kerg 

c.. energy in erg/gr
        eion    = 1.5d0 * pion*deninv 
        deiondd = (1.5d0 * dpiondd - eion)*deninv 
        deiondt = 1.5d0 * dpiondt*deninv 

c..ion degeneracy parameter (c&g 9.60)
        y         = 1.0d0/(abar*kt)
        yy        = y * sqrt(y)
        z         = xni * sifac * yy
        etaion    = log(z)
        xx        = 1.0d0/xni
        detaiondd = dxnidd*xx
        detaiondt = dxnidt*xx - 1.5d0*tempinv

c..entropy in erg/gr/kelvin
c..the last term is the usual  etaion * kerg * xni/den 
c..sometimes called the sacker-tetrode equation

        sion    = (eion + pion*deninv)*tempinv - etaion * kerg*avo*ytot1

        dsiondd = (deiondd + dpiondd*deninv - pion*deninv**2)*tempinv 
     1            - detaiondd * kerg * avo*ytot1

        dsiondt = (deiondt + dpiondt*deninv)*tempinv 
     1            - (eion + pion*deninv)*tempinv**2
     2            - detaiondt * kerg * avo*ytot1
       end if







c..electron-positron section:
       if (elemult .ne. 0) then

c..make a good guess at the the electron degeneracy parameter eta
        call etages(xni,zbar,temp,ages)
        agesav = ages


c..bracket eta
        etalo = min(0.5d0*ages,2.0d0*ages)
        etahi = max(0.5d0*ages,2.0d0*ages)
        call zbrac_eosion(geteta,etalo,etahi,success)
        if (.not.success) then
         write(6,*) 'failed to bracket eta in routine eosion'
         write(6,01) 'temp  =',temp,' den =',den
         write(6,04) (xmass(i), i=1,niso) 
         stop 'failed to bracket eta'
        end if


c..get eta
        ages = zbrent_eosion(geteta,etalo,etahi,eostol,niter)



c..with the converged value of eta, get everything else

        call getxne(ages,den,temp,aion,zion,xmass,ymass,niso,
     1              etaele,detadd,detadt,
     2              etapos,zeff,
     3              xne,dxnedd,dxnedt,
     4              xnefer,dxneferdd,dxneferdt,
     5              xnpfer,dxnpferdd,dxnpferdt,
     6              pele,dpeledd,dpeledt,
     7              ppos,dpposdd,dpposdt,
     8              pep,dpepdd,dpepdt,
     9              eele,deeledd,deeledt,
     &              epos,deposdd,deposdt,
     1              eep,deepdd,deepdt,
     2              sele,dseledd,dseledt,
     3              spos,dsposdd,dsposdt,
     4              sep,dsepdd,dsepdt,
     5              fracs,jstagemax,irowmax,
     6              potmult,
     7              eip,deipdd,deipdt,
     8              pip,dpipdd,dpipdt,
     9              sip,dsipdd,dsipdt)

       end if




c..coulomb corretions section:
       if (coulmult .ne. 0) then

        call coulomb_ion(den,temp,abar,zbar,
     1                   pion,dpiondd,dpiondt,
     2                   xne,dxnedd,dxnedt,
     3                   plasg,
     4                   pcoul,dpcouldd,dpcouldt,
     5                   ecoul,decouldd,decouldt,
     6                   scoul,dscouldd,dscouldt)


c..bomb proof the cooulomb correctins
        x   = prad + pion + pele + ppos + pcoul
        if (x .le. 0.0) then

c         write(6,*) 
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*) 

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
        end if
       end if




c..sum all the components
       pres    = prad + pion + pele + ppos + pcoul + pip
       ener    = erad + eion + eele + epos + ecoul + eip
       entr    = srad + sion + sele + spos + scoul + sip

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd + dpipdd
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt + dpipdt

       denerdd = deraddd + deiondd + deepdd + decouldd + deipdd
       denerdt = deraddt + deiondt + deepdt + decouldt + deipdt

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd + dsipdd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt + dsipdt



c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98)
c..and relativistic formula for the sound speed (c&g 14.29)

       zz    = pres/den
       chit  = temp/pres * dpresdt
       chid  = dpresdd/zz
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + clight*clight)/zz
       sound = clight * sqrt(gam1/z)




c..maxwell relations; each is zero if the consistency is perfect
c..delicate subtraction in very degenerate regions causes roundoff error

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*den**2 + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*(den**2/dpresdt) - 1.0d0




c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd

        stot_row(j)   = entr 
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad 

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion 
        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = ppos
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd

        eele_row(j)   = eele
        epos_row(j)   = epos
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd

        sele_row(j)   = sele 
        spos_row(j)   = spos 
        dsept_row(j)  = dsepdt 
        dsepd_row(j)  = dsepdd 

        xnem_row(j)   = xne
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxneferdt + dxnpferdt
        dxned_row(j)  = dxneferdd + dxnpferdd
        xnp_row(j)    = xnpfer
        zeff_row(j)   = zeff

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        etapos_row(j) = etapos

        do i=1,niso
         kend = int(zion(i)) + 1
         do k=1,kend
          frac_row(k,i,j) = fracs(k,i)
         enddo
        enddo

        eip_row(j)    = eip
        pip_row(j)    = pip
        sip_row(j)    = sip

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul 
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound


c..for debugging 
c        crap1_row(j)   = pres
c        dcrap1d_row(j) = dpresdd
c        dcrap1t_row(j) = dpresdt


c..end of pipeline loop
      enddo
      return
      end





      double precision function geteta(aa)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

c..for a root find
      double precision aa


      double precision temp,den,
     1                 xmass(irowmax),ymass(irowmax), 
     2                 aion(irowmax),
     3                 zion(irowmax)

      common /bkdoor/  temp,den,xmass,ymass,aion,zion


c..locals
      integer          i,j,jend
      double precision kt,beta,beta12,beta32,etaele,
     1                 xne,xni,frac,dfrac_deta,
     3                 dfracdd,dfracdt,f12,f12eta,f12beta,
     4                 f32,f32eta,f32beta,zz,xnefer,
     5                 etapos,xnpfer


      double precision mecc,xconst     
      parameter       (mecc    = me * clight * clight,
     1                 xconst  = 2.4883740912221807d30)


c..some common factors
       kt      = kerg * temp
       beta    = kt/mecc
       beta12  = sqrt(beta)
       beta32  = beta * beta12


c..get the number density of free electrons 
       etaele    = aa
       xne       = 0.0d0
       do i=1,niso
        if (xmass(i) .gt. 1.0e-16) then
         xni    = avo * ymass(i) * den
         jend   = int(zion(i)) + 1
         do j=2,jend

          call bigsaha(j,jend,zion(i),etaele,temp,den,
     1                 frac,dfrac_deta,dfracdd,dfracdt)

          xne       = xne + xni * float(j-1) * frac

         enddo
        endif
       enddo



c..get the electron contribution
      call dfermi(0.5d0, etaele, beta, f12, f12eta, f12beta)
      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)
      zz            = xconst * beta32
      xnefer        = zz * (f12 + beta * f32)


c..if the temperature is not too low, get the positron contributions
c..chemical equilibrium means etaele + etapos = eta_photon = 0.
      etapos        = 0.0d0
      xnpfer        = 0.0d0
      if (beta .gt. 0.02) then
       etapos       = -aa - 2.0d0/beta
       call dfermi(0.5d0, etapos, beta, f12, f12eta, f12beta)
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       xnpfer        = zz * (f12 + beta * f32)
      end if


c..charge neutrality means ne_ionizat = ne_elect - ne_posit
      geteta  = xnefer - xnpfer - xne

      return
      end





      subroutine zbrac_eosion(func,x1,x2,succes)    
      include 'implno.dek' 
c..  
c..given a function func and an initial guessed range x1 to x2, the routine  
c..expands the range geometrically until a root is bracketed by the returned 
c..values x1 and x2 (in which case succes returns as .true.) or until the  
c..range becomes unacceptably large (in which succes returns as .false.). 
c..success  guaranteed for a function which has opposite sign for sufficiently 
c..large and small arguments.    
c.. 
c..declare  
      external          func 
      logical           succes   
      integer           ntry,j   
      parameter         (ntry=50) 
      double precision  func,x1,x2,factor,f1,f2  
      parameter         (factor = 1.6d0) 

      if (x1 .eq. x2) stop ' x1 = x2 in routine zbrac_eosion' 
      f1 = func(x1)  
      f2 = func(x2)  
      succes = .true.    
      do j=1,ntry 
       if (f1*f2 .le. 0.0) return    
       if (abs(f1) .lt. abs(f2)) then    
        x1 = x1 + factor * (x1-x2) 
        f1 = func(x1)    
       else  
        x2 = x2 + factor * (x2-x1) 
        f2 = func(x2)    
       end if    
      enddo
      succes = .false.   
      return 
      end    






      double precision function zbrent_eosion(func,x1,x2,tol,niter)   
      include 'implno.dek' 

c..using brent's method this routine finds the root of a function func  
c..between the limits x1 and x2. the root is when accuracy is less than tol. 
c.. 
c..note: eps the the machine floating point precision 
c.. 
c..declare 
      external          func 
      integer           niter,itmax,iter   
      parameter         (itmax = 100)   
      double precision  func,x1,x2,tol,a,b,c,d,e,fa, 
     1                  fb,fc,xm,tol1,p,q,r,s,eps    
      parameter         (eps=3.0d-15)   


c..initialize 
      niter = 0
      a     = x1 
      b     = x2 
      fa    = func(a)   
      fb    = func(b)   

c      if ( (fa .gt. 0.0  .and. fb .gt. 0.0)  .or. 
c     1     (fa .lt. 0.0  .and. fb .lt. 0.0)       ) then 
c       write(6,100) x1,fa,x2,fb 
c 100       format(1x,' x1=',1pe11.3,' f(x1)=',1pe11.3,/, 
c     1      1x,' x2=',1pe11.3,' f(x2)=',1pe11.3) 
c       stop 'root not bracketed in routine zbrent_eosion'    
c      end if 


      c  = b 
      fc = fb    

c..rename a,b,c and adjusting bound interval d   
      do iter =1,itmax  
       niter = niter + 1  
       if ( (fb .gt. 0.0  .and. fc .gt. 0.0)  .or. 
     1      (fb .lt. 0.0  .and. fc .lt. 0.0)      ) then 
        c  = a    
        fc = fa  
        d  = b-a  
        e  = d    
       end if    
       if (abs(fc) .lt. abs(fb)) then    
        a  = b    
        b  = c    
        c  = a    
        fa = fb  
        fb = fc  
        fc = fa  
       end if    
       tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol 
       xm   = 0.5d0 * (c-b)  

c..convergence check 
       if (abs(xm) .le. tol1 .or. fb .eq. 0.0) then  
        zbrent_eosion = b   
        return   
       end if    

c..attempt quadratic interpolation   
       if (abs(e) .ge. tol1 .and. abs(fa) .gt. abs(fb)) then 
        s = fb/fa    
        if (a .eq. c) then   
         p = 2.0d0 * xm * s    
         q = 1.0d0 - s  
        else 
         q = fa/fc   
         r = fb/fc   
         p = s * (2.0d0 * xm * q *(q-r) - (b-a)*(r - 1.0d0))   
         q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0) 
        end if   

c..check if in bounds    
        if (p .gt. 0.0) q = -q    
        p = abs(p)   

c..accept interpolation  
        if (2.0d0*p .lt. min(3.0d0*xm*q - abs(tol1*q),abs(e*q))) then    
         e = d   
         d = p/q 

c..or bisect
        else 
         d = xm  
         e = d   
        end if   

c..bounds decreasing to slowly use bisection 
       else  
        d = xm   
        e = d    
       end if    

c..move best guess to a  
       a  = b 
       fa = fb   
       if (abs(d) .gt. tol1) then    
        b = b + d    
       else  
        b = b + sign(tol1,xm)    
       end if    
       fb = func(b)  
      enddo
      stop 'too many iterations in routine zbrent_eosion'   
      end    


        


       



      subroutine getxne(aa,den,temp,aion,zion,xmass,ymass,niso,
     1                  etaele,detadd,detadt,
     2                  etapos,zeff,
     3                  xne,dxnedd,dxnedt,
     4                  xnefer,dxneferdd,dxneferdt,
     5                  xnpfer,dxnpferdd,dxnpferdt,
     6                  pele,dpeledd,dpeledt,
     7                  ppos,dpposdd,dpposdt,
     8                  pep,dpepdd,dpepdt,
     9                  eele,deeledd,deeledt,
     &                  epos,deposdd,deposdt,
     1                  eep,deepdd,deepdt,
     2                  sele,dseledd,dseledt,
     3                  spos,dsposdd,dsposdt,
     4                  sep,dsepdd,dsepdt,
     5                  fracs,jstagemax,irowmax,
     6                  potmult,
     7                  eip,deipdd,deipdt,
     8                  pip,dpipdd,dpipdt,
     9                  sip,dsipdd,dsipdt)

      include 'implno.dek'
      include 'const.dek'



c..input is temperature temp, density dem, average weight abar, 
c..average charge zbar, and the free electron degeneracy parameter aa.

c..everything else is output. see the calling routine for the definitions
c..of the variables 



c..declare the pass
      integer          niso,jstagemax,irowmax,potmult
      double precision aa,den,temp,aion(niso),zion(niso),
     1                 xmass(niso),ymass(niso),
     2                 etaele,detadd,detadt,
     3                 etapos,zeff,
     4                 xne,dxnedd,dxnedt,
     5                 xnefer,dxneferdd,dxneferdt,
     6                 xnpfer,dxnpferdd,dxnpferdt,
     7                 pele,dpeledd,dpeledt,
     8                 ppos,dpposdd,dpposdt,
     9                 pep,dpepdd,dpepdt,
     &                 eele,deeledd,deeledt,
     1                 epos,deposdd,deposdt,
     2                 eep,deepdd,deepdt,
     3                 sele,dseledd,dseledt,
     4                 spos,dsposdd,dsposdt,
     5                 sep,dsepdd,dsepdt,
     6                 fracs(jstagemax,irowmax),
     7                 eip,deipdd,deipdt,
     8                 pip,dpipdd,dpipdt,
     9                 sip,dsipdd,dsipdt




c..local variables
      integer          i,j,jend
      double precision kt,beta,beta12,beta32,beta52,
     1                 xni,dxnidd,dxnidt,
     2                 f12,f12eta,f12beta,
     3                 f32,f32eta,f32beta,
     4                 f52,f52eta,f52beta,
     5                 dzeff_deta,dzeffdd,dzeffdt,
     6                 dxne_deta,dxnefer_deta,dxnefer_dbeta,
     7                 detap_deta,detap_dbeta,
     8                 dxnpfer_detap,dxnpfer_deta,dxnpfer_dbeta,
     9                 dxep_deta,dxep_dbeta,
     &                 dpele_deta,dpele_dbeta,deele_deta,deele_dbeta,
     1                 dppos_detap,dppos_deta,dppos_dbeta,
     2                 depos_detap,depos_deta,depos_dbeta,
     3                 dsfac_deta,ytot1,zz,y,yy,
     4                 frac,dfrac_deta,dfracdd,dfracdt,
     5                 dxne_dd,dxne_dt

      double precision deip_deta,deip_dd,deip_dt,
     1                 xip,ionpot



      double precision xconst,pconst,econst,mecc,dbetadt
      parameter        (xconst  = 2.4883740912221807d30,
     1                  pconst  = 1.3581730208282635d24,
     2                  econst  = 2.0372595312423953d24, 
     3                  mecc    = me * clight * clight,
     4                  dbetadt = kerg/mecc)


c..note: 
c..xconst = 8.0d0 * pi * sqrt(2.0d0) * (me/h)**3 * c**3
c..pconst = xconst * 2.0d0/3.0d0 * me * clight**2
c..econst = xconst * me * clight**2




c..some common factors
       kt      = kerg * temp
       beta    = kt/mecc
       beta12  = sqrt(beta)
       beta32  = beta * beta12
       etaele  = aa


c..get the number density of free electrons 
       xne       = 0.0d0
       dxne_deta = 0.0d0
       dxne_dd   = 0.0d0
       dxne_dt   = 0.0d0
       zeff      = 0.0d0
       eip       = 0.0d0
       deip_deta = 0.0d0
       deip_dd   = 0.0d0
       deip_dt   = 0.0d0


       do i=1,niso
        if (xmass(i) .gt. 1.0e-16) then
         xni    = avo * ymass(i) * den
         dxnidd = avo * ymass(i)
         jend   = int(zion(i)) + 1

c..fraction of neutral species, while we are here
         j = 1
          call bigsaha(j,jend,zion(i),etaele,temp,den,
     1                 frac,dfrac_deta,dfracdd,dfracdt)

         fracs(j,i) = frac


c..sweep over all ionization stages
         do j=2,jend
          call bigsaha(j,jend,zion(i),etaele,temp,den,
     1                 frac,dfrac_deta,dfracdd,dfracdt)

          xne       = xne + xni * float(j-1) * frac

          dxne_deta = dxne_deta + xni * float(j-1)*dfrac_deta

          dxne_dd   = dxne_dd   + dxnidd*float(j-1)*frac
     1                          + xni * float(j-1)*dfracdd

          dxne_dt   = dxne_dt   + xni * float(j-1)*dfracdt

          zeff      = zeff + float(j-1) * frac 

          fracs(j,i) = frac           


          xip = ionpot(zion(i),j-1) 

          eip = eip + xip * xni * float(j-1) * frac

          deip_deta = deip_deta + xip * xni * float(j-1) * dfrac_deta 

          deip_dd   = deip_dd + xip * dxnidd * float(j-1) * frac
     1                        + xip * xni * float(j-1) * dfracdd

          deip_dt   = deip_dt + xip * xni * float(j-1) * dfracdt


         enddo
        endif
       enddo




c..get the electron contribution
      call dfermi(0.5d0, etaele, beta, f12, f12eta, f12beta)
      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)

      zz            = xconst * beta32
      yy            = f12 + beta * f32
      xnefer        = zz * yy
      dxnefer_deta  = zz * (f12eta + beta * f32eta)
      dxnefer_dbeta = xconst * beta12 * (1.5d0 * yy 
     1                       +  beta * (f12beta + f32 + beta * f32beta))



c..if the temperature is not too low, get the positron contributions
c..chemical equilibrium means etaele + etapos = eta_photon = 0.
      etapos        = 0.0d0
      detap_deta    = 0.0d0
      detap_dbeta   = 0.0d0
      xnpfer        = 0.0d0
      dxnpfer_detap = 0.0d0
      dxnpfer_dbeta = 0.0d0

      if (beta .gt. 0.02) then
       etapos      = -aa - 2.0d0/beta
       detap_deta  = -1.0d0
       detap_dbeta = 2.0d0/beta**2
       call dfermi(0.5d0, etapos, beta, f12, f12eta, f12beta)
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       xnpfer        = zz * (f12 + beta * f32)
       dxnpfer_detap = zz * (f12eta + beta * f32eta)
       dxnpfer_dbeta = xconst * beta12 * (1.5d0 * (f12 + beta * f32)
     1                 +  beta * (f12beta + f32 + beta * f32beta))
      end if



c..all the derivatives are in terms of eta and beta. 
c..we want to convert to temperature, density, abar and zbar derivatives.
c..so, after the root find above on eta we have: xne = xnefer - xnpfer
c..taking the derivative of this and solving for the unknown eta derivatives
c..leads to these expressions:


      dxnpfer_deta  = dxnpfer_detap * detap_deta
      dxnpfer_dbeta = dxnpfer_dbeta + dxnpfer_detap * detap_dbeta
      dxep_deta     = dxnefer_deta  - dxnpfer_deta
      dxep_dbeta    = dxnefer_dbeta - dxnpfer_dbeta


      y      = 1.0d0/(dxep_deta - dxne_deta)

      detadd = dxne_dd  * y
      detadt = (dxne_dt - dxep_dbeta*dbetadt) * y



c..derivatives of the electron number density
      dxnedd = dxne_dd + dxne_deta*detadd
      dxnedt = dxne_dt + dxne_deta*detadt


c..derivatives of the ionization energy
      deipdd = deip_dd + deip_deta*detadd
      deipdt = deip_dt + deip_deta*detadt


c..derivatives of the fermi integral electron number densities
      dxneferdd = dxnefer_deta * detadd
      dxneferdt = dxnefer_deta * detadt + dxnefer_dbeta * dbetadt


c..derivatives of the fermi integral positron number densities
      dxnpferdd = dxnpfer_deta * detadd 
      dxnpferdt = dxnpfer_deta * detadt + dxnpfer_dbeta * dbetadt





c..now get the pressure and energy
c..for the electrons

      beta52  = beta * beta32
      yy      = pconst * beta52
      zz      = econst * beta52

      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)
      call dfermi(2.5d0, etaele, beta, f52, f52eta, f52beta)

      pele        = yy * (f32 + 0.5d0 * beta * f52)
      dpele_deta  = yy * (f32eta + 0.5d0 * beta * f52eta)
      dpele_dbeta = pconst * beta32 * (2.5d0 * (f32 + 0.5d0*beta* f52)
     1              + beta* (f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta))

       
      eele        = zz * (f32 + beta * f52)
      deele_deta  = zz * (f32eta + beta * f52eta) 
      deele_dbeta = econst * beta32 * (2.5d0 * (f32 + beta * f52)
     1              + beta * (f32beta + f52 + beta * f52beta)) 


c..for the positrons
      ppos        = 0.0d0
      dppos_detap = 0.0d0
      dppos_dbeta = 0.0d0
      epos        = 0.0d0
      depos_detap = 0.0d0
      depos_dbeta = 0.0d0

      if (beta .gt. 0.02) then
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       call dfermi(2.5d0, etapos, beta, f52, f52eta, f52beta)

       ppos        = yy * (f32 + 0.5d0 * beta * f52)
       dppos_detap = yy * (f32eta + 0.5d0*beta *f52eta)
       dppos_dbeta = pconst * beta32 * (2.5d0 * (f32 + 0.5d0*beta*f52)
     1               + beta*(f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta))

       epos        = zz * (f32 + beta * f52) 
       depos_detap = zz * (f32eta + beta * f52eta) 
       depos_dbeta = econst * beta32 * (2.5d0 * (f32 + beta * f52)
     1               + beta * (f32beta + f52 + beta * f52beta)) 
      end if


c..derivatives of the electron pressure
      dpeledd = dpele_deta * detadd
      dpeledt = dpele_deta * detadt + dpele_dbeta * dbetadt


c..derivatives of the electron energy
      deeledd = deele_deta * detadd
      deeledt = deele_deta * detadt + deele_dbeta * dbetadt



c..derivatives of the positron pressure
      dppos_deta  = dppos_detap * detap_deta
      dppos_dbeta = dppos_dbeta + dppos_detap * detap_dbeta
      dpposdd     = dppos_deta * detadd
      dpposdt     = dppos_deta * detadt + dppos_dbeta * dbetadt


c..derivatives of the positron energy
      depos_deta  = depos_detap * detap_deta
      depos_dbeta = depos_dbeta + depos_detap * detap_dbeta
      deposdd     = depos_deta * detadd
      deposdt     = depos_deta * detadt + depos_dbeta * dbetadt


       
c..electron+positron pressure and its derivatives
c..note: at high temperatures and low densities, dpepdd is very small
c..and can go negative, so limit it to be positive definite
      pep    = pele    + ppos
      dpepdd = max(dpeledd + dpposdd, 1.0d-30) 
      dpepdt = dpeledt + dpposdt 


c..electron+positron thermal energy and its derivatives
      eep    = eele    + epos
      deepdd = deeledd + deposdd 
      deepdt = deeledt + deposdt 




c..electron entropy in erg/gr/kelvin and its derivatives
      y       = kerg/den 

      sele    = ((pele + eele)/kt - etaele*xnefer) * y

      dseledd = ((dpeledd + deeledd)/kt 
     1            - detadd*xnefer)*y 
     2            - etaele*dxneferdd*y 
     3            - sele/den

      dseledt = ((dpeledt + deeledt)/kt 
     1             - detadt*xnefer 
     2             - etaele*dxneferdt  
     3             - (pele + eele)/(kt*temp))*y


c..positron entropy in erg/gr/kelvin and its derivatives
      spos    = ((ppos + epos)/kt - etapos*xnpfer) * y    

      dsposdd = ((dpposdd + deposdd)/kt 
     1           - detap_deta*detadd*xnpfer 
     2           - etapos*dxnpferdd)*y - spos/den

      dsposdt = ((dpposdt + deposdt)/kt 
     1           - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer  
     2           - etapos*dxnpferdt  
     3           - (ppos + epos)/(kt*temp))*y


c..and their sum
      sep     = sele + spos
      dsepdd  = dseledd + dsposdd
      dsepdt  = dseledt + dsposdt



c..adjust for the rest mass energy of the positrons
      y       = 2.0d0 * mecc
      epos    = epos    + y * xnpfer 
      deposdd = deposdd + y * dxnpferdd
      deposdt = deposdt + y * dxnpferdt


c..and resum
      deepdd = deeledd + deposdd 
      deepdt = deeledt + deposdt 


c..convert the electron-positron thermal energy in erg/cm**3 to 
c..a specific thermal energy in erg/gr

      eele    = eele/den
      deeledd = deeledd/den - eele/den
      deeledt = deeledt/den 

      epos    = epos/den
      deposdd = deposdd/den - epos/den
      deposdt = deposdt/den 


c..and resum
      deepdd = deeledd + deposdd 
      deepdt = deeledt + deposdt 





c..and take care of the ionization potential contributions
      if (potmult .eq. 0) then
       eip    = 0.0d0
       deipdd = 0.0d0
       deipdt = 0.0d0

       pip    = 0.0d0
       dpipdd = 0.0d0
       dpipdt = 0.0d0

       sip    = 0.0d0
       dsipdd = 0.0d0
       dsipdt = 0.0d0
      else


       pip    = 0.0d0
       dpipdd = 0.0d0
       dpipdt = 0.0d0

c..the ionization entropy in erg/gr/kelvin and its derivatives
       y      = kerg/den 

c       sip    = (eip/kt - etaele*xne) * y

c       dsipdd = (deipdd/kt 
c     1            - detadd*xne)*y 
c     2            - etaele*dxnedd*y 
c     3            - sip/den

c       dsipdt = (deipdt/kt 
c     1             - detadt*xne 
c     2             - etaele*dxnedt  
c     3             - eip/(kt*temp))*y


       sip    = eip/kt * y

       dsipdd = deipdd/kt*y - sip/den

       dsipdt = (deipdt/kt - eip/(kt*temp))*y


c..convert the ionization energy from erg/cm**3 to  erg/gr

       eip    = eip/den
       deipdd = deipdd/den - eip/den
       deipdt = deipdt/den 

      end if

      return
      end










      subroutine bigsaha(stage,maxstages,zion,eta,temp,den,
     1                   frac,dfrac_deta,dfracdd,dfracdt)

      include 'implno.dek'
      include 'const.dek'


c..this routine evaluates the fraction of atoms of species k
c..in ionization stage j relative to the total number of atoms of that species.
c..mihalas equation 5.16 and 5.17

c..declare the pass
      integer          stage,maxstages
      double precision zion,eta,temp,den,frac,dfrac_deta,
     1                 dfracdd,dfracdt


c..local variables
      integer          i,j,ijend
      double precision z0,z1,numer,denom,sp,chisum,
     1                 ww,yy,big,bigger,
     2                 kt,ktinv,stat_weight,ionpot,saha,chifac,denion,
     3                 denfac,tempinv,deninv,denioninv,
     4                 dnumer_deta,dnumer_dd,dnumer_dt,
     5                 ddenom_deta,ddenom_dd,ddenom_dt


c..the definition of what is a big number, without overflowing, depends on 
c..how many stages of ionization there are, which in turn depends on zion. 
c..here i set big such that if the atom is completely neutral, 
c..then multiplication of all the saha factors will yield a number 
c..no larger than bigger.

      big         = 200.0d0/zion
      bigger      = 10.0d0**big


c..continue initializing
      kt          = kerg * temp
      ktinv       = 1.0d0/kt 
      denion      = 0.1d0 
      denioninv   = 1.0d0/denion
      denfac      = den * denioninv
      tempinv     = 1.0d0/temp
      deninv      = 1.0d0/den 
      ijend       = maxstages - 1

      numer       = 1.0d0
      dnumer_deta = 0.0d0
      dnumer_dd   = 0.0d0
      dnumer_dt   = 0.0d0

      denom       = 1.0d0
      ddenom_deta = 0.0d0
      ddenom_dd   = 0.0d0
      ddenom_dt   = 0.0d0




c..the product numerator and summed product denominator 
      do j=1,ijend
       sp       = 1.0d0 
       chisum   = 0.0d0

       do i=j,ijend

        chifac = ionpot(zion,i) * ktinv

c..this is the proposed exponent for the exponential
        yy     = chifac + eta - denfac 


c..atom is unionized
        if (yy .gt. 200.0) then
         saha       = bigger

c..atom is completely ionized
        else if (yy .lt. -200.0) then
         saha       = 0.0d0

c..atom is well within generous bounds of being partially ionized
        else
         z0   = stat_weight(zion,i)
         z1   = stat_weight(zion,i+1)
         ww   = min(big, yy)
         saha = z0/z1 * exp(ww)
        end if 


c..uncommenting these three lines will enforce complete ionization
c        denioninv = 0.0d0
c        chifac    = 0.0d0
c        saha      = 0.0d0




c..one saha times another saha times another saha ... etc 
        sp = sp * saha


c..sum of ionization potentials for the temperature derivatives
        chisum = chisum  + chifac
   
       enddo



c..pick off the numerator, and form
c..its derivatives with eta, density and temperature

       if (j .eq. stage) then
        numer       = sp
        yy          = float(maxstages - j) * sp
        dnumer_deta = yy
        dnumer_dd   = -yy * denioninv
        dnumer_dt   = -sp * chisum * tempinv
       end if 


c..keep summing the denominator, and form
c..its derivatives with eta, density and temperature

       denom       = denom + sp
       yy          = float(maxstages - j) * sp 
       ddenom_deta = ddenom_deta + yy
       ddenom_dd   = ddenom_dd   - yy * denioninv
       ddenom_dt   = ddenom_dt   - sp * chisum * tempinv
      enddo


c..form the fraction in ionization stage j, and 
c..its derivatives with eta, density and temperature
c..write the derivatives this way to minimize roundoff error

      frac       = numer/denom
      dfrac_deta = dnumer_deta/denom - frac * ddenom_deta/denom
      dfracdd    = dnumer_dd/denom   - frac * ddenom_dd/denom
      dfracdt    = dnumer_dt/denom   - frac * ddenom_dt/denom

      return
      end




      double precision function ionpot(zion,stage)
      include 'implno.dek'
      include 'const.dek'

c..returns the ionization potential in erg

c..this routine includes all the ionization stages of the elements 
c..from hydrogen to copper. the first 21 stages of zinc are also included. 
c..all data is from the handbook of chemistry and physics 2001.

c..input is the element charge (a convenient index) and the ionization stage.


c..declare the pass
      integer          stage
      double precision zion
      
c..locals
      integer          j,k, ifirst,izion 
      double precision xip(30,30),off_table
      parameter        (off_table = 1.0e30)



c..first time flag
      data    ifirst/0/
   
c..hydrogen
      data (xip(1,k), k=1,30)/
     1      13.59844, 
     2      29*off_table/

c..helium
      data (xip(2,k), k=1,30)/
     1      24.58741, 54.41778, 
     2      28*off_table/

c..lithium
      data (xip(3,k), k=1,30)/
     1      5.39172, 75.64018, 122.45429, 
     2      27*off_table/

c..berrylium
      data (xip(4,k), k=1,30)/
     1      9.3227, 18.21116, 153.89661, 217.71865, 
     2      26*off_table/

c..boron
      data (xip(5,k), k=1,30)/
     1      8.29803, 25.15484, 37.93064, 259.37521, 340.22580,
     2      25*off_table/

c..carbon
      data (xip(6,k), k=1,30)/
     1      11.26030, 24.38332, 47.8878, 64.4939, 392.087, 489.99334,   
     2      24*off_table/

c..nitrogen
      data (xip(7,k), k=1,30)/
     1      14.53414, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 
     2      667.046,
     3      23*off_table/

c..oxygen
      data (xip(8,k), k=1,30)/
     1      13.61806, 35.11730, 54.9355, 77.41353, 113.8990, 138.1197, 
     2      739.29, 871.4101,
     3      22*off_table/

c..fluorine
      data (xip(9,k), k=1,30)/
     1      17.42282, 34.97082, 62.7084, 87.1398, 114.2428, 157.1651, 
     2      185.186, 953.9112, 1103.1176,
     3      21*off_table/


c..neon
      data (xip(10,k), k=1,30)/
     1      21.5646, 40.96328, 63.45, 97.12, 126.21, 157.93, 207.2759, 
     2      239.0989, 1195.8286, 1362.1995,
     3      20*off_table/


c..sodium
      data (xip(11,k), k=1,30)/
     1      5.13908, 47.2864, 71.6200, 98.91, 138.40, 172.18, 208.50, 
     2      264.25, 299.864, 1465.121, 1648.702,
     3      19*off_table/


c..magnesium
      data (xip(12,k), k=1,30)/
     1      7.64624, 15.03528, 80.1437, 109.2655, 141.27, 186.76, 
     2      225.02, 265.96, 328.06, 367.50, 1761.805, 1962.6650,
     3      18*off_table/


c..aluminum
      data (xip(13,k), k=1,30)/
     1      5.98577, 18.82856, 28.44765, 119.992, 153.825, 190.49, 
     2      241.76, 284.66, 330.13, 398.75, 442.00, 2085.98, 
     3      2304.1410,
     4      17*off_table/


c..silicon 
      data (xip(14,k), k=1,30)/
     1      8.15169, 16.34585, 33.49302, 45.14181, 166.767, 205.27, 
     2      246.5, 303.54, 351.12, 401.37, 476.36, 523.42, 2437.63, 
     3      2673.182,
     4      16*off_table/


c..phosphorous
      data (xip(15,k), k=1,30)/
     1      10.48669, 19.7694, 30.2027, 51.4439, 65.0251, 220.421, 
     2      263.57, 309.60, 372.13, 424.4, 479.46, 560.8, 611.74, 
     3      2816.91, 3069.842,
     4      15*off_table/


c..sulfer
      data (xip(16,k), k=1,30)/
     1      10.36001, 23.3379, 34.79, 47.222, 72.5945, 88.0530, 280.948, 
     2      328.75, 379.55, 447.5, 504.8, 564.44, 652.2, 707.01, 
     3      3223.78, 3494.1892,
     4      14*off_table/


c..chlorine
      data (xip(17,k), k=1,30)/
     1      12.96764, 23.814, 39.61, 53.4652, 67.8, 97.03, 114.1958, 
     2      348.28, 400.06, 455.63, 529.28, 591.99, 656.71, 749.76, 
     3      809.40, 3658.521, 3946.2960,
     4      13*off_table/


c..argone
      data (xip(18,k), k=1,30)/
     1      15.75962, 27.62967, 40.74, 59.81, 75.02, 91.009, 124.323, 
     2      143.460, 422.45, 478.69, 538.96, 618.26, 686.10, 755.74, 
     3      854.77, 918.03, 4120.8857, 4426.2296,
     4      12*off_table/


c..pottasium
      data (xip(19,k), k=1,30)/
     1      4.34066, 31.63, 45.806, 60.91, 82.66, 99.4, 117.56, 154.88, 
     2      175.8174, 503.8, 564.7, 629.4, 714.6, 786.6, 861.1, 968.0, 
     3      1033.4, 4610.8, 4934.046,
     4      11*off_table/


c..calcium
      data (xip(20,k), k=1,30)/
     1      6.11316, 11.87172, 50.9131, 67.27, 84.50, 108.78, 127.2, 
     2      147.24, 188.54, 211.275, 591.9, 657.2, 726.6, 817.6, 894.5, 
     3      974.0, 1087.0, 1157.8, 5128.8, 5469.864,
     4      10*off_table/


c..scandium
      data (xip(21,k), k=1,30)/
     1      6.5615, 12.79967, 24.75666, 73.4894, 91.65, 110.68, 138.0, 
     2      158.1, 180.03, 225.18, 249.798, 687.36, 756.7, 830.8, 927.5, 
     3      1009.0, 1094.0, 1213.0, 1287.97, 5674.8, 6033.712,
     4      9*off_table/


c..titanium
      data (xip(22,k), k=1,30)/
     1      6.8281, 13.5755, 27.4917, 43.2672, 99.30, 119.53, 140.8, 
     2      170.4, 192.1, 215.92, 265.07, 291.500, 787.84, 863.1, 941.9, 
     3      1044.0, 1131.0, 1221.0, 1346.0, 1425.4, 6249.0, 6625.82,
     4      8*off_table/


c..vanadium
      data (xip(23,k), k=1,30)/
     1      6.7463, 14.66, 29.311, 46.709, 65.2817, 128.13, 150.6, 
     2      173.4, 205.8, 230.5, 255.7, 308.1, 336.277, 896.0, 976.0, 
     3      1060.0, 1168.0, 1260.0, 1355.0, 1486.0, 1569.6, 6851.3, 
     4      7246.12,
     5      7*off_table/


c..chromium
      data (xip(24,k), k=1,30)/
     1      6.7665, 16.4857, 30.96, 49.16, 69.46, 90.6349, 160.18, 
     2      184.7, 209.3, 244.4, 270.8, 298.0, 354.8, 384.168, 1010.6,  
     3      1097.0, 1185.0, 1299.0, 1396.0, 1496.0, 1634.0, 1721.4, 
     4      7481.7, 7894.81,
     5      6*off_table/


c..manganese
      data (xip(25,k), k=1,30)/
     1      7.43402, 15.63999, 33.668, 51.2, 72.4, 95.6, 119.203, 194.5, 
     2      221.8, 248.3, 286.0, 314.4, 343.6, 403.0, 435.163, 1134.7, 
     3      1224.0, 1317.0, 1437.0, 1539.0, 1644.0, 1788.0, 1879.9, 
     4      8140.6, 8571.94,
     5      5*off_table/


c..iron
      data (xip(26,k), k=1,30)/
     1      7.9024, 16.1878, 30.652, 54.8, 75.0, 99.1, 124.98, 151.06, 
     2      233.6, 262.1, 290.2, 330.8, 361.0, 392.2, 457, 489.256, 
     3      1266.0, 1358.0, 1456.0, 1582.0, 1689.0, 1799.0, 1950.0, 
     4      2023.0, 8828.0, 9277.69,
     5      4*off_table/


c..cobalt
      data (xip(27,k), k=1,30)/
     1      7.8810, 17.083, 33.50, 51.3, 79.5, 102.0, 128.9, 157.8, 
     2      186.13, 275.4, 305.0, 336.0, 379.0, 411.0, 444.0, 511.96, 
     3      546.58, 1397.2, 1504.6, 1603.0, 1735.0, 1846.0, 1962.0, 
     4      2119.0, 2219.0, 9544.1, 10012.12,
     5      3*off_table/


c..nickel
      data (xip(28,k), k=1,30)/
     1      7.6398, 18.16884, 35.19, 54.9, 76.06, 108, 133.0, 162.0, 
     2      193.0, 224.6, 321.0, 352.0, 384.0, 430.0, 464.0, 499.0, 
     3      571.08, 607.06, 1541.0, 1648.0, 1756.0, 1894.0, 2011.0, 
     4      2131.0, 2295.0, 2399.2, 10288.8, 10775.40,
     5      2*off_table/


c..copper
      data (xip(29,k), k=1,30)/
     1      7.72638, 20.29240, 36.841, 57.38, 79.8, 103.0, 139.0, 166.0, 
     2      199.0, 232.0, 265.3, 369.0, 401.0, 435.0, 484.0, 520.0, 
     3      557.0, 633.0, 670.588, 1697.0, 1804.0, 1916.0, 2060.0, 
     4      2182.0, 2308.0, 2478.0, 2587.5, 11062.38, 11567.617,
     5      1*off_table/


c..zinc
c..the last nine are not in the cited reference
      data (xip(30,k), k=1,30)/
     1      9.3942, 17.96440, 39.723, 59.4, 82.6, 108.0, 134.0, 174.0, 
     2      203.0, 238.0, 274.0, 310.8, 419.7, 454.0, 490.0, 542.0, 
     3      579.0, 619.0, 698.0, 738.0, 1856.0, 1920.0, 2075.0, 
     4      2190.0, 2350.0, 2500.0, 2600.0, 11080.0, 11600.0,
     5      12000.0/       


c..convert ev to ergs
      if (ifirst .eq. 0) then
       ifirst = 1
       do k=1,30
        do j=1,30
         if (xip(j,k) .ne. off_table) xip(j,k) = xip(j,k) * ev2erg
        enddo
       enddo
      endif


c..this makes the output simple enough
      izion = int(zion)
      ionpot = xip(izion,stage)

      return
      end







      double precision function stat_weight(zion,stage)
      include 'implno.dek'
      include 'const.dek'

c..returns the statistical weight of the ionization state 

c..this routine includes most of the weights from hydrogen to nickel.
c..i am unsure of scandium, vanadium, manganese, and cobalt, so i've
c..set all those to unity. 


c..input is the element charge zion (a convenient index) and the ionization stage.


c..declare the pass
      integer          stage
      double precision zion
      
c..locals
      integer          j,k,kmax,kmaxp1,ifirst,izion 
      parameter        (kmax = 30, kmaxp1 = kmax + 1)
      double precision wgt(kmax,kmaxp1),off_table
      parameter        (off_table = 1.0e30)


c..first time flag
      data    ifirst/0/
   
c..hydrogen
      data (wgt(1,k), k=1,kmaxp1)/
     1      2.0, 1.0,      
     2      29*off_table/

c..helium
      data (wgt(2,k), k=1,kmaxp1)/
     1      1.0, 2.0, 1.0,      
     2      28*off_table/

c..lithium
      data (wgt(3,k), k=1,kmaxp1)/
     1      2.0, 1.0, 2.0, 1.0,
     2      27*off_table/

c..berrylium
      data (wgt(4,k), k=1,kmaxp1)/
     1      1.0, 2.0, 1.0, 2.0, 1.0,
     2      26*off_table/

c..boron
      data (wgt(5,k), k=1,kmaxp1)/
     1      6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     2      25*off_table/

c..carbon
      data (wgt(6,k), k=1,kmaxp1)/
     1      9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     2      24*off_table/

c..nitrogen
      data (wgt(7,k), k=1,kmaxp1)/
     1      4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      23*off_table/

c..oxygen
      data (wgt(8,k), k=1,kmaxp1)/
     1      9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      22*off_table/

c..fluorine
      data (wgt(9,k), k=1,kmaxp1)/
     1      6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      21*off_table/


c..neon
      data (wgt(10,k), k=1,kmaxp1)/
     1      1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      20*off_table/


c..sodium
      data (wgt(11,k), k=1,kmaxp1)/
     1      2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      19*off_table/


c..magnesium
      data (wgt(12,k), k=1,kmaxp1)/
     1      1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 
     2      1.0,
     3      18*off_table/


c..aluminum
      data (wgt(13,k), k=1,kmaxp1)/
     1      6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 
     2      2.0, 1.0,
     3      17*off_table/


c..silicon 
      data (wgt(14,k), k=1,kmaxp1)/
     1      9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 
     2      1.0, 2.0, 1.0,
     3      16*off_table/


c..phosphorous
      data (wgt(15,k), k=1,kmaxp1)/
     1      4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 
     2      2.0, 1.0, 2.0, 1.0,
     3      15*off_table/


c..sulfer
      data (wgt(16,k), k=1,kmaxp1)/
     1      9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 
     2      1.0, 2.0, 1.0, 2.0, 1.0,
     3      14*off_table/


c..chlorine
      data (wgt(17,k), k=1,kmaxp1)/
     1      6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 
     2      6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      13*off_table/


c..argon
      data (wgt(18,k), k=1,kmaxp1)/
     1      1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 
     2      9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      12*off_table/


c..pottasium
      data (wgt(19,k), k=1,kmaxp1)/
     1      2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 
     2      4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      11*off_table/


c..calcium
      data (wgt(20,k), k=1,kmaxp1)/
     1      1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 
     2      9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     3      10*off_table/


c..scandium
      data (wgt(21,k), k=1,kmaxp1)/
     1      10.0, 15.0, 10.0, 19*1.0,
     2      9*off_table/


c..titanium
      data (wgt(22,k), k=1,kmaxp1)/
     1      21.0, 28.0, 21.0, 10.0, 1.0, 6.0, 9.0, 4.0, 9.0, 2.0, 1.0,
     2      2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 
     3      8*off_table/


c..vanadium
      data (wgt(23,k), k=1,kmaxp1)/
     1      28.0, 25.0, 28.0, 21*1.0,
     2      7*off_table/


c..chromium
      data (wgt(24,k), k=1,kmaxp1)/
     1      7.0, 6.0, 25.0, 28.0, 21.0, 10.0, 1.0, 6.0, 9.0, 2.0, 1.0,
     2      2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 
     3      2.0, 1.0,      
     4      6*off_table/


c..manganese
      data (wgt(25,k), k=1,kmaxp1)/
     1      6.0, 7.0, 6.0, 23*1.0,
     2      5*off_table/


c..iron
      data (wgt(26,k), k=1,kmaxp1)/
     1      25.0, 30.0, 25.0, 6.0, 25.0, 28.0, 21.0, 10.0, 1.0,
     2      6.0, 9.0, 2.0, 9.0, 2.0, 1.0,        
     3      2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 
     4      4*off_table/


c..cobalt
      data (wgt(27,k), k=1,kmaxp1)/
     1      28*1.0,
     2      3*off_table/


c..nickel
      data (wgt(28,k), k=1,kmaxp1)/
     1      21.0, 10.0, 21.0, 10.0, 10.0, 10.0, 25.0, 28.0, 6.0, 6.0, 
     2      6.0, 6.0, 9.0, 6.0, 9.0, 6.0, 1.0, 
     3      2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0,
     4      2*off_table/       



c..copper
      data (wgt(29,k), k=1,kmaxp1)/
     1      30*1.0,
     2      1*off_table/            


c..zinc
      data (wgt(30,k), k=1,kmaxp1)/
     1      31*1.0/      



c..this makes the output simple enough
      izion = int(zion)

      if (izion .gt. 30 .or. stage .gt. izion+1) 
     1      stop 'bad input in stat_weight'

      stat_weight = wgt(izion,stage)

      return
      end







      subroutine coulomb_ion(den,temp,abar,zbar,
     1                       pion,dpiondd,dpiondt,
     2                       xne,dxnedd,dxnedt,
     3                       plasg,
     4                       pcoul,dpcouldd,dpcouldt,
     5                       ecoul,decouldd,decouldt,
     6                       scoul,dscouldd,dscouldt)
      include 'implno.dek'
      include 'const.dek'


c..this routine implments coulomb corrections
c..see yakovlev & shalybkov 1989, uniform background corrections 


c..declare the pass
      double precision den,temp,abar,zbar,
     1                 pion,dpiondd,dpiondt,
     2                 xne,dxnedd,dxnedt,
     3                 plasg,
     4                 pcoul,dpcouldd,dpcouldt,
     5                 ecoul,decouldd,decouldt,
     6                 scoul,dscouldd,dscouldt


c..local variables
c..for the uniform background coulomb correction
      double precision ytot1,kt,ktinv,
     1                 s,dsdd,dsdt,sinv,
     2                 aele,daeledd,daeledt,aeleinv,
     3                 eplasg,eplasgdd,eplasgdt,
     4                 aion,
     5                 lami,inv_lami,lamidd,
     6                 plasgdd,plasgdt,
     7                 x,y,z


      double precision u0,rho6,a1,b1,c1,d1,e1,a2,b2,c2
      parameter        (a1 = -0.898004d0, 
     1                  b1 =  0.96786d0, 
     2                  c1 =  0.220703d0, 
     3                  d1 = -0.86097d0,
     4                  e1 =  2.5269d0, 
     5                  a2 =  0.29561d0, 
     6                  b2 =  1.9885d0,    
     7                  c2 =  0.288675d0)


c..various derived constants
      double precision third,forth,fiveth,esqu,forthpi
      parameter        (third   = 1.0d0/3.0d0,
     1                  forth   = 4.0d0/3.0d0,
     2                  fiveth  = 5.0d0/3.0d0,
     3                  esqu    = qe*qe,
     4                  forthpi = forth * pi)



c..common variables
       ytot1   = 1.0d0/abar
       kt      = kerg * temp
       ktinv   = 1.0d0/kt



c..yakovlev & shalybkov eqs 5, 9 and 10
       s        = forthpi * xne
       dsdd     = forthpi * dxnedd
       dsdt     = forthpi * dxnedt
       sinv     = 1.0d0/s

c..electron-sphere radius aele
       aele     = sinv**third
       z        = -third * aele * sinv
       daeledd  = z * dsdd
       daeledt  = z * dsdt
       aeleinv  = 1.0d0/aele

c..electron coupling parameter eplasg
       eplasg   = esqu * ktinv * aeleinv
       z        = -eplasg * aeleinv
       eplasgdd = z * daeledd
       eplasgdt = z * daeledt - eplasg*ktinv*kerg


c..ion-sphere radius aion
       x        = zbar**third
       aion     = x * aele

c..ion coupling parameter plasg
       z        = x*x*x*x*x
       plasg    = z * eplasg
       plasgdd  = z * eplasgdd
       plasgdt  = z * eplasgdt

c       write(6,*) 
c       write(6,112) aion,aele
c       write(6,112) plasg,plasgdd,plasgdt
c       write(6,*) 



c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
       if (plasg .ge. 1.0) then
        x        = plasg**(0.25d0) 
        u0       = a1*plasg + b1*x + c1/x + d1
        ecoul    = pion/den * u0
        pcoul    = third * ecoul * den
        scoul    = -avo*ytot1*kerg * 
     1              (3.0d0*b1*x - 5.0d0*c1/x
     1             + d1 * (log(plasg) - 1.0d0) - e1)

        y        = avo/abar*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        decouldd = y * plasgdd 
        decouldt = y * plasgdt + ecoul/temp


        y        = third * den
        dpcouldd = third * ecoul + y*decouldd
        dpcouldt = y * decouldt


        y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x +1.25d0*c1/x +d1)
        dscouldd = y * plasgdd
        dscouldt = y * plasgdt



c..yakovlev & shalybkov 1989 equations 102, 103, 104
       else if (plasg .lt. 1.0) then
        x        = plasg*sqrt(plasg)
        y        = plasg**b2
        z        = c2 * x - third * a2 * y
        pcoul    = -pion * z
        ecoul    = 3.0d0 * pcoul/den
        scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

        s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        dpcouldd = -dpiondd*z - pion*s*plasgdd
        dpcouldt = -dpiondt*z - pion*s*plasgdt


        s        = 3.0d0/den
        decouldd = s * dpcouldd - ecoul/den
        decouldt = s * dpcouldt


        s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x -a2*(b2-1.0d0)*y)
        dscouldd = s * plasgdd
        dscouldt = s * plasgdt
       end if

      return
      end








      subroutine etages(xni,zbar,temp,eta)
      include 'implno.dek'
      include 'const.dek'

c..this routine makes a damn good guess for the electron degeneracy 
c..parameter eta. 
c..input is the ion number density xni, average charge zbar,
c..average atomic weigt abar, and temperature temp. 
c..output is a guess at the chemical potential eta


c..declare the pass
      double precision  xni,zbar,temp,eta

c..declare
      double precision xne,x,y,z,kt,beta,tmkt,xnefac

      double precision rt2,rt3,rtpi,cpf0,cpf1,cpf2,cpf3,
     1                 twoth,fa0,forpi,mecc
      parameter        (rt2     = 1.4142135623730951d0,
     1                  rt3     = 1.7320508075688772d0,
     2                  rtpi    = 1.7724538509055159d0,
     3                  cpf0    = h/(me*clight),
     4                  cpf1    = 3.0d0/(8.0d0*pi) * cpf0**3,
     5                  cpf2    = 4.0d0/cpf1,
     6                  cpf3    = 2.0d0*rt3*rtpi/(rt2*cpf1),
     7                  twoth   = 2.0d0/3.0d0,
     8                  fa0     = 64.0d0/(9.0d0*pi),
     9                  forpi   = 4.0d0 * pi,
     &                  mecc    = me * clight * clight)

c..notes: rt2=sqrt(2)  rt3=sqrt(3)  rtpi=sqrt(pi)


c..for the purposes of guessing eta, assume full ionization
      xne   = xni * zbar
      kt    = kerg * temp
      beta  = kt/mecc


c..number density of ionized electrons (c&g 24.354k) and number density at
c..turning point (c&g 24.354i). if either of these exceed the number density
c..as given by a saha equation, then pairs are important. set alfa = 1/2.

      if (beta .ge. 1.0) then
       x = cpf2 * beta * beta
      else
       x = cpf3 * beta * (1.0d0 + 0.75d0*beta) * exp(-1.0d0/beta)
      end if
      if (x .ge. xne) then
       eta = -0.5d0


c..get the dimensionless number density (c&g 24.313), if it is large apply the
c..formula (c&g 24.309) to get a possible alfa, if not large do a two term 
c..binomial expansion on (c&g 24.309) to estimate alfa.

      else
       z = (xne*cpf1)**twoth
       if (z .ge. 1.0e-6) then
        y = (sqrt(z + 1.0d0) - 1.0d0)/beta
       else
        y = z * (1.0d0 - z * 0.25d0) * 0.5d0/beta
       end if


c..isolate the constant in front of the number density integral. if it is
c..small enough run the divine approximation backwards with c&g 24.43. then
c..join it smoothly with the lower limit.

       x = log10(xne**0.6d0/temp)
       if (x .le. 9.5) then
        z = ((1.0d0 + fa0*beta)*sqrt(1.0d0 + fa0*beta*0.5) - 1.0d0)/fa0
        tmkt    = 2.0d0 * me/h * kt/h
        xnefac  = forpi * tmkt * sqrt(tmkt)
        eta = -log(xnefac*rtpi*(0.5d0+0.75d0*z)/xne)
        if (x .ge. 8.5) eta = eta*(9.5d0-x) + y * (1.0d0 - (9.5d0-x))
       else
        eta = y
       end if
      end if

      return
      end






c..routine dfermi gets the fermi-dirac functions and their derivaties
c..routine fdfunc1 forms the integrand of the fermi-dirac functions
c..routine fdfunc2 same as fdfunc but with the change of variable z**2=x
c..routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
c..routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
c..routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
c..routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
c..routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
c..routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
c..routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
c..routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy



      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta)
      include 'implno.dek'
c..
c..this routine computes the fermi-dirac integrals of 
c..index dk, with degeneracy parameter eta and relativity parameter theta.
c..input is dk the double precision index of the fermi-dirac function,
c..eta the degeneracy parameter, and theta the relativity parameter.
c..the output is fd is computed by applying three 10-point 
c..gauss-legendre and one 10-point gauss-laguerre rules over
c..four appropriate subintervals. the derivative with respect to eta is
c..output in fdeta, and the derivative with respct to theta is in fdtheta.
c..within each subinterval the fd kernel.
c..
c..this routine delivers at least 9 figures of accuracy
c..
c..reference: j.m. aparicio, apjs 117, 632 1998
c..
c..declare
c..declare
      external         fdfunc1,fdfunc2
      double precision dk,eta,theta,fd,fdeta,fdtheta,
     1                 d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3,
     2                 eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3),
     3                 res1,dres1,ddres1,res2,dres2,ddres2,
     4                 res3,dres3,ddres3,res4,dres4,ddres4


c   parameters defining the location of the breakpoints for the
c   subintervals of integration:
      data d   / 3.3609d 0 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d 0 /
      data b1  / 1.1418d 0 /
      data c1  / 2.9826d 0 /
      data a2  / 3.7601d 0 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d 1 /
      data e2  / 1.0056d 0 /
      data a3  / 7.5669d 0 /
      data b3  / 1.1695d 0 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d 0 /
      data e3  /-1.2819d-1 /


c   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


c   definition of xi:
      eta1=sg*(eta-d)
      if (eta1.le.5.d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

c   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2)
     +  /(1.d0+c1*xi)
      x2=(a2  +b2*xi+c2*d2*xi2)
     +  /(1.d0+e2*xi+c2*   xi2)
      x3=(a3  +b3*xi+c3*d3*xi2)
     +  /(1.d0+e3*xi+c3*   xi2)

c   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=sqrt(s1)

c   quadrature integrations: 

c 9 significant figure accuracy
c      call dqleg010(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
c      call dqleg010(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
c      call dqleg010(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
c      call dqlag010(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

c 14 significant figure accuracy
      call dqleg020(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
      call dqleg020(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
      call dqleg020(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
      call dqlag020(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

c 18 significant figure accuracy
c      call dqleg040(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
c      call dqleg040(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
c     call dqleg040(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
c     call dqlag040(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

c 32 significant figure accuracy
c      call dqleg080(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
c      call dqleg080(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
c      call dqleg080(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
c      call dqlag080(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)


c..sum the contributions
      fd      = res1 + res2 + res3 + res4
      fdeta   = dres1 + dres2 + dres3 + dres4
      fdtheta = ddres1 + ddres2 + ddres3 + ddres4
      return
      end




      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta)
      include 'implno.dek'
c..
c..forms the fermi-dirac integrand and its derivatives with eta and theta.
c..on input x is the integration variable, par(1) is the double precision 
c..index, par(2) is the degeneravy parameter, and par(3) is the relativity
c..parameter. on output fd is the integrand, fdeta is the derivative
c..with respect to eta, and fdtheta is the derivative with respect to theta.
c..
c..declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta,
     1                 factor,dxst,denom,denom2,xdk,xdkp1

c..initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xdkp1 = x * xdk
      dxst  = sqrt(1.0d0 + 0.5d0*x*theta)

c   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor  = exp(x-eta)
       denom   = factor + 1.0d0
       fd      = xdk * dxst / denom
       fdeta   = fd * factor / denom 
       denom2  = 4.0d0 * dxst * denom
       fdtheta = xdkp1 / denom2

      else
       factor   = exp(eta-x)
       fd       = xdk * dxst * factor
       fdeta    = fd
       denom2   = 4.0d0 * dxst
       fdtheta  = xdkp1/denom2 * factor
      endif

      return
      end




      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta)
      include 'implno.dek'
c..
c..forms the fermi-dirac integrand and its derivatives with eta and theta,
c..when the z**2=x variable change has been made.
c..on input x is the integration variable, par(1) is the double precision 
c..index, par(2) is the degeneravy parameter, and par(3) is the relativity
c..parameter. on output fd is the integrand, fdeta is the derivative
c..with respect to eta, and fdtheta is the derivative with respect to theta.
c..
c..declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta,
     1                 factor,dxst,denom,denom2,xdk,xdkp1,xsq

      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xdkp1 = xsq * xdk
      dxst  = sqrt(1.0d0 + 0.5d0 * xsq * theta)

c   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.d2) then
       factor  = exp(xsq - eta)
       denom   = factor + 1.0d0
       fd      = 2.0d0 * xdk * dxst/denom 
       fdeta   = fd * factor/denom
       denom2  = 4.0d0 * dxst * denom
       fdtheta = 2.0d0 * xdkp1/denom2

      else
       factor  = exp(eta - xsq)
       fd      = 2.0d0 * xdk * dxst * factor
       fdeta   = fd 
       denom2  = 4.0d0 * dxst
       fdtheta = 2.0d0 * xdkp1/denom2 * factor
      endif

      return
      end





      subroutine dqleg010(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..10 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 10-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(5),xg(5),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2

c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 20-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 20-point rule
c wg     - weights of the 20-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.4887433898 1631210884 8260011297 19984 d -1 /
      data xg (  2) /   4.3339539412 9247190799 2659431657 84162 d -1 /
      data xg (  3) /   6.7940956829 9024406234 3273651148 73575 d -1 /
      data xg (  4) /   8.6506336668 8984510732 0966884234 93048 d -1 /
      data xg (  5) /   9.7390652851 7171720077 9640120844 52053 d -1 /

      data wg (  1) /   2.9552422471 4752870173 8929946513 38329 d -1 /
      data wg (  2) /   2.6926671930 9996355091 2269215694 69352 d -1 /
      data wg (  3) /   2.1908636251 5982043995 5349342281 63192 d -1 /
      data wg (  4) /   1.4945134915 0580593145 7763396576 97332 d -1 /
      data wg (  5) /   6.6671344308 6881375935 6880989333 17928 d -2 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 10-point gauss formula

      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,5
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end




      subroutine dqleg020(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(10),xg(10),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 20-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 20-point rule
c wg     - weights of the 20-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.6526521133 4973337546 4040939883 82110 d -2 /
      data xg (  2) /   2.2778585114 1645078080 4961953685 74624 d -1 /
      data xg (  3) /   3.7370608871 5419560672 5481770249 27237 d -1 /
      data xg (  4) /   5.1086700195 0827098004 3640509552 50998 d -1 /
      data xg (  5) /   6.3605368072 6515025452 8366962262 85936 d -1 /
      data xg (  6) /   7.4633190646 0150792614 3050703556 41590 d -1 /
      data xg (  7) /   8.3911697182 2218823394 5290617015 20685 d -1 /
      data xg (  8) /   9.1223442825 1325905867 7524412032 98113 d -1 /
      data xg (  9) /   9.6397192727 7913791267 6661311972 77221 d -1 /
      data xg ( 10) /   9.9312859918 5094924786 1223884713 20278 d -1 /

      data wg (  1) /   1.5275338713 0725850698 0843319550 97593 d -1 /
      data wg (  2) /   1.4917298647 2603746787 8287370019 69436 d -1 /
      data wg (  3) /   1.4209610931 8382051329 2983250671 64933 d -1 /
      data wg (  4) /   1.3168863844 9176626898 4944997481 63134 d -1 /
      data wg (  5) /   1.1819453196 1518417312 3773777113 82287 d -1 /
      data wg (  6) /   1.0193011981 7240435036 7501354803 49876 d -1 /
      data wg (  7) /   8.3276741576 7047487247 5814322204 62061 d -2 /
      data wg (  8) /   6.2672048334 1090635695 0653518704 16063 d -2 /
      data wg (  9) /   4.0601429800 3869413310 3995227493 21098 d -2 /
      data wg ( 10) /   1.7614007139 1521183118 6196235185 28163 d -2 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end





      subroutine dqleg040(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..40 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 40-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(20),xg(20),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 40-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 40-point rule
c wg     - weights of the 40-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.8772417506 0508219331 9344402462 32946 d -2 /
      data xg (  2) /   1.1608407067 5255208483 4512844080 24113 d -1 /
      data xg (  3) /   1.9269758070 1371099715 5168520651 49894 d -1 /
      data xg (  4) /   2.6815218500 7253681141 1843448085 96183 d -1 /
      data xg (  5) /   3.4199409082 5758473007 4924811791 94310 d -1 /
      data xg (  6) /   4.1377920437 1605001524 8797458037 13682 d -1 /
      data xg (  7) /   4.8307580168 6178712908 5665742448 23004 d -1 /
      data xg (  8) /   5.4946712509 5128202075 9313055295 17970 d -1 /
      data xg (  9) /   6.1255388966 7980237952 6124502306 94877 d -1 /
      data xg ( 10) /   6.7195668461 4179548379 3545149614 94109 d -1 /
      data xg ( 11) /   7.2731825518 9927103280 9964517549 30548 d -1 /
      data xg ( 12) /   7.7830565142 6519387694 9715455064 94848 d -1 /
      data xg ( 13) /   8.2461223083 3311663196 3202306660 98773 d -1 /
      data xg ( 14) /   8.6595950321 2259503820 7818083546 19963 d -1 /
      data xg ( 15) /   9.0209880696 8874296728 2533308684 93103 d -1 /
      data xg ( 16) /   9.3281280827 8676533360 8521668452 05716 d -1 /
      data xg ( 17) /   9.5791681921 3791655804 5409994527 59285 d -1 /
      data xg ( 18) /   9.7725994998 3774262663 3702837129 03806 d -1 /
      data xg ( 19) /   9.9072623869 9457006453 0543522213 72154 d -1 /
      data xg ( 20) /   9.9823770971 0559200349 6227024205 86492 d -1 /

      data wg (  1) /   7.7505947978 4248112637 2396295832 63269 d -2 /
      data wg (  2) /   7.7039818164 2479655883 0753428381 02485 d -2 /
      data wg (  3) /   7.6110361900 6262423715 5807592249 48230 d -2 /
      data wg (  4) /   7.4723169057 9682642001 8933626132 46731 d -2 /
      data wg (  5) /   7.2886582395 8040590605 1068344251 78358 d -2 /
      data wg (  6) /   7.0611647391 2867796954 8363085528 68323 d -2 /
      data wg (  7) /   6.7912045815 2339038256 9010823192 39859 d -2 /
      data wg (  8) /   6.4804013456 6010380745 5452956675 27300 d -2 /
      data wg (  9) /   6.1306242492 9289391665 3799640839 85959 d -2 /
      data wg ( 10) /   5.7439769099 3915513666 1773091042 59856 d -2 /
      data wg ( 11) /   5.3227846983 9368243549 9647977226 05045 d -2 /
      data wg ( 12) /   4.8695807635 0722320614 3416044814 63880 d -2 /
      data wg ( 13) /   4.3870908185 6732719916 7468604171 54958 d -2 /
      data wg ( 14) /   3.8782167974 4720176399 7203129044 61622 d -2 /
      data wg ( 15) /   3.3460195282 5478473926 7818308641 08489 d -2 /
      data wg ( 16) /   2.7937006980 0234010984 8915750772 10773 d -2 /
      data wg ( 17) /   2.2245849194 1669572615 0432418420 85732 d -2 /
      data wg ( 18) /   1.6421058381 9078887128 6348488236 39272 d -2 /
      data wg ( 19) /   1.0498284531 1528136147 4217106727 96523 d -2 /
      data wg ( 20) /   4.5212770985 3319125847 1732878185 33272 d -3 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end





      subroutine dqleg080(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..80 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 80-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(40),xg(40),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 80-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 80-point rule
c wg     - weights of the 80-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   1.9511383256 7939976543 5123410745 45479 d -2 /
      data xg (  2) /   5.8504437152 4206686289 9332188341 77944 d -2 /
      data xg (  3) /   9.7408398441 5845990632 7845010493 69020 d -2 /
      data xg (  4) /   1.3616402280 9143886559 2410780007 17067 d -1 /
      data xg (  5) /   1.7471229183 2646812559 3390480112 86195 d -1 /
      data xg (  6) /   2.1299450285 7666132572 3885386663 21823 d -1 /
      data xg (  7) /   2.5095235839 2272120493 1588160350 04797 d -1 /
      data xg (  8) /   2.8852805488 4511853109 1393014347 13898 d -1 /
      data xg (  9) /   3.2566437074 7701914619 1129436273 58695 d -1 /
      data xg ( 10) /   3.6230475349 9487315619 0432863589 63588 d -1 /
      data xg ( 11) /   3.9839340588 1969227024 3796425175 33757 d -1 /
      data xg ( 12) /   4.3387537083 1756093062 3867003631 81958 d -1 /
      data xg ( 13) /   4.6869661517 0544477036 0783649358 08657 d -1 /
      data xg ( 14) /   5.0280411188 8784987593 6727503675 68003 d -1 /
      data xg ( 15) /   5.3614592089 7131932019 8572531254 00904 d -1 /
      data xg ( 16) /   5.6867126812 2709784725 4857866248 27158 d -1 /
      data xg ( 17) /   6.0033062282 9751743154 7462991640 06848 d -1 /
      data xg ( 18) /   6.3107577304 6871966247 9283872893 36863 d -1 /
      data xg ( 19) /   6.6085989898 6119801735 9671228443 17234 d -1 /
      data xg ( 20) /   6.8963764434 2027600771 2076124389 35266 d -1 /
      data xg ( 21) /   7.1736518536 2099880254 0682582938 15278 d -1 /
      data xg ( 22) /   7.4400029758 3597272316 5405279309 13673 d -1 /
      data xg ( 23) /   7.6950242013 5041373865 6160687490 26083 d -1 /
      data xg ( 24) /   7.9383271750 4605449948 6393117384 54358 d -1 /
      data xg ( 25) /   8.1695413868 1463470371 1249940122 95707 d -1 /
      data xg ( 26) /   8.3883147358 0255275616 6230439028 67064 d -1 /
      data xg ( 27) /   8.5943140666 3111096977 1921234916 56492 d -1 /
      data xg ( 28) /   8.7872256767 8213828703 7733436391 24407 d -1 /
      data xg ( 29) /   8.9667557943 8770683194 3240719673 95986 d -1 /
      data xg ( 30) /   9.1326310257 1757654164 7336561509 47478 d -1 /
      data xg ( 31) /   9.2845987717 2445795953 0459590754 53133 d -1 /
      data xg ( 32) /   9.4224276130 9872674752 2660045000 01735 d -1 /
      data xg ( 33) /   9.5459076634 3634905493 4815170210 29508 d -1 /
      data xg ( 34) /   9.6548508904 3799251452 2731556714 54998 d -1 /
      data xg ( 35) /   9.7490914058 5727793385 6452300691 36276 d -1 /
      data xg ( 36) /   9.8284857273 8629070418 2880277091 16473 d -1 /
      data xg ( 37) /   9.8929130249 9755531026 5031671366 31385 d -1 /
      data xg ( 38) /   9.9422754096 5688277892 0635036649 11698 d -1 /
      data xg ( 39) /   9.9764986439 8237688899 4942081831 22985 d -1 /
      data xg ( 40) /   9.9955382265 1630629880 0804990945 67184 d -1 /

      data wg (  1) /   3.9017813656 3066548112 8043925275 40483 d -2 /
      data wg (  2) /   3.8958395962 7695311986 2552477226 08223 d -2 /
      data wg (  3) /   3.8839651059 0519689317 7418266878 71658 d -2 /
      data wg (  4) /   3.8661759774 0764633270 7711026715 66912 d -2 /
      data wg (  5) /   3.8424993006 9594231852 1243632949 01384 d -2 /
      data wg (  6) /   3.8129711314 4776383442 0679156573 62019 d -2 /
      data wg (  7) /   3.7776364362 0013974897 7497642632 10547 d -2 /
      data wg (  8) /   3.7365490238 7304900267 0537705783 86691 d -2 /
      data wg (  9) /   3.6897714638 2760088391 5099657340 52192 d -2 /
      data wg ( 10) /   3.6373749905 8359780439 6499104652 28136 d -2 /
      data wg ( 11) /   3.5794393953 4160546028 6158881615 44542 d -2 /
      data wg ( 12) /   3.5160529044 7475934955 2659238869 68812 d -2 /
      data wg ( 13) /   3.4473120451 7539287943 6422673102 98320 d -2 /
      data wg ( 14) /   3.3733214984 6115228166 7516306423 87284 d -2 /
      data wg ( 15) /   3.2941939397 6454013828 3618090195 95361 d -2 /
      data wg ( 16) /   3.2100498673 4877731480 5649028725 06960 d -2 /
      data wg ( 17) /   3.1210174188 1147016424 4286672060 35518 d -2 /
      data wg ( 18) /   3.0272321759 5579806612 2001009090 11747 d -2 /
      data wg ( 19) /   2.9288369583 2678476927 6758601957 91396 d -2 /
      data wg ( 20) /   2.8259816057 2768623967 5319796501 45302 d -2 /
      data wg ( 21) /   2.7188227500 4863806744 1870668054 42598 d -2 /
      data wg ( 22) /   2.6075235767 5651179029 6874360026 92871 d -2 /
      data wg ( 23) /   2.4922535764 1154911051 1784700321 98023 d -2 /
      data wg ( 24) /   2.3731882865 9301012931 9252461356 84162 d -2 /
      data wg ( 25) /   2.2505090246 3324619262 2158968616 87390 d -2 /
      data wg ( 26) /   2.1244026115 7820063887 1073725061 31285 d -2 /
      data wg ( 27) /   1.9950610878 1419989288 9192871511 35633 d -2 /
      data wg ( 28) /   1.8626814208 2990314287 3541415215 72090 d -2 /
      data wg ( 29) /   1.7274652056 2693063585 8420713129 09998 d -2 /
      data wg ( 30) /   1.5896183583 7256880449 0290922917 85257 d -2 /
      data wg ( 31) /   1.4493508040 5090761169 6207458346 05500 d -2 /
      data wg ( 32) /   1.3068761592 4013392937 8682589705 63403 d -2 /
      data wg ( 33) /   1.1624114120 7978269164 6676999543 26348 d -2 /
      data wg ( 34) /   1.0161766041 1030645208 3185035240 69436 d -2 /
      data wg ( 35) /   8.6839452692 6085842640 9452204034 28135 d -3 /
      data wg ( 36) /   7.1929047681 1731275267 5570867956 50747 d -3 /
      data wg ( 37) /   5.6909224514 0319864926 9107117162 01847 d -3 /
      data wg ( 38) /   4.1803131246 9489523673 9304201681 35132 d -3 /
      data wg ( 39) /   2.6635335895 1268166929 3535831668 45546 d -3 /
      data wg ( 40) /   1.1449500031 8694153454 4171941315 63611 d -3 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end





      subroutine dqlag010(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..10 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=exp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 10-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(10),xg(10),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 10-point gauss-laguerre rule
c wg     - weights of the 10-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.3779347054 0492430830 7725056527 11188 d -1 /
      data xg (  2) /   7.2945454950 3170498160 3731216760 78781 d -1 /
      data xg (  3) /   1.8083429017 4031604823 2920075750 60883 d  0 /
      data xg (  4) /   3.4014336978 5489951448 2532221408 39067 d  0 /
      data xg (  5) /   5.5524961400 6380363241 7558486868 76285 d  0 /
      data xg (  6) /   8.3301527467 6449670023 8767197274 52218 d  0 /
      data xg (  7) /   1.1843785837 9000655649 1853891914 16139 d  1 /
      data xg (  8) /   1.6279257831 3781020995 3265393583 36223 d  1 /
      data xg (  9) /   2.1996585811 9807619512 7709019559 44939 d  1 /
      data xg ( 10) /   2.9920697012 2738915599 0879334079 91951 d  1 /

      data wg (  1) /   3.5400973860 6996308762 2268914420 67608 d -1 /
      data wg (  2) /   8.3190230104 3580738109 8296581278 49577 d -1 /
      data wg (  3) /   1.3302885617 4932817875 2792194393 99369 d  0 /
      data wg (  4) /   1.8630639031 1113098976 3988735482 46693 d  0 /
      data wg (  5) /   2.4502555580 8301016607 2693731657 52256 d  0 /
      data wg (  6) /   3.1227641551 3518249615 0818263314 55472 d  0 /
      data wg (  7) /   3.9341526955 6152109865 5812459248 23077 d  0 /
      data wg (  8) /   4.9924148721 9302310201 1485652433 15445 d  0 /
      data wg (  9) /   6.5722024851 3080297518 7668710376 11234 d  0 /
      data wg ( 10) /   9.7846958403 7463069477 0086638718 59813 d  0 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 10-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end





      subroutine dqlag020(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(20),xg(20),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.0539889691 9887533666 8900458421 50958 d -2 /
      data xg (  2) /   3.7212681800 1611443794 2413887611 46636 d -1 /
      data xg (  3) /   9.1658210248 3273564667 7162770741 83187 d -1 /
      data xg (  4) /   1.7073065310 2834388068 7689667413 05070 d  0 /
      data xg (  5) /   2.7491992553 0943212964 5030460494 81338 d  0 /
      data xg (  6) /   4.0489253138 5088692237 4953369133 33219 d  0 /
      data xg (  7) /   5.6151749708 6161651410 4539885651 89234 d  0 /
      data xg (  8) /   7.4590174536 7106330976 8860218371 81759 d  0 /
      data xg (  9) /   9.5943928695 8109677247 3672734282 79837 d  0 /
      data xg ( 10) /   1.2038802546 9643163096 2340929886 55158 d  1 /
      data xg ( 11) /   1.4814293442 6307399785 1267971004 79756 d  1 /
      data xg ( 12) /   1.7948895520 5193760173 6579099261 25096 d  1 /
      data xg ( 13) /   2.1478788240 2850109757 3517036959 46692 d  1 /
      data xg ( 14) /   2.5451702793 1869055035 1867748464 15418 d  1 /
      data xg ( 15) /   2.9932554631 7006120067 1365613516 58232 d  1 /
      data xg ( 16) /   3.5013434240 4790000062 8493590668 81395 d  1 /
      data xg ( 17) /   4.0833057056 7285710620 2956770780 75526 d  1 /
      data xg ( 18) /   4.7619994047 3465021399 4162715285 11211 d  1 /
      data xg ( 19) /   5.5810795750 0638988907 5077344449 72356 d  1 /
      data xg ( 20) /   6.6524416525 6157538186 4031879146 06659 d  1 /

      data wg (  1) /   1.8108006241 8989255451 6754059131 10644 d -1 /
      data wg (  2) /   4.2255676787 8563974520 3441725664 58197 d -1 /
      data wg (  3) /   6.6690954670 1848150373 4821149925 15927 d -1 /
      data wg (  4) /   9.1535237278 3073672670 6046847718 68067 d -1 /
      data wg (  5) /   1.1695397071 9554597380 1478222395 77476 d  0 /
      data wg (  6) /   1.4313549859 2820598636 8449948915 14331 d  0 /
      data wg (  7) /   1.7029811379 8502272402 5332616332 06720 d  0 /
      data wg (  8) /   1.9870158907 9274721410 9218392751 29020 d  0 /
      data wg (  9) /   2.2866357812 5343078546 2228546814 95651 d  0 /
      data wg ( 10) /   2.6058347275 5383333269 4989509540 33323 d  0 /
      data wg ( 11) /   2.9497837342 1395086600 2354168272 85951 d  0 /
      data wg ( 12) /   3.3253957820 0931955236 9519374217 51118 d  0 /
      data wg ( 13) /   3.7422554705 8981092111 7072932653 77811 d  0 /
      data wg ( 14) /   4.2142367102 5188041986 8080637824 78746 d  0 /
      data wg ( 15) /   4.7625184614 9020929695 2921978390 96371 d  0 /
      data wg ( 16) /   5.4217260442 4557430380 3082979899 81779 d  0 /
      data wg ( 17) /   6.2540123569 3242129289 5184903007 07542 d  0 /
      data wg ( 18) /   7.3873143890 5443455194 0300191964 64791 d  0 /
      data wg ( 19) /   9.1513287309 8747960794 3482425529 50528 d  0 /
      data wg ( 20) /   1.2893388645 9399966710 2628712874 85278 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end




      subroutine dqlag040(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(40),xg(40),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.5700394308 8883851220 8447128660 08554 d -2 /
      data xg (  2) /   1.8816228315 8698516003 5893462190 95913 d -1 /
      data xg (  3) /   4.6269428131 4576453564 9375245611 90364 d -1 /
      data xg (  4) /   8.5977296397 2934922257 2722246887 22412 d -1 /
      data xg (  5) /   1.3800108205 2733718649 8000329595 26559 d  0 /
      data xg (  6) /   2.0242091359 2282673344 2066002800 13075 d  0 /
      data xg (  7) /   2.7933693535 0681645765 3514486026 64039 d  0 /
      data xg (  8) /   3.6887026779 0827020959 1526351908 68698 d  0 /
      data xg (  9) /   4.7116411465 5497269361 8722836277 47369 d  0 /
      data xg ( 10) /   5.8638508783 4371811427 3164237995 82987 d  0 /
      data xg ( 11) /   7.1472479081 0228825068 5691951979 42362 d  0 /
      data xg ( 12) /   8.5640170175 8616376271 8522042088 13232 d  0 /
      data xg ( 13) /   1.0116634048 4519394068 4962965639 52448 d  1 /
      data xg ( 14) /   1.1807892294 0045848428 4158670436 06304 d  1 /
      data xg ( 15) /   1.3640933712 5370872283 7167636065 01202 d  1 /
      data xg ( 16) /   1.5619285893 3390738372 0196365218 80145 d  1 /
      data xg ( 17) /   1.7746905950 0956630425 7387749542 43772 d  1 /
      data xg ( 18) /   2.0028232834 5748905296 1261481017 51172 d  1 /
      data xg ( 19) /   2.2468249983 4984183513 7178622899 45366 d  1 /
      data xg ( 20) /   2.5072560772 4262037943 9608620940 09769 d  1 /
      data xg ( 21) /   2.7847480009 1688627207 5170414045 57997 d  1 /
      data xg ( 22) /   3.0800145739 4454627007 5438519619 11114 d  1 /
      data xg ( 23) /   3.3938657084 9137196090 9885858628 19990 d  1 /
      data xg ( 24) /   3.7272245880 4760043283 2076099060 74207 d  1 /
      data xg ( 25) /   4.0811492823 8869204661 5567558160 06426 d  1 /
      data xg ( 26) /   4.4568603175 3344627071 2302063449 83559 d  1 /
      data xg ( 27) /   4.8557763533 0599922809 6204880670 67936 d  1 /
      data xg ( 28) /   5.2795611187 2169329693 5202113739 17638 d  1 /
      data xg ( 29) /   5.7301863323 3936274950 3374699589 21651 d  1 /
      data xg ( 30) /   6.2100179072 7751116121 6819905789 89921 d  1 /
      data xg ( 31) /   6.7219370927 1269987990 8027755188 87054 d  1 /
      data xg ( 32) /   7.2695158847 6124621175 2192772426 19385 d  1 /
      data xg ( 33) /   7.8572802911 5713092805 4389683348 12596 d  1 /
      data xg ( 34) /   8.4911231135 7049845427 0156470966 63186 d  1 /
      data xg ( 35) /   9.1789874671 2363769923 3719348062 73153 d  1 /
      data xg ( 36) /   9.9320808717 4468082501 0905416548 68123 d  1 /
      data xg ( 37) /   1.0767244063 9388272520 7967676113 22664 d  2 /
      data xg ( 38) /   1.1712230951 2690688807 6506441235 50702 d  2 /
      data xg ( 39) /   1.2820184198 8255651192 5411043896 31263 d  2 /
      data xg ( 40) /   1.4228004446 9159997888 3488353595 41764 d  2 /

      data wg (  1) /   9.1625471157 4598973115 1169808013 74830 d -2 /
      data wg (  2) /   2.1342058490 5012080007 1933671215 12341 d -1 /
      data wg (  3) /   3.3571811668 0284673880 5107016162 92191 d -1 /
      data wg (  4) /   4.5854093503 3497560385 4323803764 52497 d -1 /
      data wg (  5) /   5.8206816577 9105168990 9963654015 43283 d -1 /
      data wg (  6) /   7.0649521636 7219392989 8300156730 16682 d -1 /
      data wg (  7) /   8.3202690300 3485238099 1129479783 49523 d -1 /
      data wg (  8) /   9.5887819879 4443111448 1226796760 28906 d -1 /
      data wg (  9) /   1.0872761620 3054971575 3869333172 02661 d  0 /
      data wg ( 10) /   1.2174623279 7778097895 4277850665 60948 d  0 /
      data wg ( 11) /   1.3496954913 5676530792 3938594423 94519 d  0 /
      data wg ( 12) /   1.4842549297 7684671120 5611786129 78719 d  0 /
      data wg ( 13) /   1.6214441628 1182197802 3168843164 54527 d  0 /
      data wg ( 14) /   1.7615953746 7676961118 4242204209 81598 d  0 /
      data wg ( 15) /   1.9050746658 9479967668 2993205972 79371 d  0 /
      data wg ( 16) /   2.0522883472 6171671760 1995822729 47454 d  0 /
      data wg ( 17) /   2.2036905532 4509588909 8283443281 40570 d  0 /
      data wg ( 18) /   2.3597925385 2320332354 0373753789 01497 d  0 /
      data wg ( 19) /   2.5211741403 7643299165 3136902874 22820 d  0 /
      data wg ( 20) /   2.6884980554 0884226415 9505447063 74659 d  0 /
      data wg ( 21) /   2.8625278132 1044881203 4763959831 04311 d  0 /
      data wg ( 22) /   3.0441506653 1151710041 0439679543 33670 d  0 /
      data wg ( 23) /   3.2344070972 6353194177 4902394288 67111 d  0 /
      data wg ( 24) /   3.4345293984 2774809220 3984818916 02464 d  0 /
      data wg ( 25) /   3.6459928249 9408907238 9656466994 90434 d  0 /
      data wg ( 26) /   3.8705845972 1651656808 4753202134 44338 d  0 /
      data wg ( 27) /   4.1104986804 3282265583 5822472639 51577 d  0 /
      data wg ( 28) /   4.3684687232 5406347450 8083382729 45025 d  0 /
      data wg ( 29) /   4.6479589840 7446688299 3033998838 83991 d  0 /
      data wg ( 30) /   4.9534461124 0989326218 6961507855 62721 d  0 /
      data wg ( 31) /   5.2908484059 0073657468 7373657188 58968 d  0 /
      data wg ( 32) /   5.6682046090 3297677000 7305290232 63795 d  0 /
      data wg ( 33) /   6.0967964147 4342030593 3760108591 98806 d  0 /
      data wg ( 34) /   6.5931088610 3999953794 4296642062 94899 d  0 /
      data wg ( 35) /   7.1824959955 3689315064 4298016266 99574 d  0 /
      data wg ( 36) /   7.9066663113 8422877369 3107423105 86595 d  0 /
      data wg ( 37) /   8.8408924928 1034652079 1255950630 26792 d  0 /
      data wg ( 38) /   1.0140899265 6211694839 0946003069 40468 d  1 /
      data wg ( 39) /   1.2210021299 2046038985 2264858758 81108 d  1 /
      data wg ( 40) /   1.6705520642 0242974052 4687743985 73553 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end





      subroutine dqlag080(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(80),xg(80),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.7960423300 6983655540 1031924740 16803 d -2 /
      data xg (  2) /   9.4639912994 3539888113 9027246521 72943 d -2 /
      data xg (  3) /   2.3262286812 5867569207 7061572163 49831 d -1 /
      data xg (  4) /   4.3199254780 2387480255 7861724977 70411 d -1 /
      data xg (  5) /   6.9282886135 2021839905 7022136354 46867 d -1 /
      data xg (  6) /   1.0152325561 8947143744 6254368599 35350 d  0 /
      data xg (  7) /   1.3993276878 4287277414 4190514309 78382 d  0 /
      data xg (  8) /   1.8452623038 3584513811 1771177695 99966 d  0 /
      data xg (  9) /   2.3532088716 0926152447 2447080161 40181 d  0 /
      data xg ( 10) /   2.9233646865 5542632483 6912342597 32862 d  0 /
      data xg ( 11) /   3.5559523140 4613405944 9673083246 38370 d  0 /
      data xg ( 12) /   4.2512200823 0987808316 4857664485 77637 d  0 /
      data xg ( 13) /   5.0094426336 2016477243 3677068182 06389 d  0 /
      data xg ( 14) /   5.8309215386 0871901982 1271132956 05083 d  0 /
      data xg ( 15) /   6.7159859778 5131711156 5500876351 99430 d  0 /
      data xg ( 16) /   7.6649934948 9177306073 4189090478 23480 d  0 /
      data xg ( 17) /   8.6783308251 6770109543 4422555426 61083 d  0 /
      data xg ( 18) /   9.7564148057 4293071316 5509446173 66591 d  0 /
      data xg ( 19) /   1.0899693371 2878553774 3610010214 89406 d  1 /
      data xg ( 20) /   1.2108646642 3656999007 0548486983 15593 d  1 /
      data xg ( 21) /   1.3383788112 7786473701 6298406038 33297 d  1 /
      data xg ( 22) /   1.4725665943 5085855393 3580762618 38437 d  1 /
      data xg ( 23) /   1.6134864371 6624665791 6585454289 90907 d  1 /
      data xg ( 24) /   1.7612005243 8144378598 6356869435 86520 d  1 /
      data xg ( 25) /   1.9157749684 2412479221 7299702056 74985 d  1 /
      data xg ( 26) /   2.0772799909 7920960924 4193790104 89579 d  1 /
      data xg ( 27) /   2.2457901204 5404583114 0959169508 77516 d  1 /
      data xg ( 28) /   2.4213844068 9586473771 9224693924 47092 d  1 /
      data xg ( 29) /   2.6041466560 1655866929 3900535654 35682 d  1 /
      data xg ( 30) /   2.7941656841 8594655558 2330692936 92111 d  1 /
      data xg ( 31) /   2.9915355964 9009855011 2704121157 37715 d  1 /
      data xg ( 32) /   3.1963560902 2089207107 7488875426 36533 d  1 /
      data xg ( 33) /   3.4087327864 7261898749 8349473428 60505 d  1 /
      data xg ( 34) /   3.6287775928 7814544588 0319884362 16948 d  1 /
      data xg ( 35) /   3.8566091009 2922104582 5630521729 08535 d  1 /
      data xg ( 36) /   4.0923530218 0312671999 0958505955 44326 d  1 /
      data xg ( 37) /   4.3361426651 7312302957 8267604682 19500 d  1 /
      data xg ( 38) /   4.5881194661 2788863456 2664899748 78378 d  1 /
      data xg ( 39) /   4.8484335660 8331891358 7372733535 63006 d  1 /
      data xg ( 40) /   5.1172444544 6070105959 8894323349 07144 d  1 /
      data xg ( 41) /   5.3947216789 5544471206 2102787875 72430 d  1 /
      data xg ( 42) /   5.6810456334 6362231341 2485032441 02122 d  1 /
      data xg ( 43) /   5.9764084342 1099549427 2959612774 71927 d  1 /
      data xg ( 44) /   6.2810148963 9264772036 2729175902 88682 d  1 /
      data xg ( 45) /   6.5950836257 4560573434 6406271607 92248 d  1 /
      data xg ( 46) /   6.9188482420 2362773741 9802886482 37373 d  1 /
      data xg ( 47) /   7.2525587544 2633453588 3896526165 68450 d  1 /
      data xg ( 48) /   7.5964831127 8641748269 4497949747 96502 d  1 /
      data xg ( 49) /   7.9509089629 0888369620 5728262599 80809 d  1 /
      data xg ( 50) /   8.3161456401 0536896630 4295068758 48705 d  1 /
      data xg ( 51) /   8.6925264419 6156234481 1659260404 48396 d  1 /
      data xg ( 52) /   9.0804112300 9407559518 4117278203 18427 d  1 /
      data xg ( 53) /   9.4801894215 9474332072 0718891387 35302 d  1 /
      data xg ( 54) /   9.8922834446 9405791648 0193727380 36790 d  1 /
      data xg ( 55) /   1.0317152750 8039130233 0470941673 45654 d  2 /
      data xg ( 56) /   1.0755298497 7539906327 6078907989 75954 d  2 /
      data xg ( 57) /   1.1207269048 4128333623 9300461662 11013 d  2 /
      data xg ( 58) /   1.1673666467 3503666318 1578881308 01099 d  2 /
      data xg ( 59) /   1.2155154249 0952625566 8638957521 10813 d  2 /
      data xg ( 60) /   1.2652466579 6515540341 5702654316 53573 d  2 /
      data xg ( 61) /   1.3166419525 2120310870 0890863080 06192 d  2 /
      data xg ( 62) /   1.3697924668 6936973947 5706372894 63788 d  2 /
      data xg ( 63) /   1.4248005891 2161601930 8265692004 55232 d  2 /
      data xg ( 64) /   1.4817820245 5004441818 6523848360 07732 d  2 /
      data xg ( 65) /   1.5408684228 1798697859 4174252655 96259 d  2 /
      data xg ( 66) /   1.6022107287 0095715935 2684168930 10646 d  2 /
      data xg ( 67) /   1.6659835193 4053918744 5211797337 12213 d  2 /
      data xg ( 68) /   1.7323907133 4249503830 9065037750 56999 d  2 /
      data xg ( 69) /   1.8016732304 9032317982 4302089977 01523 d  2 /
      data xg ( 70) /   1.8741194967 6963772390 4901345880 21771 d  2 /
      data xg ( 71) /   1.9500802244 1532991450 3904796005 99643 d  2 /
      data xg ( 72) /   2.0299898419 5074937824 8076778237 14777 d  2 /
      data xg ( 73) /   2.1143987049 4836466691 4849046955 42608 d  2 /
      data xg ( 74) /   2.2040236815 1735739654 0442066777 63168 d  2 /
      data xg ( 75) /   2.2998320607 5680004348 4109696758 44754 d  2 /
      data xg ( 76) /   2.4031908705 5841540417 5974604792 19628 d  2 /
      data xg ( 77) /   2.5161587933 0499611167 4449393109 73194 d  2 /
      data xg ( 78) /   2.6421382388 3199102097 6961086914 35553 d  2 /
      data xg ( 79) /   2.7876673304 6004563652 0141725306 11597 d  2 /
      data xg ( 80) /   2.9696651199 5651345758 8528591557 03581 d  2 /

      data wg (  1) /   4.6093103133 0609664705 2513213955 10083 d -2 /
      data wg (  2) /   1.0731300778 3932752564 1503203043 98860 d -1 /
      data wg (  3) /   1.6866442954 7948111794 2204577827 02406 d -1 /
      data wg (  4) /   2.3008808938 4940054411 2571819781 93282 d -1 /
      data wg (  5) /   2.9160130250 2437964832 1693187729 43752 d -1 /
      data wg (  6) /   3.5322675357 5408236352 7231258056 47046 d -1 /
      data wg (  7) /   4.1498817755 0940466187 1976863112 80092 d -1 /
      data wg (  8) /   4.7690979230 2936241314 7770254185 05661 d -1 /
      data wg (  9) /   5.3901621847 4955374499 5076565223 27912 d -1 /
      data wg ( 10) /   6.0133249744 7190529086 7652488407 39512 d -1 /
      data wg ( 11) /   6.6388413639 6680571849 4422407272 99214 d -1 /
      data wg ( 12) /   7.2669716361 4156688973 5672962491 40514 d -1 /
      data wg ( 13) /   7.8979818942 8428531349 7930783987 88294 d -1 /
      data wg ( 14) /   8.5321447143 8152298354 5981624313 62968 d -1 /
      data wg ( 15) /   9.1697398383 3892698590 3429000315 53302 d -1 /
      data wg ( 16) /   9.8110549100 4005747195 0601559842 18607 d -1 /
      data wg ( 17) /   1.0456386258 0654218147 5684456631 76029 d  0 /
      data wg ( 18) /   1.1106039730 0025890771 1247632597 29371 d  0 /
      data wg ( 19) /   1.1760331584 1226175056 6510765192 08666 d  0 /
      data wg ( 20) /   1.2419589444 9809359279 3517618178 71338 d  0 /
      data wg ( 21) /   1.3084153330 3134064261 1885428459 54645 d  0 /
      data wg ( 22) /   1.3754376757 4892843813 1559170934 90796 d  0 /
      data wg ( 23) /   1.4430627938 7849270398 3124172072 47308 d  0 /
      data wg ( 24) /   1.5113291075 8830693847 6550205599 17703 d  0 /
      data wg ( 25) /   1.5802767765 3099415830 2018787231 21659 d  0 /
      data wg ( 26) /   1.6499478528 0267874116 0120428193 55036 d  0 /
      data wg ( 27) /   1.7203864478 1283277182 0042814522 90770 d  0 /
      data wg ( 28) /   1.7916389147 6093832891 4426205276 88915 d  0 /
      data wg ( 29) /   1.8637540486 4909708435 9257090286 88162 d  0 /
      data wg ( 30) /   1.9367833060 3070923513 9254343278 41646 d  0 /
      data wg ( 31) /   2.0107810470 1134222912 6149881755 55546 d  0 /
      data wg ( 32) /   2.0858048023 8741046429 3039785129 89079 d  0 /
      data wg ( 33) /   2.1619155692 4159897378 3163440488 27763 d  0 /
      data wg ( 34) /   2.2391781388 2364652373 4539974474 45645 d  0 /
      data wg ( 35) /   2.3176614611 4651854068 6060480434 96370 d  0 /
      data wg ( 36) /   2.3974390514 4001430514 1172386388 49980 d  0 /
      data wg ( 37) /   2.4785894444 4973417756 3691644552 22527 d  0 /
      data wg ( 38) /   2.5611967035 7790455335 1155092225 72643 d  0 /
      data wg ( 39) /   2.6453509930 6968892850 4634410003 67534 d  0 /
      data wg ( 40) /   2.7311492228 9915138861 4102871311 69260 d  0 /
      data wg ( 41) /   2.8186957777 5934171703 1418737478 11157 d  0 /
      data wg ( 42) /   2.9081033436 8223018934 5502767774 92687 d  0 /
      data wg ( 43) /   2.9994938483 9685626832 4124518299 68724 d  0 /
      data wg ( 44) /   3.0929995346 9357468116 6951083530 33660 d  0 /
      data wg ( 45) /   3.1887641899 4712376429 3652715016 23466 d  0 /
      data wg ( 46) /   3.2869445597 5337531998 3781070122 16956 d  0 /
      data wg ( 47) /   3.3877119796 0397652334 0549087621 54571 d  0 /
      data wg ( 48) /   3.4912542659 8732012281 7324237827 64895 d  0 /
      data wg ( 49) /   3.5977779176 9613046096 2947301749 02943 d  0 /
      data wg ( 50) /   3.7075106900 1745708341 0271556592 28179 d  0 /
      data wg ( 51) /   3.8207046196 5311695152 0299594304 67622 d  0 /
      data wg ( 52) /   3.9376395977 1430720676 8005406573 30923 d  0 /
      data wg ( 53) /   4.0586276133 8354481597 4201161879 88679 d  0 /
      data wg ( 54) /   4.1840178238 1424031850 6076923345 03121 d  0 /
      data wg ( 55) /   4.3142026492 9613425820 0845732179 87912 d  0 /
      data wg ( 56) /   4.4496251505 3655906604 9828201553 77774 d  0 /
      data wg ( 57) /   4.5907880226 3617511042 9598491489 29810 d  0 /
      data wg ( 58) /   4.7382646459 8929537394 7538735058 38770 d  0 /
      data wg ( 59) /   4.8927127796 6692168696 8869367432 83567 d  0 /
      data wg ( 60) /   5.0548916853 4039512820 5725071351 75938 d  0 /
      data wg ( 61) /   5.2256837559 4272391089 2780101660 22467 d  0 /
      data wg ( 62) /   5.4061221337 9727909853 3235123407 17863 d  0 /
      data wg ( 63) /   5.5974264018 4041404016 5536941589 80053 d  0 /
      data wg ( 64) /   5.8010493213 7643943530 6261624553 94841 d  0 /
      data wg ( 65) /   6.0187389387 8099108768 0151515140 26344 d  0 /
      data wg ( 66) /   6.2526224749 1437403092 9342134800 91928 d  0 /
      data wg ( 67) /   6.5053217351 7668675787 4827196636 96133 d  0 /
      data wg ( 68) /   6.7801152120 0777294201 2873479800 59368 d  0 /
      data wg ( 69) /   7.0811712202 5414518776 1743119167 59402 d  0 /
      data wg ( 70) /   7.4138924461 5305421421 6956062266 87752 d  0 /
      data wg ( 71) /   7.7854415484 1612700386 2327403392 30532 d  0 /
      data wg ( 72) /   8.2055734781 4596472333 9050861009 17119 d  0 /
      data wg ( 73) /   8.6880138399 6161871469 4199580582 55237 d  0 /
      data wg ( 74) /   9.2528697341 5578523923 5565062019 79918 d  0 /
      data wg ( 75) /   9.9311447184 0215736008 3709865340 09772 d  0 /
      data wg ( 76) /   1.0773973641 4646829405 7508435229 90655 d  1 /
      data wg ( 77) /   1.1873891246 5097447081 9508877108 77400 d  1 /
      data wg ( 78) /   1.3422885849 7264236139 7349401540 89734 d  1 /
      data wg ( 79) /   1.5919780161 6897924449 5542522001 85978 d  1 /
      data wg ( 80) /   2.1421454296 4372259537 5210361864 15127 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,80
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end

