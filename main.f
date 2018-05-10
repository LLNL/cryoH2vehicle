c  Written by Salvador Aceves, Julio Moreno-Blanco (LLNL)
c
c  MODIFICATIONS
c  09-05-2017 , modification of vessel material volume
c  11-07-2017 , paraortho=1 to consider p-o or =0 to neglect it
c  11-08-2017 , At high pressure the maximum orthohydrogen density in REFPROP is exceeded
c               and the program does not work, then hydrogen(normal) is added form REFPROP
c               and orthohydrogen properties are calculated from normalhydrogen and parahydrogen.
c               However there is not a perfect agreement between both models.
c  12-08-2017 , Modification:  New model to calculate mass and thickness of the aluminum liner and carbon
c               fiber of the tank as a function of pressure rating, diameter, volume and safety factor.
c  03-22-2018 , Adding the option to set constant Temperature at pump outlet, in addition to constant entropy
c
c----------------------------------------------------------------------------------------------
c     Code originally written and tested using gfortran
c     How to compile : gfotran *.f -o filename.exe
c
c     Code uses real gas equations of states from REFPROP (version 9.1). Therefore, the REFPROP files written 
c     in FORTRAN must be located in the same directory 
c-----------------------------------------------------------------------------------------------
c
c     Code simulates thermodynamic states of hydrogen inside an insulated
c     pressure vessel that undergoes a defined duty cycle (pard/drive/fill)
c
c     Program setup to run multiple consecutive cycles
c     by reading the time_log.txt data file
c     time_log.txt has the information of a duty cycle as: time of the day, parking time, time of driving,
c     kilometers per hour and date 

C     Program includes para-ortho conversion and its effect on vessel thermodynamics
c     also includes the mass and volume of the tam
c     runs are for constant outer volume 

c     Units are K, kPa, mol/dm^3, mole fraction, J/mol, J/mol-K, m/s
c     uPa-s, W/m-K, and N/m

c***5****1****5****2****5****3****5****4****5****5****5****6****5****7**

      program REFPROP
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      double precision :: Pwork,V_t,Dint,t_liner,t_cf,Vliner,Vcf
      double precision :: XMalum,XMcomp,XMtank,Py,alphaP,SF,sigma_a
      double precision :: sigma_h,LD,L_cyl,Li_tank,Lext_tank,Vext,Stank
      double precision :: xalum,xcomp,Ri,PDint,sigma_y,sigma_cf 
      
      parameter(ncmax=20)   !max number of components in mixture
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),f(ncmax),xp(ncmax),
     2xot(ncmax),xpp(ncmax)
      character hrf*3, herr*255
      character*255 hf(ncmax),hfmix,htype,hmix,hcomp
      character(len=100)::timeclock,hour2,minute2,secons2
      real :: start,finish,kASME
      common /Titerc/XMcomp,XMalum,XMsteel,tt,rhoH2,U,Xmass,h1,qn,
     2               x,xliq,xvap,p,wm,ucomp,ualum,Usteel,e,eo,ep,en,s,q,
     3               k,j,i
     
c     Define input parameters in a list to be read from an external file     



c     Open files needed to run code
      open(unit=11,file='input_parameters.txt')
      open(unit=12,file='info_molecule.txt')
      open(unit=14,file='Final_conditions.txt')
      open(unit=16,file='Vessel_properties.txt')

c...If the fluid files are located in a directory that the code cannot
c.....find, make a call to SETPATH
      
      call SETPATH('C:\Vehicle\fluids')
   
c      hcomp=
c     'NBS':  NIST recommendation for specified fluid/mixture
c     some allowable choices for an equation of state:
c     'FEQ':  Helmholtz free energy model
c     'BWR':  pure fluid modified Benedict-Webb-Rubin (MBWR)
c     'ECS':  pure fluid thermo extended corresponding states

c...Call SETUP to initialize the program and set the pure fluid component name
      
      i     = 2
      hf(1) = 'parahyd.fld'
      hf(2) = 'hydrogen.fld' 
      hfmix = 'hmx.bnc'
      hrf   = 'DEF'
      htype = 'EOS'
      hmix  = 'NBS'
      hcomp = 'NBS'
      
      call SETUP (i,hf,hfmix,hrf,ierr,herr)
      if (ierr.ne.0) write (*,*) herr
      

c...Get molecular weight (wm), triple point temperature, normal boiling point
c.....temperature, critical temperature, pressure, density, and compressibility
c.....factor, dipole moment, and universal gas constant

      call INFO (i,wm,ttp,tnbp,tc,pc,dc,zc,acf,dip,rgas)
      write (12,*) 'WM,ACF,DIP,TTP,TNBP   ',wm,acf,dip,ttp,tnbp
      write (12,*) 'TC,PC,DC,RGAS         ',tc,pc,dc,rgas
      
      write (14,*)'#Iteration                   mass(kg)  density(g/L)
     2  pressure(bar) temp(K)  Fill-cycles Total-kgH2 Total-used-kgH2 
     3  Total-vented-mass(kg) %Vented-mass  Average-mass-Capacity(kg)  
     4  Total-parking-time(hr) Total-driving-time-(hr) Distance(km) 
     5   Xortho(-) Ce(-)   Qpo(W) Qheat(W)    Qap(W) hour min. sec.'
     
      write (16,*)'P_rating[bar] V_int[L] D_int[m] th_liner[m] th_cf[m]	
     2  Vliner[L] VcarbonFiber[L] Mliner[kg] Mcarbofiber[kg]	
     3  MassTank[kg] Py  Pwork/Py SF sigma_a[bar] sigma_h[bar] L/D 
     4  L_cyl[m] Lint[m] Lext[m] Dext[m] Vext[L] Surface_ext[m2]
     5  Xliner[-] Xcf[-]'

    
c     Start run cycle

      k = 0
c***5****1****5****2****5****3****5****4****5****5****5****6****5****7**

 2    continue

      k = k + 1
      
c     Enter vessel properties, Pf in Kpa, rho in g/L, V_t in L, Sfill in J/g-K,
      
      read(11,*)P0,rhoH2,Pf,pfill,V_t,Sfill,Tfill,Xdistfill,rhocomp,
     2rhoalum,sigma_y,sigma_cf,Dint,SF,xo,xpump,step1,pmin,emiss,n,Tout,
     3fuelec,fill,fuelcons10,fill0,fill10,iteration,Results,ResultsSC,
     4paraortho,model_fill
      
 
      if (Results.eq.1)then
      open(unit=13,file='Results.txt')
      endif
      
      if (pf.eq.999) stop

c     Different input parameters to start to run the program
 
      x(1)   = 1 - xo                      !para fraction concentration
      x(2)   = xo                          !ortho fraction concentration
      
      !Defining pure para-H2 
      xp(1)  = 1
      xp(2)  = 0
      
      !Defining pure normal-H2 
      xot(1) = 0
      xot(2) = 1
      
      !Defining concentration of pumped H2
      xpp(1) = 1 - xpump
      xpp(2) = xpump
            
      P0         = P0*100               !Initial pressure, kPa
      Pwork      = Pf                   !Pressure vessel rating, bar
      Pf         = Pf*100               !Pressure vessel rating, kPa
      pfill      = pfill*100            !Filling pressure, kPa
      Pvent      = Pf                   !Venting pressure set 50 bar above refueling pressure
      Xdistfill  = Xdistfill            !Remaining range in meters when refueling is done (original 80000)
      Sfill      = Sfill*wm             !Entropy of hydrogen filling the vessel, J/K mol
      Tamb       = 300                  !Ambient temperature
      vent       = 1.0                  !Vented mass step size, mol-g
      time       = 0.0                  !Initial time
      Xvent1     = 0.0
      fillCycles = 1
      ventCycles = 0.0
      Qpo        = 0.0                  !Initial energy convertion p-o,

***********************************************************************************************
c     Calculate mass and volume of composte vessel 
***********************************************************************************************
            
      PI=4.D0*DATAN(1.D0)                                      !Defining Pi value, 3.14159
      sigma_y  = sigma_y*10.                                   !Aluminum yield stress, bar
      sigma_cf = sigma_cf*10.                                  !Carbon fiber Failure stress, bar
      
C     Correlation to determine the ratio Pwork/Py as a function of the Pressure rating and diameter
c     which was obtained from the data of different tank designs manufacturated by Worthington industries
c     and Luxfer

          Ri = Dint/2.                                             
      PDint  = Pwork*Dint
      alphaP = 7.473981244766e-09*(PDint)**3 -   
     2         1.383983628075e-05*(PDint)**2 +                      !Ratio Pwork/Py
     3         8.878326882505e-03*(PDint) - 1.134328561731e-02
     
      Py     = (Pwork)/alphaP                                       !Yield pressure, bar
      
C     Calculate of liner and carbon fiber thickness. As a good aproximation it is assumed that
c     the thickness of the shell and the elipsoidal head are the same.
      
      t_liner = (3**0.5)/4*Py*Ri/sigma_y                            !Liner thickness, m
      
      sigma_a = Py*Ri/(2.*t_liner)                                  !Axial stress, bar
      
      Radical = 0.0                                                 !From yield criterion: 4.*sigma_y**2 - 3.*sigma_a**2=0
      
      sigma_h = 0.5*(sigma_a + (Radical)**0.5)                      !Hoop stress, bar
      
      t_cf    = (Pwork*SF*Ri - sigma_h*t_liner)/sigma_cf            !Carbon fiber thickness, m
         

C     Calculate mass and volume of aluminum in vessel      

      
      LD        = (4*(V_t/1000)/(PI*(Dint**3))) - (1./3)
       
      L_cyl     = Dint*LD                                            !Cylinder length, m
      
      Li_tank   = L_cyl + (Dint/2.)                                  !Tank internal length, m
      
      Lext_tank = Li_tank + 2.*(t_liner + t_cf)                      !Tank external length, m 

      Dm_liner  = Dint + t_liner                                     !Liner mean diameter of shell, [m]
 
      Dm_cf     = (Dint + 2.*(t_liner)) + t_cf                       !Carbon fiber mean diameter of shell, [m]
 
      Vliner    = (PI*Dm_liner*L_cyl*(t_liner) 
     2            + 2.*1.084*(Dm_liner**2)*(t_liner))                !Liner Volume [m3]
 
      Vcf       = (PI*Dm_cf*L_cyl*(t_cf) 
     2            + 2.*1.084*(Dm_cf**2)*(t_cf))                      !Carbon fiber Wall Volume [m3]

      Vext_l    = V_t + Vliner                                       !Liner External volume [L]

      Vext      = Vext_l + Vcf                                       !External Volume [l]

      XMalum    = Vliner*rhoAlum                                     !Aluminum liner mass, kg
      
      XMcomp    = Vcf*rhocomp                                        !Carbon fiber mass in vessel, kg
      
      Xmsteel   = 0.0                                                !stainless steel mass in vessel, kg
      
      XMtank    = XMcomp+ XMalum + Xmsteel                           !tank mass, kg
      
      xcomp     = XMcomp/XMtank                                      !mass fraction of composite
      
      xalum     = XMalum/XMtank                                      !mass fraction of aluminum
      
      D_ext     = Dint + 2.*(t_liner + t_cf)                         !External diameter in the shell, m
     
      Stank     = 2.18*(D_ext**2) + PI*D_ext*L_cyl                   !External surface area of the tank, m2 
      
      write (16,*)Pwork,V_t,Dint,t_liner,t_cf,Vliner*1000,Vcf*1000,
     2 XMalum,XMcomp,XMtank,Py,alphaP,SF,sigma_a,sigma_h,LD,L_cyl,
     3 Li_tank,Lext_tank,D_ext,Vext,Stank,xalum,xcomp
     

*************************************************************************************************

c     Calculate initial hydrogen mass in the tank
      
      rhoH2    = rhoH2/wm                                       !hydrogen density, moles/L
      
      Xmass    = rhoH2*V_t                                      !initial h2 mass, mol-g
      
      rho_fill = Xdistfill/fuelec*1000/V_t                      !density at which tank is refueled, g/L
      
      rho_fill = rho_fill/wm                                    !density when the car is filled, mol/L
      
      Xm_fill  = rho_fill*V_t                                   !mass when the car is filled, mol

c     calculate initial conditions

c...General property calculation with inputs of p,d,x
      p = p0
*****************************************************************************

c          Parahydrogen
       call PDFLSH (p,rhoh2,xp,tt,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
            if(ierr.gt.0)then
               write(*,*)'at 3',ierr
            endif
c          Hydrogen(normal)
      call PDFLSH (p,rhoh2,xot,tt,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
            if(ierr.gt.0)then
               write(*,*)'at 3',ierr
            endif
            
c     Correct properties of hydrogen            
      ep  = ep  - (347.2*wm)
      en  = en  + (180.5*wm)
      eo  = (en - (0.25*ep))/(0.75)      
      e   = x(1)*ep + x(2)*eo
      h1n = h1n + (180.5*wm)
      h1p = h1p - (347.2*wm)
      h1o = (h1n - (0.25*h1p))/(0.75) 
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      
*****************************************************************************     
      if((q.gt.0).and.(q.lt.1)) then
      
c        Parahydrogen 
         call SATP (p0,xp,2,t,dl,dv,xliq,xvap,ierr,herr)

         call PDFLSH (p0,dv,xp,ts,dls,dvs,xliq,xvap,qs,es,hs,ss,
     2                cvs,cps,ws,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 2',ierr
         endif
         
c        Hydrogen(normal) 
         call SATP (p0,xot,2,tsn,dlsn,dvsn,xliqsn,xvapsn,ierr,herr)

         call PDFLSH (p0,dvsn,xot,tsn,dlsn,dvsn,xliqsn,xvapn,qsn,esn,
     2               hsn,sso,cvs,cps,ws,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 2',ierr
         endif
         
c       In case of two phase, take enthalpy of the gas (venting gas)

c     Correct properties of hydrogen 
      hsn = hsn + (180.5*wm)
      hs  = hs  - (347.2*wm)
      hso = (hsn- (0.25*hs))/(0.75)
      h1  = x(1)*(hs) + x(2)*(hso)
        
      endif
*****************************************************************************
c     Determine initial conditions:
C     Start with thermal mass of aluminum and composite

      Ucomp  = composite(tt)*XMcomp
      Ualum  = aluminum(tt)*XMalum
      Usteel = steel(tt)*XMSteel

c     Total initial internal energy of the vessel
      
      U      = e*Xmass + Ucomp + Ualum + Usteel                         !Internal energy in joules
      
      Qheat  = (5.6704e-8)*(emiss/(n + 1))*Stank*(Tout**4 - Tt**4)      !Initial Heat transfer rate (t=0s), W
     2         + 0.001*(Tout - Tt)
      
      Qap    = Qheat
      
c     calculate initial equilibrium ortho concentration

      t100   = tt/100
      
      if (tt.lt.133.672) then
      xeq    =  0.2881*t100**6 - 3.445*t100**5 + 12.472*t100**4
     2           - 20.302*t100**3 + 15.688*t100**2 - 4.513*t100+0.4247
         else
      xeq    =  0.0241*t100**3 - 0.1877*t100**2 + 0.4901*t100+0.319
         endif
      
      p  = p0
      J  = 0
      td = 0
      Xd = 0
      time_old = 0
      
      write(14,*)'start',Xmass*wm/1000,rhoh2*wm,p/100,tt,
     2  0,0,0,0,0,0,0,0,0,xo,xeq,Qpo/Tpark,Qheat,Qap,0,0,0
     
**************************************************************************
c                 start overall driving-parking-refueling loop
**************************************************************************

 5    continue
 
      call cpu_time(start)

c     Open the driving scenario case
      open(unit=17,file='time_log.txt')
      
      J = J + 1
      time       = 0.0
      Xmvent     = 0.0
      Xvent1     = 0.0
      Xmass1     = Xmass
      Xmass11    = Xmass
      fillCycles = 1
      Tdriving   = 0.0
      Tparking   = 0.0
      Xdriving   = 0.0
      
      
      if(ResultsSC.eq.1)write(*,1001)'start',k,j,I,Xmass*wm/1000,
     2 rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,xmvent*wm/1000,
     3 Ucomp,Ualum,Usteel,U,xo
                     
c     start the time loop

      I = 0
      write(13,*)'iteration  K   J   I  time(hr)  distance(km)  mass(kg)
     2 density(g/L) pressure(bar) temp(K) u(J/g) up(J/g) uo(J/g) un(J/g)
     3 hout(J/g)  s(J/g-K) q(-) vented-mass(kg)  Ucomp(J)  Ualum(J)  
     4 Usteel(J) Utotal(J)    Xortho(-) Ce(-) Qheat(W) Qap(W)  Qpo(W)   
     5 hpo(kJ/kg)  K(1/h) dcdt'
       write(13,*)'start',k,j,I,td,Xd,Xmass*wm/1000,rhoh2*wm,p/100,Tt,
     2  e/wm,ep/wm,eo/wm,en/wm,h1/wm,s/wm,q,xmvent*wm/1000,Ucomp,Ualum,
     3  Usteel,U,xo,xeq,Qheat,Qap,Qpo
          
 10   continue
      Xmvent   = 0.0
      time_old = time
      
c     Enter conditions: time in days and Xdrive in Km      
      read(17,*)time, Xdrive
      
      if(time.eq.0) goto 10
      
      if (time.eq.9999) goto 50
      
      I = I + 1
      Tpark      = (time - time_old)*3600                      !number of second in a hour
      Xdrive     = Xdrive*1000                                 !Distance Travelled (meters) in a hour
      fuelcons00 = Xdrive/fuelec*1000/wm                       !total fuel consumption in a hour, moles         
      fuelcons0  = Xdrive/fuelec/(3.6*wm)                      !rate fuel consumption, moles/sec
         
*****************************************************************************
c        Drive
*****************************************************************************
      i1 = 0
      fuelcons   = fuelcons00
      fuelcons11 = 0.0
      Xmold      = Xmass
      Xmassold   = Xmass
         
      if (Xdrive.gt.0)then
      
      Tdriving   = Tdriving + (time - time_old)
      Xdriving   = Xdriving + Xdrive/1000
      step       = step1
      Xfed=0
          
 30   continue 
      i1 = i1 + 1
         
      Qheat     = (5.6704e-8)*(emiss/(n+1))*Stank*(Tout**4 - Tt**4)         !initial Heat transfer rate, W
     2             + 0.001*(Tout - Tt)
      Qap       = Qheat
         
      step2     = step1    
      fuelcons1 = fuelcons0*step
         
      if((abs(tt-tc).le.2).and.(fuelcons1.gt.1))then
      fuelcons1 = fuelcons10            
      if((i1/50)*50.eq.I1)print*,'Reduced fuelcons=',fuelcons1
      endif
         
      if((q.gt.0).and.(q.lt.998))then   
      fuelcons1 = fuelcons0*(step)
      if((i1/50)*50.eq.I1)print*,'Normal fuelcons=',fuelcons1
      elseif((qn.gt.0).and.(qn.lt.998))then
      fuelcons1 = fuelcons0*(step)
      if((i1/50)*50.eq.I1)print*,'Normal fuelcons=',fuelcons1
      endif 
          
**************************************************************************
      if(xmass.le.xm_fill)goto 31      
         
      if(xmass.gt.xm_fill+fuelcons1) goto 33
         
c        Limit mass extraction if there is not enough hydrogen in the tank 
        
         fuelcons1 = xmass - xm_fill
         step      = fuelcons1/fuelcons0
         Xmass     = Xmass - fuelcons1
         rhoh2     = Xmass/V_t
         Hout      = h1*fuelcons1
         
         U = U - Hout + Qheat*step
         call Titer
         
         Xfed       = Xmold - xm_fill
         Xdrive1    = fuelec*Xfed*(wm/1000)
         step       = Xfed/fuelcons0
         step2      = 3600 - step
         step       = step2/50
         fuelcons00 = fuelcons00 - Xfed
         fuelcons   = fuelcons00
         fuelcons0  = fuelcons00/(step2)
         fuelcons1  = (fuelcons0)*step
         Xdrive     = Xdrive-Xdrive1
         fuelcons11 = 0.0
      
         if(ResultsSC.eq.1)write(*,1001)'drive',k,j,i1,time,
     2   Xmass*wm/1000,rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,
     3   xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo
      
         write(13,*)'drive',k,j,i1,time,Xdrive1/1000,
     2   Xmass*wm/1000,rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,en/wm,h1/wm,
     3   s/wm,q,xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo,xeq
                  
 31   continue
 
c     star refueling   
         Told = tt
         Pold = p
       rhoold = rhoh2
           Xd = 0.0
      
***************************************************************************
c                  model vessel refueling with Linde pump
c                           start the fill loop
***************************************************************************
      I = 0
      
 32   continue
 
        Xmass0 = xmass
        fill1  = fill
        
         if((i/100)*100.eq.I)write(13,*)'refuel',k,j,I,time,Xd,
     2   Xmass*wm/1000,rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,en/wm,h1/wm,
     3   s/wm,q,xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo,xeq
      
      if(ResultsSC.eq.1)then 
         if((i/100)*100.eq.i)write(*,1001)'refuel',k,j,I,Xmass*wm/1000,
     2   rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,xmvent*wm/1000,
     3   Ucomp,Ualum,Usteel,U,xo
       endif
         
         i = i + 1

c        Find enthalpy of pumped H2 assuming an entropy of the pumped hydrogen constant
         
         if(model_fill.eq.1)then   !Refueling with constant Cryogenic vessel H2 inlet entropy
      
c         Parahydrogen
       call PSFLSH (p,Sfill,xp,tp,Rhop,dlpp,dvpp,xliq,xvapp,qpp,epp,
     2   hpp,cvp,cpp,wp,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 4',ierr
         endif
         
c         Hydrogen(normal)
       call PSFLSH (p,Sfill,xot,tp,Rhopn,dlppn,dvppn,xliqn,xvappn,qppn,
     2  eppn,hpn,cvpn,cppn,wpn,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 4',ierr
         endif
         
           elseif(model_fill.eq.0)then    !Refueling with constant Cryogenic vessel H2 inlet temperature
         
c         Parahydrogen
       call TPFLSH (Tfill,p,xp,Rhop,dlpp,dvpp,xliq,xvapp,qpp,epp,
     2   hpp,spump,cvp,cpp,wp,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 4',ierr
         endif

c         Hydrogen(normal)
       call TPFLSH (Tfill,p,xot,Rhopn,dlppn,dvppn,xliqn,xvappn,qppn,
     2  eppn,hpn,spump,cvpn,cppn,wpn,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 4',ierr
         endif         
                  endif
         
c         Correct the properties of pumped hydrogen
         
          hpn = hpn + (180.5*wm)
          hpp = hpp - (347.2*wm)
          hpo = (hpn- (0.25*hpp))/(0.75) 
          hp  = xpp(1)*(hpp) + xpp(2)*(hpo)
         
         
         if(abs(tt-tc).le.1) then
         fill1 = fill10
         if((i/100)*100.eq.i)print*,'Reduced fill=',fill1
         endif
         
         if((q.gt.0).and.(q.le.998)) then
         if((i/100)*100.eq.i)print*,'Normal fill=',fill0
         fill1 = fill0
         elseif((qn.gt.0).and.(qn.lt.998))then
         if((i/100)*100.eq.i)print*,'Normal fill=',fill0
         fill1 = fill0
         endif
         

         Xmass = Xmass + fill1
         rhoh2 = Xmass/V_t
         Hin   = hp*fill
         U     = U + Hin
         call Titer
         
c     Correct ortho fraction after refueling
c     assuming that vessel is filled with 99.7% para-H2
      
         if(paraortho.gt.0)then
         xoold = xo
         xo    = xo*Xmass0/Xmass + xpp(2)*(Xmass - Xmass0)/Xmass
         x(1)  = 1 - xo
         x(2)  = xo
         endif

      if (pfill.gt.p) goto 32
      
         Xmass1     = Xmass1 + Xmass - Xm_fill
         Xmass11    = Xmass11+ Xmass
         fillCycles = fillCycles + 1
      
      write(13,*)'refuel',k,j,I,time,Xd,Xmass*wm/1000,rhoh2*wm,p/100,Tt,
     2  e/wm,ep/wm,eo/wm,en/wm,h1/wm,s/wm,q,xmvent*wm/1000,Ucomp,Ualum,
     3  Usteel,U,xo,xeq
       
      if(ResultsSC.eq.1)write(*,1001)'refuel',k,j,I,time,Xmass*wm/1000,
     2 rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,xmvent*wm/1000,
     3 Ucomp,Ualum,Usteel,U,xo
           

         Xmassold = Xmass
***************************************************************************
  33   continue  
         fuelcons11 = fuelcons11 + fuelcons1
        
         if(fuelcons1.gt.(fuelcons-fuelcons11))then
         fuelcons1  = Xmass - (Xmassold - fuelcons) + 0.0000001
         endif
         
         Xmass = Xmass - fuelcons1
         rhoh2 = Xmass/V_t
         Hout  = h1*fuelcons1
         U     = U - Hout + Qheat*step
         call Titer
         
        if (p.le.pmin) then             !Limit pressure while is driving
            p  = pmin
*****************************************************************************

c          Para
       call PDFLSH (p,rhoh2,xp,tt,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
            if(ierr.gt.0)then
               write(*,*)'at 3',ierr
            endif

c          Hydrogen(normal)
      call PDFLSH (p,rhoh2,xot,tt,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
            if(ierr.gt.0)then
               write(*,*)'at 3',ierr
            endif

c     Correct properties of hydrogen 
      ep  = (ep - (347.2*wm))
      en  = en  + (180.5*wm)
      eo  = (en - (0.25*ep))/(0.75)      
      e   = x(1)*ep + x(2)*eo
      h1n = h1n + (180.5*wm)
      h1p = h1p - (347.2*wm)
      h1o = (h1n- (0.25*h1p))/(0.75)
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      
*****************************************************************************

         Ucomp  = composite(tt)*XMcomp
         Ualum  = Aluminum(tt)*XMalum
         Usteel = steel(tt)*XMSteel
         U      = e*Xmass + Ucomp + Ualum + Usteel
         
      if((q.gt.0).and.(q.lt.1)) then
      
c        Parahydrogen 
         call SATP (p,xp,2,t,dl,dv,xliq,xvap,ierr,herr)
         call PDFLSH (p,dv,xp,ts,dls,dvs,xliq,xvap,qs,es,hs,ss,
     2                cvs,cps,ws,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 2',ierr
         endif
         
c        Hydrogen(normal) 
         call SATP (p,xot,2,tn,dlsn,dvsn,xliqsn,xvapsn,ierr,herr)
         call PDFLSH (p,dvsn,xot,tsn,dlsn,dvsn,xliqsn,xvapn,qsn,esn,
     2               hsn,sso,cvs,cps,ws,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 2',ierr
         endif
         
c        In case of two phase, take enthalpy of the gas (venting gas)

c     Correct properties of hydrogen          
         hsn = hsn  + (180.5*wm)
         hs  = hs   - (347.2*wm)
         hso = (hsn - (0.25*hs))/(0.75) 
         h1  = x(1)*(hs) + x(2)*(hso)
      
      endif 
       endif

        if(ResultsSC.eq.1)then
         if((i1/50)*50.eq.I1)write(*,1001)'drive',k,j,i1,time,
     2   Xmass*wm/1000,rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,
     3   xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo
        endif
         if((i1/50)*50.eq.I1)write(13,*)'drive',k,j,i1,time,Xdrive/1000,
     2   Xmass*wm/1000,rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,en/wm,h1/wm,
     3   s/wm,q,xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo,xeq,Qheat,Qap
                
        if (Xmass.gt.(Xmassold-fuelcons)) goto 30
        
        write(13,*)'drive',k,j,i1,time,Xdrive/1000,
     2   Xmass*wm/1000,rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,en/wm,h1/wm,
     3   s/wm,q,xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo,xeq,Qheat,Qap
         
         endif
         
         if (Xdrive.eq.0) then
         
         Tparking = Tparking + (time - time_old)
         
*************************************************************************
c                              parking
*************************************************************************
         Qheat = (5.6704e-8)*(emiss/(n+1))*Stank*(Tout**4 - Tt**4) !Heat transfer rate, W
     2            +0.001*(Tout - Tt)
        
         Qday  = Qheat*Tpark                                       !Total heat transfer, J
        
         U     = U + Qday
         
c        Para-ortho conversion
c        Start calculating equilibrium ortho concentration

         t100  = tt/100
         if (tt.lt.133.672) then
         xeq   =   0.2881*t100**6 - 3.445*t100**5 + 12.472*t100**4
     2            -20.302*t100**3 + 15.688*t100**2 - 4.513*t100+0.4247
         else
         xeq   = 0.0241*t100**3  - 0.1877*t100**2 + 0.4901*t100+0.319
         endif

c        Calculate energy of p-o conversion as a function of T
c        use units of kJ/kg

         Qop   = 20.7902*t100**6 - 230.5418*t100**5 + 975.0159*t100**4
     2          -1893.6974*t100**3 + 1568.2454*t100**2 - 553.8621*t100
     3          + 593.8361
         Qop   = Qop/0.75
         hpo   = Qop
         
c        Calculate kinetic constant for ortho-para conversion
c        from milenko figure 3 as a function of T and rho

         slope    = 2.3816e-8
         xintc40K = 6.806159e-10*(rhoh2*wm)**2 - 1.228646e-8*(rhoh2*wm)
     2              +1.183884e-6
         xkop     = slope*(tt - 40) + xintc40k

c        Calculate new ortho concentration based on equilibrium and old concentration

         alpha = xeq/(1 - xeq)
         xold  = xo
         xo    = xeq*xold/(xold - (xold - xeq)*exp( -Tpark*alpha*xkop ))
         
         
         if(paraortho.eq.0)then !IF this statement is true not para-ortho conversion is considered
         xo    = 0.001
         endif
         
         x(1)  = 1-xo
         x(2)  = xo     
         dcdt  = -xkop*xo*(xo - xeq)/(1 - xeq)
         
c        Calculate p-o cooling in joules

         Qpo   = (Qop*1000*(xo - xold)*xmass*wm/1000)*paraortho
         U     = U - Qpo
         Qap   = Qheat - (Qpo/Tpark)
         
         write(33,*)'park',U,Qpo
         call Titer
         
******************************************************************************
c        limit vessel temperature to ambient temperature

         if (Tt.gt.Tamb) then
            Tt = Tamb
            
c          Parahydrogen
       call TDFLSH (tt,rhoh2,xp,pp,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
            if(ierr.gt.0)then
               write(*,*)'at 3',ierr
            endif
            
c           Hydrogen(normal)
      call TDFLSH (tt,rhoh2,xot,pn,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
            if(ierr.gt.0)then
               write(*,*)'at 3',ierr
            endif

c     Correct properties of hydrogen 
      ep  = ep   - (347.2*wm)
      en  = en   + (180.5*wm)
      eo  = (en  - (0.25*ep))/(0.75)     
      e   = x(1)*ep + x(2)*eo
      h1n = h1n  + (180.5*wm)
      h1p = h1p  - (347.2*wm)
      h1o = (h1n - (0.25*h1p))/(0.75)
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn  - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      po  = (pn  - (0.25*pp))/(0.75)
      p   = x(1)*(pp) + x(2)*(po)               
            
         Ucomp  = composite(tt)*XMcomp
         Ualum  = Aluminum(tt)*XMalum
         Usteel = steel(tt)*XMSteel
         U      = e*Xmass + Ucomp + Ualum + Usteel
         endif
******************************************************************************
         write(13,*)'park ',k,j,I,time,Xdrive,Xmass*wm/1000,rhoh2*wm,
     2   p/100,tt,e/wm,ep/wm,eo/wm,en/wm,h1/wm,s/wm,q,xmvent*wm/1000,
     3   ucomp,Ualum,Usteel,U,xo,xeq,Qheat,Qap,Qpo/Tpark,hpo,xkop*3600,
     4   dcdt,x(1),x(2)

       if(ResultsSC.eq.1)write(*,1001)'park3',k,j,I,time,Xmass*wm/1000,
     2 rhoh2*wm,p/100,tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,xmvent*wm/1000,
     3  ucomp,Ualum,Usteel,U,xo
       endif
       
***************************************************************** 
c                 check for evaporative losses
*****************************************************************
         do while (p.gt.pvent)
         
c        start venting hydrogen

         Xmass  = Xmass - vent
         rhoh2  = xmass/V_t
         Hout   = h1*Vent
         Xmvent = Xmvent + Vent
         U      = U - Hout
         call Titer

         write(13,*)'vent ',k,j,I,time,Xd,Xmass*wm/1000,rhoh2*wm,p/100,
     2   Tt,e/wm,ep/wm,eo/wm,en/wm,h1/wm,s/wm,q,xmvent*wm/1000,xo,xeq
         if(ResultsSC.eq.1)write(*,1001)'vent ',k,j,I,Xmass*wm/1000,
     2   rhoh2*wm,p/100,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,xmvent*wm/1000,
     3 xo
         end do
         
         xvent1     = xvent1 + xmvent                !Total vented mass
         ventCycles = ventCycles + 1                 !Total vented cycles
         
*******************************************************************
      goto 10 
        
   50 continue
      
         close(unit=17)
      
         call cpu_time(finish)
         timeiteration = finish - start
         hour   = timeiteration/3600
         if(timeiteration.lt.3600)then
         hour   = 0
         endif
         minute = (timeiteration - (3600.*hour))/60.
         secons = timeiteration  - (3600.*hour) - (60.*minute)
      
      write(timeclock,*)hour,minute,secons
      
      write(14,*)j,Xmass*wm/1000,rhoh2*wm,p/100,tt,fillCycles,
     2  Xmass1*wm/1000,(Xmass1-Xmass)*wm/1000,xvent1*wm/1000,
     3 xvent1*100/Xmass1,(Xmass11*wm/1000)/fillCycles,Tparking,Tdriving,
     4 Xdriving,xo,xeq,Qpo/Tpark,Qheat,Qap,timeclock
      
      if(ResultsSC.eq.1)write(*,1000)k,P/100,xdrive,Xmass*wm/1000,
     2  xm_fill*wm/1000,rhoh2*wm,Tt,e/wm,ep/wm,eo/wm,h1/wm,s/wm,q,
     3  xmvent*wm/1000,Ucomp,Ualum,Usteel,U,xo,xoold

           
      if (J.lt.iteration) goto 5
      
      goto 2
      
**********************************************************************
 1100 format (A10)
 1000 format (I4,30(1pe12.4))
 1001 format (A,3(I4),30(1pe12.4))
 1002 format (3(I4),30(1pe12.4))
      end


c***********************************************************************
c***5****1****5****2****5****3****5****4****5****5****5****6****5****7**

      subroutine Titer

c     iterative process to find vessel temperature 
c     as a function of density and internal energy

c     program modified to skip conditions after an error is detected

      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      character herr*255
      parameter (ncmax=20)   !max number of components in mixture
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),f(ncmax),xp(ncmax),
     2 xot(ncmax)

      common /Titerc/XMcomp,XMalum,XMsteel,tt,rhoH2,U,Xmass,h1,qn,
     2               x,xliq,xvap,p,wm,ucomp,ualum,Usteel,e,eo,ep,en,s,q,
     3               k,j,i


      rhoin  = rhoh2
      DelT   = 1.0
      xp(1)  = 1
      xp(2)  = 0
      xot(1) = 0
      xot(2) = 1 

c     start by trying to bracket the solution
      ii = 0
 5    continue
      ii = ii + 1  
      
      if((ii/1000)*1000.eq.ii)write(*,*)'waiting 1'
           
      T0     = TT
      Ucomp  = composite(tt)*XMcomp
      Ualum  = Aluminum(tt)*XMalum
      Usteel = steel(tt)*XMSteel

****************************************************************************
c          Para
      call TDFLSH (tt,rhoh2,xp,pp,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif
            
c           hydrogen(normal)
      call TDFLSH (tt,rhoh2,xot,pn,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
     
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif

c     Correct properties of hydrogen 
      ep  = ep   - (347.2*wm)
      en  = en   + (180.5*wm)
      eo  = (en  - (0.25*ep))/(0.75)      
      e   = x(1)*ep + x(2)*eo
      h1n = h1n  + (180.5*wm)
      h1p = h1p  - (347.2*wm)
      h1o = (h1n - (0.25*h1p))/(0.75) 
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn  - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      po  = (pn  - (0.25*pp))/(0.75)
      p   = x(1)*(pp) + x(2)*(po)  
      
****************************************************************************
     
      Utt  = Ucomp + Ualum + Usteel + e*Xmass
      DelU = Utt - U

      if(ierr.gt.0) then
      if((ii/1000)*1000.eq.ii)write(*,*)'Titer1',ierr
      if((ii/100)*100.eq.ii)write(*,*)ierr,tt,rhoh2*wm,p/100,ucomp,
     2  ualum,Usteel,e/wm,ep/wm,eo/wm,u,utt,T0,Tf,delT,delU,delU2,wm,
     3 xp(1),xot(1)
     
        TT    = TT - 0.02376
        rhoh2 = rhoin
        goto 5
      endif
     
      kk = 0
 10   continue
******************************************************
      kk = kk + 1
      if(((kk/1000)*1000.eq.kk).and.(abs(TT-tc).lt.2))Then
      write(*,*)'waiting 2:Not found solution 
     2 close to critical point'
      endif
      if(((kk/1000)*1000.eq.kk).and.(((q.gt.0).and.(q.lt.1))))Then
      write(*,*)'waiting 2:Not found solution at two phases region'
      endif
      if(kk.eq.1000)then
      print*,'Solution diverged!'
      stop
      endif
******************************************************
      
      if (DelU.lt.0) TT=TT+DelT
      if (DelU.gt.0) TT=TT-DelT
      Ucomp  = composite(tt)*XMcomp
      Ualum  = Aluminum(tt)*XMalum
      Usteel = steel(tt)*XMSteel
*************************************************************************************
c          Para
      call TDFLSH (tt,rhoh2,xp,pp,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif
            
c           hydrogen(normal)
      call TDFLSH (tt,rhoh2,xot,pn,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
     
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif

c     Correct properties of hydrogen 
      ep  = ep   - (347.2*wm)
      en  = en   + (180.5*wm)
      eo  = (en  - (0.25*ep))/(0.75)      
      e   = x(1)*ep + x(2)*eo
      h1n = h1n  + (180.5*wm)
      h1p = h1p  - (347.2*wm)
      h1o = (h1n - (0.25*h1p))/(0.75) 
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn  - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      po  = (pn  - (0.25*pp))/(0.75)
      p   = x(1)*(pp) + x(2)*(po)   
         
*************************************************************************************

      Utt   = Ucomp + Ualum + Usteel + e*Xmass
      DelU2 = Utt - U

      if(ierr.gt.0) then
      write(*,*)'Titer2',ierr
        TT    = TT - 0.02376
        rhoh2 = rhoin
        goto 5
      endif

      if(DelU2*DelU.gt.0) then
         T0   = TT
         delU = DelU2
         goto 10
      endif
      Tf = TT

c     Solution is now bracketed between T0 and Tf
c     make sure that T0<Tf

      if(T0.gt.Tf)then
         Tx    = T0
         T0    = Tf
         Tf    = Tx
         DelUx = DelU
         DelU  = DelU2
         DelU2 = DelUx
      endif

c     Now solve equation
c     start with a new guess for TT
      jj = 0
 20   continue
******************************************************
      !Modification
      rhoh2 = rhoin
      jj    = jj + 1
      
       if((jj/500)*500.eq.jj)write(*,*)'waiting 3: Not found solution 
     2 for Temperature'
       if((jj/500)*500.eq.jj)write(*,*)tt,p/100,rhoh2*wm,rhoin*wm,T0,Tf,
     2   e/wm,ep/wm,en/wm,Utt/1000,U/1000,delT,delU,delU2,q,qn
      if(jj.eq.500)then
       TT    = TT - 0.02376
       rhoh2 = rhoin
       goto 5
      endif
******************************************************
      TT = ( Tf + T0)/2

      Ucomp  = composite(tt)*XMcomp
      Ualum  = Aluminum(tt)*XMalum
      Usteel = steel(tt)*XMSteel
*********************************************************************************
c          Para
      call TDFLSH (tt,rhoh2,xp,pp,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif
            
c           hydrogen(normal)
      call TDFLSH (tt,rhoh2,xot,pn,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
     
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif

c     Correct properties of hydrogen 
      ep  = ep   - (347.2*wm)
      en  = en   + (180.5*wm)
      eo  = (en  - (0.25*ep))/(0.75)      
      e   = x(1)*ep + x(2)*eo
      h1n = h1n  + (180.5*wm)
      h1p = h1p  - (347.2*wm)
      h1o = (h1n - (0.25*h1p))/(0.75) 
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn  - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      po  = (pn  - (0.25*pp))/(0.75)
      p   = x(1)*(pp) + x(2)*(po)   
         
*********************************************************************************

      if(ierr.gt.0) then
      write(*,*)'Titer3',ierr

        TT    = TT - 0.02376
        rhoh2 = rhoin
        goto 5
      endif

      Utt   = Ucomp + Ualum + Usteel + e*Xmass
      DelUx = Utt - U
      
      if(DelUx.gt.0) then
      DelU2 = DelUx
      Tf    = TT
      else
      DelU  = Delux
      T0    = TT
      endif
      delT  = Tf - T0
      if ((abs(DelU).gt.10.0).or.(abs(DelT).gt.0.001)) goto 20
      
      mm = 0
 30   continue
      mm = mm + 1
      if((mm/1000)*1000.eq.mm)write(*,*)'waiting 4'
**********************************************************************
c          Para
      call TDFLSH (tt,rhoh2,xp,pp,dl,dv,xliq,xvap,q,ep,h1p,sp,cv,
     2              cp,w,ierr,herr)
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif
            
c           hydrogen(normal)
      call TDFLSH (tt,rhoh2,xot,pn,dln,dvn,xliqn,xvapn,qn,en,h1n,sn,cvn,
     2              cpn,wn,ierr,herr)
     
      if(ierr.gt.0)then
      if((ii/1000)*1000.eq.ii)write(*,*)'at 3',ierr
      endif

c     Correct properties of hydrogen 
      ep  = ep   - (347.2*wm)
      en  = en   + (180.5*wm)
      eo  = (en  - (0.25*ep))/(0.75)      
      e   = x(1)*ep + x(2)*eo
      h1n = h1n  + (180.5*wm)
      h1p = h1p  - (347.2*wm)
      h1o = (h1n - (0.25*h1p))/(0.75) 
      h1  = x(1)*(h1p) + x(2)*(h1o)
      so  = (sn  - (0.25*sp))/(0.75)
      s   = x(1)*sp + x(2)*so
      po  = (pn  - (0.25*pp))/(0.75)
      p   = x(1)*(pp) + x(2)*(po)   
         
**********************************************************************
      if(ierr.gt.0) then
      write(*,*)'Titer4',ierr
      write(*,*)ierr,tt,rhoh2*wm,p/100,ucomp,ualum,Usteel,e,u,utt,
     2     T0,Tf,delT,delU,delU2
        rhoh2 = rhoin
        goto 5
      endif

      if((q.gt.0).and.(q.lt.1)) then
      
c        Parahydrogen 
         call SATP (p,xp,2,t,dl,dv,xliq,xvap,ierr,herr)
         call PDFLSH (p,dv,xp,ts,dls,dvs,xliq,xvap,qs,es,hs,ss,
     2                cvs,cps,ws,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 2',ierr
         endif
         
c        Hydrogen(normal) 
         call SATP (p,xot,2,tsn,dlsn,dvsn,xliqsn,xvapsn,ierr,herr)
         call PDFLSH (p,dvsn,xot,tsn,dlsn,dvsn,xliqsn,xvapn,qsn,esn,
     2               hsn,sso,cvs,cps,ws,ierr,herr)
         if(ierr.gt.0)then
            write(*,*)'at 2',ierr
         endif
         
c       In case of two phase, take enthalpy of the gas (venting gas)

c     Correct properties of hydrogen 
        hsn = hsn  + (180.5*wm)
        hs  = hs   - (347.2*wm)
        hso = (hsn - (0.25*hs))/(0.75)
        h1  = x(1)*(hs) + x(2)*(hso)
      endif


 1000 format (A,I4,20(1pe12.4))
 1001 format (A,3(I4),20(1pe12.4))
 1002 format (3(I4),20(1pe12.4))
      return
      end
***************************************************************************
     
c***5****1****5****2****5****3****5****4****5****5****5****6****5****7**

      Function Composite(T)
c     calculates internal energy of the composite as a function of temperature
c     temperature is in K, U is in J/kg composite

      implicit double precision (a-h,o-z)
      U = (11.5488036*t + 3.11030698*t**2/2 - 0.0388055309*t**3/3 +
     2       0.000771849531*t**4/4 - 5.15153348e-6*t**5/5 + 
     3       1.42533974e-8*t**6/6 - 1.42078589e-11*t**7/7)
      composite = U
      return
      end



c***********************************************************************
c***5****1****5****2****5****3****5****4****5****5****5****6****5****7**

      Function Aluminum(T)
c     calculates internal energy of aluminum as a function of temperature
c     temperature is in K, U is in J/kg aluminum

      implicit double precision (a-h,o-z)
      U = (46.6553*t - 7.04488*t**2/2 + 0.296914*t**3/3 -
     2      0.00288308*t**4/4 + 1.35176e-5*t**5/5 -
     3      3.17364e-8*t**6/6 + 2.98673e-11*t**7/7)
      Aluminum = U
      return
      end



c***********************************************************************
c***5****1****5****2****5****3****5****4****5****5****5****6****5****7**

      Function Steel(T)
c     calculates internal energy of steel as a function of temperature
c     temperature is in K, U is in J/kg steel 304

c     equation obtained from data in Flynn's "cryogenics" page 307

      implicit double precision (a-h,o-z)

      U=(13.428*t-1.83278*t**2/2+0.071881*t**3/3-
     2      0.000520465*t**4/4+ 1.55678e-6*t**5/5-
     3      1.70401e-9*t**6/6)
      Steel=U
      return
      end
