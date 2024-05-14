####Terrestrial-Aquatic Model####
#By: Ceara J. Talbot
#Usage: 

#########################################################################################

### Differential equations; [units]; definition 
#dCw.dt = alCw-Lw                                               [g C m-2 day-1] ##Cw = wood carbon, alCw = allocation of NPP to wood, Lw = wood litter
#dCl.dt = L-Ll                                                  [g C m-2 day-1] ##Cl = leaf carbon, L = leaves, Ll = leaf litter
#dCs1.dt = Lf1-Ds1                                              [g C m-2 day-1] ##Cs1 = soluble fast soil C pool, Lf1= litter input from leaves and roots, Ds1= decomposition flux
#dCs2.dt = Lf2+lout-Ds2                                         [g C m-2 day-1] ##Cs2 = non-soluble fast soil C pool, Lf2= litter input from leaves and roots, lout=wood litter flux, Ds2= decomposition flux
#dCs3.dt = Bs1+Bs2-Ds3                                          [g C m-2 day-1] ##Cs3 = slow soil C pool, Bs1= C buried from Cs1, Bs2= C buried from Cs2, Ds3= decomposition flux  
#dCs4.dt = Bs3-Ds4                                              [g C m-2 day-1]  ##Cs4 = passive soil C pool, Bs3=burial from Cs3, Ds4=decomposition flux
#dCdoc1.dt = Ls1+Ls2+((P-P*pctInt)/100)*Cprecip-Bdoc-LCT1-Rhdoc [g C m-2 day-1] ##Cdoc1 = upper DOC pool, Ls1= DOC leached from Cs1, Ls2= DOC leached from Cs2, 
                                                                                ##Bdoc= vertical DOC flux from upper layer, LCT1= LCT from upper layer, Rhdoc= heterotrophic respiration of upper DOC
#dCdoc2.dt = Bdoc+Ls3-LCT2-Rhdoc2                               [g C m-2 day-1] ##Cdoc2 = middle DOC pool, Ls3= DOC leached from Cs3, LCT2= LCT from middle layer, Rhdoc2= heterotrophic respiration of middle layer DOC
#dCdoc3.dt = Bdoc2-LCT3                                         [g C m-2 day-1] ##Cdoc3 = lower DOC pool, LCT3= LCT from lower layer
#dW1.dt = (P-P*pctInt)-Q1-Q12                                   [cm water equivalence day-1] ##P = precipitation, Q1 = lateral drainage from upper water layer, D = vertical drainage
#dW2.dt = Q12-Q2-T                                              [cm water equivalence day-1] ##Q2 = lateral drainage from middle water layer, T= transpiration 
#dW3.dt = Q23-Q3                                                [cm water equivalence day-1] ##Q3 = lateral drainage from lower water layer, Q3 = export from lower water layer
#dCa.dt=LCT1*Ac + Aa*(P/100)*Cprecip-(Ca/V)*Qout-deltaA*Ca      [g C] ##Ca = aquatic carbon 
#dCdic1.dt = I1-LDIC1-Iout1-Bdic                                [g C m-2 day-1]  ##Cdic1 = upper dissolved CO2 pool, LDIC1= lateral export of dissolved CO2, Iout1= CO2 emitted to the atmosphere, Bdic1= burial of dissolved CO2 from upper to middle layerr
#dCdic2.dt = I2+Bdic-Bdic2-LDIC2-Iout2                          [g C m-2 day-1]  ##Cdic2 = middle dissolved CO2 pool, LDIC2= lateral export of dissolved CO2, Iout2= CO2 emitted to the atmosphere, Bdic2= burial of dissolved CO2 from middle to lower layerr
#dCdic3.dt = Bdic2-LDIC3                                        [g C m-2 day-1]  ##Cdic3 = lower dissolved CO2 pool, LDIC3= lateral export of dissolved CO2

# define model for simulation
tamStep<-function(t,S,p, DIC=TRUE){ #default to including dissolved CO2 export
  with(as.list(c(S,p)),{
    
    ### forcings
    PAR=PARapprox(t)  #[einsteins m-2 day-1]
    P=Papprox(t)  #[cm day-1]
    VPD=VPDapprox(t)  #[kPa]
    Tair=Tair_approx(t)  #[degrees C]
    Tsoil=Tsoil_approx(t)  #[degrees C]
    #moved leaves from this spot
    Evap=Evap_approx(t) #[cm day-1]
    Snow=Snowapprox(t) #[cm day-1]
    atmCO2=CO2_approx(t) #[ppm]
    GDDday=gdd_approx(t) #sum of growing degree days (*C)
    DOY=Doy_approx(t) #day of year
    sen=sen_approx(t) #days at end of year, where no new C is allocated to leaves
    grow=grow_approx(t)
    greenup=green_approx(t)
    
    ### supporting calculations
    #temperature, light, and vpd limitation of GPPmax
    Dtemp = max(((Tmax-Tair)*(Tair-Tmin))/(((Tmax-Tmin)/2)^2), 0) #[unitless]
    Dvpd = max(1-Kvpd*VPD^Kvpd2, 0) #[unitless]
    
    #GPPmax
    Rfo = Kf*Amax   #[nmol CO2 (g leaf)-1 s-1] 
    GPPmax = Amax*Ad+Rfo  #[nmol CO2 (g leaf)-1 s-1]
    
    LAIareal = Cl/(SLW*Cfrac) #[m2 leaves (m ground)-2]
    
    # leaves and leaf litter
    LAIi = seq(0,LAIareal,length.out=50) #[m2 leaves (m ground)-2]
    Ii = PAR*exp(-k*LAIi)  #[einsteins m-2 day-1]
    Dlighti = 1-exp(-(Ii*log(2)/PARhalf))  #[unitless]
    Dlightbar = mean(Dlighti)  #[unitless]
    
    #GPP potential
    GPPpot = GPPmax*Dtemp*Dvpd*Dlightbar*(60*60*24)*12*(1/1e9)#[g C leaf day-1] 
    GPPpotAreal = GPPpot*LAIareal*SLW   #[g C (m ground)-2 day-1]
    
    ev<-pctInt #ifelse(LAIareal > (Lmax/2), pctInt, 0) 
    #water limitation of GPP
    WUE = max(Kwue/VPD,0) #[mg CO2 (g H2O)-1] #Kwue units: mg CO2 KPa (g H2O)-1
    
    TpotAreal=ifelse(WUE > 0, GPPpotAreal/WUE*1000*(44/12)*1e-4, 0) #[cm H2O m-2 day-1]; (60*60*24)= seconds in a day, 1=g H2O to cm^3 H2O, 1e-4 m-2 to cm-2
    Wa=max(W2*f, 0)  #[cm day-1], from lower layer
    T = min(c(TpotAreal,Wa))  #[cm day-1]
    Dwater = ifelse(TpotAreal==0,0,T/TpotAreal)  #[unitless]
    ET=T+(P*ev) #cm
    
    GPP = GPPpotAreal*Dwater  #[g C (m ground)-2 day-1] 
    
    #respiration
    Rf = max(Rfo*Q10v^((Tair-Topt)/10), 0) #[nmol CO2 (g leaf)-1 s-1] 
    RfAreal = max(Rf*LAIareal*SLW*12*(1/1e9)*(60*60*24), 0)  #[g C (m ground)-2 day-1] 
    Rm = max(Ka*Cw*Q10v^(Tair/10)/365, 0) #[g C m^-2 day-1]; 365=days in year
    Ra=RfAreal+Rm  #[g C m^-2 day-1]
    
    #wood litter
    Lw= max(Cw*(Kw/365), 0)  #[g C m^-2 day-1]; 365=days in year
    #root litter
    Lr=max(Cr*(Kr/365), 0) #[g C m^-2 day-1]; 365=days in year
    #root resp
    Rr=max((Ka/365*Cr*Q10v^(Tair/10)), 0) #[g C m^-2 day-1]; 365=days in year
    
    #NPP; leaves modeled after CABLE leaf phenology (refs)
    NPP = GPP-Ra-Rr #[g C m^-2 day-1]]
    if(LAIareal < Lmax & greenup > 0){ 
      Lon = max(0.9 * NPP, 0)
      Lfall = 0
    }else if(LAIareal < Lmax & grow > 0){
      Lon= max(0.3 * NPP, 0)
      Lfall=0
    } else if(sen > 0){ #create a forcing that is 0 or 1 if GDD-GDD0 > 200 days
      Lon = 0   
      Lfall = ifelse(LAIareal > Lmin, Cl*(1/28), 0) #4 week turnover rate, 28 days
    } else if(sen==0 & grow==0 & greenup==0 | LAIareal >= Lmax){
      Lon=0
      Lfall=0
    }
    
    NPPalL = ifelse(Lon > 0, NPP-Lon, NPP)#[g C m^-2 day-1]
    
    if(NPPalL > 0){
      Lg<-max((NPPalL)*ag, 0)  #[g C m^-2 day-1] 
      #excess GPP allocation
      if(S0[[1]] != 0){ #if there is wood.. 
        alCr<-(NPPalL-Lg)*(1-aw)  #[g C m^-2 day-1]
        alCw<-(NPPalL-Lg)*aw   #[g C m^-2 day-1]
      } else{
        alCr<-(NPPalL-Lg)*(1-aw)
        alCw<-0
      }
    } else{ #if there is no excess NPP to allocate
      Lg=0
      if(S0[[1]] != 0){ #if there is wood.. 
        if(Cr > 0 & Cw > 0){
          alCr=(Cr/(Cw+Cr))*(NPPalL-Lg)
          alCw=(Cw/(Cw+Cr))*(NPPalL-Lg)
        } else if(Cr <= 0){
          alCw=NPPalL-Lg
          alCr=0
        } else if(Cw <= 0){
          alCr=NPPalL-Lg
          alCw=0
        }
      } else{
        alCr=NPPalL-Lg
        alCw=0
      }
    }
    
    #wood litter flux to soil
    lout<-max(Ccwd*Kcwd, 0) #[g C m^-2 day-1]
    
    #partition wood litter flux into soluble vs. non-soluble soil C pools
    if(Lfall > 0){
      Lf1=(Lg+Lr+Lfall)*fS1 #[g C m^-2 day-1]
      Lf2=(Lg+Lr+Lfall)*(1-fS1) #[g C m^-2 day-1]
    } else{
      Lf1=(Lg+Lr)*fS1 #[g C m^-2 day-1]
      Lf2=(Lg+Lr)*(1-fS1) #[g C m^-2 day-1]
    }
    
    
    #drainage of water with VIC implementation
    #if Tsoil >0, infiltration curve used to partition snowmelt+precip
    #into runoff vs. infiltration
    #Frozen soil does not have Q12, but does have baseflow from W2 (Q2)
    if(Tsoil > -120){
      im = Wmax1*(1+bi) #[cm]
      i0 = im-im*(((im-(1+bi)*W1)/im)^(1/(1+bi))) #updates i0 for each time step, (Liang and Lettenmaier 1994) & help from Diogo on July 15, 2020
      Q12 = max(min(Ks*((W1-r)/(Wmax1-r))^((2/Bp1)+3), (W1-0.5)), 0) #[cm H20 day-1] #need Ks, r, Bp (param values)
      Pin = P-(P*ev)+Snow
      Q1 = ifelse(((P-P*ev)+Snow) > 0, ifelse((i0+((P-P*ev)+Snow))>= im, (P-P*ev+Snow)-Wmax1+W1 , (P-P*ev+Snow)-Wmax1+W1+Wmax1*(1-((i0+(P-P*ev+Snow))/im))^(1+bi)), 0) #[cm H20 day-1] #calculate drainage 
      if((W1+Pin) < Wmax1){
        export = Q1+Q12
      }else{
        export = Pin+Q12
      }
    } else{
      im = Wmax1*(1+bi) #[cm]
      i0 = im-im*(((im-(1+bi)*W1)/im)^(1/(1+bi))) #updates i0 for each time step, (Liang and Lettenmaier 1994) 
      Q1=ifelse((P-P*ev)+Snow > 0, (P-P*ev)+Snow, 0)
      Q12=0
      Pin=0
      export=0
    }
    Q23 = max(min(Ks*((W2-r)/(Wmax2-r))^((2/Bp2)+3), (W2-0.5)),0) #[cm H20 day-1] #need Ks, r, Bp (param values)
    Q2 = ifelse(W2 > W20, (W2-W20)/Tstar, 0) #[cm H20 day-1] 
    Q3 = ifelse(W3 > W30, (W3-W30)/Tstar2, 0) #[cm H20 day-1] 
    
    #Soil DOC
    #upper soil soluble 
    Ds1=max((Cs1*deltaS1), 0)#[g C m^-2 day-1] deltaS1 = decomp rate; Ds1=decomp; rhoS1=resp frac of Ds1 
    Ls1 = max((Ds1*lambdaS1), 0) #[g C m^-2 day-1] leaching from fast sol.
    Rs1 = min(((Ds1*(1-lambdaS1))*rhoS1*Q10s^(Tsoil/10)*(W1/Wmax1)),(Ds1*(1-lambdaS1))) #[g C m^-2 day-1] , respiration in organic horizon
    Bs1 = max((Ds1-Rs1-Ls1), 0) #[g C m^-2 day-1] burial of particulate C
    Bdoc = max((Cdoc1/(W1*0.01)*(Q12*0.01)),0) #[g C m^-2 day-1] burial of DOC from organic horizon
    Rhdoc = max((Cdoc1*Kdoc*Q10s^(Tsoil/10)), 0) #[g C m^-2 day^-1]
    
    #upper soil non-soluble
    Bdoc2 = max((Cdoc2/(W2*0.01)*(Q23*0.01)),0) #[g C m^-2 day-1] burial of DOC from organic horizon
    Rhdoc2 = max((Cdoc2*Kdoc*Q10s^(Tsoil/10)), 0) #[g C m^-2 day-1]
    Ds2=max((Cs2*deltaS2), 0) #[g C m^-2 day-1] decomp from non sol
    Ls2 = max((Ds2*lambdaS2), 0) #[g C m^-2 day-1] leaching from fast non.sol.
    Rs2=min(((Ds2*(1-lambdaS2))*rhoS2*Q10s^(Tsoil/10)*(W1/Wmax1)), (Ds2*(1-lambdaS2))) #[g C m^-2 day-1] respiration from non sol.
    Bs2=max(Ds2-Rs2-Ls2, 0) #[g C m^-2 day-1] burial from non.sol.
    
    #slow soil C
    Ds3=max((Cs3*deltaS3), 0) #[g C m^-2 day-1] 
    Ls3 = max((Ds3*lambdaS3), 0) #[g C m^-2 day-1]
    Rs3 = min(((Ds3*(1-lambdaS3))*rhoS3*Q10s^(Tsoil/10)*(W2/Wmax2)), (Ds3*(1-lambdaS3))) #[g C m^-2 day^-1] heterotrophic respiration
    #passive C pool
    Bs3 = max((Ds3-Ls3-Rs3), 0) #[g C m^-2 day-1] burial from slow to passive
    Ds4=max(Cs4*deltaS4, 0) # [g C m^-2 day-1] decomposition in passive soil C pool
    
    #LCT1 is drainage DOC and is expressed as load 
    LCT2 = ifelse(Cdoc2 > 0, (Cdoc2/((W2+Q2)*0.01))*(Q2*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1] lateral DOC from middle layer
    LCT3 = ifelse(Cdoc3 > 0, (Cdoc3/((W3+Q3)*0.01))*(Q3*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1] lateral DOC from bottom layer
    
    if(Tsoil > -120){
      precipC=(((P-P*ev)+Snow)/100)*Cprecip ##g C m^-2 * (cm * 0.01 = m) = g C m3 .. /1000 m3 to L... * 1000 g to mg == mg C/L
      LCT1 = ifelse(Cdoc1 > 0, (Cdoc1/((W1+Q1)*0.01))*(Q1*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1], lateral DOC from top layer
    } else{
      LCT1=0
      precipC=0
    }
    
    #Soil DIC
    #DIC
    if(DIC==TRUE){
      #DIC pools receive DIC from heterotrophic soil respiration and root respiration
      #We assume that vertical DIC flux from the lower is lost to the atmosphere
      #and does not interact with the upper layers above it. We assume that wind speed
      #is 0 because the soil water is below the surface of the soil. 
      #Calculate the CO2 concentration in water at equilibrium with a given
      #atmospheric CO2 concentration (C*)
      A1=-160.7333
      A2=215.4152
      A3=89.8920
      A4=-1.47759
      B1=0.029941
      B2=-0.027455
      B3=0.0053407
      S=0  # ppt
      Press=1  # atm
      # atmCO2 in ppm
      K=Tair+273.15 # air T in deg C
      lnF=A1+A2*(100/K)+A3*log(K/100)+A4*(K/100)^2+S*(B1+B2*(K/100)+B3*(K/100)^2)
      Cstar=(exp(lnF)*Press*(atmCO2*1e-3))*12 #1e-6 converts mmol/m^-3 to mol l-1
      #Information and equations from Wania et al. 2010; T deg C
      SCO2 = 1911-113.7*Tair+2.967*Tair^2-0.02943*Tair^3 #to estimate Schmidt number for CO2
      kCO2 = 2.07+0.215*(SCO2/600)^(-1/2) #m/day
      I1 = Rs1+Rs2+Rhdoc #g C m^-2; Inorganic C production in upper soil layer
      I2 = Rs3+Ds4+Rhdoc2+Rr #g C m^-2; Inorganic C production in middle soil layer
      
      Bdic = max(((Cdic1/(W1*0.01))*(Q12*0.01)),0) #[g C m^-2 day-1] burial of dissolved CO2 from upper to middle
      Bdic2 = max(((Cdic2/(W2*0.01))*(Q23*0.01)),0) #[g C m^-2 day-1] burial of dissolved CO2 from middle to lower 
      
      LDIC2 = max((Cdic2/((W2+Q2)*0.01))*(Q2*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1] lateral dissolved CO2 from middle layer          
      LDIC3 = max((Cdic3/((W3+Q3)*0.01))*(Q3*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1] lateral dissolved CO2 from bottom layer          
      
      if(Tsoil > -120){
        LDIC1 = max((Cdic1/((W1+Q1)*0.01))*(Q1*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1], lateral dissolved CO2 from upper layer
      } else{
        LDIC1=0
      }
      
      Iout1 = min(-kCO2*(Cstar-(Cdic1/((W1+Q1)*0.01))), Cdic1-Bdic-LDIC1) #g C m^-2 DIC emission flux, concentrations in mols/L and then flux converted to g C
      Iout2 = min(-kCO2*(Cstar-(Cdic2/((W2+Q2)*0.01))), Cdic2+Bdic-Bdic2-LDIC2)  #g C m^-2 DIC emission flux,concentrations in mols/L and then flux converted to g C
    }
    
    ####Aquatic####
    V = Aa * zbar #aquatic volume [m^3]
    Qin = (Q1+Q2)/100 * Ac #inflow [m^3]
    Qout=ifelse(Qin+ Aa*((P+Snow)/100)-(Evap/100)*Aa < 0, 0 ,Qin+ Aa*((P+Snow)/100)-(Evap/100)*Aa)  #[m^3]; assumes constant volume 
    
    ero=er*Cs2 #erosion calculated from erosion constant and upper, non-soluble C pool
    
    ### differential equations
    dCw.dt=alCw-Rm-Lw
    dCl.dt=Lon-Lfall 
    dCs1.dt = Lf1-Ds1
    dCs2.dt = Lf2+lout-Ds2-ero
    dCs3.dt = Bs1+Bs2-Ds3
    dCs4.dt = Bs3-Ds4
    dCdoc1.dt = Ls1+Ls2+precipC-Rhdoc-Bdoc-LCT1
    dCdoc2.dt = Bdoc+Ls3-LCT2-Rhdoc2
    dCdoc3.dt = Bdoc2-LCT3 
    dW1.dt = Pin-export
    dW2.dt = Q12-Q23-Q2-T
    dW3.dt = Q23-Q3 
    dCa.dt = LCT1*Ac + LCT2*Ac + Aa*((P+Snow)/100)*Cprecip-(Ca/V)*Qout-deltaA*Ca
    dCr.dt = alCr-Rr-Lr
    dCcwd.dt = Lw-lout
    
    if(DIC==TRUE){
      dCdic1.dt = I1-LDIC1-Iout1-Bdic
      dCdic2.dt = I2+Bdic-Bdic2-LDIC2-Iout2
      dCdic3.dt = Bdic2-LDIC3 
      
      return(list(c(dCw.dt,dCl.dt,dCs1.dt, dCs2.dt, dCs3.dt, dCs4.dt, dCdoc1.dt, 
                    dCdoc2.dt, dW1.dt, dW2.dt, dCa.dt, dCr.dt, dCcwd.dt, dCdic1.dt, 
                    dCdic2.dt, dW3.dt, dCdic3.dt, dCdoc3.dt),
                  c(GPP=GPP,Q1=Q1, Q2=Q2, Rf = Rf, Ra = Ra, 
                    NPP = NPP, LCT1 = LCT1, Rs3=Rs3, Dwater=Dwater, 
                    Wa=Wa, TpotAreal=TpotAreal, i0=i0, T=T, 
                    Dvpd=Dvpd, Dtemp=Dtemp, Dlightbar=Dlightbar, 
                    GPPmax=GPPmax, LAIareal=LAIareal,
                    GPPpotAreal=GPPpotAreal, WUE=WUE, 
                    V=V, Lw=Lw, Lfall=Lfall, Lon=Lon, Rs1=Rs1, 
                    Qout=Qout, RfAreal=RfAreal, Bdoc=Bdoc, Q12=Q12, 
                    LCT2=LCT2, Lr=Lr, alCr=alCr, lout=lout, 
                    Lg=Lg, Bs1=Bs1, Ds1=Ds1, Ds3=Ds3, Ds4=Ds4, Ls3=Ls3, 
                    Ls1=Ls1, Rhdoc=Rhdoc, Rhdoc2=Rhdoc2, Bs3=Bs3, 
                    Rs2=Rs2, Ls2=Ls2, Ds2=Ds2, Lf1=Lf1, Rm=Rm, 
                    Bs2=Bs2, ET=ET, LDIC1=LDIC1, LDIC2=LDIC2, 
                    Iout1=Iout1,Iout2=Iout2, I1=I1, I2=I2, Bdic=Bdic, 
                    ero=ero, SCO2=SCO2,kCO2=kCO2, Cstar=Cstar, Q3=Q3, 
                    Q23=Q23, LCT3=LCT3, LDIC3=LDIC3, Bdic2=Bdic2)))
    } else{
      return(list(c(dCw.dt,dCl.dt,dCs1.dt, dCs2.dt, dCs3.dt, dCs4.dt, dCdoc1.dt, 
                    dCdoc2.dt, dW1.dt, dW2.dt, dCa.dt, dCr.dt, dCcwd.dt, dW3.dt),
                  c(GPP=GPP,Q1=Q1, Q2=Q2, Rf = Rf, Ra = Ra, 
                    NPP = NPP, LCT1 = LCT1, Rs3=Rs3, Dwater=Dwater, 
                    Wa=Wa, TpotAreal=TpotAreal, i0=i0, T=T, 
                    Dvpd=Dvpd, Dtemp=Dtemp, Dlightbar=Dlightbar, 
                    GPPmax=GPPmax, LAIareal=LAIareal,
                    GPPpotAreal=GPPpotAreal, WUE=WUE, 
                    V=V, Lw=Lw, Lflux=Lflux, Rs1=Rs1, Qout=Qout, 
                    RfAreal=RfAreal, Bdoc=Bdoc, Q12=Q12, LCT2=LCT2, 
                    Lr=Lr, alCr=alCr, lout=lout, Lg=Lg, Bs1=Bs1, 
                    Ds1=Ds1, Ds3=Ds3, Ds4=Ds4, Ls3=Ls3, Ls1=Ls1, 
                    Rhdoc=Rhdoc, Rhdoc2=Rhdoc2, Bs3=Bs3, Rs2=Rs2, 
                    Ls2=Ls2, Ds2=Ds2, Lf1=Lf1, Rm=Rm, 
                    Bs2=Bs2, ET=ET, ero=ero,  Q3=Q3, Q23=Q23, LCT3=LCT3)))
    }
  }) 
}
