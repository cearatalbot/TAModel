#########Terrestrial-aquatic C and hydrology model####################################################################
#########################################################################################

### Differential equations
#dCw.dt = alCw-Lw                                 [g C m-2 day-1] ##Cw = wood carbon, alCw = allocation of NPP to wood, Lw = wood litter
#dCl.dt = L-Ll                                    [g C m-2 day-1] ##Cl = leaf carbon, L = leaves, Ll = leaf litter
#dCs1.dt = Lf1-Ds1                                 [g C m-2 day-1] ##Cs1 = soluble fast soil C pool, Lf1= litter input from leaves and roots, Ds1= decomposition flux
#dCs2.dt = Lf2+lout-Ds2                           [g C m-2 day-1] ##Cs2 = non-soluble fast soil C pool, Lf2= litter input from leaves and roots, lout=wood litter flux, Ds2= decomposition flux
#dCs3.dt = Bs1+Bs2-Ds3                               [g C m-2 day-1] ##Cs3 = slow soil C pool, Bs1= C buried from Cs1, Bs2= C buried from Cs2, Ds3= decomposition flux  
#dCs4.dt = Bs3-Ds4                                  [g C m-2 day-1]  ##Cs4 = passive soil C pool, Bs3=burial from Cs3, Ds4=decomposition flux
#dCdoc1.dt = Ls1+Ls2+((P-P*pctInt)/100)*Cprecip-Bdoc-LCT1-Rhdoc [g C m-2 day-1] ##Cdoc1 = upper DOC pool, Ls1= DOC leached from Cs1, Ls2= DOC leached from Cs2, 
                                                                  ##Bdoc= vertical DOC flux from upper layer, LCT1= LCT from upper layer, Rhdoc= heterotrophic respiration of upper DOC
#dCdoc2.dt = Bdoc+Ls3-LCT2-Rhdoc2                  [g C m-2 day-1] ##Cdoc2 = lower DOC pool, Ls3= DOC leached from Cs3, LCT2= LCT from lower layer, Rhdoc2= heterotrophic respiration of lower DOC
#dW1.dt = (P-P*pctInt)-Q1-Q12                               [cm water equivalence day-1] ##P = precipitation, Q1 = lateral drainage from upper water layer, D = vertical drainage
#dW2.dt = Q12-Q2-T                               [cm water equivalence day-1] ##Q2 = lateral drainage from lower water layer, T= transpiration 
#dCa.dt=LCT1*Ac + Aa*(P/100)*Cprecip-(Ca/V)*Qout-deltaA*Ca [g C]##Ca = aquatic carbon 
#dCdic1

# define model for simulation
tamStep<-function(t,S,p, DIC=TRUE, PFT="DE"){
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
    
    LAIareal=(Cl)/(SLW*Cfrac) #[m2 leaves (m ground)-2] #leaves originally set to 4 g C (ground m)-2 day-1
    
    ### supporting calculations
    #temperature, light, and vpd limitation of GPPmax
    Dtemp = max(((Tmax-Tair)*(Tair-Tmin))/(((Tmax-Tmin)/2)^2), 0) #[unitless]
    Dvpd = max(1-Kvpd*VPD^2, 0) #[unitless]
    
    LAIi = seq(0,Cl,length.out=50)/(SLW*Cfrac)  #[m2 leaves (m ground)-2]
    Ii = PAR*exp(-k*LAIi)  #[einsteins m-2 day-1]
    Dlighti = 1-exp(-(Ii*log(2)/PARhalf))  #[unitless]
    Dlightbar = mean(Dlighti)  #[unitless]
    
    #GPPmax
    Rfo = Kf*Amax   #[nmol CO2 (g leaf)-1 s-1] 
    GPPmax = Amax*Ad+Rfo  #[nmol CO2 (g leaf)-1 s-1]
    
    #GPP potential
    GPPpot = GPPmax*Dtemp*Dvpd*Dlightbar*(60*60*24)*12*(1/1e9)#[g C leaf day-1] ####changed to per day
    GPPpotAreal = GPPpot*LAIareal*SLW   #[g C (m ground)-2 day-1]
    
    #water limitation of GPP
    WUE = max(Kwue/VPD,0) #[mg CO2 (g H2O)-1]
    #Tpot = GPPpot/WUE*44*1e-6 #[g H2O (g leaf)-1 s-1]; 44=MW CO2, 1e-6=nmol to mmol
    TpotAreal=GPPpotAreal/WUE*1000*(44/12)*1e-4 #[cm H2O m-2 day-1]; (60*60*24)= seconds in a day, 1=g H2O to cm^3 H2O, 1e-4 m-2 to cm-2
    Wa=max(W2*f, 0)  #[cm day-1], from lower layer
    T = min(c(TpotAreal,Wa))  #[cm day-1]
    Dwater = ifelse(TpotAreal==0,0,T/TpotAreal)  #[unitless]
    ET=T+(P*pctInt) #cm
    
    #GPP
    GPP = GPPpotAreal*Dwater  #[g C (m ground)-2 day-1] *GPP is either equal to GPPpot or 0
    
    #respiration
    Rf = max(Rfo*Q10v^((Tair-Topt)/10), 0) #[nmol CO2 (g leaf)-1 s-1] 
    RfAreal = max(Rf*LAIareal*SLW*12*(1/1e9)*(60*60*24), 0)  #[g C (m ground)-2 day-1] 
    Rm = max(Ka*Cw*Q10v^(Tair/10)/365, 0) #[g C m^-2 day-1]; 365=days in year
    #Rh = Cs*Kh*Q10s^(Tsoil/10)/365*(W/Wc) #[g C m^-2 day-1]; 365=days in year #****W/Wc could be made more complex to account for O2 limitation of soil respiration..
    Ra=RfAreal+Rm  #[g C m^-2 day-1]
    
    #wood litter
    Lw= max(Cw*Kw/365, 0)  #[g C m^-2 day-1]; 365=days in year
    #root litter
    Lr=max(Cr*Kr/365, 0) #[g C m^-2 day-1]; 365=days in year
    #root resp
    Rr=max(Ka*Cr*Q10v^(Tair/10)/365, 0) #[g C m^-2 day-1]; 365=days in year
    
    #NPP
    NPP = GPP-Ra-Rr #[g C m^-2 day-1]
    
    # leaves and leaf litter
    if(PFT=="EG"){
      Clmax=Lmax*SLW*Cfrac #converts Lmax (LAI units) to g C m^-2
      L=max(lt*(1-(Cl/Clmax))*NPP, 0) #[g C m^-2 day-1]
      Ll=l*Cl #[g C m^-2 day-1]
    } else {
      L=Lapprox(t)  #[g C m^-2 day-1]
      Ll=Ll_approx(t)  #[g C m^-2 day-1]
    }
    
    Lg<-max((NPP-L)*ag, 0)  #[g C m^-2 day-1] proportion of NPP to reproduction
    
    #excess GPP allocation
      alCw<-(NPP-L-Lg)*aw     #[g C m^-2 day-1]
      alCr<-(NPP-L-Lg)*(1-aw)  #[g C m^-2 day-1]
      lout<-max(Ccwd*Kcwd, 0)     #[g C m^-2 day-1]
    
    #drainage of water with VIC implementation
    #implementing frozen soil, may cause weird dynamics
    #if Tsoil >0, infiltration curve used to partition snowmelt+precip
      #into runoff vs. infiltration
      #Frozen soil does not have Q12, but does have baseflow from W2 (Q2)
    if(Tsoil > 0){
      im = Wmax1*(1+bi) #[cm]
      i0 = im-im*(((im-(1+bi)*W1)/im)^(1/(1+bi))) #updates i0 for each time step, (Liang and Lettenmaier 1994) & help from Diogo on July 15, 2020
      Q12 = Ks*((W1-r)/(Wmax1-r))^((2/Bp)+3)  #[cm H20 day-1] #need Ks, r, Bp (param values)
      Pin = (P-P*pctInt+Snow)
      Q1 = ifelse(((P-P*pctInt)+Snow) > 0, ifelse((i0+((P-P*pctInt)+Snow))>= im, (P-P*pctInt+Snow)-Wmax1+W1 , (P-P*pctInt+Snow)-Wmax1+W1+Wmax1*(1-((i0+(P-P*pctInt+Snow))/im))^(1+bi)), 0) #[cm H20 day-1] #calculate drainage 
        if(W1 < Wmax1){
            export = Q1+Q12
        }else{
           export = Pin
         }
    } else{
      im = Wmax1*(1+bi) #[cm]
      i0 = im-im*(((im-(1+bi)*W1)/im)^(1/(1+bi))) #updates i0 for each time step, (Liang and Lettenmaier 1994) & help from Diogo on July 15, 2020
      Q1=ifelse((P-P*pctInt)+Snow > 0, (P-P*pctInt)+Snow, 0)
      Q12=0
      Pin=0
      export=0
    }
    Q2 = ifelse(W2 > W20, (W2-W20)/Tstar, 0) #[cm H20 day-1] W20 = reference height for W2 pool (lower pool)
    #Soil DOC
    #organic horizon
    Ds1=max((Cs1*deltaS1), 0)#[g C m^-2 day-1] deltaS1 = decomp rate; Ds1=decomp; rhoS1=resp frac of Ds1 
    Ls1 = max((Ds1*lambdaS1), 0) #[g C m^-2 day-1] leaching from fast sol.
    Rs1 = min(((Ds1*(1-lambdaS1))*rhoS1*Q10s^(Tsoil/10)*(W1/Wmax1)),(Ds1*(1-lambdaS1))) #[g C m^-2 day-1] , respiration in organic horizon
    Bs1 = max((Ds1-Rs1-Ls1), 0) #[g C m^-2 day-1] burial of particulate C
    Bdoc = max((Cdoc1/(W1*0.01)*(Q12*0.01)),0) #[g C m^-2 day-1] burial of DOC from organic horizon
    Rhdoc = max((Cdoc1*Kh*Q10s^(Tsoil/10)), 0) #[g C m^-2 day^-1]
    Rhdoc2 = max((Cdoc2*Kh*Q10s^(Tsoil/10)), 0) #[g C m^-2 day-1]
    Ds2=max((Cs2*deltaS2), 0) #[g C m^-2 day-1] decomp from non sol
    Ls2 = max((Ds2*lambdaS2), 0) #[g C m^-2 day-1] leaching from fast non.sol.
    Rs2=min(((Ds2*(1-lambdaS2))*rhoS2*Q10s^(Tsoil/10)*(W1/Wmax1)), (Ds2*(1-lambdaS2))) #[g C m^-2 day-1] respiration from non sol.
    Bs2=max(Ds2-Rs2-Ls2, 0) #[g C m^-2 day-1] burial from non.sol.
  
    
    #slow soil C
    Ds3=max((Cs3*deltaS3), 0) #[g C m^-2 day-1] deltaS3=decomp rate; Ds3=decomp; Km=resp frac of Ds3
    Ls3 = max((Ds3*lambdaS3), 0) #[g C m^-2 day-1]
    Rs3 = min(((Ds3*(1-lambdaS3))*rhoS3*Q10s^(Tsoil/10)*(W2/Wmax2)), (Ds3*(1-lambdaS3))) #[g C m^-2 day^-1] heterotrophic respiration
    #passive C pool
    Bs3 = max((Ds3-Ls3-Rs3), 0) #[g C m^-2 day-1] burial from slow to passive
    Ds4=max(Cs4*deltaS4, 0) # [g C m^-2 day-1] decomposition in passive soil C pool
    #LCT1 is drainage DOC and is expressed as load 
    if(Tsoil > 0){
      precipC=(((P-P*pctInt)+Snow)/100)*Cprecip ##g C m^-2 * (cm * 0.01 = m) = g C m3 .. /1000 m3 to L... * 1000 g to mg == mg C/L
      #so LCT1 is too high and taking too much C  
      LCT1 = ifelse(Cdoc1 > 0, (Cdoc1/((W1+Q1)*0.01))*(Q1*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1], lateral DOC from top layer
           # if(LCT1 >= Cdoc1){
           #   LCT1=Cdoc1+Ls1+Ls2+precipC-Rhdoc-Bdoc 
          #  }
      LCT2 = ifelse(Cdoc2 > 0, (Cdoc2/((W2+Q2)*0.01))*(Q2*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1] lateral DOC from bottom layer

    } else{
      LCT1=0
      LCT2=0
      precipC=0
    }
    
    #separate litter flux into soluble vs. non-soluble
    Lf1=(Lg+Lr+Ll)*fS1 #[g C m^-2 day-1]
    Lf2=(Lg+Lr+Ll)*(1-fS1) #[g C m^-2 day-1]
   
    #Soil DIC
    #DIC
    if(DIC==TRUE){
        #DIC pools receive DIC from heterotrophic soil respiration and root respiration
        #We assume that DIC flux from the lower is just lost to the atm
        #and does not interact with the upper layer. We assume that wind speed
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
          P=1  # atm
          # atmCO2 in ppm
          # air T in deg C
          K=Tair+273.15
          lnF=A1+A2*(100/K)+A3*log(K/100)+A4*(K/100)^2+S*(B1+B2*(K/100)+B3*(K/100)^2)
          Cstar=exp(lnF)*P*(atmCO2*1e-6)  # mol l-1
        #Information and equations from Wania et al. 2010; T deg C
        SCO2 = 1911-113.7*Tair+2.967*Tair^2-0.02943*Tair^3 #to estimate Schmidt number for CO2
        kCO2 = 2.07+0.215*(SCO2/600)^(-1/2) #m/day
        I1 = Rs1+Rs2+Rhdoc #g C m^-2; Inorganic C production in upper soil layer
        I2 = Rs3+Ds4+Rhdoc2+Rr #g C m^-2; Inorganic C production in lower soil layer
        Bdic = max(((Cdic1/(W1*0.01))*(Q12*0.01)),0) #[g C m^-2 day-1] burial of DOC from organic horizon
        #need to change?? vol--no longer 1 m^3
        #double check units
        if(Tsoil > 0){
          Iout1 = max((-kCO2*(Cstar-((Cdic1/12)/(W1*0.01*1000)))*12), 0) #g C m^-2 DIC emission flux, concentrations in mols/L and then flux converted to g C
          Iout2 = max((-kCO2*(Cstar-((Cdic2/12)/(W2*0.01*1000)))*12), 0) #g C m^-2 DIC emission flux,concentrations in mols/L and then flux converted to g C
          LDIC1 = max((Cdic1/((W1+Q1)*0.01))*(Q1*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1], lateral DIC from top layer
            #if(LDIC1 >= Cdic1){
          #    LDIC1=Cdic1+I1-Iout1-Bdic
          #  }
          LDIC2 = max((Cdic2/((W2+Q2)*0.01))*(Q2*0.01), 0) # g C m^-3 * m^-3 = [g C day^-1] lateral DIC from bottom layer
          } else{
          Iout2=0
          Iout1=0
          LDIC1=0
          LDIC2=0
        }
       }
    ####Aquatic####
    V = Aa * zbar #aquatic volume [m^3]
    Qin = (Q1+Q2)/100 * Ac #inflow [m^3]
    Qout=ifelse(Qin+ Aa*((P+Snow)/100)-(Evap/100)*Aa < 0, 0 ,Qin+ Aa*((P+Snow)/100)-(Evap/100)*Aa)  # [m^3] 
    
    ### differential equations
    dCw.dt=alCw-Lw
    dCl.dt=L-Ll
    dCs1.dt = Lf1-Ds1
    dCs2.dt = Lf2+lout-Ds2
    dCs3.dt = Bs1+Bs2-Ds3
    dCs4.dt = Bs3-Ds4
    dCdoc1.dt = Ls1+Ls2+precipC-Rhdoc-Bdoc-LCT1
    dCdoc2.dt = Bdoc+Ls3-LCT2-Rhdoc2
    dW1.dt = Pin-export
    dW2.dt = Q12-Q2-T
    dCa.dt = LCT1*Ac + LCT2*Ac + Aa*((P+Snow)/100)*Cprecip-(Ca/V)*Qout-deltaA*Ca
    dCr.dt = alCr-Lr
    dCcwd.dt = Lw-lout
    
  if(DIC==TRUE){
      dCdic1.dt = I1-LDIC1-Iout1-Bdic
      dCdic2.dt = I2+Bdic-LDIC2-Iout2
      
      return(list(c(dCw.dt,dCl.dt,dCs1.dt, dCs2.dt, dCs3.dt, dCs4.dt, dCdoc1.dt, dCdoc2.dt, dW1.dt, dW2.dt, dCa.dt, dCr.dt, dCcwd.dt, dCdic1.dt, dCdic2.dt),
                  c(GPP=GPP,Q1=Q1, Q2=Q2, Rf = Rf, 
                    Ra = Ra, NPP = NPP, LCT1 = LCT1, 
                    Rs3=Rs3, Dwater=Dwater, Wa=Wa, 
                    TpotAreal=TpotAreal, i0=i0, T=T, 
                    Dvpd=Dvpd, Dtemp=Dtemp, Dlightbar=Dlightbar, 
                    GPPmax=GPPmax, LAIareal=LAIareal,
                    GPPpotAreal=GPPpotAreal, WUE=WUE, 
                    V=V, Lw=Lw, Ll=Ll, L=L, Rs1=Rs1, Qout=Qout, 
                    RfAreal=RfAreal, Bdoc=Bdoc, Q12=Q12, 
                    LCT2=LCT2, Lr=Lr, alCr=alCr, lout=lout, 
                    Lg=Lg, Bs1=Bs1, Ds1=Ds1, Ds3=Ds3, Ds4=Ds4, Ls3=Ls3, 
                    Ls1=Ls1, Rhdoc=Rhdoc, Rhdoc2=Rhdoc2, Bs3=Bs3, 
                    Rs2=Rs2, Ls2=Ls2, Ds2=Ds2, Lf1=Lf1, Rm=Rm, 
                    Bs2=Bs2, ET=ET, LDIC1=LDIC1, LDIC2=LDIC2, Iout1=Iout1,
                    Iout2=Iout2, I1=I1, I2=I2)))
  } else{
      return(list(c(dCw.dt,dCl.dt,dCs1.dt, dCs2.dt, dCs3.dt, dCs4.dt, dCdoc1.dt, dCdoc2.dt, dW1.dt, dW2.dt, dCa.dt, dCr.dt, dCcwd.dt),
                  c(GPP=GPP,Q1=Q1, Q2=Q2, Rf = Rf, 
                    Ra = Ra, NPP = NPP, LCT1 = LCT1, 
                    Rs3=Rs3, Dwater=Dwater, Wa=Wa, 
                    TpotAreal=TpotAreal, i0=i0, T=T, 
                    Dvpd=Dvpd, Dtemp=Dtemp, Dlightbar=Dlightbar, 
                    GPPmax=GPPmax, LAIareal=LAIareal,
                    GPPpotAreal=GPPpotAreal, WUE=WUE, 
                    V=V, Lw=Lw, Ll=Ll, L=L, Rs1=Rs1, Qout=Qout, 
                    RfAreal=RfAreal, Bdoc=Bdoc, Q12=Q12, 
                    LCT2=LCT2, Lr=Lr, alCr=alCr, lout=lout, 
                    Lg=Lg, Bs1=Bs1, Ds1=Ds1, Ds3=Ds3, Ds4=Ds4, Ls3=Ls3, 
                    Ls1=Ls1, Rhdoc=Rhdoc, Rhdoc2=Rhdoc2, Bs3=Bs3, 
                    Rs2=Rs2, Ls2=Ls2, Ds2=Ds2, Lf1=Lf1, Rm=Rm, 
                    Bs2=Bs2, ET=ET)))
      }
      }) 
}
