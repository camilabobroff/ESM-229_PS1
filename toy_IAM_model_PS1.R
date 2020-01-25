# ESM 229: Economics and policy of climate change
# Problem Set 1
# This script builds a toy IAM model with no abatement and welfare. Based loosely on David Anthoff's assignments
# Author: Kyle Meng
# Date: 10/7/14. Updated: 10/1/15.


###################### Part 0: Program file directory, initial conditions, parameters ###################### 

####### Part 0a: Set file directory, clear workspace ######

    # #set your own directory here 
    # setwd('~/Documents/Dropbox/work/teaching/UCSB/2016_2017/ESM_229/problem_sets/Problem_set1')
    # rm(list=ls()) # clears workspace

#######  Part 0b: Model parameters ###### 

    # Fixed Economic parameters (LEAVE THIS ALONE FOR NOW)
    p_saving<-c(.22) # output saved, default at 0.22
    p_capShare<-c(.3) # Cobb-Douglas capital share, defaut at 0.3
    p_capDep<-c(.1) # capital depreciation, default at 0.1
    p_TFP_g<-(.005) # TFP (total factor productivity) growth rate, default at 0.005
    p_damage<-c(.009) # default at 0.009

    # Variable Economic parameters (TRY DIFFERENT VALUES)
    p_emInt_g<-c(-.005) # default at -.005 (emissions intensity)
    p_enInt_g<-c(-.01) # default at -.01 (energy intensity)
    p_pop_g<-c(0.01) # default at .01
    
    # Fixed climate parameters (LEAVE THIS ALONE)
    p_carDecay<-400 # carbon decay, default 400
    p_climDelay<-c(.015) # climate decay, default 0.015

    #Variable climate parameters (TRY DIFFERENT VALUES)
    p_climSens<-c(2.5) # climate sensitivity default 2.5

###### Part 0c: Initial conditions (in 2010) ###### 
    
    # years
    p_years<-c(291)

    # Economic
    ic_pop<-c(6900) # in millions
    ic_enInt<-c(5.98) # in EJ/$trill
    ic_emInt<-c(68.27) # in Mt CO2 /EJ            
    ic_K<-c(139.65) # $trill   
    ic_TFP<-c(.03)

    # Climate
    ic_preCO2<-c(275) # ppm
    ic_nowCO2<-c(380) # ppm
    ic_T<-c(0.8) # degree C above preindustrial

    # creating empty matrices for storing data
    year<-matrix(0,p_years,1)    
    pop<-matrix(0, p_years,1)
    enInt<-matrix(0, p_years,1)
    emInt<-matrix(0, p_years,1)    
    TFP<-matrix(0, p_years,1)
    Y_gross<-matrix(0, p_years,1)
    K<-matrix(0, p_years,1)
    gdppc<-matrix(0, p_years,1)
    CO2ppm<-matrix(0, p_years,1)
    em_MtCO2<-matrix(0,p_years,1)
    em_ppm<-matrix(0,p_years,1)
    Y_net<-matrix(0,p_years,1)
    C<-matrix(0,p_years,1)
    C_pc<-matrix(0,p_years,1)
    Teq<-matrix(0,p_years,1)
    climateDamage<-matrix(0,p_years,1)

############################ Part 1: economic model ################################## 
####### Part 1a: Getting exogenous changes in population, energy intensity, emissions intensity, TFP####### 


    # initializing for first period    
    year[1]<-2010
    pop[1]<-ic_pop
    enInt[1]<-ic_enInt    
    emInt[1]<-ic_emInt
    TFP[1]<-ic_TFP

    for (j in 2:p_years) {
      year[j]<-year[j-1]+1
      pop[j]=pop[j-1]*(1+p_pop_g)
      enInt[j]=enInt[j-1]*(1+p_enInt_g)
      emInt[j]=emInt[j-1]*(1+p_emInt_g)
      TFP[j]=TFP[j-1]*(1+p_TFP_g)
    }


######## Part 1b: Getting gross output, capital, emissions ############## 

    # initializing for first period
    #Part 1: Economic model w/o climate damages
    K[1]<-ic_K
    Y_gross[1]<-TFP[1]*K[1]^p_capShare*pop[1]^(1-p_capShare)
      
    #Part 2: climate model
    #CO2 concentration
    CO2ppm[1]<-ic_nowCO2                      
    T[1]=ic_T

    #Part 3: Economic model with climate damages
    #climate damages
    climateShare<-(p_damage*T[1]^2)/(1+p_damage*T[1]^2) # damage function
    climateDamage[1]<-climateShare*Y_gross[1] #climate damages
    Y_net[1]<-Y_gross[1]-climateDamage[1] # net of damages output          
    C[1]<-Y_net[1]-p_saving*Y_net[1]
    C_pc[1]<-(C[1]/pop[1])*1000      


    # looping over rest of the years
    for (j in 2:p_years) {
      
      #Part 1: Economic model w/o climate damages     
      #Economic output 
      K[j]<-(1-p_capDep)*K[j-1]+p_saving*Y_net[j-1]  #capital accumulation
      Y_gross[j]<-TFP[j]*K[j]^p_capShare*pop[j]^(1-p_capShare) #Cobb-Douglas production function
      gdppc[j]=Y_gross[j]/pop[j]

      #emissions
      em_MtCO2[j]<-pop[j]*gdppc[j]*enInt[j]*emInt[j] #Kaya identity
      em_ppm[j]<-em_MtCO2[j]/7810 # convert to atmospheric CO2 concentrations
      
      #Part 2: climate model     
      #CO2 concentration
      CO2ppm[j]<-CO2ppm[j-1]+0.5*((em_ppm[j])- (CO2ppm[j-1]-(ic_preCO2))/(p_carDecay)) #CO2 concentration accumulation                             
    
      #Radiative forcing
      rad_CO2_Wpm2<-5.35*log(CO2ppm[j]/ic_preCO2)

      # Converting to temperature  
      Teq[j]<-rad_CO2_Wpm2*(p_climSens/(5.3*log(2))) #Equilibrium temperature
      T[j]=T[j-1]+p_climDelay*(Teq[j]-T[j-1]) #Transient temperature

      #Part 3: Economic model with climate damages
      climateShare<-(p_damage*T[j]^2)/(1+p_damage*T[j]^2) # damage function
      climateDamage[j]<-climateShare*Y_gross[j] # climate damages
      Y_net[j]<-Y_gross[j]-climateDamage[j] # output net of damages          
      C[j]<-Y_net[j]-p_saving*Y_net[j] # consumption 
      C_pc[j]<-(C[j]/pop[j])*1000 #per capita consumption     
    }
  


###################### Part 4: plots ###################### 

par(mfrow=c(3,2), mar=c(4,4,2,2))
plot(year[2:j], em_MtCO2[2:j], ylab="Emissions (MtCO2")
plot(year[2:j], CO2ppm[2:j], ylab="CO2 concentrations (ppm)")
plot(year[2:j], T[2:j], ylab="Transient temperature (C)")
plot(year[2:j], Teq[2:j], ylab="Equilibrium temperature (C)")
plot(year[2:j], Y_net[2:j], ylab="Net output ($ trillion)")
plot(year[2:j], C_pc[2:j], ylab="Consumption per cap ($ thousand)")









