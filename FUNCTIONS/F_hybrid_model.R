# source('./help_functions_analysis.R')
source('./FUNCTIONS/F_EDM_analysis.R')

require('tidyverse')
require('lubridate')
require("rEDM",lib.loc = "./lib")

oxysat = 13 # oxygen saturation at 5∞C [mg/L]
cr_dm = 250 # critical depth of mixed layer to consider lake as undergoing complete mixing
oxyR=10
g = 9.81
alpha = 100e-6
Vbot = 1.4e10
H = 310 - cr_dm
# time_hybrid[1,] = c(1981,5,1,0,0,0);
time_hybrid_0 <- ymd("1981-5-1")
# oxygen[1] = 9.8958



Simstrat_gen_warming_scenario <- function(delta_T = 0,
                                          run=F,
                                          # initial_output_path = './Outputs/LakeGeneva/',
                                          initial_par="INPUTS/simstrat_LakeGeneva_Deep_Mixing.par",
                                          initial_forcing="INPUTS/LAKEGENEVA/ForcingGeneva.dat"){
  

  
  # INPUTS:
  # delta_T - the warming scenario in degrees Celsius
  # initial_par - file path to Simstrat parameter files used for base simluation
  # initial_forcing - file path to the meteorological forcing file used for base simulation
  
  # Create a new forcing file, parameter file path, and output path
  scenario_forcing = str_replace(initial_forcing,".dat",paste0("_",delta_T,"degC.dat"))
  scenario_par = str_replace(initial_par,".par",paste0("_",delta_T,"degC.dat"))

  # The lines from the initial par file will be used again later if "run=T".
  lines_initial_par <- readLines(initial_par)
  initial_output_path = str_subset(lines_initial_par,"Path") %>% str_split("\"") %>% {.[[1]][4]}
  scenario_output_path = paste0(initial_output_path,"Scenario_",delta_T,"degC_warming/")
  
  
  if(run){
  # The forcing data has symbols in the column names that can cause trouble
  # for machines in some regions. Thus we first copy over the original header.
  forcing_header <- read.table(file = initial_forcing,header = F,nrows = 1,sep="\t")
  write.table(forcing_header,file=scenario_forcing,quote=F,col.names=F,row.names=F,sep="\t")
  
  # Now we import the forcing data, keeping only the initial part of colum
  # names.
  
  data_forcing <- read.table(file=initial_forcing,sep="\t",header=T,check.names = F)
  names(data_forcing) <- str_extract(names(data_forcing),"^[:alpha:]+")
  
  # We detrend the air temperature, preserving the mean, and replace the forcing
  # with the "delta_T" offset of this, writing the result to the new forcing
  # file.
  
  Tair_detrend <- pracma::detrend(data_forcing$Tair) + mean(data_forcing$Tair)
  
  data_forcing_scenario <- data_forcing %>%
    mutate(Tair = round(Tair_detrend + delta_T,4))
  
  write.table(data_forcing_scenario,
              file=scenario_forcing,append=T,
              row.names = F,
              col.names = F,
              sep="\t")

  
  
  # Create a par file to match, including possibly creating a new directory
  # to hold the output file.

  if(!dir.exists(scenario_output_path)){
    dir.create(scenario_output_path)
  }
  
  lines_scneario_par <- gsub(pattern=initial_forcing, replacement=scenario_forcing, lines_initial_par )
  lines_scneario_par <- gsub(pattern=initial_output_path, replacement=scenario_output_path, lines_scneario_par )
  
  cat(lines_scneario_par, file=scenario_par, sep="\n")
  
  
  # Run Simstrat
  
  # if(run){
    system(paste("simstrat_windows_301.exe",scenario_par))
  }
  
  Output_Simstrat_T <- read.csv(paste0(scenario_output_path,"T_out.dat"),header = FALSE)
  
  time_mo <- ymd_h("1981-1-1-0") + hours(Output_Simstrat_T[-1,1]*24)
  time_model_1d = seq(round(first(time_mo),"days"),round(last(time_mo),"days"),by="day");
  # time_mo <- as.Date(Output_Simstrat_T[-1,1],origin="1980-1-1-0-0-0",tz="UTC")
  
  z <- -as.numeric(Output_Simstrat_T[1, -1])
  dz <- -diff(z);
  Temp <- as.matrix(Output_Simstrat_T[-1,-1])
  
  zdeep <- numeric(length(time_mo))
  Tsurf <-numeric(length(time_mo))
  
  ## Extraction of surface temperature and thermocline depth
  for(k in 1:length(time_mo)){
    dT <- diff(Temp[k,])
    n2 <- g*alpha*dT/dz 
    # [val,ind] = max(n2(1:end));
    ind <- which.max(n2)
    zdeep[k]=z[last(ind)] # thermocline depth estimated based on the max of stratification
    Tsurf[k] = Temp[k,NCOL(Temp)]
  }
  
  
  # Down-sample (with linear interpolation if slight mismatches) to 1 day.
  zdeep <- approx(x=time_mo,y=zdeep,xout = time_model_1d)$y
  T_surf <- approx(x=time_mo,y=Tsurf,xout = time_model_1d)$y
  
  h_mix = round(zdeep); # daily mixed-layer depth
  
  return(data.frame(time=time_model_1d,h_mix=h_mix,T_surf=T_surf))
  
}


Simstrat_physics_inital <- function(run = FALSE,
                                    file.output = './Outputs/LakeGeneva/T_out.dat',
                                    file.params = "./Inputs/simstrat_LakeGeneva_Deep_Mixing.par"){
  
  
  
  # INPUTS: Simstrat parameter file name
  # OUTPUT: list, [1] "mixing"
  if(run){
    system(paste("simstrat_windows_301.exe",file.params))
  }
  
  Output_Simstrat_T <- read.csv(file.output,header = FALSE)
  
  time_mo <- ymd_h("1981-1-1-0") + hours(Output_Simstrat_T[-1,1]*24)
  time_model_1d = seq(round(first(time_mo),"days"),round(last(time_mo),"days"),by="day");
  # time_mo <- as.Date(Output_Simstrat_T[-1,1],origin="1980-1-1-0-0-0",tz="UTC")
  
  z <- -as.numeric(Output_Simstrat_T[1, -1])
  dz <- -diff(z);
  Temp <- as.matrix(Output_Simstrat_T[-1,-1])
  
  zdeep <- numeric(length(time_mo))
  Tsurf <-numeric(length(time_mo))
  
  ## Extraction of surface temperature and thermocline depth
  for(k in 1:length(time_mo)){
    dT <- diff(Temp[k,])
    n2 <- g*alpha*dT/dz 
    # [val,ind] = max(n2(1:end));
    ind <- which.max(n2)
    zdeep[k]=z[last(ind)] # thermocline depth estimated based on the max of stratification
    Tsurf[k] = Temp[k,NCOL(Temp)]
  }
  
  
  # Down-sample (with linear interpolation if slight mismatches) to 1 day.
  zdeep <- approx(x=time_mo,y=zdeep,xout = time_model_1d)$y
  T_surf <- approx(x=time_mo,y=Tsurf,xout = time_model_1d)$y
  
  h_mix = round(zdeep); # daily mixed-layer depth
  
  return(data.frame(time=time_model_1d,h_mix=h_mix,T_surf=T_surf))
  
}


EDM_DO_module <- function(block_observed,block_scenario,time_0,DO_0,dt=6){
  # INPUT: time, TP, initial DO
  # OUTPUT: DO
  
  # L_model <- c('PO4_lake','PO4_epi','h_mix_model','T_surf_model','oxydeep','chl')
  L_model <- c('PO4_lake','PO4_epi','Z_thermo_model','T_surf_model','oxydeep','chl')
  
  tp_max <- dt
  month_start <- month(time_0)
  
  i_start <- which.min(abs(block_scenario$date - time_0))
  
  # i_start <- which(year(block_scenario$date) == year(time_0) & month(block_scenario$date) == month(time_0))
  
  block_sim_i <- block_scenario[i_start + 0:tp_max,]
  block_sim_i[1,'oxydeep'] <- as.numeric(DO_0)
  
  out_sim_i <-  sim_EDM_simplex_diff(block_train = block_observed %>% rename(time=date),
                                     block_sim = block_sim_i %>% rename(time=date),
                                     t_sim = tp_max,
                                     # num_neighbors = 4,
                                     tp = 1,
                                     lib_train = c(1,NROW(block_observed)),
                                     pred_sim = c(1,NROW(block_sim_i)),
                                     sim_col = 'oxydeep', # could set to which columns are NAs
                                     predictor_col=L_model,
                                     normalize = TRUE,
                                     # first_column_time = TRUE,
                                     # bootstrap = F,
                                     exclusion_radius = 0) %>%
    filter(type=="sim") %>%
    tail(1)
  
  # out_sim_i <-  sim_EDM_smap_diff(block_train = block_observed %>% rename(time=date),
  #                                    block_sim = block_sim_i %>% rename(time=date),
  #                                    t_sim = tp_max,
  #                                    num_neighbors = 0,
  #                                    tp = 1,
  #                                 theta=4,
  #                                    lib_train = c(1,NROW(block_observed)),
  #                                    pred_sim = c(1,NROW(block_sim_i)),
  #                                    sim_col = 'oxydeep', # could set to which columns are NAs
  #                                    predictor_col=L_model,
  #                                    normalize = TRUE,
  #                                    first_column_time = TRUE,
  #                                    bootstrap = F,
  #                                    exclusion_radius = 0) %>%
  #   filter(type=="sim") %>%
  #   tail(1)
  
  
  DO_1 <- out_sim_i$oxydeep 
  
  return(DO_1)
}

# EDM_DO_module(block_observed=lake_geneva_interp,block_scenario=lake_geneva_interp,time=ymd("2010-05-01"),DO_0=4.5)

hybrid_model_historic <- function(
  block_observed,
  block_scenario = block_observed, # for historical the scenario is the same as observed
  # t0=c(1981,5,1,0,0,0),
  t0="1981-5-15",    # a "yyyy-mm-dd" string or POSIX date
  DO_0=9.8958,
  n_steps=70,
  dt=6              # in months
){
  if(is.character(t0)){
    t0 = ymd(t0)
  }
  # 
  # n_steps = NROW(block_scenario)
  # df_hybrid_out <- data.frame(time = block_scenario$date[1:n_steps],
  #                             oxydeep = NA)
  
  # df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
  #                             oxydeep = NA)

  df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
                              oxydeep = NA)
  df_hybrid_out$oxydeep[1] = DO_0
  
  # Perform initializing run of Simstrat to get lake physics over model duration
  Simstrat_out <- Simstrat_physics_inital()
  
  #cut
  
  # df_cut_hybrid <- block_observed %>%
  #   select(date) %>%
  #   mutate(date_0 = lag(date,1))
  # 
  # df_max_zdeep <- pmap_dfr(df_cut_hybrid,function(date,date_0,..) {
  #   h_mix = Simstrat_out %>% filter(time > date_0 & time <= date) %>% pull(h_mix) %>% max(na.rm=T)
  #   return(data.frame(date=date,h_mix=h_mix))
  # })
  # 
  # block_train <- block_observed %>%
  #   left_join(df_max_zdeep)
  # 
  # df_cut_hybrid <- Simstrat_out %>%
  #   mutate(time_chunk)
  
  
  # Exclude observations from training data in months with deep h_mix
  
  # block_train
  # block_observed
  
  # read in high frequency Rhone data
  # Q_rhone <- R.matlab::readMat("./INPUTS/Q_rhone.mat")
  # 
  # Rhone <- as.data.frame(Q_rhone$tps) %>% setNames(c("Y","M","D","h","m","s")) %>%
  #   mutate(time = make_datetime(Y,M,D,hour=h)) %>%
  #   select(time) %>%
  #   mutate(Q = t(Q_rhone$Q))
  # 
  # save(Rhone,file="./INPUTS/Q_rhone.Rdata")
  
  load(file= "./INPUTS/Q_rhone.Rdata")
  
  Rhone <- approx(x=Rhone$time,y=Rhone$Q,xout=Simstrat_out$time,method="linear",rule=2) %>%
    data.frame() %>% rename(time=x,Q_rhone=y)
  
  k=1
  
  while(k<n_steps){
    k=k+1
    
    time_k = df_hybrid_out$time[k-1] # time at start of interval k
    month = month(time_k)
    yy = year(time_k)
    
    # time_hybrid(k,:) = ([1981,5+6*(k-1),1,0,0,0]);  
    # month = time_hybrid(k,2); 
    
    # DO_k_0 = round(df_hybrid_out$oxydeep[k-1],1); #initial DO in kth interval
    DO_k_0 = df_hybrid_out$oxydeep[k-1]
    
    # # check if model interval is over summer
    # if(month %in% 5:10){
    #   # generate block_scenario
    #   # load(paste0('scenario.C_T=0.PO4_lake=',num2str(TP),'.Rdata'))
    #   
    #   # run EDM module
    #   DO_k <- EDM_DO_module(time_0=time_k,
    #                         block_observed=block_observed,
    #                         block_scenario=block_scenario,
    #                         DO_0=DO_k_0,
    #                         dt=dt)
    #   
    #   df_hybrid_out$oxydeep[k] = DO_k
    # 
    # }else{
    #  
    {
      # if the model is over winter, invoke parametric relationships described in Schwefel et al. for 2-box oxygen model
      
      # check maximum h_mix over interval
      # I_year <- which(year(Simstrat_out$time) == yy);
      
      deep_mixing = Simstrat_out %>%
        filter(time >= time_k, time <= time_k + months(dt)) %>%
        pull(h_mix) %>%
        max()
      
      # deep_mixing = max(Simstrat_out$h_mix[I_year]);
      
      if(deep_mixing < cr_dm){
        # run EDM module over winter period without complete mixing

        DO_k <- EDM_DO_module(time_0=time_k,
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k_0,
                              dt=dt)
        
        # account for Rhone river inputs when there is no mixing mixing using parametric structure
        
        # intQ = Rhone %>%
        #   filter(time >= time_k, time <= time_k + months(dt)) %>%
        #   pull(Q_rhone)
        # Qtot = mean(intQ)*length(intQ)*86400; # multiply integrated flux by elapsed time in seconds
        # DO_k = Qtot/Vbot*oxyR +(Vbot-Qtot)/Vbot*DO_k;
        
        df_hybrid_out$oxydeep[k]=DO_k;

      }else{
        h1 = deep_mixing-cr_dm;
        h2 = 310-deep_mixing;
        # val = round(h1/H*oxysat + h2/H*DO_k_0,1); 
        DO_k = h1/H*oxysat + h2/H*DO_k_0
        
        DO_k =  EDM_DO_module(time_0=time_k,
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k,
                              dt=dt)
        
        df_hybrid_out$oxydeep[k] = DO_k
      }
      
    }
    
    
  }
  
  return(df_hybrid_out)
}



hybrid_model_scenario <- function(
  block_observed,
  block_scenario = block_observed, # for historical the scenario is the same as observed
  Simstrat_file_output = './Outputs/LakeGeneva/T_out.dat',
  # t0=c(1981,5,1,0,0,0),
  t0="1981-5-15",    # a "yyyy-mm-dd" string or POSIX date
  DO_0=9.8958,
  n_steps=70,
  dt=6              # in months
){
  if(is.character(t0)){
    t0 = ymd(t0)
  }
  # 
  # n_steps = NROW(block_scenario)
  # df_hybrid_out <- data.frame(time = block_scenario$date[1:n_steps],
  #                             oxydeep = NA)
  
  # df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
  #                             oxydeep = NA)
  
  df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
                              oxydeep = NA)
  df_hybrid_out$oxydeep[1] = DO_0
  
  # Perform initializing run of Simstrat to get lake physics over model duration
  # Simstrat_out <- Simstrat_physics_inital(file.output = Simstrat_file_output)
  # Simstrat_out <- Simstrat_physics_inital()
  
  #cut
  
  # df_cut_hybrid <- block_observed %>%
  #   select(date) %>%
  #   mutate(date_0 = lag(date,1))
  # 
  # df_max_zdeep <- pmap_dfr(df_cut_hybrid,function(date,date_0,..) {
  #   h_mix = Simstrat_out %>% filter(time > date_0 & time <= date) %>% pull(h_mix) %>% max(na.rm=T)
  #   return(data.frame(date=date,h_mix=h_mix))
  # })
  # 
  # block_train <- block_observed %>%
  #   left_join(df_max_zdeep)
  # 
  # df_cut_hybrid <- Simstrat_out %>%
  #   mutate(time_chunk)
  
  
  # Exclude observations from training data in months with deep h_mix
  
  # block_train
  # block_observed
  
  # read in high frequency Rhone data
  # Q_rhone <- R.matlab::readMat("./INPUTS/Q_rhone.mat")
  # 
  # Rhone <- as.data.frame(Q_rhone$tps) %>% setNames(c("Y","M","D","h","m","s")) %>%
  #   mutate(time = make_datetime(Y,M,D,hour=h)) %>%
  #   select(time) %>%
  #   mutate(Q = t(Q_rhone$Q))
  # 
  # save(Rhone,file="./INPUTS/Q_rhone.Rdata")
  
  # load(file= "./INPUTS/Q_rhone.Rdata")
  
  # Rhone <- approx(x=Rhone$time,y=Rhone$Q,xout=Simstrat_out$time,method="linear",rule=2) %>%
  #   data.frame() %>% rename(time=x,Q_rhone=y)
  
  k=1
  
  while(k<n_steps){
    k=k+1
    
    time_k = df_hybrid_out$time[k-1] # time at start of interval k
    month = month(time_k)
    yy = year(time_k)
    
    # time_hybrid(k,:) = ([1981,5+6*(k-1),1,0,0,0]);  
    # month = time_hybrid(k,2); 
    
    # DO_k_0 = round(df_hybrid_out$oxydeep[k-1],1); #initial DO in kth interval
    DO_k_0 = df_hybrid_out$oxydeep[k-1]
    
    # # check if model interval is over summer
    # if(month %in% 5:10){
    #   # generate block_scenario
    #   # load(paste0('scenario.C_T=0.PO4_lake=',num2str(TP),'.Rdata'))
    #   
    #   # run EDM module
    #   DO_k <- EDM_DO_module(time_0=time_k,
    #                         block_observed=block_observed,
    #                         block_scenario=block_scenario,
    #                         DO_0=DO_k_0,
    #                         dt=dt)
    #   
    #   df_hybrid_out$oxydeep[k] = DO_k
    # 
    # }else{
    #  
    {
      # if the model is over winter, invoke parametric relationships described in Schwefel et al. for 2-box oxygen model
      
      # check maximum h_mix over interval
      # I_year <- which(year(Simstrat_out$time) == yy);
      
      # deep_mixing = Simstrat_out %>%
      #   filter(time >= time_k, time <= time_k + months(dt)) %>%
      #   pull(h_mix) %>%
      #   max()
      
      
      deep_mixing <- block_scenario %>%
        filter(date >= time_k, date <= time_k + months(dt)) %>%
        pull(h_mix) %>%
        max()
      
      
      # deep_mixing = max(Simstrat_out$h_mix[I_year]);
      
      if(deep_mixing < cr_dm){
        # run EDM module over winter period without complete mixing
        
        DO_k <- EDM_DO_module(time_0=time_k,
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k_0,
                              dt=dt)
        
        # account for Rhone river inputs when there is no mixing mixing using parametric structure
        
        # intQ = Rhone %>%
        #   filter(time >= time_k, time <= time_k + months(dt)) %>%
        #   pull(Q_rhone)
        # Qtot = mean(intQ)*length(intQ)*86400; # multiply integrated flux by elapsed time in seconds
        # DO_k = Qtot/Vbot*oxyR +(Vbot-Qtot)/Vbot*DO_k;
        
        df_hybrid_out$oxydeep[k]=DO_k;
        
      }else{
        
        h1 = deep_mixing-cr_dm;
        h2 = 310-deep_mixing;
        # val = round(h1/H*oxysat + h2/H*DO_k_0,1); 
        DO_k = h1/H*oxysat + h2/H*DO_k_0
        
        DO_k =  EDM_DO_module(time_0=time_k,
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k,
                              dt=dt)
        
        df_hybrid_out$oxydeep[k] = DO_k
      }
      
    }
    
    
  }
  
  return(df_hybrid_out)
}

# hybrid_model(block_observed = lake_geneva_interp,block_scenario = lake_geneva_interp)



hybrid_model_v1 <- function(
  block_observed,
  block_scenario = block_observed, # for historical the scenario is the same as observed
  # t0=c(1981,5,1,0,0,0),
  t0="1981-5-15",    # a "yyyy-mm-dd" string or POSIX date
  DO_0=9.8958,
  n_steps=70,
  dt=6              # in months
){
  if(is.character(t0)){
    t0 = ymd(t0)
  }
  # 
  # n_steps = NROW(block_scenario)
  # df_hybrid_out <- data.frame(time = block_scenario$date[1:n_steps],
  #                             oxydeep = NA)
  
  # df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
  #                             oxydeep = NA)
  
  df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
                              oxydeep = NA)
  df_hybrid_out$oxydeep[1] = DO_0
  
  # Perform initializing run of Simstrat to get lake physics over model duration
  Simstrat_out <- Simstrat_physics_inital()
  
  #cut
  
  # df_cut_hybrid <- block_observed %>%
  #   select(date) %>%
  #   mutate(date_0 = lag(date,1))
  # 
  # df_max_zdeep <- pmap_dfr(df_cut_hybrid,function(date,date_0,..) {
  #   h_mix = Simstrat_out %>% filter(time > date_0 & time <= date) %>% pull(h_mix) %>% max(na.rm=T)
  #   return(data.frame(date=date,h_mix=h_mix))
  # })
  # 
  # block_train <- block_observed %>%
  #   left_join(df_max_zdeep)
  # 
  # df_cut_hybrid <- Simstrat_out %>%
  #   mutate(time_chunk)
  
  
  # Exclude observations from training data in months with deep h_mix
  
  # block_train
  # block_observed
  
  # read in high frequency Rhone data
  # Q_rhone <- R.matlab::readMat("./INPUTS/Q_rhone.mat")
  # 
  # Rhone <- as.data.frame(Q_rhone$tps) %>% setNames(c("Y","M","D","h","m","s")) %>%
  #   mutate(time = make_datetime(Y,M,D,hour=h)) %>%
  #   select(time) %>%
  #   mutate(Q = t(Q_rhone$Q))
  # 
  # save(Rhone,file="./INPUTS/Q_rhone.Rdata")
  
  load(file= "./INPUTS/Q_rhone.Rdata")
  
  Rhone <- approx(x=Rhone$time,y=Rhone$Q,xout=Simstrat_out$time,method="linear",rule=2) %>%
    data.frame() %>% rename(time=x,Q_rhone=y)
  
  k=1
  
  while(k<n_steps){
    k=k+1
    
    time_k = df_hybrid_out$time[k-1] # time at start of interval k
    month = month(time_k)
    yy = year(time_k)
    
    # time_hybrid(k,:) = ([1981,5+6*(k-1),1,0,0,0]);  
    # month = time_hybrid(k,2); 
    
    # DO_k_0 = round(df_hybrid_out$oxydeep[k-1],1); #initial DO in kth interval
    DO_k_0 = df_hybrid_out$oxydeep[k-1]
    
    # # check if model interval is over summer
    if(month %in% 5:10){
      # generate block_scenario
      # load(paste0('scenario.C_T=0.PO4_lake=',num2str(TP),'.Rdata'))

      # run EDM module
      DO_k <- EDM_DO_module(time_0=time_k,
                            block_observed=block_observed,
                            block_scenario=block_scenario,
                            DO_0=DO_k_0,
                            dt=dt)

      df_hybrid_out$oxydeep[k] = DO_k

    }else{

    # {
      # if the model is over winter, invoke parametric relationships described in Schwefel et al. for 2-box oxygen model
      
      # check maximum h_mix over interval
      # I_year <- which(year(Simstrat_out$time) == yy);
      
      # deep_mixing = Simstrat_out %>%
      #   filter(time >= time_k, time <= time_k + months(dt)) %>%
      #   pull(h_mix) %>%
      #   max()
      
      deep_mixing <- block_scenario %>%
        filter(date >= time_k, date <= time_k + months(dt)) %>%
        pull(h_mix_model) %>%
        max()
      
      # deep_mixing = max(Simstrat_out$h_mix[I_year]);
      
      if(deep_mixing < cr_dm){
        # run EDM module over winter period without complete mixing

        DO_k <- EDM_DO_module(time_0=time_k-months(dt),
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k_0,
                              dt=dt)
        
        # account for Rhone river inputs when there is no mixing mixing using parametric structure
        
        # intQ = Rhone %>%
        #   filter(time >= time_k, time <= time_k + months(dt)) %>%
        #   pull(Q_rhone)
        # Qtot = mean(intQ)*length(intQ)*86400; # multiply integrated flux by elapsed time in seconds
        # DO_k = Qtot/Vbot*oxyR +(Vbot-Qtot)/Vbot*DO_k;
        
        df_hybrid_out$oxydeep[k]=DO_k;
        
      }else{
        h1 = deep_mixing-cr_dm;
        h2 = 310-deep_mixing;
        # val = round(h1/H*oxysat + h2/H*DO_k_0,1); 
        DO_k = h1/H*oxysat + h2/H*DO_k_0

        DO_k =  EDM_DO_module(time_0=time_k-months(dt),
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k,
                              dt=dt)

        df_hybrid_out$oxydeep[k] = DO_k
      }
      
    }
    
    if(df_hybrid_out$oxydeep[k] < 0){
      df_hybrid_out$oxydeep[k] = 0
    }
    
  } # while(k)
  
  return(df_hybrid_out)
}


hybrid_model <- function(
  block_observed,
  block_scenario = block_observed, # for historical the scenario is the same as observed
  # t0=c(1981,5,1,0,0,0),
  t0="1981-5-15",    # a "yyyy-mm-dd" string or POSIX date
  DO_0=9.8958,
  n_steps=70,
  dt=6              # in months
){
  if(is.character(t0)){
    t0 = ymd(t0)
  }
 
  df_hybrid_out <- data.frame(time = seq(t0,by=paste(dt,"months"),length.out = n_steps),
                              oxydeep = NA)
  df_hybrid_out$oxydeep[1] = DO_0
  
  # Perform initializing run of Simstrat to get lake physics over model duration
  Simstrat_out <- Simstrat_physics_inital()
  
 
  # read in high frequency Rhone data
  load(file= "./INPUTS/Q_rhone.Rdata")
  
  Rhone <- approx(x=Rhone$time,y=Rhone$Q,xout=Simstrat_out$time,method="linear",rule=2) %>%
    data.frame() %>% rename(time=x,Q_rhone=y)
  
  k=1
  
  while(k<n_steps){
    k=k+1
    
    time_k = df_hybrid_out$time[k-1] # time at start of interval k
    month = month(time_k)
    yy = year(time_k)
    
    # time_hybrid(k,:) = ([1981,5+6*(k-1),1,0,0,0]);  
    # month = time_hybrid(k,2); 
    
    # DO_k_0 = round(df_hybrid_out$oxydeep[k-1],1); #initial DO in kth interval
    DO_k_0 = df_hybrid_out$oxydeep[k-1]
    
    # # check if model interval is over summer
 
    # deep_mixing = Simstrat_out %>%
    #   filter(time >= time_k, time <= time_k + months(dt)) %>%
    #   pull(h_mix) %>%
    #   max()
      
      deep_mixing <- block_scenario %>%
        filter(date >= time_k, date <= time_k + months(dt)) %>%
        pull(h_mix_model) %>%
        max()
      
      # deep_mixing = max(Simstrat_out$h_mix[I_year]);
      
      if(deep_mixing < cr_dm){
        # run EDM module over winter period without complete mixing
        
        DO_k <- EDM_DO_module(time_0=time_k,
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k_0,
                              dt=dt)
        
        
        # DO_k <- EDM_DO_module(time_0=time_k-months(dt),
        #                       block_observed=block_observed,
        #                       block_scenario=block_scenario,
        #                       DO_0=DO_k_0,
        #                       dt=dt)
        
        # account for Rhone river inputs when there is no mixing mixing using parametric structure
        
        # intQ = Rhone %>%
        #   filter(time >= time_k, time <= time_k + months(dt)) %>%
        #   pull(Q_rhone)
        # Qtot = mean(intQ)*length(intQ)*86400; # multiply integrated flux by elapsed time in seconds
        # DO_k = Qtot/Vbot*oxyR +(Vbot-Qtot)/Vbot*DO_k;
        
        df_hybrid_out$oxydeep[k]=DO_k;
        
      }else{
        h1 = deep_mixing-cr_dm;
        h2 = 310-deep_mixing;
        # val = round(h1/H*oxysat + h2/H*DO_k_0,1); 
        DO_k = h1/H*oxysat + h2/H*DO_k_0
        
        # DO_k =  EDM_DO_module(time_0=time_k-months(dt),
        #                       block_observed=block_observed,
        #                       block_scenario=block_scenario,
        #                       DO_0=DO_k,
        #                       dt=dt)
        
        DO_k =  EDM_DO_module(time_0=time_k,
                              block_observed=block_observed,
                              block_scenario=block_scenario,
                              DO_0=DO_k,
                              dt=dt)
        
        df_hybrid_out$oxydeep[k] = DO_k
      }
      
    if(df_hybrid_out$oxydeep[k] < 0){
      df_hybrid_out$oxydeep[k] = 0
    }
    
  } # while(k)
  
  return(df_hybrid_out)
}

