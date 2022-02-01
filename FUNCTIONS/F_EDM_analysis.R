require(rEDM,lib.loc="./lib")
require(tidyverse)

do_mEDM_models <- function(block=lake_geneva_interp_norm,L_models,target_var='oxydeep',tp=0,exclusion_radius=6){
  
  theta_list <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  
  results_mEDM <- map(L_models,function(L){
    
    multi_stats <- map_df(theta_list, function(theta.i)
      block_lnlp(block = block,
                 target_column = match(target_var,names(block)),
                 columns = match(L,names(block)),
                 tp = tp,
                 lib = c(1,NROW(block)),
                 pred = c(1,NROW(block)),
                 method = 's-map',
                 theta = theta.i,
                 num_neighbors = 0,
                 exclusion_radius = exclusion_radius,
                 stats_only = TRUE) )
    
    theta.i <- multi_stats$theta[which.max(multi_stats$rho)]
    
    multi_preds <- block_lnlp(block = block,
                              target_column = match(target_var,names(block)),
                              columns = match(L,names(block)),
                              tp = tp,
                              method = 's-map',
                              theta = theta.i,
                              num_neighbors = 0,
                              exclusion_radius = exclusion_radius,
                              stats_only = FALSE)$model_output[[1]]
    
    
    
    return(list(multi_stats=multi_stats,
                multi_preds=multi_preds))
  } # function(L)
  ) # map(L_models)
  
  
  results_mEDM <- transpose(results_mEDM)
  
} #do_mEDM_models


do_mEDM_greedy <- function(block,embed0,L_variables,target_col=1,tp=1,max_E=length(L_variables),...){
  
  result_embed0 <- do_mEDM_models(block=block, list(embed0),target_var=target_col,tp=tp,...) %>%
    {do.call(bind_rows,.$multi_stats)} %>%
    mutate(embedding_label=label_embeddings(embedding,
                                            names(block),
                                            dictionary = dictionary_extended) %>% as.factor()) %>%
    mutate(embedding=embedding_int_to_chr(embedding,names(block)))
  
  g_steps <- vector(mode = 'list')
  
  while(T){
    
    L_models_greedy <- map(setdiff(L_variables,embed0),~c(embed0,.))
    
    df_RESULT_i <- do_mEDM_models(block=block, L_models_greedy,target_var=target_col,tp=tp,...) %>%
      { do.call(bind_rows,.$multi_stats)} %>%
      mutate(embedding_label=label_embeddings(embedding,
                                              names(block),
                                              dictionary = dictionary_extended) %>% as.factor()) %>%
      mutate(embedding=embedding_int_to_chr(embedding,names(block)))
    
    g_step_i <- df_RESULT_i %>%
      ggplot(aes(x=theta,y=rho,color=embedding_label)) + 
      geom_line(lwd=.75) + 
      geom_line(data=result_embed0,lwd=1) +
      theme_bw()
    
    print(g_step_i)
    
    g_steps <- c(g_steps,g_step_i)
    
    if(max(result_embed0$rho) > max(df_RESULT_i$rho)){
      break
    } # if (max rho)
    
    embed0 <- df_RESULT_i %>% top_n(1,rho) %>% pull(embedding) %>% unlist
    
    result_embed0 <- df_RESULT_i %>% filter(paste(embedding)==paste(list(embed0)))
    
    if(length(embed0) >= max_E){
      break
    } # if length(embed0)
    
  } #  while(T)
  
  
  return(result_embed0)
  
  
} # function(do_mEDM_greedy)

make_crossval_lib <- function(nrow,n_fold=4){
  map(1:n_fold,~c(floor((.-1)*nrow/n_fold)+1, floor((.)*nrow/n_fold)  )) %>% do.call(rbind,.)
} # function(make_crossval_lib)

sim_EDM_simplex_diff <- function(block_train,
                                 block_sim,
                                 t_sim,
                                 tp = 1,
                                 num_neighbors="e+1",
                                 lib_train = c(1,NROW(block_train)),
                                 pred_sim = c(1,NROW(block_sim)),
                                 sim_col = NULL, # could set to which columns are NAs
                                 normalize = TRUE,
                                 keep_trajectories = TRUE,
                                 first_column_time = TRUE,
                                 predictor_col = 1:( NCOL(block_train) - first_column_time ),
                                 bootstrap = FALSE,
                                 ...){
  
  # NOTES:
  # - sim_col are the columns (names or indices) to iteratively predict with s-maps.
  #    > the predictions in this version of the code are done on the first-differences.
  # - The set of points for testing L_test must be at least t_sim points long. This ensures that
  # when fix_col is specified there is information for the entire simulation.
  # - bootstrap specified as either a logical or an integer. If TRUE, will be set to floor(NROW(block_train)/2).
  # - We are predicting differences in the sim_col variables, but not using differences in the embeddings.
  
  library('dplyr')
  
  if(!identical(names(block_sim),names(block_train))){
    warning("FATAL ERROR: columns of block_train and block_sim do not match")
    return()
  }
  
  if(NROW(block_sim) < t_sim){
    warning("FATAL ERROR: t_sim set longer than block_sim")
    return()
  }
  
  if(bootstrap){
    I_lib_train <- apply(rbind(lib_train),1,function(v) seq(v[1],v[2]-tp)) %>% unlist()
    lib_train <- do.call(rbind,map(I_lib_train,~c(.,.+tp)))
    
    if(identical(bootstrap,TRUE)){
      bootstrap <- floor(NROW(lib_train)/2)
    }
  }
  
  if(!is.character(sim_col)){
    sim_col <- names(block_sim)[sim_col]
  }
  
  # add time column if not specified
  if(!first_column_time){
    block_train <- block_train %>%
      mutate(time = 1:n()) %>%
      select(time,everything())
    block_sim <- block_sim %>%
      mutate(time = NROW(block_train) + (1:n()) ) %>%
      select(time,everything())
  }
  
  n_var <- NCOL(block_train) - 1
  
  # Normalize and take differences
  v_norms <- block_train %>%
    summarise_all(list(~sd(.,na.rm=T)))
  
  block.t_norm <- bind_rows( block_train %>%
                               mutate_at( set_names(sim_col,sim_col), list(delta = ~ c(NA,diff(.)))  ) %>%
                               mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]/col_norm,.)) ,
                             block_sim %>%
                               mutate_at( set_names(sim_col,sim_col), list(delta = ~ c(NA,diff(.)))  ) %>%
                               mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]/col_norm,.))
  )

  #### ITERATE THROUGH TIME
  t_pred <- 1
  while(t_pred < t_sim + 1){
    t_pred <- t_pred+1
    delY <- as.data.frame(array(dim=c(1,n_var),dimnames=list(NULL,names(block_train)[-1])))
    
    if(bootstrap){
      lib_temp <- lib_train[sample(1:NROW(lib_train),size = bootstrap,replace=F),]
    }else{
      lib_temp <- lib_train
    }
    
    #### EVALUATE SIMPLEX
    for(i.col in sim_col){
      
      out.temp <- block_lnlp(block.t_norm,
                             method='simplex',
                             lib=lib_temp,
                             pred=c((NROW(block_train))+(t_pred-1),(NROW(block_train))+t_pred),
                             columns=predictor_col,
                             target_column = paste0(i.col,"_delta"),
                             num_neighbors = num_neighbors,
                             stats_only = FALSE,
                             silent = TRUE,
                             first_column_time = TRUE,
                             ...)
      
      delY[i.col] <- out.temp$model_output[[1]]$pred[1]
      
      
      #### REPLACE IN BLOCK BY ADDING TO x_i.col(t_pred-1)
      block.t_norm[NROW(block_train)+t_pred,paste0(i.col,"_delta")] <- delY[i.col]
      block.t_norm[NROW(block_train)+t_pred,i.col] <- block.t_norm[NROW(block_train)+t_pred-1,i.col] + delY[i.col] 
      
    } # for(i.col)
  }
  
  # Undo normalize
  block_sim_out <- block.t_norm[(NROW(block_train))+1+(1:t_sim),] %>%
    mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]*col_norm,.))
  
  return(
    bind_rows( block_sim_out %>% mutate(type = 'sim'), #time=lib[2]+ 1:t_sim,
               block_sim[1+(1:t_sim),] %>% mutate(type = 'true') ) #time=lib[2]+ 1:t_sim,
  )
} # sim_EDM_simplex_diff




sim_EDM_smap_diff <- function(block_train,
                              block_sim,
                              t_sim,
                              tp = 1,
                              theta = 0,
                              num_neighbors="e+1",
                              lib_train = c(1,NROW(block_train)),
                              pred_sim = c(1,NROW(block_sim)),
                              sim_col = NULL, # could set to which columns are NAs
                              normalize = TRUE,
                              keep_trajectories = TRUE,
                              first_column_time = TRUE,
                              predictor_col = 1:( NCOL(block_train) - first_column_time ),
                              bootstrap = FALSE,
                              ...){
  
  # NOTES:
  # - sim_col are the columns (names or indices) to iteratively predict with s-maps.
  #    > the predictions in this version of the code are done on the first-differences.
  # - The set of points for testing L_test must be at least t_sim points long. This ensures that
  # when fix_col is specified there is information for the entire simulation.
  # - bootstrap specified as either a logical or an integer. If TRUE, will be set to floor(NROW(block_train)/2).
  # - We are predicting differences in the sim_col variables, but not using differences in the embeddings.
  
  library('dplyr')
  
  if(!identical(names(block_sim),names(block_train))){
    warning("FATAL ERROR: columns of block_train and block_sim do not match")
    return()
  }
  
  if(NROW(block_sim) < t_sim){
    warning("FATAL ERROR: t_sim set longer than block_sim")
    return()
  }
  
  if(bootstrap){
    I_lib_train <- apply(rbind(lib_train),1,function(v) seq(v[1],v[2]-tp)) %>% unlist()
    lib_train <- do.call(rbind,map(I_lib_train,~c(.,.+tp)))
    
    if(identical(bootstrap,TRUE)){
      bootstrap <- floor(NROW(lib_train)/2)
    }
  }
  
  if(!is.character(sim_col)){
    sim_col <- names(block_sim)[sim_col]
  }
  
  # add time column if not specified
  if(!first_column_time){
    block_train <- block_train %>%
      mutate(time = 1:n()) %>%
      select(time,everything())
    block_sim <- block_sim %>%
      mutate(time = NROW(block_train) + (1:n()) ) %>%
      select(time,everything())
  }
  
  n_var <- NCOL(block_train) - 1
  
  # Normalize and take differences
  v_norms <- block_train %>%
    summarise_all(list(~sd(.,na.rm=T)))
  
  block.t_norm <- bind_rows( block_train %>%
                               mutate_at( set_names(sim_col,sim_col), list(delta = ~ c(NA,diff(.)))  ) %>%
                               mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]/col_norm,.)) ,
                             block_sim %>%
                               mutate_at( set_names(sim_col,sim_col), list(delta = ~ c(NA,diff(.)))  ) %>%
                               mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]/col_norm,.))
  )
  
  #### ITERATE THROUGH TIME
  t_pred <- 1
  while(t_pred < t_sim + 1){
    t_pred <- t_pred+1
    delY <- as.data.frame(array(dim=c(1,n_var),dimnames=list(NULL,names(block_train)[-1])))
    
    if(bootstrap){
      lib_temp <- lib_train[sample(1:NROW(lib_train),size = bootstrap,replace=F),]
    }else{
      lib_temp <- lib_train
    }
    
    #### EVALUATE SIMPLEX
    for(i.col in sim_col){
      
      out.temp <- block_lnlp(block.t_norm,
                             method='s-map',
                             lib=lib_temp,
                             pred=c((NROW(block_train))+(t_pred-1),(NROW(block_train))+t_pred),
                             columns=predictor_col,
                             theta=theta,
                             target_column = paste0(i.col,"_delta"),
                             num_neighbors = num_neighbors,
                             stats_only = FALSE,
                             silent = TRUE,
                             first_column_time = TRUE,
                             ...)
      
      delY[i.col] <- out.temp$model_output[[1]]$pred[1]
      
      
      #### REPLACE IN BLOCK BY ADDING TO x_i.col(t_pred-1)
      block.t_norm[NROW(block_train)+t_pred,paste0(i.col,"_delta")] <- delY[i.col]
      block.t_norm[NROW(block_train)+t_pred,i.col] <- block.t_norm[NROW(block_train)+t_pred-1,i.col] + delY[i.col] 
      
    } # for(i.col)
  }
  
  # Undo normalize
  block_sim_out <- block.t_norm[(NROW(block_train))+1+(1:t_sim),] %>%
    mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]*col_norm,.))
  
  return(
    bind_rows( block_sim_out %>% mutate(type = 'sim'), #time=lib[2]+ 1:t_sim,
               block_sim[1+(1:t_sim),] %>% mutate(type = 'true') ) #time=lib[2]+ 1:t_sim,
  )
} # sim_EDM_simplex_diff




sim_EDM_smap <- function(block_train,
                         block_sim,
                         t_sim,
                         tp = 1,
                         lib_train = c(1,NROW(block_train)),
                         pred_sim = c(1,NROW(block_sim)),
                         sim_col = NULL, # could set to which columns are NAs
                         normalize = TRUE,
                         theta = 1,
                         keep_trajectories = TRUE,
                         first_column_time = TRUE,
                         predictor_col = 1:( NCOL(block_train) - first_column_time ),
                         bootstrap = FALSE,
                         ...){
  
  # NOTES:
  # - sim_col are the columns (names or indices) to iteratively predict with s-maps.
  #    > the predictions in this version of the code are done on the first-differences.
  # - The set of points for testing L_test must be at least t_sim points long. This ensures that
  # when fix_col is specified there is information for the entire simulation.
  # - bootstrap specified as either a logical or an integer. If TRUE, will be set to floor(NROW(block_train)/2).
  
  library('dplyr')
  
  if(!identical(names(block_sim),names(block_train))){
    warning("FATAL ERROR: columns of block_train and block_sim do not match")
    return()
  }
  
  if(NROW(block_sim) < t_sim){
    warning("FATAL ERROR: t_sim set longer than block_sim")
    return()
  }
  
  if(bootstrap){
    I_lib_train <- apply(rbind(lib_train),1,function(v) seq(v[1],v[2]-tp)) %>% unlist()
    lib_train <- do.call(rbind,map(I_lib_train,~c(.,.+tp)))
    
    if(identical(bootstrap,TRUE)){
      bootstrap <- floor(NROW(lib_train)/2)
    }
  }
  
  if(!is.character(sim_col)){
    sim_col <- names(block_sim)[sim_col]
  }
  
  # add time column if not specified
  if(!first_column_time){
    block_train <- block_train %>%
      mutate(time = 1:n()) %>%
      select(time,everything())
    block_sim <- block_sim %>%
      mutate(time = NROW(block_train) + (1:n()) ) %>%
      select(time,everything())
  }
  
  
  
  # #### Normalize
  # if(normalize){
  #   block <- block %>%
  #     mutate_at(vars(-1),funs(normalize_sqrt))
  # }
  
  n_var <- NCOL(block_train) - 1
  
  block.t <- bind_rows(block_train,
                       block_sim)
  
  
  v_norms <- block_train %>%
    summarise_all(list(~sd(.,na.rm=T)))
  
  block.t_norm <- bind_rows( block_train %>%
                               # mutate_at( set_names(sim_col,sim_col), list(delta = ~ c(NA,diff(.)))  ) %>%
                               mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]/col_norm,.)) ,
                             block_sim %>%
                               # mutate_at( set_names(sim_col,sim_col), list(delta = ~ c(NA,diff(.)))  ) %>%
                               mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]/col_norm,.))
  )
  
  # block.t_norm <- block.t %>%
    # mutate_at(-1,list( ~./sd(block_train$.,na.rm=TRUE)))
  
  # print(names(block.t_norm))
  
  #### fit S-map thetas to prediction of each simulation column
  if(length(theta)>1){
    theta.star <- rep(NA,n_var)
    for(i.col in sim_col){
      
      out.temp <- do.call(bind_rows,lapply(theta,function(theta.i){
        block_lnlp(block.t_norm,lib=lib_train,pred=lib_train,
                   method = 's-map',
                   theta = theta.i,num_neighbors = 0,
                   columns=predictor_col,
                   target_column = i.col,
                   first_column_time = TRUE,
                   ...)
      }))
      
      theta.star[i.col] <- theta[which.max(out.temp$rho)]
    }}else{
      theta.star=rep(theta,n_var)
    }
  ####
  
  
  t_pred <- 1
  
  #### ITERATE THROUGH TIME
  while(t_pred < t_sim + 1){
    t_pred <- t_pred+1
    Y <- as.data.frame(array(dim=c(1,n_var),dimnames=list(NULL,names(block_train)[-1])))
    
    if(bootstrap){
      lib_temp <- lib_train[sample(1:NROW(lib_train),size = bootstrap,replace=F),]
    }else{
      lib_temp <- lib_train
    }
    
    #### EVALUATE SMAP
    for(i.col in sim_col){
      
      out.temp <- block_lnlp(block.t_norm,
                             method='s-map',
                             lib=lib_temp,
                             pred=c((NROW(block_train))+(t_pred-1),(NROW(block_train))+t_pred),
                             columns=predictor_col,
                             target_column = i.col,
                             theta=theta.star[i.col],
                             num_neighbors = 0,
                             stats_only = FALSE,
                             silent = TRUE,
                             first_column_time = TRUE,
                             ...)
      
      Y[i.col] <- out.temp$model_output[[1]]$pred[1]
      
      
      #### REPLACE IN BLOCK
      block.t_norm[NROW(block_train)+t_pred,i.col] <- Y[i.col] 
      
    } # for(i.col)
    
    
  }
  
  # Undo normalize
  # block_sim_out <- block.t_norm[(NROW(block_train))+1+(1:t_sim),] %>%
  block_sim_out <- block.t_norm[(NROW(block_train))+1+(1:t_sim),] %>%
    mutate(!!!imap(v_norms[-1],function(col_norm, col_name, data) data[[col_name]]*col_norm,.))
  
  return(
    bind_rows( block_sim_out %>% mutate(type = 'sim'), #time=lib[2]+ 1:t_sim,
               block_sim[1:t_sim,] %>% mutate(type = 'true') ) #time=lib[2]+ 1:t_sim,
  )
} # sim_EDM_smap

