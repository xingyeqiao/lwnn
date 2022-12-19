library(modelr)
library(dplyr)
library(tidyr)
library(FNN)      # fast nearest neighbor search
library(purrr)    # required by map_dfr. Efficiently apply a function to each row of a dataframe

wnn <- function (train, test, cl, param = NULL, alpha = 0) 
{
  train <- as.matrix(train)
  if (is.null(dim(test))) 
    dim(test) <- c(1, length(test))
  test <- as.matrix(test)
  if (any(is.na(train)) || any(is.na(test)) || any(is.na(cl))) 
    stop("no missing values are allowed")
  p <- ncol(train)
  ntr <- nrow(train)
  if (length(cl) != ntr) 
    stop("'train' and 'class' have different lengths")
  nte <- nrow(test)
  if (ncol(test) != p) 
    stop("dims of 'test' and 'train' differ")
  
  cl = as.factor(as.character(cl))
  cl_label = levels(cl)
  cl_num = as.numeric(cl) - 1
  
  
  # train = scale(train)
  # test = scale(test)

  ### slice kNN (Gabet et al.)
  ## density estimate
  if(p < 6){
    den_te = ks::kde(x = train, eval.points = test)$estimate
  }else{
    KDEval <- function(xx){  # use the estimated joint density
     dx <- density(xx[1:ntr])
     denFun = approxfun(dx$x,dx$y,)
     denFun(xx[(ntr+1):(ntr+nte)])
    }
    den_te <- apply(apply(rbind(train,test),2,KDEval),1,prod)
    den_te[is.na(den_te)] <- 1e-50
    den_te[den_te == 0] <- 1e-50
  }
  
  if((!is.null(param$sliceknn)) & (length(param$sliceknn) > 1)){
    k0_slice = param$sliceknn
    k_sliceknn = matrix(0, nte, length(k0_slice))
    # k_sliceknn[den_te >= ntr^(-alpha/(2+alpha+p)), ] = matrix(1,numeric(nte),1) %*% floor(k0_slice)
    k_sliceknn[den_te >= ntr^(-alpha/(2+alpha+p)), ] = matrix(1,sum(den_te >= ntr^(-alpha/(2+alpha+p))),1) %*% floor(k0_slice)
    for(j in 0:min(ceiling(log(ntr^(-alpha/(2+alpha+p))/min(den_te))/log(2)), 100)){
      index = (den_te >= ntr^(-alpha/(2+alpha+p)) / 2^(j+1)) & (den_te < ntr^(-alpha/(2+alpha+p)) / 2^(j))
      k_sliceknn[index, ] = matrix(1,sum(index),1) %*% pmax(floor(k0_slice * 2^(-2*(j+1)/(2+p))), 1)
    }
  }else{
    if(is.null(param$sliceknn)){k0_slice = ntr^(2/(2+alpha+p))*log(ntr)}else{k0_slice= param$sliceknn}
    k_sliceknn = numeric(nte)
    k_sliceknn[den_te >= ntr^(-alpha/(2+alpha+p))] = floor(k0_slice)
    for(j in 0:ceiling(log(ntr^(-alpha/(2+alpha+p))/min(den_te))/log(2))){
      k_sliceknn[(den_te >= ntr^(-alpha/(2+alpha+p)) / 2^(j+1)) & (den_te < ntr^(-alpha/(2+alpha+p)) / 2^(j))] = pmax(floor(k0_slice * 2^(-2*(j+1)/(2+p))), 1)
    }
  }
  k_sliceknn = pmin(k_sliceknn, ntr-1)
  k_sliceknn = as.data.frame(k_sliceknn)
  names(k_sliceknn) <- paste('sliceknn_',1:length(k_sliceknn),sep='')
  
  ### OWNN weights
  ownn_w_fcn <- function(k_ownn_t){
    ii <- 1:k_ownn_t
    alpha_ownn <- ii^(1 + 2/p) - (ii - 1)^(1 + 2/p)
    return((1 + p/2 - p/(2 * k_ownn_t^(2/p)) * alpha_ownn)/k_ownn_t)
  }
  
  if((!is.null(param$ownn)) & (length(param$ownn) > 1)){
    k_ownn = param$ownn
    ownn_weights = k_ownn %>% map(ownn_w_fcn)
  }else{
    if(is.null(param$ownn)){
      k_ownn <- min(floor(ntr^(4/(4+p))), ntr-1)  
    }else{
      k_ownn = param$ownn
    }
    ownn_weights = list()
    ownn_weights[[1]] = ownn_w_fcn(k_ownn)
  }
  names(ownn_weights) <- paste('ownn_',1:length(ownn_weights),sep='')
    
  ### knn
  if(!is.null(param$knn)){
    k_knn = param$knn
  }else{
    k_knn = floor(ntr^(2/(2+d)))
  }
  k_knn = pmin(k_knn, ntr-1)  
  names(k_knn) <- paste('knn_',1:length(k_knn),sep='')
  
  
  ### LWNN
  if(!is.null(param$lwnn)){
    k0_lwnn = param$lwnn
  }else{
    k0_lwnn = ceiling(ntr^(2/(2+d)))           ## theorem 1
  }
  kc_lwnn = pmin(k0_lwnn, ntr-1)
  names(kc_lwnn) <- paste('lwnn_',1:length(kc_lwnn),sep='')

  
  ### kstarnn
  if(!is.null(param$kstarnn)){
    LC_kstarnn = param$kstarnn
  }else{
    LC_kstarnn = 1
  }
  names(LC_kstarnn) <- paste('kstarnn_',1:length(LC_kstarnn),sep='')
  
  kstarnn_w_fcn <- function(beta.unsorted)
  {
    beta=sort(beta.unsorted)
    n=length(beta.unsorted)
    lambda=beta[1]+1
    k=0
    sumbeta=0
    sumbeta2=0
    while(lambda>beta[k+1]&k<=(n-1))
    {
      k=k+1
      sumbeta=sumbeta+beta[k]
      sumbeta2=sumbeta2+beta[k]^2
      lambda=(1/k)*(sumbeta+sqrt(k+sumbeta^2-k*sumbeta2))
    }
    alpha=(lambda-beta.unsorted)*as.integer(beta.unsorted<lambda)
    alpha/sum(alpha)
  }
  
  num_2_label = function(x){
    as.numeric(cl_label[x + 1])
  }
  
  # get the nearest neighbor index and distance
  nn = get.knnx(train, test, k = ntr, algorithm=c("kd_tree"))
  
  predict.nn <- function(i){ # return the predictions from knn, sliceNN, ownn, kstarnn, and lwnn
    data.frame(
      LC_kstarnn %>% map_dfc(function(x){
          kstarnn_weights = kstarnn_w_fcn(nn$nn.dist[i,] * x)
          kstarnn = as.numeric(sum(cl_num[nn$nn.index[i,]] * kstarnn_weights)    >= 0.5)
          kstarnn = num_2_label(kstarnn)
          return(kstarnn)}),
      kc_lwnn %>% map_dfc(function(x){
          nndist = nn$nn.dist[i,] + (1e-8)*(1:length(nn$nn.dist[i,]))
          lwnn_weights = nndist[x + 1] - nndist[1:x]
          lwnn_weights = lwnn_weights/sum(lwnn_weights)
          lwnn = as.numeric(sum(cl_num[nn$nn.index[i,1:x]] * lwnn_weights)    >= 0.5)
          lwnn = num_2_label(lwnn)
          return(lwnn)}),
      k_knn %>% map_dfc(function(x){num_2_label(getmode(cl_num[nn$nn.index[i,1:x]]))}),
      k_sliceknn %>% map_dfc(function(x){num_2_label(getmode(cl_num[nn$nn.index[i,1:(x[i])]]))}),
      ownn_weights %>% map_dfc(function(x){num_2_label(as.numeric(sum(cl_num[nn$nn.index[i,1:length(x)]] * x) >= 0.5))})
    )
  }
  pred <- (1:nte) %>% map_dfr(predict.nn)
  return(pred)
}


getmode <- function(v) {  # get the mode of a sample (i.e. cast majority vote)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

get0margin<-function(x,delta)
{
  if(x<0.5-delta){x/(1-2*delta)}
  else if(x>0.5+delta){(x-2*delta)/(1-2*delta)}
  else {0.5}
}

get0marginpar<-function(x,delta)
{
  n=length(x)
  result=rep(NA,n)
  for(i in 1:n){result[i]=get0margin(x[i],delta)}
  result
}


adjust_n_opt <- function(param,n,d,fold){
  
  if ((!exists("distn"))|!(distn %in% c("Uniform", "Normal", "Cauchy", "Pareto", "Truncated_Normal", "Laplace"))){
    return(param)  # no adjustment
  }
  else{
      if (distn %in% c("Uniform", "Normal", "Cauchy", "Pareto", "Truncated_Normal", "Laplace")){
        param$sliceknn = param$sliceknn * (fold/(fold-1))^(2/(2+d)) * (log(n)/log(n*(fold-1)/fold))
        param$ownn = param$ownn * (fold/(fold-1))^(4/(4+d))
      }
      
      if((distn == "Uniform")|(distn == "Truncated_Normal")){
        
        param$lwnn = param$lwnn * (fold/(fold-1))^(2/(2+d))
        param$knn = param$knn * (fold/(fold-1))^(2/(2+d))
      
      }else if((distn == "Normal")|(distn == "Laplace")){
        
        param$lwnn = param$lwnn * (fold/(fold-1))^(2/(2+d)) * (log(n)/log(n*(fold-1)/fold))^(2*d/(d+2))
        param$knn = param$knn * (fold/(fold-1))^(2/(3+d))
    
      }else if(distn == "Cauchy"){
        
        param$lwnn = param$lwnn * (fold/(fold-1))^(2/(2+0.5*d))
        param$knn = param$knn * (fold/(fold-1))^(2/(4+d))
      }else if(distn == "Pareto"){
        
        param$lwnn = param$lwnn * (fold/(fold-1))^(2/(2+2*d/3))
        param$knn = param$knn * (fold/(fold-1))^(2/(3.5+d))
      }
      return(param)
  }
  
}

cv.wnn <- function (train, cl, fold = 10, adjust_n = FALSE, param = list(lwnn = c(1),
                                            knn = c(1),
                                            kstarnn = c(1), 
                                            ownn = c(1), 
                                            sliceknn = c(1))){
  # ptm <- proc.time()
  
  cv_modelr <- data.frame(cl,train) %>% 
    crossv_kfold(k = fold) %>% # split to training set and validation set
    mutate(prediction_y = map2(train, test, function(.x,.y) wnn(data.frame(.x)[,-1], data.frame(.y)[,-1], data.frame(.x)[,1], param = param))) %>% 
    mutate(misclass_error = map2(prediction_y, test, function(.x,.y) colMeans(as.matrix(.x) != data.frame(.y)[,1])))
  
  cv_misclass <- cv_modelr$misclass_error %>% 
    as.data.frame(.) %>% 
    apply(., 1, mean)
  cv_misclass_sd <- cv_modelr$misclass_error %>% 
    as.data.frame(.) %>% 
    apply(., 1, sd)/sqrt(fold)
  
  param2 = param
  for(x in names(param2)){
    param2[[x]] <- data.frame(param2[[x]])
    names(param2[[x]]) <- c('par_val')
    rownames(param2[[x]]) <- paste(x, 1:nrow(param2[[x]]), sep='_')
  }
  
  par_df = do.call(rbind,param2)
  par_df$method2 = rownames(par_df)
  par_df <- par_df %>% separate(method2, c("method","method_par"),sep="\\.") %>% select(par_val,method_par)
  
  cv_error_df = data.frame(t(rbind(cv_misclass, cv_misclass_sd)))
  cv_error_df$method_par = row.names(cv_error_df)
  
  optimal_par <- cv_error_df %>% left_join(par_df, by = "method_par") %>% 
    separate(method_par,c('method','par'),sep = '_') %>% 
    group_by(method) %>% 
    summarise(opt_par_1se = ifelse(method[1] == 'kstarnn',
                               min(par_val[which(cv_misclass <= min(cv_misclass) +  cv_misclass_sd[which.min(cv_misclass)])]),
                               max(par_val[which(cv_misclass <= min(cv_misclass) +  cv_misclass_sd[which.min(cv_misclass)])])), 
              opt_par_min = par_val[which.min(cv_misclass)], 
              .groups = 'drop')
  

  optimal_par_1se <- split(as.numeric(optimal_par$opt_par_1se), optimal_par$method)
  optimal_par_min <- split(as.numeric(optimal_par$opt_par_min), optimal_par$method)
  
  if(adjust_n == TRUE){
    optimal_par_1se = adjust_n_opt(optimal_par_1se,length(cl),ncol(train),fold)
    optimal_par_min = adjust_n_opt(optimal_par_min,length(cl),ncol(train),fold)
  }
  
  
  # print(paste(round((proc.time() - ptm)[3], 2), " seconds elapsed for cv.wnn", sep=''))
  
  return(list(error = cv_error_df, param_1se = optimal_par_1se, param_min = optimal_par_min, param = param))
}



runRealData <- function(nrep, fold_in, fold_out){
  n = nrow(features)
  # n_tr_cv = n/fold_out*(fold_out-1)/fold_in*(fold_in-1)
  n_tr_cv = n/fold_out*(fold_out-1)
  param = list(lwnn = c(floor(seq(1, n_tr_cv/2-1, length.out = 20))),
               knn = c(floor(seq(1, n_tr_cv/2-1, length.out = 20))),
               sliceknn = c(floor(seq(1, n_tr_cv/2-1, length.out = 20))),
               ownn = c(floor(seq(1, n_tr_cv/2-1, length.out = 20))),
               kstarnn = c(exp(seq(log(0.001), log(10), length.out = 20)))
  )
  
  
  cv_modelr <- map(1:nrep,function(x){print(paste('Replication #', x, sep=''));
    data.frame(cl,features) %>% 
      crossv_kfold(k = fold_out) %>%
      mutate(prediction_y = map2(train, test, function(.x, .y){
        temp  = cv.wnn(data.frame(.x)[,-1], data.frame(.x)[,1], param = param, fold = fold_in)
        temp_opt_param = Map(c, temp$param_1se, temp$param_min)
        wnn(data.frame(.x)[,-1], data.frame(.y)[,-1], data.frame(.x)[,1], param = temp_opt_param)
      })) %>% 
      mutate(misclass_error = map2(prediction_y, test, function(.x,.y) colMeans(as.matrix(.x) != data.frame(.y)[,1])))})
  
  cv_modelr <- do.call(rbind,cv_modelr)
  
  cv_misclass <- cv_modelr$misclass_error %>% as.data.frame(.) %>% apply(., 1, mean)
  
  cv_misclass_sd <- cv_modelr$misclass_error %>% as.data.frame(.) %>% apply(., 1, sd)/sqrt(fold_out*nrep)
  
  return(list(name = name_data,
              dim = dim(features),
              report = rbind(cv_misclass, cv_misclass_sd)[,c(1,3,5,7,9,2,4,6,8,10)]))
}
