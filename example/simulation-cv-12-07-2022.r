## same setting as in the paper, with optimal k adjusted and fine grid.

setting = "I"
filename = paste("LWNN_12_07_2022_simulation_cv_",setting,".Rdata",sep = '')
filename2 = paste("result_12_07_2022_",setting,".csv",sep = '')

source("../R/wnn-r.R")

library(EnvStats)  # required by rpareto
library(purrr)    # required by map_dfr. Efficiently apply a function to each row of a dataframe
library(rmutil)   # required by rlaplace
library(truncnorm) # required by rtruncnorm


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

fold_in = 10

############
tmax = 100 # number of repetitions
rho_simulation = 0.001 # the parameter rho in the regression function in simulation examples

dist.list = c("Uniform","Truncated_Normal","Normal","Laplace","Pareto","Cauchy")

if(setting == "I"){
  dist.list = dist.list[1:2]
}else if(setting == "II"){
  dist.list = dist.list[3:4]
}else if(setting == "III"){
  dist.list = dist.list[5:6]
}

d.list = c(4,8)
n.list = seq(50,500,75)
nt = 5000  # test data sample size
alpha = 0  # margin index is set to 0
set.seed(100)

errorM = NULL

for(distn in dist.list){
  for(d in d.list){
    for(n in n.list){
      for(t in 1:tmax){
        print(c(distn,d,n,t))
        
        if(distn == "Uniform"){
          lowerlimit=-pi/2
          upperlimit=pi/2
          X = matrix(runif(n*d,lowerlimit,upperlimit), nrow = n)
          Xt = matrix(runif(nt*d,lowerlimit,upperlimit), nrow = nt) 
          
          dens <- function(x){prod(dunif(x,lowerlimit,upperlimit))}
          
          k0_lwnn = ceiling(n^(2/(2+d)))           ## theorem 1
          k_knn = floor(n^(2/(2+d)))
          
          eta = function(x){  # regression function for each test data point x
            f =get0margin(abs(sin(x[1]+x[2])),rho_simulation)
            rbinom(1,1,f)
          }
          
          Ybayes =get0marginpar((abs(sin(Xt[,1]+Xt[,2]))),rho_simulation) >= 0.5
          
        }else if(distn == "Truncated_Normal"){
          
          lowerlimit=-pi/2
          upperlimit=pi/2
          X=matrix(rtruncnorm(n*d,lowerlimit,upperlimit,0,1),nrow=n)
          Xt=matrix(rtruncnorm(nt*d,lowerlimit,upperlimit,0,1),nrow=nt)  
          
          dens <- function(x){prod(dtruncnorm(x,lowerlimit,upperlimit,0,1))}
          
          k0_lwnn = ceiling(n^(2/(2+d)))  ## theorem 1
          k_knn = floor(n^(2/(2+d)))
          
          eta = function(x){  # regression function for each test data point x
            f =get0margin(abs(sin(x[1]+x[2])),rho_simulation)
            rbinom(1,1,f)
          }
          
          Ybayes =get0marginpar((abs(sin(Xt[,1]+Xt[,2]))),rho_simulation)>= 0.5
          
        }else if(distn == "Normal"){
          
          X = matrix(rnorm(n*d), nrow = n)
          Xt = matrix(rnorm(nt*d), nrow = nt)
          X[,1:2] = X[,1:2] * 2
          Xt[,1:2] = Xt[,1:2] * 2
          
          dens <- function(x){0.25*prod(dnorm((x)/c(2,2,rep(1,d-2))))}
          
          k0_lwnn = ceiling(n^(2/(2+d))*log(n)^(2*d/(d+2)))   ## theorem 2
          k_knn = floor(n^(2/(3+alpha+d)))
          
          eta = function(x){  # regression function for each test data point x
            f =sum(x)
            rbinom(1,1,get0margin(1/(1+exp(-f)),rho_simulation))
          }
          ft=1/(1+exp(-1*apply(Xt,1,sum)))
          Ybayes=get0marginpar(ft,rho_simulation) >= 0.5
          
        }else if(distn == "Laplace"){
          
          X = matrix(rmutil::rlaplace(n*d), nrow = n)
          Xt = matrix(rmutil::rlaplace(nt*d), nrow = nt)
          
          dens <- function(x){prod(rmutil::dlaplace(x))}
          
          k0_lwnn = ceiling(n^(2/(2+d))*log(n)^(2*d/(d+2)))   ## theorem 2
          k_knn = floor(n^(2/(3+alpha+d)))
          
          eta = function(x){  # regression function for each test data point x
            f =sum(x)
            rbinom(1,1,get0margin(1/(1+exp(-f)),rho_simulation))
          }
          ft=1/(1+exp(-1*apply(Xt,1,sum)))
          Ybayes=get0marginpar(ft,rho_simulation)>= 0.5
          
        }else if(distn == "Cauchy"){
          
          X = matrix(rcauchy(n*d), nrow = n)
          Xt = matrix(rcauchy(nt*d), nrow = nt)
          
          dens <- function(x){prod(dcauchy(x))}
          
          k0_lwnn = ceiling(n^(2/(2+0.5*d)))  ## theorem 3/(a)
          k_knn = floor(n^(2/(4+alpha+d)))  ## Gabat them 4.5
          
          eta = function(x){  # regression function for each test data point x
            f =sum(x)
            rbinom(1,1,get0margin(1/(1+exp(-f)),rho_simulation))
          }
          ft=1/(1+exp(-1*apply(Xt,1,sum)))
          Ybayes=get0marginpar(ft,rho_simulation)>= 0.5
          
        }else if(distn == "Pareto"){
          loca_para=2
          scale_para=2
          X = matrix(EnvStats::rpareto(n*d,loca_para,scale_para)-loca_para*scale_para/(scale_para-1), 
                     nrow = n)
          Xt = matrix(EnvStats::rpareto(nt*d,loca_para,scale_para)-loca_para*scale_para/(scale_para-1), 
                      nrow = nt)
          
          dens <- function(x){prod(EnvStats::dpareto(x+loca_para*scale_para/(scale_para-1),
                                                     loca_para,scale_para))}
          
          k0_lwnn = ceiling(n^(2/(2+scale_para*d/(scale_para+1))))  ## theorem 3/(a)
          k_knn = floor(n^(2/(3.5+alpha+d)))  ## Gabat them 4.5
          
          eta = function(x){  # regression function for each test data point x
            f =sum(x)
            rbinom(1,1,get0margin(1/(1+exp(-f)),rho_simulation))
          }
          ft=1/(1+exp(-1*apply(Xt,1,sum)))
          Ybayes=get0marginpar(ft,rho_simulation)>= 0.5
        }
        
        k0_sliceknn = floor(n^(2/(2+alpha+d))*log(n))
        k0_ownn = min(floor(n^(4/(4+d))), n-1)
        
        # generate Y and Yt
        Y = apply(X, 1, eta)
        Yt = apply(Xt, 1, eta)
        
        
        ## set up cross-valiation grid
        
       param = list(lwnn = seq(1, n/2-1, 2),
                    knn = seq(1, n/2-1, 2),
                    sliceknn = seq(1, n/2-1, 2),
                    ownn = seq(1, n/2-1, 2),
                    kstarnn = c(exp(seq(log(0.001), log(1000), length.out = 30))))   

        # cross-valiate on training, then predict testing
        cv_wnn_test <- function(dat){
          temp  = cv.wnn(dat[,-1], dat[,1], adjust_n = TRUE, param = param, fold = fold_in)
          temp_opt_param = Map(c, temp$param_1se, temp$param_min)
          wnn(dat[,-1], Xt, dat[,1], param = temp_opt_param)
        } 

        error = data.frame(Y,X) %>% cv_wnn_test %>% {colMeans(. != Yt)*100}
        error_row = c(distn, n, d, mean(Ybayes != Yt)*100, error)
        names(error_row) <- c("dist","n","d","Bayes",
                          "kstarnn_1se","kstarnn_min",
                          "lwnn_1se","lwnn_min",
                          "knn_1se","knn_min",
                          "sliceknn_1se","sliceknn_min",
                          "ownn_1se","ownn_min")
        
        errorM = rbind(errorM, error_row)
        save.image(filename)
      }
      
    }
  }
}


library(tidyverse)
library(dplyr)
errorM_I <- data.frame(errorM) %>% mutate(across(c(n, Bayes:ownn_min), as.numeric))

risk = errorM_I %>% group_by(dist,d,n) %>% summarise(across(kstarnn_1se:ownn_min, mean))  # average over repetitions
bayes = errorM_I %>% group_by(dist, d) %>% summarise(Bayes = mean(Bayes))
excess_risk = risk %>% left_join(bayes, by = c("dist", "d")) %>%
  mutate(across(kstarnn_1se:ownn_min, ~ .x - Bayes)) %>% 
  select(- Bayes) %>% 
  select(- ends_with("1se")) %>% 
  rename_at(vars(ends_with("_min")), ~str_replace(., "_min", "")) %>% 
  relocate(knn, .before = kstarnn) %>% 
  relocate(lwnn, .after = ownn) %>% 
  mutate(across(kstarnn:lwnn, .fns = ~ .x / knn, .names = "{col}_RR"))

write_csv(excess_risk, filename2)