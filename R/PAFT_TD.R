if(!require("R.utils")){
  install.packages("R.utils")
}
R.utils::use("nloptr")
R.utils::use("survival")
R.utils::use("dplyr")
R.utils::use("nnls")
#------------ Function ------------#

#' Parametric AFT model with time dependent covariates
#'
#' @param dat Input data
#' @param X List of vectors indicating the parameters of any covariate that the user needs to introduce in this model
#' @param initi The initial values of betas
#' @param beta List of vectors indicating the effect of the corresponding covariate
#' @param gab The initial support(=SE)
#' @param dist Distribution assumed for log(T_0)
#' @param tol Allowable margin of error
#' @param maxiter The maximum number of iteration
#' @param algorithm Algorithm used for iteration
#'
#' @return $Estimated Betas $itration $status $message
#' @export
#'
#' @examples PAFT_TD(dat=dat1,X=c("Z1","X"),dist="log_normal") # True data: log-normal with b0=1, b1=1, b2=1)
#' @examples PAFT_TD(dat=dat4,X=c("Z1","X"),dist="log_normal") # True data: log-t(3) with b0=1, b1=1, b2=1
#'
PAFT_TD <- function(formula,data,dist,initial=F,beta0=NA,scale0=NA,tol=1.0e-5,maxiter=2000,algorithm="NLOPT_LN_COBYLA"){
  
  ## Basic options
  varnames <- all.vars(formula)
  last_data <- as.data.frame(data%>%dplyr::group_by(data[,1])%>%dplyr::arrange(varnames[2])%>%dplyr::slice(n()))[,c(names(data)[1],varnames[2:length(varnames)])]
  
  if(length(varnames)==3){
    # when no covariates
    covr_names <- NULL
    names(last_data) <- c("ID","Time","delta")
  }else{
    covr_names <- varnames[c(4:length(varnames))]
    names(last_data) <- c("ID","Time","delta",covr_names)
  }
  
  ##  Check the types of covariate(s)
  if(nrow(last_data)==nrow(unique(as.data.frame(data[,c(names(data)[1],covr_names)])))){
    
    # 1. Time-independent covariate(s) type : survreg modeling
    
    if(length(covr_names)==0){ # when no covariate
      model_result <- survival::survreg(formula=Surv(Time,delta)~1,data=last_data,dist="lognormal")
      Betas_results <- round(as.numeric(c(model_result$coefficients,model_result$scale)),4)
      names(Betas_results) <- c("Beta0","Scale")
    }else{
      formula_indep <- stats::as.formula(paste0("Surv(Time,delta)~",paste(covr_names, collapse="+")))
      model_result <- survival::survreg(formula_indep,data=last_data,dist="lognormal")
      Betas_results <- round(as.numeric(c(model_result$coefficients,model_result$scale)),4)
      names(Betas_results) <- c(paste("Beta",c(0:length(covr_names)),sep=""),"Scale")
    }
    final <- list(EstBetas=Betas_results,Model="Time-independent AFT")
    
  }else{
    
    # 2. Time-dependent covariate(s) type : using optimization method
    
    
    # initial values: beta0 and scale0
    if(initial==T){
      beta <- beta0; scale <- scale0
    }else{
      formula_indep <- stats::as.formula(paste0("Surv(Time,delta)~",paste(covr_names, collapse="+")))
      PAFT <- survival::survreg(formula_indep,data=last_data,dist="lognormal")
      beta <- PAFT$coefficients;scale <- PAFT$scale
    }
    
    
    # dist1: log-normal distribution
    if(dist=="lognormal"){
      
      # Optimization function
      PAFT_opt <- function(theta){
        
        # 1. Calculation of the observation time
        dat <- data[,c(names(data)[1],varnames)]
        names(dat) <- c("ID",varnames)
        obs_time <- last_data[,c("ID","Time")];names(obs_time) <- c("ID","obs_time")
        dat <- merge(dat,obs_time,by="ID",all.x=T)
        dat <- dat[order(dat$ID,dat$start),] # ordering
        
        # 2. Likelihood function of the PAFT with time-dependent covariates
        
        ## Term1 : Calculation of the Capital psi 
        AA <- 0
        for(i in 1:length(covr_names)){
          AA <- -theta[i+1]*dat[,covr_names[i]]+AA
        }
        dat$psi <- (dat$stop-dat$start)*exp(AA)
        final <- as.data.frame(dat%>%group_by(ID)%>%summarise(psi=sum(psi)))
        
        ## Term2 : Calculation of psi : exp(-beta*z(yi))
        BB <- 0
        for(i in 1:length(covr_names)){
          BB <- -theta[i+1]*last_data[,covr_names[i]]+BB
        }
        final$psi_d <- exp(BB)
        
        ## Joint delta
        delta <- last_data[,c("ID","delta")]
        final <- merge(final,delta,by="ID",all.x=T)
        final <- final[order(final$ID),] # ordering
        
        ## likelihood function ## 
        scale <- theta[length(theta)]
        logPsi <- log(final$psi)
        dnorm_2 <- dnorm((logPsi-theta[1])/scale)*(final$psi_d/(scale*final$psi))
        S <- final$delta*dnorm_2 + (1-final$delta)*(1 - pnorm((logPsi-theta[1])/scale))
        S <- ifelse(S<=1.0e-20,1.0e-20,S)
        
        lik <- sum(log(S))
        lik_inv <- -lik
        return(lik_inv)
      }
      
      b <- as.numeric(c(beta,scale))
      
      # bounded values
      bound <- c(max(abs(b))*10,rep(max(abs(b))*10,length(covr_names)),max(abs(b))*10)
      
      # ----  nloptr function ----#
      res2 <- nloptr(x0=b,eval_f=PAFT_opt,lb = -bound,ub = bound,
                     opts = list("algorithm"=algorithm,"xtol_rel"=tol,maxeval=maxiter))
      
      # Results
      Betas_results <- res2$solution[1:(1+length(covr_names)+1)] # Estimated betas
      names(Betas_results) <- c(paste("Beta",c(0:length(covr_names)),sep=""),"Scale")
      itr <- res2$iterations # number of iterations of the algorithm
      status <- res2$status # integer value with the status of the optimization
      message <- res2$message # more informative message with the status of the optimization
      
      final <- list(EstBetas=Betas_results,itr=itr,status=status,message=message)
      
    }else{
      
      final <- NULL
      
    }
  }
  return(final)
}



