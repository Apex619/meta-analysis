#### Custom functions
#  A set of custom functions named: `make_VCV_matrix`, `I2()`, `R2()`, `get_est()`, `get_pred()`, `cont_gen()`, `compare_levels()`



#' @title Variance-covariance and correlation matrix function basing on shared level ID
#' @description Function for generating simple variance-covariance and correlation matrices 
#' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#' @param V Name of the variable (vector) containing effect size variances variances
#' @param cluster Name of the variable (vector) indicating which effects belong to the same cluster. Same value of 'cluster' are assumed to be nonindependent (correlated).
#' @param obs Name of the variable (vector) containing individual IDs for each value in the V (Vector of variances). If this parameter is missing, label will be labelled with consecutive integers starting from 1.
#' @param rho Known or assumed correlation value among effect sizes sharing same 'cluster' value. Default value is 0.5.
#' @param type Optional logical parameter indicating whether a full variance-covariance matrix (default or "vcv") is needed or a correlation matrix ("cor") for the non-independent blocks of variance values.
#' @export
#Value: Labelled full variance-covariance or correlation matrice of the size and labels matching initial dataframe will be returned 

make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  
  if (missing(data)) 
    stop("Must specify dataframe via 'data' argument.")
  if (missing(V)) 
    stop("Must specify name of the variance variable via 'V' argument.")
  if (missing(cluster)) 
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  if (missing(obs)) 
    obs <- 1:length(V)   
  if (missing(type)) 
    type <- "vcv" 
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  # find start and end coordinates for the subsets
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  return(new_matrix)
}



# General modeling functions 
# Functions for I2

#' Title Function to obtain total and separate I2 from multilevel-meta-analytic model
#'
#' @param model 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if(method == "Wolfgang"){
    W <- solve(model$V) 
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each  <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
    
    # or my way
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k-1) / (sum(1/model$vi)^2 - sum((1/model$vi)^2)) 
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) #s^2_t = total variance
    I2_each  <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
  }
  return(I2s)
}

I2_Hamza <- function(model){
  W <- diag(1/model$vi)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total <- 100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  return(I2_total)
}


# test <- dataset$fit4.1[[3]]
# I2(test, method = "Wolfgang")
# I2(test, method = "Shinichi")


#' Title: R2 based on Nakagawa & Schielzeth 2013
#'
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
R2 <- function(model){
  warning("Conditional R2 is not meaningful and the same as marginal R2\n")
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  R2
  #Rm <- round(100*R2m, 3)
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}


#' Title: the function to get estimates from rma objects (metafor)
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_est <- function (model, mod = " ") {
  
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub 
  
  table <- tibble(name = name, estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)
}


#' Title: the function to get prediction intervals (crediblity intervals) from rma objects (metafor)
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_pred <- function (model, mod = " ") {
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  len <- length(name)
  
  if(len != 1){
    newdata <- matrix(NA, ncol = len, nrow = len)
    for(i in 1:len) {
      # getting the position of unique case from X (design matrix)
      pos <- which(model$X[,i] == 1)[[1]]
      newdata[, i] <- model$X[pos,]
    }
    pred <- predict.rma(model, newmods = newdata)
  }
  else {
    pred <- predict.rma(model)
  }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub 
  
  table <- tibble(name = name, lowerPR = lowerPR, upperPR = upperPR)
}

#Here are links for how to do confidence regions for rma.mv regression lines
#https://www.rdocumentation.org/packages/metafor/versions/1.9-9/topics/predict.rma
#https://stackoverflow.com/questions/50804464/out-of-sample-prediction-for-rma-object-in-metafor


#' Title: Contrast name geneator
#'
#' @param name: a vector of character strings
cont_gen <- function (name) {
  combination <- combn(name,2)
  name_dat <- t(combination)
  names <- paste(name_dat[ ,1], name_dat[, 2], sep = "-")
  return(names)
}


#' Title: Contrast models geneator
#'
#' @param model: a model to be updated (ideally, in itercept model for a single factor moderator)
#' @param factor: a factor to be repeatedly releveled
#' 
compare_levels <- function(model, factor) {
  sym <- substitute(factor)
  results <- list() # empty list to store results
  for (i in levels(factor)) {
    change <- eval(bquote(~ relevel(.(sym), ref=.(i))))
    new_model <- update(model, change)
    results[[i]] <- tibble(names = dimnames((new_model)$b)[[1]], 
                           estimate = summary(new_model)$b[,1], 
                           se = summary(new_model)$se,
                           #zval = summary(new_model)$zval,
                           #pval = summary(new_model)$pval),
                           ci.lb = summary(new_model)$ci.lb,
                           ci.ub = summary(new_model)$ci.ub
    )
  }
  return(results)
}

#levels(data_OF$Trait)
#cont_gen(levels(data_OF$Trait))
#data_OF$Trait <- ordered(data_OF$Trait, levels = c("Body_Weight", "Adiposity", "Triglycerides", "Glucose_TT", "Glucose_FBG", "Insulin_FI", "Insulin_TT", "Leptin"))



