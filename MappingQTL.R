# R code for mapping resistance QTL
#----------------------------------
# parameters:
#   - data: data.frame for dataset
#   - LG: list containing IBD matrices
#   - type.resp: an integer 1 if the fisrt symptom, 2 for the death
#   - spatial: a boolean: TRUE if spatial effect in models (3) and (4)
#               FALSE otherwise
#   - m: an integer for the number of parents

MappingQTL <- function(data, LG, type.resp = 1, spatial=TRUE, m){

# p the number of position at testing for mapping QTL
#----------------------------------------------------
p <- length(LG)

if (type.resp == 1){
  if (spatial == TRUE){
  form.null <- as.formula(paste("Surv(time=Y_T1S,event=EVENT_T1S,type=\"right\") ~ ",paste(c("1","SPA","PROGENY"), collapse = "+"), sep="" ))
  form.1 <- as.formula(paste("Surv(time=Y_T1S,event=EVENT_T1S,type=\"right\") ~ ",paste(c("1","SPA","PROGENY","(1|iden)"), collapse = "+"), sep="" ))
  }else{
    form.null <- as.formula(paste("Surv(time=Y_T1S,event=EVENT_T1S,type=\"right\") ~ ",paste(c("1","PROGENY"), collapse = "+"), sep="" ))
    form.1 <- as.formula(paste("Surv(time=Y_T1S,event=EVENT_T1S,type=\"right\") ~ ",paste(c("1","PROGENY","(1|iden)"), collapse = "+"), sep="" ))
  }
}else{
  if (spatial == TRUE){
  form.null <- as.formula(paste("Surv(time=EVENT_TD,event=Y_TD,type=\"right\") ~ ",paste(c("1","SPA","PROGENY"), collapse = "+"), sep="" ))
  form.1 <- as.formula(paste("Surv(time=EVENT_TD,event=Y_TD,type=\"right\") ~ ",paste(c("1","SPA","PROGENY","(1|iden)"), collapse = "+"), sep="" ))
  }else{
    form.null <- as.formula(paste("Surv(time=EVENT_TD,event=Y_TD,type=\"right\") ~ ",paste(c("1","PROGENY"), collapse = "+"), sep="" ))
    form.1 <- as.formula(paste("Surv(time=EVENT_TD,event=Y_TD,type=\"right\") ~ ",paste(c("1","PROGENY","(1|iden)"), collapse = "+"), sep="" ))
  }
}

res.coxme <- data.frame()

# The null model defined by model (3)
#------------------------------------
mod.cox.null <- coxph(form.null,data,ties="efron")

# A loop for testing QTLs at p positions with model (4)
#------------------------------------------------------

for (i in 1:p){
    
    M <- LG[[i]]
    colnames(M)[-(1:m)] <- data$iden
    rownames(M)[-(1:m)] <- data$iden
    
    M <- M+diag(0.01,nrow(M))
    
    mod.cox.1 <- coxme(form.1, data,ties="efron",varlist =list(coxmeMlist(M)))
    LRT <- 2*(mod.cox.1$loglik[2]-mod.cox.null$loglik[2])
      
    res.coxme <- rbind(res.coxme,LRT)
    rownames(res.coxme)[nrow(res.coxme)] <- paste("pos",i,sep="")
  }
print(res.coxme)
  return(res.coxme)
}

