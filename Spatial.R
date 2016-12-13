# R code for predicting spatial effects
#--------------------------------------
# parameters:
#   - data: data.frame for dataset
#   - Ped: matrix containing pedigree information
#   - nb: an integer for the distance to consider for defining neighborhood
#   - m: an integer for the number of parents

SpatialFunction <- function(data, Ped, nb , m ){

library(BayesX)
library(spdep)
 
# For creating the variance-covariance matrix based on the predigree
#-------------------------------------------------------------------
mped <- with(Ped, 2*kinship(ID,FATHER,MOTHER))
colnames(mped)[-(1:m)] <- data$iden
rownames(mped)[-(1:m)] <- data$iden

# For constructing the structure matrix at the 'PLOT' level
#----------------------------------------------------------
data.mean <- with(data,aggregate(data.frame(x=X_POSITION,y=Y_POSITION),
                                 by=list(PLOT=PLOT),mean))
W <- NeighMatrix(data.mean, type = "euclidian", nb = nb, level = "PLOT")$C

# For coputing the CAR structure
#-------------------------------
res.GMRF <- GMRF(W , alpha = 0.80, fit, level="PLOT")
Binv <- res.GMRF$Binv
Ws <- res.GMRF$Ws
M <- res.GMRF$M
Id <- res.GMRF$Id


# Analyses for predicting spatial effects with the coxme function
# with model (2)
#----------------------------------------------------------------

# For the first symptom
#----------------------
form.T1S <- as.formula(paste("Surv(time=Y_T1S,event=EVENT_T1S,type=\"right\") ~ ",paste(c("1","PROGENY","(1|iden)", "(1|PLOT)"), collapse = "+"), sep="" )) 

mod.Y.T1S <- coxme(form.T1S,data,ties="efron",
                   varlist =list(coxmeMlist(mped),gexchange(Binv)))

blup.plot.Y.T1S <- ranef(mod.Y.T1S)$PLOT

# For the death
#--------------
# form.TD <- as.formula(paste("Surv(time=Y_TD,event=EVENT_TD,type=\"right\") ~ ",paste(c("1","PROGENY","(1|iden)", "(1|PLOT)"), collapse = "+"), sep="" ))   
# 
# mod.Y.TD <- coxme(form.TD,data,ties="efron",
#                 varlist =list(coxmeMlist(mped),gexchange(Binv)))
# 
# blup.plot.Y.TD <- ranef(mod.Y.TD)$PAREXP

return(list(T1S = blup.plot.Y.T1S))
}