#####################Installation################################
devtools::install_github("royarkaprava/OPM/OPMpackage")
library(OPMpackage)
#########################################################

################################################################################################
########################  USE OF THE fit.ALASSOOPM FUNCTION ####################################
################################################################################################

##################################################
##################Simulated data##################
##################################################

load("~/GitHub/OPM/data.rda")
phenotypes2 <- Macro.micro_2_test
data <- phenotypes2
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.ALASSOOPM(Ymat[, 1], X, LAMBDA = 0.5, ncores = 1, test, train, IndMa)

##################################################
##################Realdata data##################
##################################################
load("~/GitHub/OPM/Real data/Macro.micro_1.rda")
data <- Macro.micro
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.ALASSOOPM(Ymat[, 1], X, LAMBDA = 0.5, ncores = 1, test, train, IndMa)

#############################################################################################
#############################################################################################



##############################################################################################
########################  USE OF THE fit.SCADOPM FUNCTION ####################################
##############################################################################################

##################################################
##################Simulated data##################
##################################################

load("data.rda")
phenotypes2 <- Macro.micro_2_test
data <- phenotypes2
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.SCADOPM(Ymat[, 1], X, LAMBDA = 0.5, 1, test, train, IndMa)

##################################################
##################Realdata data##################
##################################################
load("~/GitHub/OPM/Real data/Macro.micro_1.rda")
data <- Macro.micro
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.SCADOPM(Ymat[, 1], X, LAMBDA = 0.5, ncores = 1, test, train, IndMa)


#############################################################################################
#############################################################################################




#############################################################################################
########################  USE OF THE fit.RAWOPM FUNCTION ####################################
#############################################################################################

##################################################
##################Simulated data##################
##################################################


load("data.rda")
phenotypes2 <- Macro.micro_2_test
data <- phenotypes2
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.RAWOPM(Ymat[, 1], X, LAMBDA = 0.5, 1, test, train, IndMa)


##################################################
##################Realdata data##################
##################################################
load("~/GitHub/OPM/Real data/Macro.micro_1.rda")
data <- Macro.micro
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.RAWOPM(Ymat[, 1], X, LAMBDA = 0.5, ncores = 1, test, train, IndMa)

#############################################################################################
#############################################################################################






#############################################################################################
########################  USE OF THE fit.AENOPM FUNCTION ####################################
#############################################################################################


##################################################
##################Simulated data##################
##################################################


load("data.rda")
phenotypes2 <- Macro.micro_2_test
data <- phenotypes2
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.AENOPM(Ymat[, 1], X, LAMBDA = 0.5, 1, test, train, IndMa)

##################################################
##################Realdata data##################
##################################################
load("~/GitHub/OPM/Real data/Macro.micro_1.rda")
data <- Macro.micro
Ymat <- data[, 2:6] #The 5 different class label vectors with different proportions of 1
X    <- data[, -c(1:6)] #The design matrix with first 20 columns are macrovariables and rest are microvariables
test <- 121:277 #Indices of the test set
train <- 1:120 #Indices of the training set
IndMa <- 20 #The number of macrovariables

out <- fit.AENOPM(Ymat[, 1], X, LAMBDA = 0.5, ncores = 1, test, train, IndMa)


#############################################################################################
#############################################################################################