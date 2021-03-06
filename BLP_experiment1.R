rm(list = ls())
#Observing the performance of BLPestimatoR of Dan Brunner (Dec 2019)
#install.packages("BLPestimatoR")
library("BLPestimatoR")
library(dplyr)

#Simulate BLP data set
K<-2 #number of random coefficients
data <- simulate_BLP_dataset(nmkt = 25, nbrn = 20,
                        Xlin = c("price", "x1", "x2", "x3", "x4", "x5"),
                        Xexo = c("x1", "x2", "x3", "x4", "x5"),
                        Xrandom = paste0("x",1:K),
                        instruments = paste0("iv",1:10),
                        true.parameters = list(Xlin.true.except.price = rep(0.2,5),
                                               Xlin.true.price = -0.2,
                                               Xrandom.true = rep(2,K),
                                               instrument.effects = rep(2,10),
                                               instrument.Xexo.effects = rep(1,5)),
                        price.endogeneity = list( mean.xi = -2,
                                                  mean.eita = 0,
                                                  cov = cbind( c(1,0.7), c(0.7,1))),
                        printlevel = 0, seed = 234234 )
head(data)

#Specify the model
model <- as.formula("shares ~ price + x1 + x2 + x3 + x4 + x5 | 
                    x1 + x2 + x3 + x4 + x5 | 
                    0+ x1 + x2 | 
                    iv1 + iv2 + iv3 + iv4 + iv5 + iv6 + iv7 + iv8 +iv9 +iv10" )

#prepare the blp data
blp_data <- BLP_data(model = model, market_identifier="cdid",
                     product_id = "prod_id",
                     productData = data,
                     integration_method = "MLHS" ,
                     integration_accuracy = 40,
                     integration_seed = 1)
theta_guesses <- matrix(c(0.5,2), nrow=2)
rownames(theta_guesses) <- c("x1","x2")
colnames(theta_guesses) <- "unobs_sd"

#estimate blp
blp_est <- estimateBLP(blp_data =blp_data,
                       par_theta2 = theta_guesses,
                       extremumCheck = FALSE ,
                       printLevel = 1 )
summary(blp_est)

#Get elasticities
shareObj <- getShareInfo( blp_data=blp_data,
                          par_theta2 = theta_guesses,
                          printLevel = 1)
get_elasticities(blp_data=blp_data,
                 share_info = shareObj ,
                 theta_lin = 1,
                 variable = "price",
                 products = c("4","20"),
                 market = 1)