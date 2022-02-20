Original_Data <- read_delim("Skin_NonSkin.txt", 
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE, trim_ws = TRUE)
Original_Data<-Original_Data[,c(4,1:3)]
Original_Data[,1]<-ifelse(Original_Data[,1]==1,0,1)
no_of_Variables<-ncol(Original_Data[,-1])
Original_Data[,-1]<-apply(Original_Data[,-1],2,scale)
Original_Data<-cbind(Original_Data[,1],1,Original_Data[,-1],Original_Data[,-1]^2)
colnames(Original_Data)<-c("Y","X0",paste0("X",1:no_of_Variables),
                           paste0("X",1:no_of_Variables,"^2"))

glm(Y~.-1,data=Original_Data,family = "binomial")->Full_Model

stepAIC(Full_Model)

glm(Y~X0+X1+X2+X3,data=Original_Data,family = "binomial")->Covariate_Model

stepAIC(Covariate_Model)

glm(Y~.-1,data=Original_Data[,-c(4,7)],family = "binomial")->Alternate_Model_2

stepAIC(Alternate_Model_2)

Full_Model$aic
Alternate_Model$aic
Alternate_Model_2$aic

corrr::correlate(Original_Data[,c(3:5)])

plot(Original_Data[,c(3,4)])
plot(Original_Data[,c(3,5)])
plot(Original_Data[,c(4,5)])

Squared_Terms<-paste0("X",1:no_of_Variables,"^2")
term_no <- 2

All_Models <- list(c("X0",paste0("X",1:no_of_Variables)))

for (i in 1:no_of_Variables)
  {
  x <- as.vector(combn(Squared_Terms,i,simplify = FALSE))
  for(j in 1:length(x))
    {
    All_Models[[term_no]] <- c("X0",paste0("X",1:no_of_Variables),x[[j]])
    term_no <- term_no+1
    }
  }

Subsample_Size<-2000; 
r0<-Nc_size<-200; 
Replicates <- 1000; 
Choice <-100*seq(4,20,1)

save(list = c("Original_Data","Subsample_Size","Replicates","Nc_size","r0","All_Models","Choice"),
     file=here("Identical_r0","Generate_Big_Data","Scaled.RData"))

save(list = c("Original_Data","Subsample_Size","Replicates","Nc_size","r0","All_Models","Choice"),
     file=here("Non_Identical_r0","Generate_Big_Data","Scaled.RData"))

rm(list = ls())

