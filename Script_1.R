library(MASS)
library(fic)
studentnumber =  0787524
fulldata = read.table("dataHW.txt",header=T)
set.seed(studentnumber)
rownumbers = sample(1:nrow(fulldata),600,replace=F)
mydata = fulldata[rownumbers,]


hist(mydata$Y,
main = "Empirical and Poisson distributions comparison", xlab="Y values",
 )

par(new=TRUE)
plot(0:20, dpois( x=0:20, lambda=3.5), xlim=c(-2,20),xaxt='n', ann=FALSE, 
     yaxt='n',main = paste("Empirical and Poisson distributions comparison"))
normden <- function(x){dnorm(x, mean=3.5, sd=sqrt(3.5))}
curve(normden, from=-2, to=20, add=TRUE, col="red")

 
hw.full=glm(Y~.^2, data = mydata, family = poisson())
summary(hw.full)
#use stepAIC to find the best models based on AIC
hw.AIC=stepAIC(hw.full, k=2, scope =list(upper=~.^2,lower=~1) ,
               direction = "backward")


#fit the best three models selected by AIC(the ones with the LOWEST AIC)
fit1.aic=glm(Y~X1+X2+X3+X4+X1*X2+X3*X4, data = mydata, family = poisson())
fit2.aic=update(fit1.aic,.~.+ X5
)
fit3.aic=update(fit2.aic,.~.+ X4*X5)


summary(fit1.aic)
summary(fit2.aic)
summary(fit3.aic)

#calculate aic
AIC1=AIC(fit1.aic,k=2)
AIC2=AIC(fit2.aic,k=2)
AIC3=AIC(fit3.aic,k=2)

#Build dataframe with 3 best models based on aic


#use stepAIC to find the best models based on BIC
hw.BIC=stepAIC(hw.full, k=log(nrow(mydata)), scope =list(upper=~.^2,lower=~1) ,
               direction = "backward")


### fit 3 best models according to 3 final steps of stepAIC based on BIC
###
###
fit1.bic=glm(Y~X1+X2+X3+X4+X3*X4+X1*X2, data = mydata, family = poisson())
fit2.bic=update(fit1.bic,.~.+X5)
fit3.bic=update(fit2.bic,.~.+X5*X4)


summary(fit1.bic)
summary(fit2.bic)
summary(fit3.bic)

#calculate bic : in both BIC and AIC the lower the IC value the best is the model
BIC1=AIC(fit1.bic,k=log(nrow(mydata)))
BIC2=AIC(fit2.bic,k=log(nrow(mydata)))
BIC3=AIC(fit3.bic,k=log(nrow(mydata)))



#(b) FIC

library(fic)

#(i) specify the focus function
## we are using the log function as the link function and we focus on the mean
## response value
focus=function(par,X) exp(X%*%par)

variable=model.matrix(hw.full)
matrix.cor=cor(variable)

matrix.cor=matrix.cor[2:nrow(matrix.cor),]
matrix.cor=matrix.cor[,2:ncol(matrix.cor)]
matrix.cor1=matrix.cor[,1:3]
#plot of correlation matrix
library('corrplot') 
corrplot(matrix.cor, method = "circle")





## remove interactions that cause singularity problems
## remove all interactions: X1&X4, X1&X5 , X2&X4, X3&X4, X3&X5,X4&X5,
hw.full.singular=glm(Y~.^2-X1:X4-X1:X5-X2:X4-X3:X4-X3:X5-X4:X5-X2:X5, data = mydata, 
                 family = poisson())
# hw.full.singular=glm(Y~.^2-X1:X4-X1:X5-X2:X4-X3:X4-X3:X5-X4:X5, data = mydata,
#                      family = poisson())
summary(hw.full.singular)



inds0=c(1,rep(0,length(hw.full.singular$coefficients)-1))
combs=all_inds(wide = hw.full.singular, inds0 = inds0)

#exclude models with interactions that do not include both main effects
combs.excl=with(combs,combs[!(
                            combs[,2]==0 & (combs[,8]==1|combs[,7]==1)|
                            combs[,3]==0 & (combs[,7]==1|combs[,9]==1)|
                            combs[,4]==0 & (combs[,8]==1|combs[,9]==1 ))
                                ,])
# combs.excl=with(combs,combs[!(
#     combs[,2]==0 & (combs[,8]==1|combs[,7]==1)|
#         combs[,3]==0 & (combs[,7]==1|combs[,9]==1 |combs[,10]==1)|
#         combs[,4]==0 & (combs[,8]==1|combs[,9]==1 )|
#         combs[,6]==0 & (combs[,10]==1))
#     ,])


# specify X matrix, which contains the focus we'll use for the evaluation
X_eval1=round(c(1,mean(mydata$X1),
                mean(mydata$X2),
                mean(mydata$X3),
                mean(mydata$X4),
                mean(mydata$X5),
                mean(mydata$X1*mydata$X2),
                mean(mydata$X1*mydata$X3),
                # mean(mydata$X2*mydata$X5),
                mean(mydata$X2*mydata$X3)),
              digits = 2)
X_eval2=round(c(1,var(mydata$X1),
                var(mydata$X2),
                var(mydata$X3),
                var(mydata$X4),
                var(mydata$X5),
                var(mydata$X1*mydata$X2),
                var(mydata$X1*mydata$X3),
                # mean(mydata$X2*mydata$X5),
                var(mydata$X2*mydata$X3)),
              digits = 2)


# Implement fic for focus 1
fic1=fic(wide = hw.full.singular, inds = combs.excl, inds0 = inds0, 
         focus = focus, X=X_eval1)

index1=order(fic1$rmse.adj)[1:3]
fic1[index1,]
#Y~Intercept + X4
#Y~Intercept + X3
#Y~Intercept + X5




# Implement fic for focus 2
fic2=fic(wide = hw.full.singular, inds = combs.excl, inds0 = inds0, 
         focus = focus, X=X_eval2)

index2=order(fic2$rmse.adj)[1:3]
fic2[index2,]
#Y~Intercept + X4
#Y~Intercept + X3
#Y~Intercept + X3 + X4




