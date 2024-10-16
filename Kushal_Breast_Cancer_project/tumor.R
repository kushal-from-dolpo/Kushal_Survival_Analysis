#install.packages("survminer","survival","car","glmnet","cmprsk")
library(survminer)
library(survival)
library(car)
library(glmnet)
library(cmprsk)

# Data
bio<-read.table("data/SurvData6.txt", sep="",header=T)
names(bio)<-c("hormon","age","menopause","size","cstage","posnods","prog","oestr","time","censor")
attach(bio)

#transform some variable values
cstagef<-as.factor(cstage)
menopausef<-as.factor(menopause)
hormonf<-as.factor(hormon)
time<-time/356

split_variable<-function(variable ,mid1, mid2){
  for(i in 1:length(variable)){
    if(variable[i] < mid1){
      variable[i]<-1
    }
    else if(variable[i] >= mid1 & variable[i] < mid2){
      variable[i]<-2
    } else {
      variable[i]<-3
    }
  }
  return(variable)
}

#split age variable in 3 
agef<-as.factor(split_variable(age,40,60))

#split size variable in 3 
sizef<-as.factor(split_variable(size,40,80))

#split posnods variable in 3 
fpnods<-as.factor(split_variable(fpnods,10,20))

#split prog variable in 3 
prog3<-as.factor(split_variable(prog,15,100))

#split oestr variable in 3 
oestr3<-as.factor(split_variable(oestr3,15,100))

# create a new dataframe with the transformed variables
bio1<-data.frame(hormonf, age, menopausef, size, cstagef, posnods, prog, oestr, time, censor)
attach(bio1)

## descriptives
barplot(table(age3), xlab = "Age", ylim = c(0,350), main="Patient Age", xaxt = "n", col ="lightpink2")
axis(1, at=1:3, labels=c("<40",">40 & <60",">60"))
range(age)
hist(size,col="lightpink3", main = "Tumor size", xlab = "Size (mm)")
barplot(table(size3), xlab = "Size(mm)", ylim = c(0,470), main="Tumor size categorized", col ="slategray3", xaxt = "n")
axis(1, at=1:3, labels=c("<40",">40 & <80",">80"))

hist(posnods, col = "purple")
hist(oestr, col="purple")
barplot(table(cstage), xlab = "Stage", main = "Cancer stage", col="slategray4")

## Survival - Kaplan Mayer
s<-Surv(time, censor)

tumor.surv<-survfit( s~ 1, conf.type="none") # The K-M
summary(tumor.surv)

# Plot the K-M to see the survival probability over years
plot(tumor.surv, xlab="Time", ylab="Survival Probability", col="deepskyblue4", lwd=4, main="Survival Propability" )
abline(h=0.5, col=2, lwd=2)

## Kaplan Mayer with age
# Plot the K-M by age
KMage <- survfit(s ~ age3, data=bio)
ggsurvplot(KMage, surv.median.line = "hv", legend.labs=c("20-40","40-60","60-80"), title="Survival probablity depending on age")

## Kaplan Mayer with menopause
# Plot the K-M by menopause
KMmeno <- survfit(s ~ menopause, data=bio)
ggsurvplot(KMmeno, surv.median.line = "hv", legend.labs=c("no","yes"), title="Survival probablity depending on menopause")

## Kaplan Mayer with Hormone
# Plot the K-M by hormone
KMhorm <- survfit(s ~ hormon, data=bio)
ggsurvplot(KMhorm, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Placebo", "Hormone"), legend.title="Therapy",  
           palette=c("orchid2","darkcyan"), 
           title="Overall survival by the choice of the treatment", 
           risk.table.height=.15)


## Kaplan Mayer with Stage
# Plot the K-M by stage
KMstage<- survfit(s ~ cstage, data=bio)
ggsurvplot(KMstage, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("1", "2","3"), legend.title="Stage",  
           palette=c("orchid2","darkcyan","dodgerblue2"), 
           title="Overall survival by the cancer stage", 
           risk.table.height=.15)

## Highly statisticaly significant the differences by cancer stage

## Kaplan Mayer with size
# Plot the K-M by size
KMsize <- survfit(s ~ size3, data=bio)
ggsurvplot(KMsize, surv.median.line = "hv", legend.labs=c("<40mm","40-80mm",">80mm"), title="Overall survival by the tumor size")


## Kaplan Mayer with POSITIVE NODES
# Plot the K-M by pnods
KMnods <- survfit(s ~ fpnods, data=bio)
ggsurvplot(KMnods, surv.median.line = "hv", legend.labs=c("<10",">10 & <34.3",">34.5"), title="Overall survival by the number of positive nods")

## Highly significant the number of positive nods

## Kaplan Mayer with PROGESTERON
# Plot the K-M by prog
KMprog <- survfit(s ~ prog3, data=bio)
ggsurvplot(KMprog, surv.median.line = "hv", legend.labs=c("<15", ">15 & <100",">100"), title="Overall survival by the levels of progesteron")

## Highly significant the levels of progesteron

## Kaplan Mayer with OESTROGEN
# Plot the K-M by prog
KMoestr <- survfit(s ~ oestr3, data=bio)
ggsurvplot(KMoestr, surv.median.line = "hv", legend.labs=c("<15", ">15 & <100",">100"), title="Overall survival by the levels of Oestrogen")

## Highly significant the levels of oestrogen

## COX MODEL
cstage<-as.factor(cstage)

model_full<-coxph(s ~ hormonf+age+menopausef+size+cstagef+posnods+prog+oestr, data=bio1)
summary(model_full)

vif(model_full)

# STEP
step_model<-step(model_full, direction = "both", data=bio1)
summary(step_model)

## LASSO
X<-model.matrix(model_full)
cv.fit_full<-cv.glmnet(X, Surv(time,censor), family="cox")
plot(cv.fit_full)
cv.fit_full$lambda.min
coeffici_lasso<-coef(cv.fit_full, s=cv.fit_full$lambda.min)
coeffici_lasso

## DIAGNOSTICS

## COX_SNELL
rc<-abs(censor-step_model$residuals)
km.rc<-survfit(Surv(rc, censor)~1)
summary.km.rc<-summary(km.rc)
rcu<-summary.km.rc$time
surv.rc<-summary.km.rc$surv
plot(rcu, -log(surv.rc),type="p", pch=".", xlab="Cox-Snell residuals", ylab="Cumulative hazard", lwd=2, main="Cox-Snell Residuals plot")
abline(a=0, b=1, col=2, ld=3);

## Schoenfeld residuals(Testing proportional Hazards assumption)
test.ph <- cox.zph(step_model)
test.ph
ggcoxzph(test.ph)

## Deviance residuals
ggcoxdiagnostics(step_model, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

