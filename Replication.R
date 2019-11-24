library(survival)
library(dplyr)
library(survminer)
library(stargazer)
library(ggplot2)
library(SurvRegCensCov)

###### OPEN DATA FILES ######

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
io_data <- read.csv("io_data.csv")

#membership form (KT: converting numbers to text)
io_data$memb_formc <- ifelse(io_data$memb_form == 1, 'open', NA)
io_data$memb_formc <- ifelse(io_data$memb_form == 2, 'geo_resticted', io_data$memb_formc)
io_data$memb_formc <- ifelse(io_data$memb_form == 3, 'purp_resticted', io_data$memb_formc)
io_data$memb_formc <- as.factor(io_data$memb_formc)



######### MODELS ############
##Cox model
#model

#KT:the following line is wrong and causes an error because the scientiest caled them as a scale and not as funktion
#io_data[8,13] <- 'economic'
io_data$function_c <- droplevels(io_data$function_c)


io_data$membership_growth <- io_data$full_membership_end - io_data$full_membership_start

#KT: io_data1 <- io_data[-83,] 
io_data$memb_form <- as.factor(io_data$memb_form)

io_data$Global <- ifelse(io_data$region_c == "Global", 0,1 )
io_data$region_c = factor(io_data$region_c,levels(io_data$region_c)[c(5,1:4, 6:7)])
#KT: The author removes MiddleEast OceaniaAntartica and changes the order


write.csv(io_data, "final_dataset.csv")


coxph <- coxph(Surv(lifespan,event) ~ full_membership_end + region_c + function_c + scope_c, data = io_data, method="breslow")

summary(coxph)



#KT: Technical organisations
#org_technical <- subset(io_data, function_c == "technical")
#org_technical$membership_growth <- org_technical$full_membership_end - org_technical$full_membership_start
#Cox nur für Technical Organisations
#coxph_tech <- coxph(Surv(lifespan,event) ~ function_c + memb_formc + scope_c  , data = org_technical, method="breslow")
#summary(coxph_tech)


### Table 1 model specification ####
coxphT1_succ1 <- coxph(Surv(lifespan,event) ~ full_membership_end  + region_c + function_c  + memb_formc + scope_c  , data = io_data, method="breslow")
coxphT1_succ2 <- coxph(Surv(life_1,event_1) ~ full_membership_end  + region_c + function_c  + memb_formc + scope_c  , data = io_data, method="breslow")

summary(coxphT1_succ1)
summary(coxphT1_succ2)

#KT:
ggsurvplot(survfit(coxphT1_succ1), data = io_data, color = "#2E9FDF",
           risk.table = TRUE, risk.table.col="strata",
           ggtheme = theme_minimal())

ggsurvplot(survfit(coxphT1_succ2), data = io_data, color = "#2E9FDF",
           risk.table = TRUE, risk.table.col="strata",
           ggtheme = theme_minimal())

#KT: Grambsch and Therneau test, including Harrel's rho for individual covariates: 
test.ph1 <- cox.zph(coxphT1_succ1)
test.ph1
ggcoxzph(test.ph1, newpage = TRUE, font.main = 4, font.x  = c(9, "plain", "grey"), 
         font.y  = c(9, "plain", "grey"), font.legend = c(9, "plain", "grey"))
write.csv(data.frame(test.ph1$table), "test.ph1.csv")

test.ph2 <- cox.zph(coxphT1_succ2)
test.ph2
ggcoxzph(test.ph2, newpage = TRUE, font.main = 4, font.x  = c(9, "plain", "grey"), 
         font.y  = c(9, "plain", "grey"), font.legend = c(9, "plain", "grey"))
write.csv(data.frame(test.ph2$table), "test.ph2.csv")

#KT: Plot Schoenfeld residuals 
resid <- residuals(coxphT1_succ1, type = "schoenfeld")
class(resid)
head(resid)
plot(resid)
#KT:
scatter.smooth(resid[, "full_membership_end"], lpars = list(lwd=3, col="red"))


coxphT1_succ_01 <- coxph(Surv(lifespan,event) ~  geo_resticted + purp_resticted + open, data = io_data, method="breslow")
summary(coxphT1_succ_01)

#KT: levels(io_data$memb_formc)
#io_data$memb_formc


#coxph1 <- coxph(Surv(life_1,event_1) ~  full_membership_end + full_membership_start  + region_c + function_c  + scope_c  + WW , data = io_data, method="breslow")
#summary(coxph1)
#coxph4 <- coxph(Surv(life_1,event_1) ~  full_membership_end  + region_c + function_c  + scope_c, data = io_data, method="breslow")
#summary(coxph4)
#coxph5 <- coxph(Surv(lifespan,event) ~ full_membership_end + membership_growth + region_c + function_c  + scope_c, data = io_data, method="breslow")
#summary(coxph5)
#coxph6 <- coxph(Surv(life_1,event_1) ~  full_membership_end  + membership_growth  +region_c + function_c  + scope_c, data = io_data, method="breslow")
#summary(coxph6)

options(scipen=999)
write.csv(data.frame(summary(coxphT1_succ1)$coefficients), "succ=death.csv")
write.csv(data.frame(summary(coxphT1_succ2)$coefficients), "succ=living.csv")


#write.csv(data.frame(summary(coxph1)$coefficients), "succ=death_end.csv")
#write.csv(data.frame(summary(coxph4)$coefficients), "succ=living_end.csv")
#write.csv(data.frame(summary(coxph5)$coefficients), "succ=death_end_growth.csv")
#write.csv(data.frame(summary(coxph6)$coefficients), "succ=living_end_growth.csv")


### Baseline hazard plot ####

baseline <- basehaz(coxphT1_succ1)
baseline1 <- basehaz(coxphT1_succ2)
plot(baseline$time, baseline$hazard, type='l',main="Cummulative Hazard rates") 

#names(baseline) <- c("hazard_1", "time")
#names(baseline1) <- c("hazard_2", "time")

baseline$model <- 1
baseline1$model <- 2

baselines_test <- rbind(baseline, baseline1)
summary(baselines_test)

baselines <- left_join(baseline, baseline1)

base_line_hazard <- ggplot(baselines_test, aes(time, hazard, colour = factor(model))) + 
      geom_line() + 
      labs(x="Time",y="Baseline hazard") + 
      theme(plot.title = element_text( color="#666666", face="bold", size=22, hjust=0)) +
      theme(axis.title = element_text( color="#666666", face="bold", size=22)) 
 

pdf('base_line_hazard.pdf',  onefile=T, paper='A4r', width = 10)
base_line_hazard
dev.off()

##### alternative plots ###
ggplot(baseline, aes(time, hazard)) + geom_line() 


cox.fit <- survfit(coxphT1_succ1, data = io_data)
p <- ggsurvplot(cox.fit, risk.table = TRUE)
p

cox.fit <- survfit(coxphT1_succ2, data = io_data)
p <- ggsurvplot(cox.fit, risk.table = TRUE)
p

p <- ggsurvplot(cox.fit)
p



##Weibull model
#model
weibull <- survreg(Surv(lifespan,event) ~ members_censor +  function_c + region + WW + memb_form + scope_c, data = io_data, dist = "weibull")

weibull <- survreg(Surv(lifespan,event) ~ members_censor +  function_c + region + WW + memb_form + scope_c , data = io_data)
summary(weibull)


weibull <- survreg(Surv(lifespan,event) ~ members_censor +  function_c + region  + memb_form + scope_c , data = io_data, dist="weibull")
summary(weibull)

weibull1 <- survreg(Surv(lifespan,event) ~ full_membership_end  + region_c + function_c  + memb_formc + scope_c  , data = io_data, dist="weibull")
weibull2 <- coxph(Surv(life_1,event_1) ~ full_membership_end  + region_c + function_c  + memb_formc + scope_c  , data = io_data, dist="weibull")

summary(weibull1)




# Plot
#KT: applying formula from previous line
weibull.fit <- survfit(Surv(lifespan,event) ~ members_censor +  function_c + region + WW + memb_form + scope_c, data = io_data)


intercept<- summary(weibull)$coefficients[['(Intercept)']]
scale<- weibull$scale

par(mfrow=c(1,2),mar=c(5.1,5.1,4.1,2.1)) # Make room for the hat.
# S(t), the survival function
curve(pweibull(x, scale=exp(intercept), shape=1/scale, lower.tail=FALSE), 
      from=0, to=200, col='red', lwd=2, ylab=expression(hat(S)(t)), xlab='t',bty='n',ylim=c(0,1))
# h(t), the hazard function
curve(dweibull(x, scale=exp(intercept), shape=1/scale)
      /pweibull(x, scale=exp(intercept), shape=1/scale, lower.tail=FALSE), 
      from=0, to=200, col='blue', lwd=2, ylab=expression(hat(h)(t)), xlab='t',bty='n')
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

WeibullDiag(Surv(lifespan, event) ~ factor(scope_c), data = io_data)
#WeibullDiag needed pacjkage library(SurvRegCensCov)

library(SurvRegCensCov)

# Kaplan-Meier non-parametric analysis
kmsurvival <- survfit(Surv(lifespan,event) ~ 1, data=io_data)
kmsurvival1 <- survfit(Surv(life_1,event_1) ~ 1, data=io_data)


p <- ggsurvplot(kmsurvival, risk.table = TRUE)
p1 <- ggsurvplot(kmsurvival1, risk.table = TRUE)

library(survminer)

pdf('KM_0.pdf',  onefile=T, paper='A4r', width = 10)
p
dev.off()

pdf('KM_1.pdf',  onefile=T, paper='A4r', width=10)
p1
dev.off()

# Kaplan-Meier non-parametric analysis by group

kmsurvival <- survfit(Surv(lifespan, event) ~ memb_formc, data=io_data)
kmsurvival1 <- survfit(Surv(life_1, event_1) ~ region_c, data=io_data)

p <- ggsurvplot(kmsurvival, risk.table = TRUE) 
p
p1 <- ggsurvplot(kmsurvival1, risk.table = TRUE)

pdf('KM_region_0.pdf',  onefile=T, paper='A4r', width = 10)
p
dev.off()

pdf('KM_region_1.pdf',  onefile=T, paper='A4r', width = 10)
p1
dev.off()

ggsurvplot(
  kmsurvival1,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
 # conf.int = TRUE,         # show confidence intervals for 
 # point estimaes of survival curves.
  xlim = c(0,200),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 10,     # break X axis in time intervals by 500.
   # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)





