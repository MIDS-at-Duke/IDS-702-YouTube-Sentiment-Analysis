mg = read.csv("C:/Users/Vishaal/Desktop/IDS Final/df_average_compound_score_bycategory(many Nans)1.csv")
library(tseries)
library(car)
library(ggplot2)
library(lme4) #for the lmer function
library(rcompanion)
library(jtools)
library(interactions)

#Removing some outliers and leverage points
mg = mg[-c(11323, 5800, 17933, 479), ]

mg_complete = mg[complete.cases(mg), ]

c = table(mg_complete$compound_score)
c[names(c)==0.0]

#performing the Tukey Ladder thing. To get power transformation for normality...
T_tuk =
  transformTukey(mg_complete$L.D,
                 plotit=FALSE)


#plot of normality
hist((mg_complete$L.D)^(0.225), main = 'Histogram Checking Normality')

#converting category_id to factor
mg_complete$category_id = factor(mg_complete$category_id)

#testing a bunch of models...
regwage <- lm((mg_complete$L.D)^(0.225)~views+likes_x+dislikes+comment_total+compound_score+category_id,data= mg_complete)
summary(regwage)
plot(regwage)
plot((mg_complete$L.D)^(0.225)~mg_complete$likes_x), main = 'L/D Ratio Vs Compound Score')


regwage1 = lm((mg_complete$L.D)^(0.225)~compound_score , data = mg_complete)
summary(regwage1)
plot(regwage1)
crPlots(regwage)


regwage2 = lm((mg_complete$L.D)^(0.225)~compound_score + poitive_score + likes_x+ category_id:compound_score, data = mg_complete)
summary(regwage2)
plot(regwage2)


regwage3 = lm((mg_complete$L.D)^(0.225)~compound_score + likes_x + poitive_score + category_id +  category_id*likes_x, data = mg_complete)
summary(regwage3)
plot(regwage3)

regwage4 = lm((mg_complete$L.D)^(0.225)~compound_score  + views + category_id, data = mg_complete)
summary(regwage4)
plot(regwage4)

regwage5 = lm(((mg_complete$L.D)^(0.225))~compound_score  + views  + category_id + category_id*compound_score, data = mg_complete)
summary(regwage5)
plot(regwage5)

#plot of residual vs a continuous predictor
res = resid(regwage5)
plot(res~mg_complete$views, xlim = c(0, 100000), xlab = 'Views', ylab = 'Residual', main = 'Residual Vs Views (Continuous Predictor)')
abline(0,0)

#nested f-test to check for the significance of interaction term
anova(regwage4,regwage5)

#VIF to check for multicollinearity
vif(regwage4)

#step-wise regression backwards...
base_model = lm((mg_complete$L.D)^(0.225)~compound_score  + views + category_id + likes_x +dislikes + poitive_score+negative_score+neutral_score + compound_score*category_id, data = mg_complete)
summary(base_model)
Model_backward = step(base_model, direction='backward',trace=0)
Model_backward$call

#plot to check percent contribution to R^2 by predcitor
regwagelibrary(relaimpo)
calc.relimp(regwage,type=c("lmg","last","first","pratt"),
            rela=TRUE)
# Bootstrap Measures of Relative Importance (1000 samples)
boot <- boot.relimp(regwage, b = 1000, type = c("lmg",
                                            "last", "first", "pratt"), rank = TRUE,
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result


#dataframe containing averages by category... to get overview of all categories...
bp = read.csv("C:/Users/Vishaal/Desktop/IDS Final/boxplot_dataframe_by_categry_id1.csv")
bp$category_id = factor(bp$category_id)
boxplot(bp$category_id, bp$compound_score)

plot(L.D~compound_score,data=mg_complete,ylab="Amount of arsenic in well",
        xlab="Switched to safe well?")

#subset for a couple of categories to check for plotting interaction  plot
v = mg_complete[mg_complete$category_id == 10 | mg_complete$category_id == 2 |mg_complete$category_id == 1 |mg_complete$category_id == 43 |mg_complete$category_id ==24,]
fiti <- lm(L.D ~ category_id * compound_score + compound_score, data = v)
interact_plot(fiti, pred = compound_score, modx = category_id, main = 'Interaction Plot')

#calculating in-model and out-of-model MSE for contreversial youtube channels
c1 = read.csv("C:/Users/Vishaal/Desktop/IDS Final/controversial1.csv")
c1$category_id = factor(c1$category_id)
c = c1[complete.cases(c1), ]
plot((c$L.D)^0.225,predict(regwage5,c))
abline(0,1)
training_res = mg_complete$L.D - predict(regwage5,mg_complete)
holdout_res = c$L.D - predict(regwage5,c)
t.test(holdout_res,mu=0)
t.test(training_res,holdout_res)
sqrt(sum((predict(regwage5,c) - (c$L.D)^(0.225))^2)/9)

mse.train <- summary(regwage5)$sigma^2
mse.test  <- sum((predict(regwage5,c) - (c$L.D)^0.225)^2)/(9-2)

mse.train
mse.test

#comments data

comm = read.csv("C:/Users/Vishaal/Desktop/IDS Final/comment_corr.csv")

##important - correlation b/w title and top comment
cor(comm$compound_score_comm, comm$compound_score)

a = table(comm$compound_score)
a[names(a)==0.0]

b = table(comm$compound_score_comm)
b[names(b) == 0.0]

commo = read.csv("C:/Users/Vishaal/Desktop/IDS Final/other_comm.csv")


length(commo$compound_score)
comemm = head(comm$compound_score, -13)
length(comemm)

##important - correlation b/w top comment and all other comments. 
cor(commo$compound_score, comemm)
mg_complete

## CI for estimates of slopes and intercepts. 
confint(regwage5)
