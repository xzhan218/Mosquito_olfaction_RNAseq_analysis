#Blood-feeding success

library(lme4)
library(lmerTest)
library(ggplot2)
library(readxl)

bm <- read_excel("fed_female_number.xlsx")

str(bm)
summary(bm)
head(bm)

bm$Treatment <- factor(bm$Treatment,
                       levels=c("Clean air",
                                "Odor stimulation",
                                "Passive exposure"))


###from here use generalized linear mixed model 
bm$Rep <- factor(bm$Rep)

model <- glmer(
  cbind(Blood_fed,Unfed) ~ Treatment + (1|Rep),
  data=bm,
  family = binomial
)

summary(model)
#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: cbind(Blood_fed, Unfed) ~ Treatment + (1 | Rep)
#Data: bm

#AIC       BIC    logLik -2*log(L)  df.resid 
#139.1     143.3     -65.5     131.1        17 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.5741 -1.4782 -0.3055  1.2290  2.8469 

#Random effects:
#  Groups Name        Variance Std.Dev.
#Rep    (Intercept) 0        0       
#Number of obs: 21, groups:  Rep, 7

#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)   
#(Intercept)                 0.3808     0.1433   2.658  0.00787 **
#TreatmentOdor stimulation   0.1357     0.2074   0.654  0.51302   
#TreatmentPassive exposure  -0.5943     0.2051  -2.898  0.00375 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#  (Intr) TrtmOs
#TrtmntOdrst -0.691       
#TrtmntPssve -0.699  0.483
#optimizer (Nelder_Mead) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')

exp(0.1357)
#1.145338
exp(-0.5943)
#0.5519488

confint(model, method = "profile")
#2.5 %     97.5 %
#.sig01                     0.0000000  0.3167389
#(Intercept)                0.1011911  0.6645366
#TreatmentOdor stimulation -0.2705130  0.5433033
#TreatmentPassive exposure -0.9985538 -0.1940577

exp(confint(model, method = "profile"))
#2.5 %    97.5 %
# sig01                    1.0000000 1.3726441
#(Intercept)               1.1064880 1.9435896
#TreatmentOdor stimulation 0.7629880 1.7216848
#TreatmentPassive exposure 0.3684118 0.8236104

plot(model)

qqnorm(resid(model))
qqline(resid(model))

#Blood meal mass
library(lme4)
library(lmerTest)
library(ggplot2)
library(readxl)

bm <- read_excel("bloodmass_007022026.xlsx")

str(bm)
summary(bm)
head(bm)

bm$Treatment <- factor(bm$Treatment,
                       levels=c("Clean air",
                                "Odor stimulation",
                                "Passive exposure"))


model_log <- lmer(log1p(meal_per_female) ~ Treatment +
                    (1|Rep),
                  data=bm)

summary(model_log)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(meal_per_female) ~ Treatment + (1 | Rep)
#   Data: bm

#REML criterion at convergence: -22.6

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-2.0149 -0.4785  0.0092  0.5118  1.5766 

#Random effects:
# Groups   Name        Variance Std.Dev.
# Rep      (Intercept) 0.011451 0.10701 
# Residual             0.007491 0.08655 
#Number of obs: 22, groups:  Rep, 7

#Fixed effects:
#                          Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)                0.54747    0.05202 10.46350  10.524 6.79e-07 ***
#TreatmentOdor stimulation  0.14031    0.04505 12.97770   3.115  0.00823 ** 
#TreatmentPassive exposure  0.12488    0.04626 12.92482   2.699  0.01830 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) TrtmOs
#TrtmntOdrst -0.457       
#TrtmntPssve -0.445  0.513

confint(model_log, method = "profile")
#2.5 %    97.5 %
#.sig01                    0.04729572 0.2007202
#.sigma                    0.05848594 0.1212153
#(Intercept)               0.44412690 0.6508093
#TreatmentOdor stimulation 0.05259505 0.2281524
#TreatmentPassive exposure 0.03470246 0.2150533

plot(model_log)

qqnorm(resid(model_log))
qqline(resid(model_log))


#Diuresis

library(lme4)
library(lmerTest)
library(ggplot2)
library(readxl)

bm <- read_excel("diuresis.xlsx")

str(bm)
summary(bm)
head(bm)

bm$Treatment <- factor(bm$Treatment,
                       levels=c("Clean air",
                                "Odor stimulation",
                                "Passive exposure"))


###from here use linear mix model 
bm$Rep <- factor(bm$Rep)

model <- lmer(log1p(diuresis) ~ Treatment +
                (1|Rep),
              data=bm)

summary(model)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(diuresis) ~ Treatment + (1 | Rep)
#   Data: bm

#REML criterion at convergence: 52.4

#Scaled residuals: 
#     Min       1Q   Median       3Q      Max 
#-1.74037 -0.36391 -0.08278  0.25908  2.00676 

#Random effects:
# Groups   Name        Variance Std.Dev.
# Rep      (Intercept) 1.9279   1.3885  
# Residual             0.1768   0.4204  
#Number of obs: 23, groups:  Rep, 8

#Fixed effects:
#                          Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)                13.8491     0.5129  7.8274  27.001 5.23e-09 ***
#TreatmentOdor stimulation  -0.1365     0.2102 13.0181  -0.649  0.52732    
#TreatmentPassive exposure  -0.8359     0.2210 13.0544  -3.782  0.00227 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) TrtmOs
#TrtmntOdrst -0.205       
#TrtmntPssve -0.195  0.476


confint(model, method = "profile")
#2.5 %     97.5 %
#.sig01                     0.8336448  2.3559492
#.sigma                     0.2843208  0.5873096
#(Intercept)               12.7950486 14.9030632
#TreatmentOdor stimulation -0.5459982  0.2729457
#TreatmentPassive exposure -1.2650058 -0.4039581

plot(model)

qqnorm(resid(model))
qqline(resid(model))

#Egg production
library(lme4)
library(lmerTest)
library(ggplot2)
library(readxl)

bm <- read_excel("egg_production_07022026_use this date.xlsx")

str(bm)
summary(bm)
head(bm)

bm$Treatment <- factor(bm$Treatment,
                       levels=c("Clean air",
                                "Odor stimulation",
                                "Passive exposure"))

###from here use linear mix model 
bm$Rep <- factor(bm$Rep)

model <- lmer(log1p(egg_per_female) ~ Treatment +
                (1|Rep),
              data=bm)

summary(model)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(egg_per_female) ~ Treatment + (1 | Rep)
#   Data: bm

#REML criterion at convergence: 68.8

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-2.8776 -0.3298  0.1396  0.5173  1.7412 

#Random effects:
# Groups   Name        Variance Std.Dev.
# Rep      (Intercept) 0.04828  0.2197  
# Residual             0.75273  0.8676  
#Number of obs: 27, groups:  Rep, 3

#Fixed effects:
#                          Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)                2.58314    0.31580  8.97609   8.180 1.88e-05 ***
#TreatmentOdor stimulation  0.06651    0.40899 22.00000   0.163    0.872    
#TreatmentPassive exposure  0.04209    0.40899 22.00000   0.103    0.919    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) TrtmOs
#TrtmntOdrst -0.648       
#TrtmntPssve -0.648  0.500


confint(model, method = "profile")
#2.5 %    97.5 %
#.sig01                     0.0000000 0.8310954
#.sigma                     0.6414133 1.1209497
#(Intercept)                1.9791681 3.1871068
#TreatmentOdor stimulation -0.7326452 0.8656609
#TreatmentPassive exposure -0.7570618 0.8412443

plot(model)

qqnorm(resid(model))
qqline(resid(model))

#Enzyme activity
library(readxl)

#=====================================trypsin activity==========================
df <- read_excel("trypsin_R.xlsx")

df$condition <- factor(
  df$condition,
  levels = c("clean", "human", "normal")
)

library(lme4)
library(lmerTest)

m_enzyme_log <- lmer(log1p(enzyme) ~ condition + (1 | replicate), data = df)
summary(m_enzyme_log)

#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(enzyme) ~ condition + (1 | replicate)
#   Data: df

#REML criterion at convergence: -12.9

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-2.6303 -0.5857  0.2316  0.6343  1.9236 

#Random effects:
# Groups    Name        Variance Std.Dev.
# replicate (Intercept) 0.004666 0.06831 
# Residual              0.037368 0.19331 
#Number of obs: 58, groups:  replicate, 4

#Fixed effects:
#                Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      2.58759    0.05537  7.33554  46.733 2.38e-10 ***
#conditionhuman  -0.06406    0.05888 52.67524  -1.088   0.2816    
#conditionnormal  0.16190    0.06794 53.03381   2.383   0.0208 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) cndtnh
#conditinhmn -0.582       
#conditnnrml -0.486  0.454

#calculate CI_regular
confint(m_enzyme_log, method = "Wald")

#lmer preferred CI_use this one (very similar)
confint(m_enzyme_log, method = "profile")
#2.5 %     97.5 %
#  .sig01           0.00000000 0.17613259
#.sigma           0.15895975 0.23247380
#(Intercept)      2.47937409 2.69925389
#conditionhuman  -0.17695464 0.05578292
#conditionnormal  0.02287283 0.29188040

#=======================chymotrypsin activity===================================
df <- read_excel("chymo_R.xlsx")

df$condition <- factor(
  df$condition,
  levels = c("clean", "human", "normal")
)

m_enzyme_log <- lmer(log1p(enzyme) ~ condition + (1 | replicate), data = df)
summary(m_enzyme_log)

#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(enzyme) ~ condition + (1 | replicate)
#   Data: df

#REML criterion at convergence: -82.3

#Scaled residuals: 
#     Min       1Q   Median       3Q      Max 
#-2.68644 -0.33481  0.09538  0.70743  1.98443 

#Random effects:
# Groups    Name        Variance  Std.Dev.
# replicate (Intercept) 0.0006305 0.02511 
# Residual              0.0101551 0.10077 
#Number of obs: 56, groups:  replicate, 4

#Fixed effects:
#                Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      2.52537    0.02591 10.93465  97.474   <2e-16 ***
#conditionhuman   0.02044    0.03132 51.56834   0.653    0.517    
#conditionnormal  0.05601    0.03532 51.59740   1.586    0.119    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) cndtnh
#conditinhmn -0.634       
#conditnnrml -0.547  0.449

#lmer preferred CI_use this one (very similar)
confint(m_enzyme_log, method = "profile")
#2.5 %     97.5 %
#  .sig01           0.00000000 0.07262701
#.sigma           0.08253082 0.12126834
#(Intercept)      2.47566373 2.57660590
#conditionhuman  -0.04000442 0.08258090
#conditionnormal -0.01291145 0.12503784

#P.falciparum infection
library(readxl)
df <- read_excel("infection_for_R.xlsx")

library(dplyr)

df_intensity <- df %>%
  dplyr::filter(oocyst_count > 0) %>%          # intensity only
  mutate(
    condition = factor(condition, levels = c("clean", "human", "normal")),
    replicate = factor(replicate)
  )

library(lme4)
library(lmerTest)

model_lmm <- lmer(
  log1p(oocyst_count) ~ condition + (1 | replicate),
  data = df_intensity
)

summary(model_lmm)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(oocyst_count) ~ condition + (1 | replicate)
#   Data: df_intensity

#REML criterion at convergence: 368.8

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-2.3109 -0.5170  0.1472  0.6867  1.9225 

#Random effects:
# Groups    Name        Variance Std.Dev.
# replicate (Intercept) 0.08952  0.2992  
# Residual              0.90589  0.9518  
#Number of obs: 132, groups:  replicate, 4

#Fixed effects:
#                 Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)       2.86164    0.20695   5.89383  13.828 1.03e-05 ***
#conditionhuman   -0.39310    0.20822 128.97627  -1.888   0.0613 .  
#conditionnormal   0.01743    0.20278 127.51105   0.086   0.9317    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) cndtnh
#conditinhmn -0.458       
#conditnnrml -0.479  0.453

#calculate CI
confint(model_lmm, method = "profile")
#2.5 %     97.5 %
#  .sig01           0.0000000 0.71000227
#.sigma           0.8396349 1.07323020
#(Intercept)      2.4345262 3.27608687
#conditionhuman  -0.7931717 0.02702193
#conditionnormal -0.3836081 0.41104529

#Microbiota
library(readxl)

df <-  read_excel("microbiota_24h_R.xlsx")

df$condition <- factor(
  df$condition,
  levels = c("clean", "human", "normal")
)

library(lme4)
library(lmerTest)

m_microbiota_log <- lmer(log1p(microbiota) ~ condition + (1 | replicate), data = df)
summary(m_microbiota_log)

#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: log1p(microbiota) ~ condition + (1 | replicate)
#   Data: df

#REML criterion at convergence: -10

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-0.4021 -0.2956 -0.2516 -0.1753  6.3723 

#Random effects:
# Groups    Name        Variance  Std.Dev.
# replicate (Intercept) 0.0006219 0.02494 
# Residual              0.0479546 0.21899 
#Number of obs: 112, groups:  replicate, 5

#Fixed effects:
#                  Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept)       0.076837   0.034864  17.748603   2.204    0.041 *
#conditionhuman   -0.024770   0.046688 105.182770  -0.531    0.597  
#conditionnormal  -0.005265   0.056344  97.958647  -0.093    0.926  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#            (Intr) cndtnh
#conditinhmn -0.670       
#conditnnrml -0.560  0.414

confint(m_microbiota_log, method = "profile")
#2.5 %     97.5 %
#  .sig01           0.00000000 0.09020466
#.sigma           0.19093365 0.24900587
#(Intercept)      0.00964839 0.14418807
#conditionhuman  -0.11628303 0.06674392
#conditionnormal -0.11372118 0.10482862
