library(openxlsx) 

Body_Weight_F0 <- Body_Weight <- read_excel("data/Archive/Body_Weight.xlsx", 
                                            sheet = "Body_Weight F0")

library(metafor)
#Calculate Effect Size (lnRR)

Effect_Size_F0 <- escalc(measure="ROM", m1i=F0_Treatment_Mean, sd1i=F0_Treatment_SD, n1i=F0_Treatment_n, m2i=F0_Control_Mean, sd2i=F0_Control_SD, n2i=F0_Control_n, data=Body_Weight_F0)

modelf0 <- rma.mv(yi, vi, random = ~ 1|Paper_ID, data=Effect_Size_F0)
summary(modelf0)
forest(modelf0)
funnel(modelf0)

modelf1 <- rma.mv(yi, vi, mod = ~Age_Days, random = ~ 1|Paper_ID, data=Effect_Size_F0)
summary(modelf1)

prediction <- predict(modelf1, newmods = c(50, 100, 150))
prediction

model_offspring <- rma.mv(yi, vi, mod = ~Age_Days, random = list(~1|Cohort_ID,~1|ES_ID), data=Body_Weight_lnRR_MG)

prediction_f2 <- predict(model_offspring, newmods = c(50, 100, 150))

prediction_f2


