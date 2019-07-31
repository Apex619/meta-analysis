library(openxlsx) 

Body_Weight_F0 <- read.xlsx("C:\\Users\\hamanw\\Dropbox\\Meta Analysis Project\\Extraction\\Transgenerational Meta-Analysis\\data\\Body_Weight.xlsx", sheet = 2)

Effect_Size_F0$ES_ID <-ES_ID
library(metafor)
#Calculate Effect Size (lnRR)
Effect_Size_F0 <- escalc(measure="ROM", m1i=F0_Treatment_Mean, sd1i=F0_Treatment_SD, n1i=F0_Treatment_n, m2i=F0_Control_Mean, sd2i=F0_Control_SD, n2i=F0_Control_n, data=Body_Weight_F0)

modelf0 <- rma.mv(yi, vi, mod = ~ Study, random = ~ (1|Paper_ID), data=Effect_Size_F0)
summary(modelf0)
