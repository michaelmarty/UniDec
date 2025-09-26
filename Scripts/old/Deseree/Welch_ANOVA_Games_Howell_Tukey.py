install.packages(c("ggplot2", "ggpubr", "tidyverse", "broom", "AICcmodavg"))
install.packages('ufs')
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)


DaptoOxRaw.data <- read.csv("Z:/Group Share/Deseree/FPOP Paper 2/4R_Raw_June28.csv", header = TRUE, colClasses = c("numeric", "factor", "factor", "numeric"))

summary(DaptoOxRaw.data)
ggplot(DaptoOxRaw.data, aes(x=LIPID, y=OX, fill=LIPID)) + geom_boxplot()

one.way <- aov(OX ~ LIPID, data = DaptoOxRaw.data)

one.way.ratio <- aov(OX ~ RATIO, data = DaptoOxRaw.data)

two.way <- aov(OX ~ RATIO + LIPID, data = DaptoOxRaw.data)

interaction <- aov(OX ~ RATIO * LIPID, data = DaptoOxRaw.data)

oneway.test <- aov(OX ~ LIPID, data = DaptoOxRaw.data)

summary(one.way)
summary(two.way)
summary(interaction)
summary(oneway.test)

model.set <- list(one.way, two.way, interaction, oneway.test, one.way.ratio)
model.names <- c("one.way", "two.way", "Interaction","Welch", "ratio")

aictab(model.set, modnames = model.names)

#Welch<-oneway.test(OX ~ LIPID, data = DaptoOxRaw.data)

one.way <- oneway(DaptoOxRaw.data$LIPID, y = DaptoOxRaw.data$OX, posthoc = 'games-howell') one.way

tukey.two.way<-TukeyHSD(interaction)

tukey.two.way