#install.packages(c("ggplot2", "ggpubr", "tidyverse", "broom", "AICcmodavg"))

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)

DaptoOxRaw.data <- read.csv("Z:\\Group Share\\Deseree\\FPOP Paper 2\\4R_Raw_June14.csv", header = TRUE, colClasses = c("numeric", "factor", "factor", "numeric"))

summary(DaptoOxRaw.data)

one.way <- aov(OX ~ LIPID, data = DaptoOxRaw.data)

one.way.ratio <- aov(OX ~ RATIO, data = DaptoOxRaw.data)

two.way <- aov(OX ~ RATIO + LIPID, data = DaptoOxRaw.data)

interaction <- aov(OX ~ RATIO * LIPID, data = DaptoOxRaw.data)

interactionnew <- aov(OX ~ LIPID:RATIO, data = DaptoOxRaw.data)

summary(interaction)

#model.set <- list(one.way, two.way, interaction, one.way.ratio, interactionnew)
#model.names <- c("one.way", "two.way", "interaction", "ratio", "interactionnew")
#aictab(model.set, modnames = model.names)

tres<-TukeyHSD(interaction, conf.level = 0.95)

#tukey.two.way

dim(tres[[1]])

plot(tres)