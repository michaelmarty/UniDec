library(readr)
library(car)
library(emmeans)
FRET <- read_csv("r_formatted_FRET.csv")

my_anova <- lm(FRET_corrected ~ Group, data = FRET)
Anova(my_anova)
emmeans(my_anova, pairwise~Group)

aggregate(FRET_corrected ~ Group, data=FRET ,FUN="mean")
aggregate(FRET_corrected ~ Group, data = FRET, FUN="sd")

summary_df <- data.frame(
  FRET_corrected = c(879, 1485.667, 853, 3630), # From Aggregate
  sds = c(290.0879, 207.4231, 272.4408, 218.8241), # From Aggregate
  Group = c("DMPC Annexin No Ca", "DMPC Annexin with Ca", "DMPG Annexin No Ca", "DMPG Annexin with Ca"),
  Calcium = c("No", "Yes", "No", "Yes"))

ggplot(summary_df, aes(x=Group, y=FRET_corrected, fill = Calcium)) + 
  geom_col() + geom_errorbar(aes(ymin = FRET_corrected - sds, ymax=  FRET_corrected + sds), width=0.2) + 
  theme_classic() + scale_fill_manual(values=c("#E5CDF6", "#Fb6b6b")) + 
  ylab("Corrected FRET Intensity") + xlab("Sample Name")

ggsave("FRET_barplot.pdf")
dev.off()
