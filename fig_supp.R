## fig_supp.m
# Simon Frew | NNL | BCCHRI
# supplementary figure generation and data processing in R

library(tidyverse)
library(ggplot2)
theme_set(theme_classic())
theme_update(axis.line = element_line(size=2), axis.ticks = element_line(size=2), legend.position = c(0.8,0.9), 
             text = element_text(size = 20), aspect.ratio = 0.6, axis.ticks.length = unit(0.25, "cm"))

setwd("path/to/wd")

dat.hm = as.tibble(read.csv("hm_table.csv"))
head(dat.hm)


# Interleaved bar chart of head motion  ---------------------------------------------------

dat2 = dat.hm[c("age", "age_rounded", "r1_mean_fd", "mDM_mean_fd")] %>%
  pivot_longer(cols = c("r1_mean_fd", "mDM_mean_fd"), names_to = "condition", values_to = "mean_FD")

dat2$condition = recode(dat2$condition, r1_mean_fd = "Rest1", mDM_mean_fd = "Movie-H")


dat.summary <- dat2 %>%
  group_by(age_rounded, condition) %>%
  summarize(
    Mean = mean(mean_FD),
    SD = sd(mean_FD)
  ) %>%
  arrange( age_rounded, desc(condition) )

ggplot(dat.summary, aes(x = age_rounded, y = Mean, fill = condition, order = condition)) + 
  geom_col(width = 0.8, position=position_dodge2(reverse = TRUE, padding = 0.1, width = 0.5)) + 
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), width = 0.8, position = position_dodge2(width = 0.5, reverse = TRUE)) + 
  labs(fill = "", title = "HBN-1388: Mean FD by Age, Condition", x = "Age (y)", y = "Mean FD (mm)") + 
  scale_x_continuous(breaks = 5:21) + 
  scale_y_continuous(expand = c(0,0), lim = c(0, 5)) +
  scale_fill_manual(values = c("#e68b2c","#225f80"))
ggsave("figureS2.png", device = "png")
ggsave("figureS2.tiff", device = "tiff")
ggsave("figureS2.pdf", device = "pdf")



# CONNERS -----------------------------------------------------------------


dat.c3sr = read.csv(file.path("..", "data", "folder", "filename.csv"))
dat.c3sr = as.tibble(dat.c3sr[-1,c(1,5,51:62)] )
dat.c3sr[3:14] = lapply( dat.c3sr[3:14], as.numeric)
head(dat.c3sr)

dat.all = left_join(dat.hm, dat.c3sr, by = c("id" = "EID"))
dat.all = distinct(dat.all, id, .keep_all = TRUE)

# export for PRISM
write.csv(data.matrix(dat.all[, c(6,8,12,16,18)]), file = "dat-conners.csv")


# export for histogram figure, Aggression
dat.all %>% 
  mutate( sex = factor(sex, levels = c(0,1), labels = c("Male", "Female")), sex = fct_relevel(sex, "Female", "Male")) %>%
  ggplot(aes(x=C3SR_AG_T, fill=as.factor(sex))) + geom_histogram(binwidth = 5, boundary = 37.5, col="white") + 
    geom_vline(xintercept = 65, linetype = "dashed", color = "red", size = 2) + 
    scale_x_continuous(lim = c(), breaks = c(40, 65, 90) ) + 
    scale_y_continuous(lim = c(0, 200), expand = c(0,0)) + 
    scale_fill_manual(values = c("#000000", "#A0A0A4")) +
    labs(fill = "Mean = 56.9 ± 14.5", title="C.i. Conner's 3 Aggression Total T distribution (n=1048)", x = "Aggression Total T", y = "# of subjects")
ggsave("figureS1_C3-AG.png", device = "png")
ggsave("figureS1_C3-AG.pdf", device = "pdf")
mean(dat.all$C3SR_AG_T, na.rm=TRUE)
sd(dat.all$C3SR_AG_T, na.rm=TRUE)

# export for histogram figure, Hyperactivity
dat.all %>% 
  mutate( sex = factor(sex, levels = c(0,1), labels = c("Male", "Female")), sex = fct_relevel(sex, "Female", "Male")) %>%
  ggplot( aes(x=C3SR_HY_T, fill=as.factor(sex))) + geom_histogram(binwidth = 5, boundary = 37.5, col="white") + 
    geom_vline(xintercept = 65, linetype = "dashed", color = "red", size = 2) + 
    scale_x_continuous(lim = c(), breaks = c(40, 65, 90) ) + 
    scale_y_continuous(lim = c(0, 200), expand = c(0,0)) + 
    scale_fill_manual(values = c("#000000", "#A0A0A4")) +
    labs(fill = "Mean = 60.3 ± 13.0", title="C.ii. Conner's 3 Hyperactivity Total T distribution (n=1048)", x = "Hyperactivity Total T", y = "# of subjects")
ggsave("figureS1_C3-HY.png", device = "png")
ggsave("figureS1_C3-HY.pdf", device = "pdf")
mean(dat.all$C3SR_HY_T, na.rm=TRUE)
sd(dat.all$C3SR_HY_T, na.rm=TRUE)

# export for histogram figure, Inattention
dat.all %>% 
  mutate( sex = factor(sex, levels = c(0,1), labels = c("Male", "Female")), sex = fct_relevel(sex, "Female", "Male")) %>%
  ggplot( aes(x=C3SR_IN_T, fill=as.factor(sex))) + geom_histogram(binwidth = 5, boundary = 37.5, col="white") + 
    geom_vline(xintercept = 65, linetype = "dashed", color = "red", size = 2) + 
    scale_x_continuous(lim = c(), breaks = c(40, 65, 90) ) + 
    scale_y_continuous(lim = c(0, 200), expand = c(0,0)) + 
    scale_fill_manual(values = c("#000000", "#A0A0A4")) +
    labs(fill = "Mean = 64.6 ± 14.4", title="C.iii. Conner's 3 Inattentive Total T distribution (n=1048)", x = "Inattentive Total T", y = "# of subjects")
ggsave("figureS1_C3-IN.png", device = "png")
ggsave("figureS1_C3-IN.pdf", device = "pdf")
mean(dat.all$C3SR_IN_T, na.rm=TRUE)
sd(dat.all$C3SR_IN_T, na.rm=TRUE)



# SCARED P  ---------------------------------------------------------------

dat.scaredP = as.tibble(read.csv(file.path("..", "data", "folder", "filename.csv")))


dat.scaredP = as.tibble(dat.scaredP[-1,c(1,5,53:58)] )
dat.scaredP$SCARED_P_Total = as.numeric(dat.scaredP$SCARED_P_Total)
head(dat.scaredP)

dat.all = left_join(dat.hm, dat.scaredP, by = c("id" = "EID"))
dat.all = distinct(dat.all, id, .keep_all = TRUE)


# export for PRISM
write.csv(data.matrix(dat.all[,c(6,8,11:16)]), file = "dat-scared_p.csv")

# export for histogram figure, 
dat.all %>% 
  mutate( sex = factor(sex, levels = c(0,1), labels = c("Male", "Female")), sex = fct_relevel(sex, "Female", "Male")) %>%
  ggplot(aes(x=SCARED_P_Total, fill=as.factor(sex))) + geom_histogram(binwidth = 5, boundary = -2.5, col="white") + 
    geom_vline(xintercept = 25, linetype = "dashed", color = "red", size = 2) + 
    scale_x_continuous(lim = c(), breaks = c(0,25,60)) + 
    scale_y_continuous(lim = c(0,300), expand = c(0,0)) + 
    scale_fill_manual(values = c("#000000", "#A0A0A4")) +
    labs(fill = "Mean = 14.3 ± 11.6", title="B. SCARED P distribution (n=1257)", x = "SCARED-P Total", y = "# of subjects")
ggsave("figureS1_scared.png", device = "png")
ggsave("figureS1_scared.pdf", device = "pdf")

mean(dat.all$SCARED_P_Total, na.rm=TRUE)
sd(dat.all$SCARED_P_Total, na.rm=TRUE)
sum(!is.na(dat.all$SCARED_P_Total))






# WISC subscales ---------------------------------------------------------------

dat.WISC = as.tibble(read.csv(file.path("..", "data", "folder", "filename.csv")))


dat.WISC = as.tibble(dat.WISC[-1,c(1,5,50)] )
dat.WISC$WISC_FSIQ = as.numeric(dat.WISC$WISC_FSIQ)

head(dat.WISC)

dat.all = left_join(dat.hm, dat.WISC, by = c("id" = "EID"))
dat.all = distinct(dat.all, id, .keep_all = TRUE)

# export for PRISM
write.csv(data.matrix(dat.all[,c(3,6,8,12)]), file = "dat-WISC.csv")

# export for histogram figure 
dat.all %>% 
  mutate( sex = factor(sex, levels = c(0,1), labels = c("Male", "Female")), sex = fct_relevel(sex, "Female", "Male")) %>%
  ggplot(aes(x=WISC_FSIQ, fill=as.factor(sex))) + geom_histogram(binwidth = 5, boundary = 2.5, col="white") + 
    scale_x_continuous(breaks = seq(50, 150, 50), lim = c(40, 160)) + 
    scale_y_continuous(lim = c(0, 150), expand = c(0,0)) + 
    scale_fill_manual(values = c("#000000", "#A0A0A4")) +
    labs(fill = "Mean = 99.6 ± 16.7", title="A. WISC-V FSIQ distribution (n=1193)", x = "WISC-V FSIQ", y = "# of subjects")
ggsave("figureS1_wisc.png", device = "png")
ggsave("figureS1_wisc.pdf", device = "pdf")


mean(dat.all$WISC_FSIQ, na.rm=TRUE)
sd(dat.all$WISC_FSIQ, na.rm=TRUE)
sum(!is.na(dat.all$WISC_FSIQ))

