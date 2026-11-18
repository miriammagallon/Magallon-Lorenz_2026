###################################
#     Kapplan meyer Analysis      #
###################################
# Packages needed
library(lubridate)
library(ggsurvfit)
library(gtsummary)
# install.packages(c("tidycmprsk"))
library(tidycmprsk)
# devtools::install_github("zabore/condsurv")
library(condsurv)
library(dplyr)
library(survival)
library(scales)
library(survminer)
library(forcats)
library(gtsummary)
library(webshot2)
library(kableExtra)
 
#Functions
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}



# Data generation
execution.dir <- "."
kapplan.df <- read.table(file.path(execution.dir,"results","clinical_data_Miller.csv"), header = T, sep = "\t")

# Table 1
(t1 <- kapplan.df %>%
    mutate(Decesed = factor(Decesed, levels = c(FALSE, TRUE), labels = c("censored", "died"))) %>%
    tbl_summary(by = "Groups", include = c(Overall.survival.months, Decesed, Groups)) %>%
    add_overall()) %>%
  as_flex_table() %>%
  flextable::set_caption("Summary statistics.") %>%
  flextable::theme_zebra()
t1 <- as_gt(t1)
gt::gtsave(data = t1, filename =file.path("./images",paste0("SummaryStatisticsTablebyGroups_v5.png")))
kapplan.df$Age.months.onset.MPNST<- as.numeric(kapplan.df$Age.months.onset.MPNST)

# Table 2
(t1 <- kapplan.df %>%
    mutate(Decesed = factor(Decesed, levels = c(FALSE, TRUE), labels = c("censored", "died"))) %>%
    tbl_summary(by = "Groups", include = c(Age.months.onset.MPNST, Decesed, Groups)) %>%
    add_overall()) %>%
  as_flex_table() %>%
  flextable::set_caption("Summary statistics.") %>%
  flextable::theme_zebra()
t1 <- as_gt(t1)
gt::gtsave(data = t1, filename =file.path("./images",paste0("SummaryStatisticsTablebyGroupsAgeMPNSTdiagnosis_v2.png")))

#Participant censoring
kapplan.df %>%
  filter(Decesed == FALSE ) %>%
  ggplot(aes(x = Overall.survival.months, y = fct_rev(Groups))) +
  geom_point() +
  theme_light() +
  labs(title = "Participant Censoring", x = "Time (days)", y = NULL)

################ Kaplan Meier curves 1 group vs the others #######################

# workaround functions to save kaplan meier curves function
ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}


#representation
groups <- c("1","2","3","4")
names(groups) <-c("#7C2C00", "tomato","goldenrod", "#E09B9B")

plot.data <- kapplan.df
plot.data$Diagnosed <- TRUE
rownames(plot.data) <- kapplan.df$Sample.names

barplot(table(plot.data[plot.data$Groups ==4, "Sex"]))
i=1

for(i in seq_len(length(groups))){
  col.g <- names(groups)[i]
  
  #### Overall survival
  #We create a variable containing the group selected vs the others
  plot.data <- plot.data %>% mutate(Group.vs.others = ifelse(Groups == i, paste0("Group ",i),"Others"))
  plot.data$Group.vs.others <- factor(plot.data$Group.vs.others, levels = c(paste0("Group ",i), "Others"))
  fit <- survfit(Surv(Overall.survival.months, Decesed)~Group.vs.others, data =plot.data)
  pp <- ggsurvplot(fit, 
              data = plot.data, 
              pval =T,
              pval.size =35,
              lwd = 10,
              font.main = c(65, "plain"),
              font.x = c(80),
              font.y = c(80),
              font.tickslab = c(65),
              font.legend = c(50),
              censor.size=40,
              # risk.table = TRUE,
              # risk.table.fontsize = 6,
              # risk.table.height = 0.25,
              ggtheme = theme_bw(base_size = 80),
              palette = c(col.g,"grey65"),
              legend.labs = levels(plot.data$Group.vs.others),
              legend.title = "",
              # title = "Overall Survival"
              )#, facet.by = "Groups", title = "Test")
  g_to_save <- ggsave_workaround(pp)
  
  ggsave(plot = g_to_save,
         filename =paste0("kaplanmeier_Survival_g",i,"vsOthers_v5.svg"), path =file.path("./images"),
         device = "svg",
         width = 18, height =18, dpi = 300)
  
  # OS NF1vsSP per group
  plot.data.g <- plot.data[plot.data$Groups ==i,]
  fit <- survfit(Surv(Overall.survival.months, Decesed)~NF.status, data =plot.data.g) 
  pp <-  ggsurvplot(fit, 
                    data = plot.data.g, 
                    pval =T,
                    pval.size =35,
                    lwd = 10,
                    font.main = c(65, "plain"),
                    font.x = c(80),
                    font.y = c(80),
                    font.tickslab = c(65),
                    font.legend = c(50),
                    censor.size=40,
                    # risk.table = TRUE,
                    # risk.table.fontsize = 6,
                    # risk.table.height = 0.25,
                    ggtheme = theme_bw(base_size = 80),
                    palette = c("#815AD1","#5AC1D1"),
                    legend.labs = levels(plot.data.g$NF.status),
                    legend.title = "",
                    # title = "Age (months) of onset"
  )#, facet.by = "Groups", title = "Test")
  g_to_save <- ggsave_workaround(pp)
  ggsave(plot = g_to_save, filename =paste0("kaplanmeier_OS_g",i,"_NF1vsSP_v6.svg"), path =file.path("./images"),
         device = "svg",
         width = 18, height =18, dpi = 300)
  
  # OS NF1 group vs NF1other groups
  plot.data.g <- plot.data[plot.data$NF.status ==TRUE,]
  plot.data.g <- plot.data.g[plot.data.g$Groups!=4,]
  fit <- survfit(Surv(Overall.survival.months, Decesed)~Group.vs.others, data =plot.data.g) 
  pp <-  ggsurvplot(fit, 
                    data = plot.data.g, 
                    pval =T,
                    pval.size =35,
                    lwd = 10,
                    font.main = c(65, "plain"),
                    font.x = c(80),
                    font.y = c(80),
                    font.tickslab = c(65),
                    font.legend = c(50),
                    censor.size=40,
                    # risk.table = TRUE,
                    # risk.table.fontsize = 6,
                    # risk.table.height = 0.25,
                    ggtheme = theme_bw(base_size = 80),
                    palette = c("#815AD1","#5AC1D1"),
                    legend.labs = levels(plot.data.g$NF.status),
                    legend.title = "",
                    # title = "Age (months) of onset"
  )#, facet.by = "Groups", title = "Test")
  g_to_save <- ggsave_workaround(pp)
  ggsave(plot = g_to_save, filename =paste0("kaplanmeier_Survival_g",i,"_NF1vsSP_v6.svg"), path =file.path("./images"),
         device = "svg",
         width = 18, height =18, dpi = 300)
  #### Age of onset
  
  fit <- survfit(Surv(Age.months.onset.MPNST, Diagnosed)~Group.vs.others, data =plot.data) #crear una columna que sean diagnosed todos true, porque se censura en el momento que se diagnostica
  pp <-  ggsurvplot(fit, 
                    data = plot.data, 
                    pval =T,
                    pval.size =35,
                    lwd = 10,
                    font.main = c(65, "plain"),
                    font.x = c(80),
                    font.y = c(80),
                    font.tickslab = c(65),
                    font.legend = c(50),
                    censor.size=40,
                    # risk.table = TRUE,
                    # risk.table.fontsize = 6,
                    # risk.table.height = 0.25,
                    ggtheme = theme_bw(base_size = 80),
                    palette = c(col.g,"grey65"),
                    legend.labs = levels(plot.data$Group.vs.others),
                    legend.title = "",
                    # title = "Age (months) of onset"
                    )#, facet.by = "Groups", title = "Test")
  g_to_save <- ggsave_workaround(pp)
  ggsave(plot = g_to_save, filename =paste0("kaplanmeier_AgeOnset_g",i,"vsOthers_v6.svg"), path =file.path("./images"),
         device = "svg",
         width = 19, height =18, dpi = 300)
  
  # AO NF1vsSP per group
  plot.data.g <- plot.data[plot.data$Groups ==i,]
  fit <- survfit(Surv(Age.months.onset.MPNST, Diagnosed)~NF.status, data =plot.data.g) #crear una columna que sean diagnosed todos true, porque se censura en el momento que se diagnostica
  pp <-  ggsurvplot(fit, 
                    data = plot.data.g, 
                    pval =T,
                    pval.size =35,
                    lwd = 10,
                    font.main = c(65, "plain"),
                    font.x = c(80),
                    font.y = c(80),
                    font.tickslab = c(65),
                    font.legend = c(50),
                    censor.size=40,
                    # risk.table = TRUE,
                    # risk.table.fontsize = 6,
                    # risk.table.height = 0.25,
                    ggtheme = theme_bw(base_size = 80),
                    palette = c("#815AD1","#5AC1D1"),
                    legend.labs = levels(plot.data.g$NF.status),
                    legend.title = "",
                    # title = "Age (months) of onset"
  )#, facet.by = "Groups", title = "Test")
  g_to_save <- ggsave_workaround(pp)
  ggsave(plot = g_to_save, filename =paste0("kaplanmeier_AgeOnset_g",i,"_NF1vsSP_v6.svg"), path =file.path("./images"),
         device = "svg",
         width = 19, height =18, dpi = 300)
}




