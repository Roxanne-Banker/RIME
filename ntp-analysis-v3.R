# rm(list = ls())



save.image("~/Desktop/RIIME/RIIME/final-recompile.RData")
#save.image("~/Desktop/RIIME/RIIME/final-analyses-full-dat-backup.RData")

setwd("/Users/rockdoc/Desktop/RIIME/")

load("~/Desktop/RIIME/RIIME/final-recompile.RData")

# load packages
library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
#library(ggpubr)
library(grid)

# library(rjags)
# library(BEST)

library(brms)
library(rstan)

library(bayesplot)
library(tidybayes)
library(modelr)



set.seed(3645)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



# I. Data Import and Manipulation ----
# importing main dataset
all <- read.csv("~/Desktop/RIIME/Datasets/RIME_all_collapsed_data-R-final.csv")

an <- data.frame(all[,c(1:3,4:9)], rep("anisian",189) )
ca <- data.frame(all[,c(1:3,10:15)], rep("carnian",189) )
ba <- data.frame(all[,c(1:3,16:21)], rep("bathonian",189) )
ap <- data.frame(all[,c(1:3,22:27)], rep("aptian",189) )
names(an) <- c("guild", "functional_group","guild_no", "no_sp", "ntp", "mean_size", "N", "sd", "lc_mean", "stage")
names(ca) <- c("guild", "functional_group","guild_no", "no_sp", "ntp", "mean_size", "N", "sd", "lc_mean", "stage")
names(ba) <- c("guild", "functional_group","guild_no", "no_sp", "ntp", "mean_size", "N", "sd", "lc_mean", "stage")
names(ap) <- c("guild", "functional_group","guild_no", "no_sp", "ntp", "mean_size", "N", "sd", "lc_mean", "stage")

# Instead of using absolute richness I want to create a column of percent richness
an <- an %>% mutate(no_sp_pa = no_sp /sum(no_sp) * 100 )
ca <- ca %>% mutate(no_sp_pa = no_sp /sum(no_sp) * 100 )
ba <- ba %>% mutate(no_sp_pa = no_sp /sum(no_sp) * 100 )
ap <- ap %>% mutate(no_sp_pa = no_sp /sum(no_sp) * 100 )

data <- rbind(an, ca, ba, ap)

# going to remove a body size entry that is a clear outlier for all communities
data$mean_size[396] <- NA

# modifying a few columns
# log normalize body size
data$log_mean_size <- log10(data$mean_size)

# now center body size data around zero in case we want that later
data$body_scale <- as.numeric(scale(data$log_mean_size, center = TRUE, scale = TRUE))
data$sp_pa_scale <- as.numeric(scale(data$no_sp_pa, center = TRUE, scale = TRUE))
data$sp_scale <- as.numeric(scale(data$no_sp, center = TRUE, scale = TRUE))

# make stage an ordered factor
data$stage <- factor(data$stage,
                      levels = c("anisian", "carnian", "bathonian", "aptian"))

# remove entries without body size data
data.bs <- na.omit(data)

# keep only guilds with ntp > 2
data.bs.2 <- filter(data.bs, ntp > 2)




# II. Figure 1: General plots ----

# ntp distribution of each Stage/community

# ntp
p1 <- ggplot(data.bs.2, aes(x=stage, y=ntp, fill=stage)) +
  geom_boxplot(alpha=0.70, outlier.color="white") +
  geom_jitter(width = 0.25) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_fill_manual(values=viridis(4),
                    name = "Stage",
                    breaks = c("anisian", "carnian", "bathonian","aptian"),
                    labels = c("Anisian", "Carnian", "Bathonian","Aptian")) +
  xlab("Stage") + ylab(expression(italic(bar(ntp)) ) ) +
  scale_x_discrete(labels=c("anisian" = "Anisian", "carnian" = "Carnian",
                            "bathonian" = "Bathonian", "aptian" = "Aptian")) +
  geom_text(x=0.75, y=4.1, label="A", size = 8, fontface=1, family = "Helvetica")
p1

# mcl
p2 <- ggplot(data.bs.2, aes(x=stage, y=lc_mean, fill=stage)) +
  geom_boxplot(alpha=0.70, outlier.color="white") +
  geom_jitter(width = 0.25) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "none") +
  scale_fill_manual(values=viridis(4),
                    name = "Stage",
                    breaks = c("anisian", "carnian", "bathonian","aptian"),
                    labels = c("Anisian", "Carnian", "Bathonian","Aptian")) +
  xlab("Stage") + ylab(expression(italic(bar(mcl)) ) ) +
  scale_x_discrete(labels=c("anisian" = "Anisian", "carnian" = "Carnian",
                            "bathonian" = "Bathonian", "aptian" = "Aptian")) +
  geom_text(x=0.75, y=4.3, label="B", size = 8, fontface=1, family = "Helvetica")
p2


# log mean size
p3 <- ggplot(data.bs.2, aes(x=stage, y=log_mean_size, fill=stage)) +
  geom_boxplot(alpha=0.70, outlier.color="white") +
  geom_jitter(width = 0.25) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10)) +
  scale_fill_manual(values=viridis(4),
                    name = "Stage",
                    breaks = c("anisian", "carnian", "bathonian","aptian"),
                    labels = c("Anisian", "Carnian", "Bathonian","Aptian")) +
  xlab("Stage") + ylab(expression(italic("log"~bar(bs)) )) +
  scale_x_discrete(labels=c("anisian" = "Anisian", "carnian" = "Carnian",
                            "bathonian" = "Bathonian", "aptian" = "Aptian")) +
  geom_text(x=0.75, y=3.4, label="C", size = 8, fontface=1, family = "Helvetica")
p3

# combine into one figure
fig1a <- ggplotGrob(p1)
fig1b <- ggplotGrob(p2)
fig1c <- ggplotGrob(p3)

fig1.final <- cbind(fig1a, fig1b, fig1c, size = "first")

grid.newpage()
grid.draw(fig1.final)




# III. Figure 2: all compare ----


# ANOVA-like bayesian analysis between all stages
# comparing ntp and longest chain distributions at guild level

# trying a Kruschke-like approach with the brms package:
# https://rpruim.github.io/Kruschke-Notes/nominal-predictors.html


# using this model, which exlucdes ntp < 2
# improves fit
all_ntp_brm_2 <- brm(ntp ~ stage, data = data.bs.2)

report_table(all_ntp_brm_2)

summary(all_ntp_brm_2)
contrasts(data.bs.2$stage) # Shows the coding of the categorical variable
plot(all_ntp_brm_2)
hist(resid(all_ntp_brm_2))
pp_check(all_ntp_brm_2) # posterior parameter check. Checks if density plot for observed (y) is similar to that predicted (yrep)


### Z. overall ntp distributions in each comm ----
hyp_all_dist = hypothesis(
  all_ntp_brm_2,
  hypothesis = c(
    "Intercept = 0",
    # Carnian
    "Intercept + stagecarnian = 0",
    # Bathonian
    "Intercept + stagebathonian = 0",
    # Aptian
    "Intercept + stageaptian =0" ))

toplot = hyp_all_dist$samples # extracts the distributions for the hypotheses
# change labels to the columns which defaulted to "H1", "H2", "H3", "H4" for the four hypotheses
summary(hyp_all_dist)
names(toplot) = c("Anisian", "Carnian", "Bathonian", "Aptian")
mcmc_dens(toplot)
mcmc_areas(hyp_all_dist$samples)


### A. ntp comparison ----

# set up hypothesis tests
hyp_2.1 <- hypothesis(
  all_ntp_brm_2,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "stagecarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "stagebathonian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagebathonian = 0"
  ))
plot(hyp_2.1)
hyp_2.1

mcmc_dens(hyp_2.1$samples) 





# final ntp plot
color_scheme_set("gray")
plot2.1 <- mcmc_areas(hyp_2.1$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
  )) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = -0.6)) +
  ggtitle(expression("A:"~ bar(ntp)))

plot2.1




### B. mcl comparison ----

all_lc_brm <- brm(lc_mean ~ stage, data = data.bs.2)
summary(all_lc_brm)  
plot(all_lc_brm)
hist(resid(all_lc_brm))
pp_check(all_lc_brm) # posterior parameter check. Checks if density plot for observed (y) is similar to that predicted (yrep)
# this is okay

hyp_3.1 <- hypothesis(
  all_lc_brm,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "stagecarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "stagebathonian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagebathonian = 0"
  ))
plot(hyp_3.1)
hyp_3.1

# final mcl plot
color_scheme_set("gray")
plot3.1 <- mcmc_areas(hyp_3.1$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
  )) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = -0.6)) +
  ggtitle(expression("B:"~ bar(mcl) ))
plot3.1




### C. log bs  comparison ----

all_logbs_brm <- brm(log_mean_size ~ stage, data = data.bs.2)

summary(all_logbs_brm)  
plot(all_logbs_brm)
hist(resid(all_logbs_brm))
pp_check(all_logbs_brm)

hyp_4.1 <- hypothesis(
  all_logbs_brm,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "stagecarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "stagebathonian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagebathonian = 0"
  ))
plot(hyp_4.1)
hyp_4.1

# final log bs plot
color_scheme_set("gray")
plot4.1 <- mcmc_areas(hyp_4.1$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
  )) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = -0.6)) +
  ggtitle(expression("C:"~ italic(log) ~ bar(bs)))
plot4.1


### D. percent richness comparison ----

all_rich_brm <- brm(no_sp_pa ~ stage, data = data.bs.2)

summary(all_rich_brm)  
# plot(all_logbs_brm)
# hist(resid(all_logbs_brm))
# pp_check(all_logbs_brm)

hyp_5.1 <- hypothesis(
  all_rich_brm,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "stagecarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "stagebathonian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagebathonian = 0"
  ))
plot(hyp_5.1)
hyp_5.1

# final rich plot
color_scheme_set("gray")
plot5.1 <- mcmc_areas(hyp_5.1$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  #geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
  )) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = -0.6)) +
  ggtitle("A")
plot5.1


### E. in-degree ----

# load in data
degree.data <- read.csv("degData.csv")
# make stage an ordered factor
degree.data$Stage <- factor(degree.data$Stage,
                            levels = c("Anisian", "Carnian", "Bathonian", "Aptian"))



all_indeg_brm <- brm(inDegree ~ Stage, data = degree.data)
# summary(all_indeg_brm)  
# plot(all_logbs_brm)
# hist(resid(all_logbs_brm))
# pp_check(all_logbs_brm)

hyp_6.1 <- hypothesis(
  all_indeg_brm,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "StageCarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "StageBathonian - StageCarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "StageAptian - StageBathonian = 0"
  ))
plot(hyp_6.1)
hyp_6.1

# final log bs plot
color_scheme_set("gray")
plot6.1 <- mcmc_areas(hyp_6.1$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  #geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
  )) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = -0.6)) +
  ggtitle("B")
plot6.1



### F. out-degree ----
# 
# all_outdeg_brm <- brm(outDegree ~ Stage, data = degree.data)
# summary(all_outdeg_brm)  
# # plot(all_logbs_brm)
# # hist(resid(all_logbs_brm))
# # pp_check(all_logbs_brm)
# 
# hyp_7.1 <- hypothesis(
#   all_outdeg_brm,
#   hypothesis = c(
#     # compare Carnian to Anisian is just the stagecarnian effect
#     "StageCarnian = 0",
#     # compare Carnian to Bathonian [bathonian - carnian]
#     "StageBathonian - StageCarnian = 0",
#     # compare Bathonian to Aptian [aptian - bathonian]
#     "StageAptian - StageBathonian = 0"
#   ))
# plot(hyp_7.1)
# hyp_7.1

# grobify
fig2.1 <- ggplotGrob(plot2.1) # A
fig3.1 <- ggplotGrob(plot3.1) # B
fig4.1 <- ggplotGrob(plot4.1) # C

# fig5.1 <- ggplotGrob(plot5.1) # D
# fig6.1 <- ggplotGrob(plot6.1) # E
# fig7.1 <- ggplotGrob(plot7.1) # F

# final main text fig with ntp, mcl, bs
fig3.final <- rbind(fig2.1, fig3.1, fig4.1, size = "first")
grid.newpage()
grid.draw(fig3.final)

# supp fig final with rich, in deg
# supp.fig.final <- rbind(fig5.1, fig6.1, size = "first")
# grid.newpage()
# grid.draw(supp.fig.final)




# IV. Figure 3: Degree distributions ----

# reorder the column for this figure
degree.data$Stage <- factor(degree.data$Stage,
                            levels = c("Aptian", "Bathonian", "Carnian", "Anisian"))


# Replicating figure from Carrie's code:
# function for the figure 
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

custom_breaks <- seq(0, 75, 5)

in.deg <- ggplot(degree.data, aes(x= inDegree, fill = Stage)) +
  geom_histogram(binwidth = 5, color = "black") +
  facet_grid(Stage~.) +
  scale_fill_manual(values=rev(viridis(4)) ) +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=10),
        strip.background = element_blank(),
        strip.text.y = element_blank()) +
  ylab("Frequency") + xlab("In Degree") +
  scale_x_continuous(breaks = custom_breaks,
                     labels = every_nth(custom_breaks, 4, inverse = TRUE))
in.deg

# out.deg <- ggplot(degree.data, aes(x= outDegree, fill = Stage)) +
#   geom_histogram(binwidth = 5, color = "black") +
#   facet_grid(Stage~.) +
#   scale_fill_manual(values=rev(viridis(4)) ) +
#   theme_classic() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size=10),
#         strip.background = element_blank(),
#         strip.text.y = element_blank()) +
#   ylab("") + xlab("Out Degree")

# n1 <- ggplotGrob(in.deg)
# n2 <- ggplotGrob(out.deg)
# degree <- cbind(n1, n2, size = "first")
# grid.newpage()
# grid.draw(degree)






# VI. Figure 4-5: EPN ----

# Find in separate .R file. 




# V. Figure 6: Bayesian linear regression ----

## Classic generalized linear model ----

# # this was the best linear model I came up with from the glm package
# fin1 <- glm(ntp ~ body_scale + lc_mean + stage + 
#               body_scale:lc_mean + lc_mean:stage, 
#             family = Gamma(link = log), 
#             na.action = "na.fail",
#             data = data3)
# summary(fin1)


# okay, now going to use what I learned with Carrie and apply a bayesian approach
# I was going to use the BAS package basedon on this tutorial:
# see: https://statswithr.github.io/book/bayesian-model-choice.html#calculating-posterior-probability-in-r
# but I am going to use the brms package instead, based on use in these sources:
# https://tem11010.github.io/regression_brms/
# https://www.rensvandeschoot.com/tutorials/r-linear-regression-bayesian-using-brms/

# figures based on this code:
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html


data.2 <- filter(data, ntp >2)
# names(data.2) <- c("guild", "functional_group", "guild_no", "no_sp", "ntp", "mean_size",       
#      "N", "sd", "lc_mean", "stage", "no_sp_pa", "log_mean_size",   
#      "body_scale", "rich_pa_scale", "rich_scale")



## A. NTP ~ Stage * lc_mean ----

# the model
lc_model <- brm(ntp ~ lc_mean*stage, data = data.2)
#lc_model2 <- brm(ntp ~ lc_mean*stage, data = data.bs.2)

summary(lc_model)
posterior_summary(lc_model)




data.2 %>%
  group_by(stage) %>%
  data_grid(lc_mean = seq_range(lc_mean, n = 51)) %>%
  add_epred_draws(lc_model) %>%
  ggplot(aes(x = lc_mean, y = ntp, color = ordered(stage))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = data.2) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(face="bold", size = 10))


plot7.1 <- data.2 %>%
  group_by(stage) %>%
  modelr::data_grid(lc_mean = modelr::seq_range(lc_mean, n = 101)) %>%
  # NOTE: this shows the use of ndraws to subsample within add_epred_draws()
  # ONLY do this IF you are planning to make spaghetti plots, etc.
  # NEVER subsample to a small sample to plot intervals, densities, etc.
  tidybayes::add_epred_draws(lc_model, ndraws = 100) %>%
  ggplot(aes(x = lc_mean, y = ntp, color = ordered(stage))) +
  geom_line(aes(y = .epred, group = paste(stage, .draw)), alpha = .1) +
  geom_point(data = data.2) +
  scale_color_manual(values = viridis(4),
                     name = "Stage",
                     breaks = c("anisian", "carnian", "bathonian","aptian"),
                     labels = c("Anisian", "Carnian", "Bathonian","Aptian")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        legend.position = "none") +
  xlab(expression(italic(bar(mcl)) )) +
  ylab(expression(italic(bar(ntp)) ) )
plot7.1






## B. NTP ~ Stage * log_mean_size ----

# the model
body_model <- brm(ntp ~ mean_size*stage, data = data.2)
body_model2 <- brm(ntp ~ log_mean_size*stage, data = data.2)

summary(body_model)
#posterior_summary(body_model)

plot7.2 <- data.2 %>%
  group_by(stage) %>%
  data_grid(mean_size = seq_range(mean_size, n = 101)) %>%
  # NOTE: this shows the use of ndraws to subsample within add_epred_draws()
  # ONLY do this IF you are planning to make spaghetti plots, etc.
  # NEVER subsample to a small sample to plot intervals, densities, etc.
  add_epred_draws(body_model, ndraws = 100) %>%
  ggplot(aes(x = mean_size, y = ntp, color = ordered(stage))) +
  geom_line(aes(y = .epred, group = paste(stage, .draw)), alpha = .1) +
  geom_point(data = data.2) +
  scale_color_manual(values = viridis(4)) +
  scale_color_manual(values = viridis(4),
                     name = "Stage",
                     breaks = c("anisian", "carnian", "bathonian","aptian"),
                     labels = c("Anisian", "Carnian", "Bathonian","Aptian")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        legend.position = "none") +
  xlab("Guild mean body size (cm)") +
  ylab(expression(italic(bar(ntp)) ) )
plot7.2



# posterior distributions of slopes for different stages
hyp_body_compare = hypothesis(
  body_model,
  hypothesis = c(
    # ntp-mcl slope in the Anisian
    "lc_mean = 0",
    # ntp-mcl slope in the Carnian
    "lc_mean + lc_mean:stagecarnian = 0",
    # ntp-mcl slope in the Bathonian
    "lc_mean + lc_mean:stagebathonian = 0",
    # ntp-mcl slope in the Aptian
    "lc_mean + lc_mean:stageaptian = 0"))
hyp_body_compare




## C. NTP ~ Stage * richness ----

# the model
rich_model <- brm(ntp ~ no_sp*stage, data = data.2)
rich_model_2.1 <- brm(ntp ~ no_sp_pa*stage, data = data.2)

plot(rich_model_2.1)
summary(rich_model_2.1)
posterior_summary(rich_model)

plot7.3 <- data.2 %>%
  group_by(stage) %>%
  data_grid(no_sp_pa = seq_range(no_sp_pa, n = 101)) %>%
  # NOTE: this shows the use of ndraws to subsample within add_epred_draws()
  # ONLY do this IF you are planning to make spaghetti plots, etc.
  # NEVER subsample to a small sample to plot intervals, densities, etc.
  add_epred_draws(rich_model_2.1, ndraws = 100) %>%
  ggplot(aes(x = no_sp_pa, y = ntp, color = ordered(stage))) +
  geom_line(aes(y = .epred, group = paste(stage, .draw)), alpha = .1) +
  geom_point(data = data.2) +
  scale_color_manual(values = viridis(4),
                     name = "Stage:",
                     breaks = c("anisian", "carnian", "bathonian","aptian"),
                     labels = c("Anisian", "Carnian", "Bathonian","Aptian")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        legend.position = "bottom") +
  xlab("Guild percent richness") +
  ylab(expression(italic(bar(ntp)) ) )
plot7.3


# grobify
fig7.1 <- ggplotGrob(plot7.1) # A
fig7.2 <- ggplotGrob(plot7.2) # B
fig7.3 <- ggplotGrob(plot7.3) # C

# final main text fig with ntp, mcl, bs
fig7.final <- rbind(fig7.1, fig7.2, fig7.3, size = "first")
grid.newpage()
grid.draw(fig7.final)






# VII. Figure 7: Changing NTP  ----


# for guilds that persisted all four stages using the epn dataset
all.persist <- all %>% filter( An.C == "persist" & C.B == "persist" & B.Ap == "persist")
# now just the guild names
all.persist.guilds <- all.persist$guild

# now pull out the guilds from that list from the main dataset
all.persist.dat <- filter(data, guild %in% all.persist.guilds)
all.persist.dat.pp <- filter(all.persist.dat, ntp <= 2) # ntps <= 2
all.persist.dat.2 <- filter(all.persist.dat, ntp > 2) # ntps > 2

all_persist_ntp_brm <- brm(ntp ~ stage, data = all.persist.dat.2)

summary(all_persist_ntp_brm)  
# plot(all_logbs_brm)
# hist(resid(all_logbs_brm))
# pp_check(all_logbs_brm)

hyp_persist <- hypothesis(
  all_persist_ntp_brm,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "stagecarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "stagebathonian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagebathonian = 0"
  ))
plot(hyp_persist)
hyp_persist


color_scheme_set("blue")
plot8.3 <- mcmc_areas(hyp_persist$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  #geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
    )) +
  theme(axis.text = element_text(size = 15) ) 
plot8.3


# Changing NTP figure 
ntp.change <- all[,c(1:3, 5, 11, 17, 23, 28:30)]
ntp.all.change <- ntp.change %>% filter(An.C == "persist" & C.B == "persist" & B.Ap == "persist") %>% filter(Ap_ntp > 2)
# just used the Aptian ntps to filter to ntp >2 for all guilds

# ntp differences
var1 <- data.frame(ntp.all.change$guild, ntp.all.change$C_ntp - ntp.all.change$An_ntp, rep(1, 37))
names(var1) <- c("guild","diff","txn")
var2 <- data.frame(ntp.all.change$guild, ntp.all.change$B_ntp - ntp.all.change$C_ntp, rep(2, 37))
names(var2) <- c("guild", "diff","txn")
var3 <- data.frame(ntp.all.change$guild, ntp.all.change$Ap_ntp - ntp.all.change$B_ntp, rep(3, 37))
names(var3) <- c("guild", "diff","txn")

ntp.diffs <- rbind(var1,var2,var3)
ntp.diffs$txn <- as.factor(ntp.diffs$txn)


# ntp diff figure
plot8.2 <- ggplot(ntp.diffs, aes(x=txn, y=diff, color = diff, group=guild)) +
  geom_hline(yintercept = 0) +
  geom_path(color = "gray") +
  geom_point() +
  xlab("Comparison") + ylab(expression(Delta(italic(bar(ntp)) ) ) ) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_color_continuous(name = expression(Delta(italic(bar(ntp)) ) )) +
  scale_x_discrete(breaks=c("1", "2", "3"), 
                   labels=c("Anisian to \nCarnian", "Carnian to \nBathonian", "Bathonian to \nAptian")) 
plot8.2


# ntps
v1 <- data.frame(ntp.all.change$guild, ntp.all.change$An_ntp, rep("anisian", 37))
names(v1) <- c("guild","ntp", "stage")
v2 <- data.frame(ntp.all.change$guild, ntp.all.change$C_ntp, rep("carnian", 37))
names(v2) <- c("guild","ntp", "stage")
v3 <- data.frame(ntp.all.change$guild, ntp.all.change$B_ntp, rep("bathonian", 37))
names(v3) <- c("guild","ntp", "stage")
v4 <- data.frame(ntp.all.change$guild, ntp.all.change$Ap_ntp, rep("aptian", 37))
names(v4) <- c("guild","ntp", "stage")

ntp.change.fig.dat <- rbind(v1,v2,v3,v4)
ntp.change.fig.dat$stage <- factor(ntp.change.fig.dat$stage, levels = c("anisian", "carnian", "bathonian", "aptian"))


# ntp fig
plot8.1 <- ggplot(ntp.change.fig.dat, aes(x=stage, y=ntp, group=guild, color = ntp)) +
  geom_path(color = "gray") +
  geom_point() +
  xlab("Stage") + ylab(expression(italic(bar(ntp)) )) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_color_continuous(name = expression(italic(bar(ntp)) ) ) +
  scale_x_discrete(breaks=c("anisian", "carnian", "bathonian", "aptian"), 
                   labels=c("Anisian", "Carnian", "Bathonian", "Aptian")) 
plot8.1


# grobify
fig8.1 <- ggplotGrob(plot8.1) # A
fig8.2 <- ggplotGrob(plot8.2) # B
fig8.3 <- ggplotGrob(plot8.3) # C

# final main text fig with ntp, mcl, bs
fig8.final <- rbind(fig8.1, fig8.2, size = "first")
grid.newpage()
grid.draw(fig8.final)


# VIII. Fig 8: Compare Consumer Size ----

pred.data <- read.csv("Datasets/RIME_predator_data_R-v1.csv")
pred.size <- pred.data[,c(1:3,6,9,12,15)]
names(pred.size)[4:7] <- c("anisian", "carnian",  "bathonian",  "aptian")

pred.size.long <- pivot_longer(pred.size, cols = c(4:7), names_to ="stage",
             values_to = "mean_size")
pred.size.long.na <- na.omit(pred.size.long)

pred_size_brm <- brm(mean_size ~ stage, data = pred.size.long.na)

summary(pred_size_brm)
pred_size_brm


hyp_pred_diffs <- hypothesis(
  pred_size_brm,
  hypothesis = c(
    # compare Carnian to Anisian is just the stagecarnian effect
    "stagecarnian = 0",
    # compare Carnian to Bathonian [bathonian - carnian]
    "stagebathonian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagebathonian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian - stagecarnian = 0",
    # compare Bathonian to Aptian [aptian - bathonian]
    "stageaptian = 0"
  ))
plot(hyp_pred_diffs)
hyp_pred_diffs

toplot1 = hyp_pred_diffs$samples # extracts the distributions for the hypotheses
# change labels to the columns which defaulted to "H1", "H2", "H3", "H4" for the four hypotheses
summary(hyp_pred_diffs)


# USE THIS FIGURE
color_scheme_set("gray")
plot_pred <- mcmc_areas(hyp_pred_diffs$samples, prob=.95, pars=c("H1", "H2", "H3")) +
  geom_vline(xintercept = 0, color = "darkgray") +
  scale_y_discrete(labels = c(
    "H1" = expression(paste(mu[Carnian],"-",mu[Anisian])),
    "H2" = expression(paste(mu[Bathonian],"-",mu[Carnian])),
    "H3" = expression(paste(mu[Aptian],"-",mu[Bathonian]))
  )) +
  theme(axis.text = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = -0.6)) 
plot_pred



hyp_pred = hypothesis(
  pred_size_brm,
  hypothesis = c(
    # Aptian
    "Intercept + stageaptian = 0",
    # Bathonian
    "Intercept + stagebathonian = 0",
    # Carnian
    "Intercept + stagecarnian = 0",
    # Anisian
    "Intercept = 0" )) 

toplot = hyp_pred$samples # extracts the distributions for the hypotheses
# change labels to the columns which defaulted to "H1", "H2", "H3", "H4" for the four hypotheses
summary(hyp_pred)
names(toplot) = c("Aptian", "Bathonian", "Carnian", "Anisian")
mcmc_dens(toplot)

# USE THIS FIGURE
mcmc_areas(toplot, prob=.95)





# X. Caribbean Comparison ----

jamaica.full <- read.csv("~/Desktop/RIIME/Jamaica/jamaica_modern_data.csv")

names(jamaica.full)

guilds <- unique(jamaica.full$guild)
Dataset <- jamaica.full

ntp.mean <- NULL
lc.mean <- NULL
group.names <- NULL

# loop returns mean, sd, n, of ntp and group name of each functional group
for (k in guilds) {
  sub <- filter(Dataset, guild == k)
  ntp.mean <- c(mean(sub$ntp), ntp.mean)
  lc.mean <- c(mean(sub$long_chain), lc.mean)
  group.names <- c(sub$guild[1], group.names)
}

jamaica.reduced <- data.frame(group.names, ntp.mean, lc.mean)

jamaica.2 <- filter(jamaica.reduced, ntp.mean > 2)
jamaica.2$stage <- rep("jamaica", nrow(jamaica.2))
names(jamaica.2) <- c("guild", "ntp", "lc_mean", "stage")


# combine modern jamaica data with the Mesozoic so that plotting is easier

combined.1 <- data.2 %>% select(guild, ntp, lc_mean, stage)
combined.2 <- rbind(jamaica.2, combined.1)

# set factors in stage column
combined.2$stage <- factor(combined.2$stage,
                     levels = c("anisian", "carnian", "bathonian", "aptian", "jamaica"))


lapply(combined.2, class)


# the model
jamaica_mesozoic_model <- brm(ntp ~ lc_mean*stage, data = combined.2)

summary(jamaica_mesozoic_model)
posterior_summary(jamaica_mesozoic_model)



# posterior distributions of slopes for different stages
hyp_modern_compare = hypothesis(
  jamaica_mesozoic_model,
  hypothesis = c(
    # ntp-mcl slope in the Anisian
    "lc_mean = 0",
    # ntp-mcl slope in the Carnian
    "lc_mean + lc_mean:stagecarnian = 0",
    # ntp-mcl slope in the Bathonian
    "lc_mean + lc_mean:stagebathonian = 0",
    # ntp-mcl slope in the Aptian
    "lc_mean + lc_mean:stageaptian = 0",
    # ntp-mcl slope in the modern
    "lc_mean + lc_mean:stagejamaica = 0"))
hyp_modern_compare




combined.2 %>%
  group_by(stage) %>%
  data_grid(lc_mean = seq_range(lc_mean, n = 51)) %>%
  add_epred_draws(jamaica_mesozoic_model) %>%
  ggplot(aes(x = lc_mean, y = ntp, color = ordered(stage))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = combined.2) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(face="bold", size = 10))


plot.modern <- combined.2 %>%
  group_by(stage) %>%
  modelr::data_grid(lc_mean = modelr::seq_range(lc_mean, n = 101)) %>%
  # NOTE: this shows the use of ndraws to subsample within add_epred_draws()
  # ONLY do this IF you are planning to make spaghetti plots, etc.
  # NEVER subsample to a small sample to plot intervals, densities, etc.
  tidybayes::add_epred_draws(jamaica_mesozoic_model, ndraws = 100) %>%
  ggplot(aes(x = lc_mean, y = ntp, color = ordered(stage))) +
  geom_line(aes(y = .epred, group = paste(stage, .draw)), alpha = .1) +
  geom_point(data = combined.2) +
  scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF", "#F594DC"),
                     name = "Stage",
                     breaks = c("anisian", "carnian", "bathonian","aptian","jamaica"),
                     labels = c("Anisian", "Carnian", "Bathonian","Aptian","Modern")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(face="bold", size = 10)) +
  xlab(expression(italic(bar(mcl)) )) +
  ylab(expression(italic(bar(ntp)) ) )
plot.modern




