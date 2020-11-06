###############################################################################
# Analysis code for "Metabolic depression in sea urchin barrens associated    #
# with food deprivation."                                                     #
#                                                                             #
# Authors: N_B_ Spindel & D_K_ Okamoto                                        #
###############################################################################

# Load packages ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(Hmisc)
library(car)
library(cowplot)
library(sjPlot)
library(xtable)

# Load functions --------------------------------------------------------------
source("./Functions.R")

# Get calculated data ---------------------------------------------------------
dat <- read_excel("./Derived_Data/Respirometry_Data_M_franciscanus.xlsx", sheet = 1)
dat$date <- as.character(dat$date)
dat <- dat %>% 
  mutate(
    date = if_else(date %in% c("2019-08-11", "2019-08-12"), "2019-08-11", date)
  ) 
dat$date <- as.Date(dat$date)
dat <- dat %>% 
  mutate(
    total_gonad_approx_g = gonad_wet_mass_g*5,
    volumetric_rate = I(ifelse(volumetric_rate>0,0.00001,-volumetric_rate)),
    habitat = as.factor(habitat),
    species = as.factor(species)
  ) %>%
  mutate(
    site = if_else(site %in% c("Faraday N", "Faraday South"), "Faraday", site)
  ) %>%
  mutate(
    site = if_else(site %in% c("Murchison NorthWest", "Murchison NE"), "Murchison", site)
  ) %>%
  mutate(
    site = as.factor(site)
  ) %>%
  filter(
    !is.na(total_gonad_approx_g) 
  )

# Wrangle wild individual data: -----------------------------------------------
# Data frame of wild M. franciscanus with wet weight measurements
dat_filtered <- dat %>%
  filter(
    species == "Mesocentrotus franciscanus" &
    is_background == "N" &
    is_trial_subject == "N" &
    !is.na(animal_wet_mass_g)
  ) 

# Data frame of wild M. franciscanus with AFDM measurements
dat_filtered_AFDM <- dat_filtered %>%
  filter(
      !is.na(whole.afdm.g.approx)
  ) %>%
  mutate(
    habitat = as.factor(habitat),
    site = as.factor(site),
    chamber_id = as.factor(chamber_id)
  )

# Data frame of wild M. franciscanus with wet weight measurements from Surge Narrows
dat_filtered_SN <- dat_filtered %>%
  filter(site == "Surge Narrows")

# RMR Analysis (Total Wet Mass): effects of site and habitat  --------------------
# Model selection: Fit several reasonable GLMMs for fixed effects on volumetric
# oxygen consumption rate using Restricted Maximum Likelihood. 
f1 <- glmmTMB(volumetric_rate ~ log(total_gonad_approx_g) + (1 | date/chamber_id), 
              data = dat_filtered,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(total_gonad_approx_g) + habitat + (1 | date/chamber_id))
f3 <- update(f1, formula = . ~ log(total_gonad_approx_g) + habitat + site + (1 | date/chamber_id))
f4 <- update(f1, formula = . ~ log(total_gonad_approx_g) * habitat + site + (1 | date/chamber_id))
f5 <- update(f1, formula = . ~ log(total_gonad_approx_g) + habitat * site + (1 | date/chamber_id))
f6 <- update(f1, formula = . ~ log(total_gonad_approx_g) * habitat * site + (1 | date/chamber_id))
f7 <- update(f1, formula = . ~ log(total_gonad_approx_g) : habitat + site + (1 | date/chamber_id))
f8 <- update(f1, formula = . ~ log(total_gonad_approx_g) : habitat : site + (1 | date/chamber_id))
f9 <- update(f1, formula = . ~ log(total_gonad_approx_g) : habitat * site + (1 | date/chamber_id))

# f1 <- glmmTMB(volumetric_rate ~ log(animal_wet_mass_g) + (1 | date/chamber_id), 
#               data = dat_filtered,
#               family = gaussian(link = "log"),
#               REML = TRUE)
# f2 <- update(f1, formula = . ~ log(animal_wet_mass_g) + habitat + (1 | date/chamber_id))
# f3 <- update(f1, formula = . ~ log(animal_wet_mass_g) + habitat + site + (1 | date/chamber_id))
# f4 <- update(f1, formula = . ~ log(animal_wet_mass_g) * habitat + site + (1 | date/chamber_id))
# f5 <- update(f1, formula = . ~ log(animal_wet_mass_g) + habitat * site + (1 | date/chamber_id))
# f6 <- update(f1, formula = . ~ log(animal_wet_mass_g) * habitat * site + (1 | date/chamber_id))
# f7 <- update(f1, formula = . ~ log(animal_wet_mass_g) : habitat + site + (1 | date/chamber_id))
# f8 <- update(f1, formula = . ~ log(animal_wet_mass_g) : habitat : site + (1 | date/chamber_id))
# f9 <- update(f1, formula = . ~ log(animal_wet_mass_g) : habitat * site + (1 | date/chamber_id))
# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9)

# Although f3 has the lowest AICc, we select f6 to examine habitat:body size  and habitat:site interactions. Re-fit selected model using 
# Maximum Likelihood for inference.
f6 <- update(f6, REML = FALSE)
fit_glmmTMB_vo2_wetmass_site <- f6

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_vo2_wetmass_site, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_vo2_wetmass_site_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_vo2_wetmass_site,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_vo2_wetmass_site_anova
anova_table <- autoformat(xtable(red_vo2_wetmass_site_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/VO2_wetmass_habitat_site_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f6, file = "./Model_Summaries/VO2_wetmass_habitat_site")

# Visualize predictions and data 
fixed_predictions_wetmass_site <- predict(object = fit_glmmTMB_vo2_wetmass_site, 
                                      newdata = dat_filtered,
                                      re.form = ~ 0,
                                      se.fit = TRUE)
dat_filtered$pred_fixed_wetmass_site <- exp(fixed_predictions_wetmass_site$fit)
dat_filtered$pred_fixed_se_wetmass_site_lci <- exp(fixed_predictions_wetmass_site$fit-
                                                 fixed_predictions_wetmass_site$se.fit)
dat_filtered$pred_fixed_se_wetmass_site_uci <- exp(fixed_predictions_wetmass_site$fit+
                                                 fixed_predictions_wetmass_site$se.fit)

p_vo2_wetmass_site_red <- ggplot(aes(x = total_gonad_approx_g, 
                                 
                                 # p_vo2_wetmass_site_red <- ggplot(aes(x = spheroid.volume.ml, 
                                 y = volumetric_rate), data = dat_filtered)+
  # geom_point(aes(color = habitat, shape = habitat, size = total_gonad_approx_g), alpha = 0.5)+
  geom_point(aes(color = habitat, shape = habitat, size = spheroid.volume.ml),alpha = 0.5)+
  geom_line(aes(y = pred_fixed_wetmass_site, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_wetmass_site_lci, 
                  ymax = pred_fixed_se_wetmass_site_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  # scale_size_discrete(name = "gonad mass (g)")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Gonadal Wet Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = c(0.15,0.7));

p_vo2_wetmass_site_red

ggsave("./Figures/test.pdf", plot = p_vo2_wetmass_site_red, units = "in", width = 9, height = 5)
ggsave("./Figures/test.jpg", plot = p_vo2_wetmass_site_red, units = "in", width = 9, height = 5)

# RMR Analysis (Test Volume): effects of site and habitat  --------------------
# Model selection: Fit several reasonable GLMMs for fixed effects on volumetric
# oxygen consumption rate using Restricted Maximum Likelihood. 
f1 <- glmmTMB(volumetric_rate ~ log(spheroid.volume.ml) + (1 | date/chamber_id), 
              data = dat_filtered,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + (1 | date/chamber_id))
f3 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + site + (1 | date/chamber_id))
f4 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + site + (1 | date/chamber_id))
f5 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * site + (1 | date/chamber_id))
f6 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat * site + (1 | date/chamber_id))
f7 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat + site + (1 | date/chamber_id))
f8 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat : site + (1 | date/chamber_id))
f9 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat * site + (1 | date/chamber_id))

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9)

# Although f3 has the lowest AICc, we select f6 to examine habitat:body size  and habitat:site interactions. Re-fit selected model using 
# Maximum Likelihood for inference.
f6 <- update(f6, REML = FALSE)
fit_glmmTMB_vo2_vol_site <- f6

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_vo2_vol_site, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_vo2_vol_site_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_vo2_vol_site,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_vo2_vol_site_anova
anova_table <- autoformat(xtable(red_vo2_vol_site_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/VO2_volume_habitat_site_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f6, file = "./Model_Summaries/VO2_volume_habitat_site")

# Visualize predictions and data 
fixed_predictions_vol_site <- predict(object = fit_glmmTMB_vo2_vol_site, 
                                 newdata = dat_filtered,
                                 re.form = ~ 0,
                                 se.fit = TRUE)
dat_filtered$pred_fixed_vol_site <- exp(fixed_predictions_vol_site$fit)
dat_filtered$pred_fixed_se_vol_site_lci <- exp(fixed_predictions_vol_site$fit-
                                                          fixed_predictions_vol_site$se.fit)
dat_filtered$pred_fixed_se_vol_site_uci <- exp(fixed_predictions_vol_site$fit+
                                                          fixed_predictions_vol_site$se.fit)

p_vo2_vol_site_red <- ggplot(aes(x = spheroid.volume.ml,
                            y = volumetric_rate), data = dat_filtered)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_vol_site, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_vol_site_lci, 
                  ymax = pred_fixed_se_vol_site_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Test Volume (mL)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = c(0.1,0.8));

p_vo2_vol_site_red

# RMR Analysis (Test Volume): effect of habitat at Surge Narrows --------------
# Model selection: Fit several reasonable GLMMs for fixed effects on volumetric
# oxygen consumption rate using Restricted Maximum Likelihood. 
f1 <- glmmTMB(volumetric_rate ~ log(spheroid.volume.ml) + (1 | date/chamber_id), 
              data = dat_filtered_SN,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + (1 | date/chamber_id))
f3 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + (1 | date/chamber_id))
f4 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat + (1 | date/chamber_id))

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4)

# Although f2 has the lowest AICc, we select f3 to examine the habitat:body size interaction. Re-fit selected model using 
# Maximum Likelihood for inference.
f3 <- update(f3, REML = FALSE)
fit_glmmTMB_vo2_vol_SN <- f3

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_vo2_vol_SN, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type II Wald chisquare tests)
red_vo2_vol_SN_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_vo2_vol_SN,
  type = c("II"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_vo2_vol_SN_anova
anova_table <- autoformat(xtable(red_vo2_vol_SN_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/VO2_volume_habitat_SN_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f3, file = "./Model_Summaries/VO2_volume_habitat_SN")

# Visualize predictions and data 
fixed_predictions_vol_SN <- predict(object = fit_glmmTMB_vo2_vol_SN, 
                                      newdata = dat_filtered_SN,
                                      re.form = ~ 0,
                                      se.fit = TRUE)
dat_filtered_SN$pred_fixed_vol_SN <- exp(fixed_predictions_vol_SN$fit)
dat_filtered_SN$pred_fixed_se_vol_SN_lci <- exp(fixed_predictions_vol_SN$fit-
                                                 fixed_predictions_vol_SN$se.fit)
dat_filtered_SN$pred_fixed_se_vol_SN_uci <- exp(fixed_predictions_vol_SN$fit+
                                                 fixed_predictions_vol_SN$se.fit)

p_vo2_vol_SN_red <- ggplot(aes(x = spheroid.volume.ml, 
                                 y = volumetric_rate), data = dat_filtered_SN)+
  geom_point(aes(color = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_vol_SN, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_vol_SN_lci, 
                  ymax = pred_fixed_se_vol_SN_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  xlab("Test Volume (mL)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = c(0.2,0.8));

p_vo2_vol_SN_red

# RMR Analysis (Total Ash-Free Dry Mass): effect of habitat  ------------------------
# Model selection: Fit several reasonable GLMMs for fixed effects on volumetric
# oxygen consumption rate using Restricted Maximum Likelihood. 
f1 <- glmmTMB(volumetric_rate ~ log(whole.afdm.g.approx) + (1 | date/chamber_id),
              dispformula = ~ log(whole.afdm.g.approx)*habitat, # account for heteroskedasticity
              data = dat_filtered_AFDM,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(whole.afdm.g.approx) + habitat + (1 | date/chamber_id))
f3 <- update(f1, formula = . ~ log(whole.afdm.g.approx) * habitat + (1 | date/chamber_id))
f4 <- update(f1, formula = . ~ log(whole.afdm.g.approx) : habitat + (1 | date/chamber_id))

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4)

# Although f2 has the lowest AICc, we select f3 to examine habitat:body size interaction. Re-fit selected model using 
# Maximum Likelihood for inference.
f3 <- update(f3, REML = FALSE)
fit_glmmTMB_vo2_afdm <- f3

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_vo2_afdm, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_vo2_afdm_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_vo2_afdm,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_vo2_afdm_anova
anova_table <- autoformat(xtable(red_vo2_afdm_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/VO2_AFDM_habitat_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f3, file = "./Model_Summaries/VO2_AFDM_habitat")

# Visualize predictions and data 
fixed_predictions_afdm <- predict(object = fit_glmmTMB_vo2_afdm, 
                                      newdata = dat_filtered_AFDM,
                                      re.form = ~ 0,
                                      se.fit = TRUE)
dat_filtered_AFDM$pred_fixed_afdm <- exp(fixed_predictions_afdm$fit)
dat_filtered_AFDM$pred_fixed_se_afdm_lci <- exp(fixed_predictions_afdm$fit-
                                                 fixed_predictions_afdm$se.fit)
dat_filtered_AFDM$pred_fixed_se_afdm_uci <- exp(fixed_predictions_afdm$fit+
                                                 fixed_predictions_afdm$se.fit)

p_vo2_afdm <- ggplot(aes(x = whole.afdm.g.approx, 
                                 y = volumetric_rate), data = dat_filtered_AFDM)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_afdm, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_afdm_lci, 
                  ymax = pred_fixed_se_afdm_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  xlab("Total Ash-Free Dry Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = "");

p_vo2_afdm

# RMR Analysis (Gonadal Ash-Free Dry Mass): effect of habitat  ----------------
# Model selection: Fit several reasonable GLMMs for fixed effects on volumetric
# oxygen consumption rate using Restricted Maximum Likelihood. 
f1 <- glmmTMB(volumetric_rate ~ log(gonad.afdm.g.total.approx) + (1 | date/chamber_id), 
              data = dat_filtered_AFDM,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(gonad.afdm.g.total.approx) + habitat + (1 | date/chamber_id))
f3 <- update(f1, formula = . ~ log(gonad.afdm.g.total.approx) * habitat + (1 | date/chamber_id))
f4 <- update(f1, formula = . ~ log(gonad.afdm.g.total.approx) : habitat + (1 | date/chamber_id))

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4)

# Select model f3 based on having lowest AICc. Re-fit selected model using 
# Maximum Likelihood for inference.
f3 <- update(f3, REML = FALSE)
fit_glmmTMB_vo2_afdm_gonad <- f3

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_vo2_afdm_gonad, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_vo2_afdm_gonad_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_vo2_afdm_gonad,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_vo2_afdm_gonad_anova
anova_table <- autoformat(xtable(red_vo2_afdm_gonad_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/VO2_gonadAFDM_habitat_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f3, file = "./Model_Summaries/VO2_gonadAFDM_habitat")

# Visualize predictions and data 
fixed_predictions_afdm_gonad <- predict(object = fit_glmmTMB_vo2_afdm_gonad, 
                                  newdata = dat_filtered_AFDM,
                                  re.form = ~ 0,
                                  se.fit = TRUE)
dat_filtered_AFDM$pred_fixed_afdm_gonad <- exp(fixed_predictions_afdm_gonad$fit)
dat_filtered_AFDM$pred_fixed_se_afdm_gonad_lci <- exp(fixed_predictions_afdm_gonad$fit-
                                                  fixed_predictions_afdm_gonad$se.fit)
dat_filtered_AFDM$pred_fixed_se_afdm_gonad_uci <- exp(fixed_predictions_afdm_gonad$fit+
                                                  fixed_predictions_afdm_gonad$se.fit)

p_vo2_afdm_gonad <- ggplot(aes(x = gonad.afdm.g.total.approx, 
                         y = volumetric_rate), data = dat_filtered_AFDM)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_afdm_gonad, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_afdm_gonad_lci, 
                  ymax = pred_fixed_se_afdm_gonad_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  xlab("Gonadal Ash-Free Dry Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = "");

p_vo2_afdm_gonad

# Combine RMR regression visualizations ---------------------------------------
# 3-panel figure with RMR regressed against test volume, total AFDM, and 
# gonadal AFDM by habitat for urchins from Surge Narrows
RMR_multiplot_SN <- cowplot::plot_grid(p_vo2_afdm, p_vo2_afdm_gonad,
                                     ncol = 2,
                                     labels = c('A', 'B'), 
                                     label_size = 25)
RMR_multiplot_SN

# 3-panel figure with RMR regressed against test volume by habitat and site for 
# urchins from Faraday Island, Murchison Island, and Surge Narrows
RMR_multiplot_hab_x_site <- cowplot::plot_grid(RMR_multiplot_SN, p_vo2_vol_site_red,
                              ncol = 1,
                              # labels = c('', 'C'), 
                              label_size = 25)
RMR_multiplot_hab_x_site

ggsave("RMR_multiplot.pdf", path = "./Figures", plot = RMR_multiplot_hab_x_site, 
       units = "in", width = 12, height = 10)
ggsave("RMR_multiplot.jpg", path = "./Figures", plot = RMR_multiplot_hab_x_site, 
       units = "in", width = 12, height = 10)

# Gonad Content Analysis: effects of site and habitat -------------------------
# Model selection: Fit several reasonable GLMMs for fixed effects on gonad
# content using Restricted Maximum Likelihood. 
f1 <- glmmTMB(formula = total_gonad_approx_g ~ log(spheroid.volume.ml),
              dispformula = ~ habitat * site, # account for heteroskedasticity
              data = dat_filtered,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat)
f3 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + site)
f4 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + site)
f5 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * site)
f6 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat * site)
f7 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat + site)
f8 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat : site)
f9 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat * site)

# f1 <- glmmTMB(formula = gonad_wet_mass_g ~ log(test_diam_mm),
#               dispformula = ~ site * habitat, # account for heteroskedasticity 
#               data = dat_filtered,
#               family = gaussian(link = "log"),
#               REML = TRUE)
# f2 <- update(f1, formula = . ~ log(test_diam_mm) + habitat)
# f3 <- update(f1, formula = . ~ log(test_diam_mm) + habitat + site)
# f4 <- update(f1, formula = . ~ log(test_diam_mm) * habitat + site)
# f5 <- update(f1, formula = . ~ log(test_diam_mm) + habitat * site)
# f6 <- update(f1, formula = . ~ log(test_diam_mm) * habitat * site)
# f7 <- update(f1, formula = . ~ log(test_diam_mm) : habitat + site)
# f8 <- update(f1, formula = . ~ log(test_diam_mm) : habitat : site)
# f9 <- update(f1, formula = . ~ log(test_diam_mm) : habitat * site)

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9)

# Select model f4 based on having lowest AICc. Re-fit selected model using 
# Maximum Likelihood for inference.
f4 <- update(f4, REML = FALSE)
fit_glmmTMB_gonad_vol_site <- f4

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_gonad_vol_site, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_gonad_vol_site_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_gonad_vol_site,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_gonad_vol_site_anova
anova_table <- autoformat(xtable(red_gonad_vol_site_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/gonadwetmass_volume_habitat_site_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f4, file = "./Model_Summaries/gonadwetmass_volume_habitat_site")

# Visualize predictions and data 
fixed_predictions_gonad_vol_site <- predict(object = fit_glmmTMB_gonad_vol_site, 
                                      newdata = dat_filtered,
                                      re.form = ~ 0,
                                      se.fit = TRUE)
dat_filtered$pred_fixed_gonad_vol_site <- exp(fixed_predictions_gonad_vol_site$fit)
dat_filtered$pred_fixed_se_gonad_vol_site_lci <- exp(fixed_predictions_gonad_vol_site$fit-
                                                 fixed_predictions_gonad_vol_site$se.fit)
dat_filtered$pred_fixed_se_gonad_vol_site_uci <- exp(fixed_predictions_gonad_vol_site$fit+
                                                 fixed_predictions_gonad_vol_site$se.fit)

p_gonad_vol_site_red <- ggplot(aes(x = test_diam_mm, 
                                 y = total_gonad_approx_g), data = dat_filtered)+
  # geom_point(aes(color = habitat, shape = habitat), size = 3,  alpha = 0.5)+
  geom_point(aes(color = habitat, shape = habitat),  alpha = 0.5)+

  geom_line(aes(y = pred_fixed_gonad_vol_site, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_gonad_vol_site_lci,
                  ymax = pred_fixed_se_gonad_vol_site_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Test Diameter (mm)")+
  ylab("Gonadal Wet Mass (g)")+
  ggtitle("Natural History Note Biometrics")+
  gg_options()+
  theme(legend.position = c(0.15,0.7));

p_gonad_vol_site_red

ggsave("./Figures/test.pdf", plot = p_gonad_vol_site_red, units = "in", width = 9, height = 5)

# Gonad Content Analysis: effect of habitat at Surge Narrows ------------------
# Model selection: Fit several reasonable GLMMs for fixed effects on gonad
# content using Restricted Maximum Likelihood. 
f1 <- glmmTMB(formula = total_gonad_approx_g ~ log(spheroid.volume.ml),
              dispformula = ~ habitat*log(spheroid.volume.ml), # account for heteroskedasticity
              data = dat_filtered_SN,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat)
f3 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat)
f4 <- update(f1, formula = . ~ log(spheroid.volume.ml) : habitat)

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4)

# Select model f3 based on having lowest AICc. Re-fit selected model using 
# Maximum Likelihood for inference.
f3 <- update(f3, REML = FALSE)
fit_glmmTMB_gonad_vol_SN <- f3

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB_gonad_vol_SN, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_gonad_vol_SN_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB_gonad_vol_SN,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_gonad_vol_SN_anova
anova_table <- autoformat(xtable(red_gonad_vol_SN_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/gonadwetmass_volume_habitat_SN_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f3, file = "./Model_Summaries/gonadwetmass_volume_habitat_SN")

# Visualize predictions and data 
fixed_predictions_gonad_vol_SN <- predict(object = fit_glmmTMB_gonad_vol_SN, 
                                            newdata = dat_filtered_SN,
                                            re.form = ~ 0,
                                            se.fit = TRUE)
dat_filtered_SN$pred_fixed_gonad_vol_SN <- exp(fixed_predictions_gonad_vol_SN$fit)
dat_filtered_SN$pred_fixed_se_gonad_vol_SN_lci <- exp(fixed_predictions_gonad_vol_SN$fit-
                                                       fixed_predictions_gonad_vol_SN$se.fit)
dat_filtered_SN$pred_fixed_se_gonad_vol_SN_uci <- exp(fixed_predictions_gonad_vol_SN$fit+
                                                       fixed_predictions_gonad_vol_SN$se.fit)

p_gonad_vol_SN_red <- ggplot(aes(x = spheroid.volume.ml, 
                                   y = total_gonad_approx_g), data = dat_filtered_SN)+
  geom_point(aes(color = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_gonad_vol_SN, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_gonad_vol_SN_lci, 
                  ymax = pred_fixed_se_gonad_vol_SN_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  xlab("Test Volume (mL)")+
  ylab("Gonadal Wet Mass (g)")+
  gg_options()+
  theme(legend.position = "");

p_gonad_vol_SN_red

# Gonad Content Analysis: Haida Gwaii -----------------------------------------
dat <- read_excel("./Derived_Data/Urchin_Biometrics.xlsx", sheet = 1)
dat$date <- as.Date(dat$date)
dat$taxon <- paste(dat$genus, dat$species, sep = " ")
dat <- dat %>%
  mutate(
    semi.major.axis.mm = 0.5 * test_diam_mm,
    semi.minor.axis.mm = 0.5 * test_ht_mm
  ) %>%
  mutate(
    spheroid.volume.ml = (4/3 * pi * semi.major.axis.mm ^2 * semi.minor.axis.mm)/1000, # Equation for volume of oblate spheroid = 4/3pia^2c in mL
    site = if_else(site %in% c("Faraday North", "Faraday South"), "Faraday", site)
  ) %>%
  mutate(
    site = if_else(site %in% c("Murchison NorthWest", "Murchison NorthEast"), "Murchison", site)
  ) %>%
  mutate(
    site = as.factor(site)
  ) %>%
  mutate(
    year_habitat = paste(year, habitat, sep = " ")
  )

dat_filtered <- dat %>%
  filter(
    taxon == "Mesocentrotus franciscanus"
    & urch_no %nin% c(137)
    )

dat_2018 <- dat_filtered %>% filter(year == 2018)

dat_2019 <- dat_filtered %>% filter(year == 2019)

dat_murchison <- dat_filtered %>% filter(site == "Murchison")
dat_faraday <- dat_filtered %>% filter(site == "Faraday")

# Model selection: Fit several reasonable GLMMs for fixed effects on gonad
# content using Restricted Maximum Likelihood. 
f1 <- glmmTMB(formula = total_gonad_wet_wt_g ~ log(spheroid.volume.ml),
              # dispformula = ~ site * habitat, # account for heteroskedasticity
              dispformula = ~habitat,
              data = dat_filtered,
              # family = gaussian(link = "log"),
              # start = 0,
              family = tweedie(),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat)
f3 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat)
f4 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + year)
f5 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + year)
f6 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * year)
f7 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + year + site)
f8 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + year + site)
f9 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * year + site)
f10 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + year * site)
f11 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * year * site)
f12 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + year * site)
f13 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat * year * site)


AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13)
# AICc(f1, f2, f3)

f4 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + site)
f5 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + site)
f6 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * site)
f7 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat * site)
AICc(f1, f2, f3, f4, f5, f6, f7)
# f7 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + site + year)
# f8 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat + site + year + taxon)
# f9 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * site + year + taxon)
# f10 <- update(f1, formula = . ~ log(spheroid.volume.ml) * habitat + site + year + taxon)
# f11 <- update(f1, formula = . ~ log(spheroid.volume.ml) + habitat * site * year + taxon)





# f1 <- glmmTMB(formula = gonad_wet_mass_g ~ log(test_diam_mm),
#               dispformula = ~ site * habitat, # account for heteroskedasticity 
#               data = dat_filtered,
#               family = gaussian(link = "log"),
#               REML = TRUE)
# f2 <- update(f1, formula = . ~ log(test_diam_mm) + habitat)
# f3 <- update(f1, formula = . ~ log(test_diam_mm) + habitat + site)
# f4 <- update(f1, formula = . ~ log(test_diam_mm) * habitat + site)
# f5 <- update(f1, formula = . ~ log(test_diam_mm) + habitat * site)
# f6 <- update(f1, formula = . ~ log(test_diam_mm) * habitat * site)
# f7 <- update(f1, formula = . ~ log(test_diam_mm) : habitat + site)
# f8 <- update(f1, formula = . ~ log(test_diam_mm) : habitat : site)
# f9 <- update(f1, formula = . ~ log(test_diam_mm) : habitat * site)

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9)

# Select model f4 based on having lowest AICc. Re-fit selected model using 
# Maximum Likelihood for inference.
# f11 <- update(f11, REML = FALSE)
final_model <- f9

# Model diagnostics
simulation_output <- simulateResiduals(final_model, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
red_gonad_vol_site_anova <- glmmTMB:::Anova.glmmTMB(
  final_model,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
red_gonad_vol_site_anova
anova_table <- autoformat(xtable(red_gonad_vol_site_anova), zap = 3)
print.xtable(anova_table, type = "html", file = "./ANOVA_Results/gonadwetmass_volume_habitat_site_ANOVA")

# Tabulate mixed model summary in HTML format
tab_model(f4, file = "./Model_Summaries/gonadwetmass_volume_habitat_site")

# Visualize predictions and data 
preds <- predict(object = final_model, 
          newdata = dat_filtered,
          re.form = NA,
          se.fit = TRUE, 
          type = "link")
dat_filtered$pred_fixed_gonad_vol_site <- exp(preds$fit)
dat_filtered$pred_fixed_se_gonad_vol_site_lci <- exp(preds$fit-
                                                   preds$se.fit)
dat_filtered$pred_fixed_se_gonad_vol_site_uci <- exp(preds$fit+
                                                   preds$se.fit)

p <- ggplot(aes(x = test_diam_mm, 
                                   y = total_gonad_wet_wt_g), data = dat_filtered)+
  geom_point(aes(color = habitat, shape = year_habitat), size = 3)+
  # geom_line(aes(y = pred_fixed_gonad_vol_site, color = habitat, alpha = year))+
  # geom_ribbon(aes(ymin = pred_fixed_se_gonad_vol_site_lci,
  #                 ymax = pred_fixed_se_gonad_vol_site_uci, fill = habitat, alpha = year))+
  
  #   geom_point(aes(color = habitat, shape = habitat), size = 3,  alpha = 0.5)+
  # geom_line(aes(y = pred_fixed_gonad_vol_site, color = habitat))+
  # geom_ribbon(aes(ymin = pred_fixed_se_gonad_vol_site_lci,
  #                 ymax = pred_fixed_se_gonad_vol_site_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_shape_manual(values = c(21,24,16,17))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Test Diameter (mm)")+
  ylab("Gonadal Wet Mass (g)")+
  ggtitle("")+
  gg_options()+
  # xlim(0,180)+
  ylim(0,90)+
  theme(legend.position = c(0.15,0.7),
        legend.title = element_blank());

p

ggsave("./Figures/test.pdf", plot = p, units = "in", width = 6, height = 6)
ggsave("./Figures/test.jpg", plot = p, units = "in", width = 10, height = 6)


p_faraday <- ggplot(aes(x = test_diam_mm, 
                                   y = total_gonad_wet_wt_g), data = dat_faraday)+
  geom_point(aes(color = habitat, shape = year_habitat), size = 3)+
  scale_color_manual(values = c("darkgray", "chartreuse3"))+
  scale_fill_manual(values = c("darkgray", "chartreuse3"))+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Test Diameter (mm)")+
  ylab("Gonadal Wet Mass (g)")+
  gg_options()+
  ylim(0,90)+
  theme(legend.position = c(0.15,0.7),
        legend.title = element_blank());

p_faraday

p_multi <- cowplot::plot_grid(p_faraday, p_murchison,
                              ncol = 2)

ggsave("./Figures/test.pdf", plot = p_gonad_vol_site_red, units = "in", width = 6, height = 6)
ggsave("./Figures/test.jpg", plot = p_multi, units = "in", width = 10, height = 6)
