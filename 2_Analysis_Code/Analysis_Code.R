###############################################################################
# Analysis code for "Zombies of the Nearshore: metabolic depression in sea    #
# urchin barrens associated with food deprivation."                           #                          #
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

# Load functions --------------------------------------------------------------
source("./2_Analysis_Code/Functions.R")

# Set general plot aesthetics
plot_colors <- c("red", "blue")

# Load data -------------------------------------------------------------------
dat <- read_excel("./1_Data/Respirometry_Data_M_franciscanus.xlsx", 
                  sheet = 1)
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

# RMR Analysis (Wet Mass): effects of site and habitat  --------------------
# Model selection: Fit several reasonable GLMMs for fixed effects on volumetric
# oxygen consumption rate using Restricted Maximum Likelihood. 
f1 <- glmmTMB(volumetric_rate ~ log(animal_wet_mass_g) + (1 | date/chamber_id), 
              data = dat_filtered,
              family = gaussian(link = "log"),
              REML = TRUE)
f2 <- update(f1, formula = . ~ log(animal_wet_mass_g) + habitat + (1 | date/chamber_id))
f3 <- update(f1, formula = . ~ log(animal_wet_mass_g) + habitat + site + (1 | date/chamber_id))
f4 <- update(f1, formula = . ~ log(animal_wet_mass_g) * habitat + site + (1 | date/chamber_id))
f5 <- update(f1, formula = . ~ log(animal_wet_mass_g) + habitat * site + (1 | date/chamber_id))
f6 <- update(f1, formula = . ~ log(animal_wet_mass_g) * habitat * site + (1 | date/chamber_id))
f7 <- update(f1, formula = . ~ log(animal_wet_mass_g) : habitat + site + (1 | date/chamber_id))
f8 <- update(f1, formula = . ~ log(animal_wet_mass_g) : habitat : site + (1 | date/chamber_id))
f9 <- update(f1, formula = . ~ log(animal_wet_mass_g) : habitat * site + (1 | date/chamber_id))

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9)

f5 <- update(f5, REML = FALSE)
fit_glmmTMB <- f5
# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
wet_mass_x_site_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
wet_mass_x_site_anova

# Visualize predictions and data 
fixed_predictions <- predict(object = fit_glmmTMB, 
                                 newdata = dat_filtered,
                                 re.form = ~ 0,
                                 se.fit = TRUE)
dat_filtered$pred_fixed_fit <- exp(fixed_predictions$fit)
dat_filtered$pred_fixed_se_lci <- exp(fixed_predictions$fit-
                                                 fixed_predictions$se.fit)
dat_filtered$pred_fixed_se_uci <- exp(fixed_predictions$fit+
                                                 fixed_predictions$se.fit)

p_wet_mass_x_site <- ggplot(aes(x = animal_wet_mass_g,
                            y = volumetric_rate), data = dat_filtered)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_fit, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_lci, 
                  ymax = pred_fixed_se_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = plot_colors)+
  scale_fill_manual(values = plot_colors)+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Wet Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = c(0.1,0.8));

p_wet_mass_x_site

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

f5 <- update(f5, REML = FALSE)
fit_glmmTMB <- f5
# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
volume_x_site_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
volume_x_site_anova

# Visualize predictions and data 
fixed_predictions <- predict(object = fit_glmmTMB, 
                             newdata = dat_filtered,
                             re.form = ~ 0,
                             se.fit = TRUE)
dat_filtered$pred_fixed_fit <- exp(fixed_predictions$fit)
dat_filtered$pred_fixed_se_lci <- exp(fixed_predictions$fit-
                                        fixed_predictions$se.fit)
dat_filtered$pred_fixed_se_uci <- exp(fixed_predictions$fit+
                                        fixed_predictions$se.fit)

p_volume_x_site <- ggplot(aes(x = spheroid.volume.ml,
                                y = volumetric_rate), data = dat_filtered)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_fit, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_lci, 
                  ymax = pred_fixed_se_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = plot_colors)+
  scale_fill_manual(values = plot_colors)+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Wet Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = c(0.1,0.8));

p_volume_x_site

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

# Although f2 has the lowest AICc, we select f3 to examine habitat:body size interaction. 
# Re-fit selected model using Maximum Likelihood for inference.
f3 <- update(f3, REML = FALSE)
fit_glmmTMB <- f3

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
afdm_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
afdm_anova

# Visualize predictions and data 
fixed_predictions <- predict(object = fit_glmmTMB, 
                                      newdata = dat_filtered_AFDM,
                                      re.form = ~ 0,
                                      se.fit = TRUE)
dat_filtered_AFDM$pred_fixed_fit <- exp(fixed_predictions$fit)
dat_filtered_AFDM$pred_fixed_se_lci <- exp(fixed_predictions$fit-
                                                 fixed_predictions$se.fit)
dat_filtered_AFDM$pred_fixed_se_uci <- exp(fixed_predictions$fit+
                                                 fixed_predictions$se.fit)

p_afdm <- ggplot(aes(x = whole.afdm.g.approx, 
                                 y = volumetric_rate), data = dat_filtered_AFDM)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_fit, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_lci, 
                  ymax = pred_fixed_se_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = plot_colors)+
  scale_fill_manual(values = plot_colors)+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  xlab("Total Ash-Free Dry Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = "");

p_afdm

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
fit_glmmTMB <- f3

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
gonad_afdm_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
gonad_afdm_anova

# Visualize predictions and data 
fixed_predictions <- predict(object = fit_glmmTMB, 
                                  newdata = dat_filtered_AFDM,
                                  re.form = ~ 0,
                                  se.fit = TRUE)
dat_filtered_AFDM$pred_fixed_fit <- exp(fixed_predictions$fit)
dat_filtered_AFDM$pred_fixed_se_lci <- exp(fixed_predictions$fit-
                                                        fixed_predictions$se.fit)
dat_filtered_AFDM$pred_fixed_se_uci <- exp(fixed_predictions$fit+
                                                        fixed_predictions$se.fit)

p_gonad_afdm <- ggplot(aes(x = gonad.afdm.g.total.approx, 
                         y = volumetric_rate), data = dat_filtered_AFDM)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_fit, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_lci, 
                  ymax = pred_fixed_se_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = plot_colors)+
  scale_fill_manual(values = plot_colors)+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  xlab("Gonadal Ash-Free Dry Mass (g)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = "");

p_gonad_afdm

# Combine RMR regression visualizations ---------------------------------------
# 3-panel figure with RMR regressed against test volume, total AFDM, and 
# gonadal AFDM by habitat for urchins from Surge Narrows
RMR_multiplot_SN <- cowplot::plot_grid(p_afdm, p_gonad_afdm,
                                     ncol = 2,
                                     labels = c('A', 'B'), 
                                     label_size = 25)
RMR_multiplot_SN

# 3-panel figure with RMR regressed against test volume by habitat and site for 
# urchins from Faraday Island, Murchison Island, and Surge Narrows
RMR_multiplot_hab_x_site <- cowplot::plot_grid(RMR_multiplot_SN, p_wet_mass_x_site,
                              ncol = 1,
                              # labels = c('', 'C'), 
                              label_size = 25)
RMR_multiplot_hab_x_site

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

# Calculate Second-order Akaike Information Criterion for various model fits.
AICc(f1, f2, f3, f4, f5, f6, f7, f8, f9)

# Select model f4 based on having lowest AICc. Re-fit selected model using 
# Maximum Likelihood for inference.
f4 <- update(f4, REML = FALSE)
fit_glmmTMB <- f4

# Model diagnostics
simulation_output <- simulateResiduals(fit_glmmTMB, n = 250)

# Plot scaled residuals:
plot(simulation_output)

# Test for uniformity, outliers, and dispersion:
testResiduals(simulation_output)

# Run Analysis of Deviance (Type III Wald chisquare tests)
gonad_anova <- glmmTMB:::Anova.glmmTMB(
  fit_glmmTMB,
  type = c("III"),
  test.statistic = c("Chisq"),
  component = "cond"
)

# Render test results
gonad_anova

# Visualize predictions and data 
fixed_predictions <- predict(object = fit_glmmTMB, 
                             newdata = dat_filtered,
                             re.form = ~ 0,
                             se.fit = TRUE)
dat_filtered$pred_fixed_fit <- exp(fixed_predictions$fit)
dat_filtered$pred_fixed_se_lci <- exp(fixed_predictions$fit-
                                             fixed_predictions$se.fit)
dat_filtered$pred_fixed_se_uci <- exp(fixed_predictions$fit+
                                             fixed_predictions$se.fit)

p_gonad <- ggplot(aes(x = spheroid.volume.ml, 
                           y = total_gonad_approx_g), data = dat_filtered)+
  geom_point(aes(color = habitat, shape = habitat), size = 3, alpha = 0.5)+
  geom_line(aes(y = pred_fixed_fit, color = habitat))+
  geom_ribbon(aes(ymin = pred_fixed_se_lci,
                  ymax = pred_fixed_se_uci, fill = habitat), alpha = 0.5)+
  scale_color_manual(values = plot_colors)+
  scale_fill_manual(values = plot_colors)+
  scale_x_continuous(trans = "log10")+
  annotation_logticks(sides = "b")+
  facet_wrap(vars(site))+
  xlab("Test Volume (mL)")+
  labs(y=expression(VO["2"]*" (mgO"["2"]*"h"^-1*")"))+
  gg_options()+
  theme(legend.position = "");

p_gonad