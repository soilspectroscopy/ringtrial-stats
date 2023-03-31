
## Loading packages
library("tidyverse")
library("lubridate")
library("readr")
library("ggbeeswarm")
library("corrr")
library("multcompView")
library("cowplot")

options(scipen = 999)

## Mounted disk for storing big files
mnt.dir <- "~/projects/mnt-ringtrial/"

## Files
list.files("outputs")

perf.plsr <- read_csv("outputs/tab_CT-KSSL_PLSR_test_performance.csv") %>%
  mutate(model_type = "plsr", .before = 1)

perf.mbl <- read_csv("outputs/tab_CT-KSSL_MBL_test_performance.csv") %>%
  mutate(model_type = "mbl", .before = 1)

perf.cubist <- read_csv("outputs/tab_CT-KSSL_Cubist_test_performance.csv") %>%
  mutate(model_type = "cubist", .before = 1)

performance <- bind_rows(perf.plsr, perf.mbl, perf.cubist)

unique(performance$prep_spectra)

performance <- performance %>%
  mutate(prep_spectra = recode(prep_spectra, "SNVplusSG1stDer" = "SNV+SG1stDer")) %>%
  mutate(prep_spectra = factor(prep_spectra,
                               levels = c("raw",
                                          "BOC",
                                          "SG1stDer",
                                          "SNV",
                                          "SNV+SG1stDer",
                                          "wavelet",
                                          "SST"))) %>%
  mutate(model_type = recode(model_type,
                             "cubist" = "Cubist",
                             "plsr" = "PLSR",
                             "mbl" = "MBL")) %>%
  mutate(model_type = factor(model_type,
                             levels = c("PLSR",
                                        "MBL",
                                        "Cubist"))) %>%
  mutate(soil_property = recode(soil_property,
                                "carbon_org_perc" = "OC",
                                "clay_perc" = "Clay",
                                "pH_H20" = "pH",
                                "potassium_cmolkg" = "K")) %>%
  mutate(soil_property = factor(soil_property,
                                levels = c("OC",
                                           "Clay",
                                           "pH",
                                           "K"))) %>%
  filter(!(prep_spectra == "wavelet"))

f <- function(x) {
  p05 <- quantile(x, probs = 0.10)
  p50 <- quantile(x, probs = 0.50)
  p95 <- quantile(x, probs = 0.90)
  data.frame(ymin = p05, y = p50, ymax = p50)
}

## Plot preprocessing

p.dispersion.preps <- performance %>%
  ggplot(aes(x = prep_spectra, y = ccc)) +
  facet_wrap(~soil_property, ncol = 1) +
  geom_beeswarm(size = 0.35, cex = 0.85, color = "gray30", method = "hex", show.legend = F) +
  stat_summary(fun.data = f, geom = "crossbar", fill = NA, color = "gray40",
               linewidth = 0.25, width = 0.75, show.legend = F) +
  scale_y_continuous(limits = c(-0.2, 1.2), breaks = c(0,0.25,0.50,0.75,1.00)) +
  labs(tittle = "",
       y = "Lin's CCC", x = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.dispersion.preps

plot.labels.median <- read_csv("outputs/tab_statistical_test_preprocessings_ccc_all.csv") %>%
  mutate(soil_property = recode(soil_property,
                                "carbon_org_perc" = "OC",
                                "clay_perc" = "Clay",
                                "pH_H20" = "pH",
                                "potassium_cmolkg" = "K")) %>%
  mutate(soil_property = factor(soil_property,
                                levels = c("OC",
                                           "Clay",
                                           "pH",
                                           "K")))

plot.labels.median <- plot.labels.median %>%
  filter(!(prep_spectra == "wavelet")) %>%
  mutate(ccc = 1.1)

p.cld.preps <- p.dispersion.preps +
  geom_text(data = plot.labels.median, aes(label = letter),
            color = "gray30", size = 3) +
  theme(plot.caption = element_text(size = 8, face = "italic")); p.cld.preps

## Plot model types

p.dispersion.models <- performance %>%
  ggplot(aes(x = model_type, y = ccc)) +
  facet_wrap(~soil_property, ncol = 1) +
  geom_beeswarm(size = 0.35, cex = 0.85, color = "gray30", method = "hex", show.legend = F) +
  stat_summary(fun.data = f, geom = "crossbar", fill = NA, color = "gray40",
               linewidth = 0.25, width = 0.75, show.legend = F) +
  scale_y_continuous(limits = c(-0.2, 1.2), breaks = c(0,0.25,0.50,0.75,1.00)) +
  labs(tittle = "",
       y = "Lin's CCC", x = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.dispersion.models

plot.labels.median.models <- read_csv("outputs/tab_statistical_test_model_types_ccc_all.csv") %>%
  mutate(soil_property = recode(soil_property,
                                "carbon_org_perc" = "OC",
                                "clay_perc" = "Clay",
                                "pH_H20" = "pH",
                                "potassium_cmolkg" = "K")) %>%
  mutate(soil_property = factor(soil_property,
                                levels = c("OC",
                                           "Clay",
                                           "pH",
                                           "K")))

plot.labels.median.models <- plot.labels.median.models %>%
  mutate(ccc = 1.1)

p.cld.models <- p.dispersion.models +
  geom_text(data = plot.labels.median.models, aes(label = letter),
            color = "gray30", size = 3) +
  theme(plot.caption = element_text(size = 8, face = "italic")); p.cld.models

## Together

p.together <- plot_grid(p.cld.preps, p.cld.models, ncol = 2, labels = "auto", rel_widths = c(1.25,0.75))
p.together

ggsave("outputs/plot_paper_preps_models.png", p.together,
       width = 8, height = 6, dpi = 300, scale = 1, units = "in")
