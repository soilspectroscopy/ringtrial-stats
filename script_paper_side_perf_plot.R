
## Loading packages
library("tidyverse")
library("lubridate")
library("readr")
library("corrr")
library("qs")
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

options(scipen = 999)

## Mounted disk for storing big files
mnt.dir <- "~/projects/mnt-ringtrial/"
dir.dissimilarity <- paste0(mnt.dir, "dissimilarity/")

## Files
test.ids <- qread("outputs/RT_test_ids.qs")

# Predictions
list.files("outputs")

perf.plsr <- read_csv("outputs/tab_CT-KSSL_PLSR_test_performance.csv") %>%
  mutate(model_type = "plsr", .before = 1)

perf.mbl <- read_csv("outputs/tab_CT-KSSL_MBL_test_performance.csv") %>%
  mutate(model_type = "mbl", .before = 1)

perf.cubist <- read_csv("outputs/tab_CT-KSSL_Cubist_test_performance.csv") %>%
  mutate(model_type = "cubist", .before = 1)

performance <- bind_rows(perf.plsr, perf.mbl, perf.cubist)

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
  filter(!(prep_spectra == "wavelet")) %>%
  mutate(soil_property = factor(soil_property,
                                levels = c("OC",
                                           "Clay",
                                           "pH",
                                           "K")))

performance <- performance %>%
  filter(prep_spectra %in% c("SNV", "SST")) %>%
  select(organization, soil_property, model_type, prep_spectra, ccc)

performance

## Load 10CV data
performance.10CV <- read_csv(paste0("outputs/tab_int10CVrep10_PLSR_performance_metrics.csv"))

performance.10CV <- performance.10CV %>%
  mutate(prep_spectra = recode(prep_spectra, "SNVplusSG1stDer" = "SNV+SG1stDer")) %>%
  mutate(prep_spectra = factor(prep_spectra,
                               levels = c("raw",
                                          "BOC",
                                          "SG1stDer",
                                          "SNV",
                                          "SNV+SG1stDer",
                                          "wavelet",
                                          "SST"))) %>%
  mutate(soil_property = recode(soil_property,
                                "carbon_org_perc" = "OC",
                                "clay_perc" = "Clay",
                                "pH_H20" = "pH",
                                "potassium_cmolkg" = "K")) %>%
  filter(!(prep_spectra == "wavelet")) %>%
  mutate(soil_property = factor(soil_property,
                                levels = c("OC",
                                           "Clay",
                                           "pH",
                                           "K")))

performance.10CV <- performance.10CV %>%
  filter(prep_spectra %in% c("SNV")) %>%
  mutate(model_type = "10CVrep10") %>%
  select(organization, soil_property, model_type, prep_spectra, ccc)

performance.plot <- bind_rows(performance.10CV, performance) %>%
  mutate(label = ifelse(model_type == "10CVrep10",
                        "10CVrep10-PLSR",
                        paste0("KSSL-", model_type)),
         preprocessing = ifelse(prep_spectra == "SNV",
                               "beforeSST",
                               "afterSST")) %>%
  mutate(label = factor(label, levels = c("10CVrep10-PLSR",
                                          "KSSL-PLSR",
                                          "KSSL-MBL",
                                          "KSSL-Cubist")),
         preprocessing = factor(preprocessing, levels = c("beforeSST",
                                                          "afterSST")))
## Visualization

p.metrics <- performance.plot %>%
  mutate(organization = as.factor(organization)) %>%
  # filter(soil_property == "OC") %>%
  ggplot() +
  geom_col_pattern(aes(x = organization, y = ccc, fill = label, pattern = preprocessing),
                   position = position_dodge(preserve = "single"),
                   color = "gray10",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   show.legend = TRUE) +
  labs(x = "Instrument", y = "Lin's CCC", fill = "", pattern = "") +
  facet_wrap(~soil_property, ncol = 1) +
  scale_fill_manual(values = c("gray20", "gray40", "gray60", "gray80")) +
  scale_pattern_manual(values = c(beforeSST = "none", afterSST = "stripe")) +
  theme_light() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))

ggsave(paste0("outputs/plot_paper_final_performance.png"),
       p.metrics, dpi = 300, width = 8, height = 8,
       units = "in", scale = 1)
