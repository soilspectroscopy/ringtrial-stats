
## Loading packages
library("tidyverse")
library("lubridate")
library("readr")
library("corrr")
library("qs")

options(scipen = 999)

## Mounted disk for storing big files
mnt.dir <- "~/projects/mnt-ringtrial/"

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
                                        "Cubist")))

## Median SST
performance %>%
  group_by(soil_property, prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T)) %>%
  filter(prep_spectra == "SST")

## Median SNV
performance %>%
  group_by(soil_property, prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(prep_spectra == "SNV")

## Median PLSR
performance %>%
  group_by(model_type, soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  View()
