
## Loading packages
library("tidyverse")
library("lubridate")
library("readr")
library("corrr")
library("qs")

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
  filter(!(prep_spectra == "wavelet"))

## Median soil properties
performance %>%
  group_by(soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T))

## KSSL overall

performance %>%
  filter(!(prep_spectra == "wavelet")) %>%
  filter(organization == 16) %>%
  group_by(soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T))

## KSSL supplement

# performance %>%
#   filter(!(prep_spectra == "wavelet")) %>%
#   arrange(soil_property, model_type, prep_spectra) %>%
#   filter(organization == 16) %>%
#   select(soil_property, model_type, prep_spectra, rmse, bias, rsq, ccc, rpiq) %>%
#   clipr::write_clip()

## Median SNV
performance %>%
  group_by(soil_property, prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(prep_spectra == "SNV")

## Median SST
performance %>%
  group_by(soil_property, prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T)) %>%
  filter(prep_spectra == "SST")

## Median PLSR
performance %>%
  group_by(model_type, soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T))

## Lowest performances in preprocessings
performance %>%
  group_by(soil_property, prep_spectra) %>%
  mutate(p10_ccc = quantile(ccc, p = 0.10, na.rm = T)) %>%
  filter(ccc <= p10_ccc) %>%
  count(soil_property, prep_spectra, organization)

performance %>%
  group_by(soil_property, prep_spectra) %>%
  mutate(p10_ccc = quantile(ccc, p = 0.10, na.rm = T)) %>%
  filter(ccc <= p10_ccc) %>%
  count(soil_property, prep_spectra, organization) %>%
  group_by(organization) %>%
  summarise(sum = sum(n))

## Lowest performances in model types
performance %>%
  group_by(soil_property, model_type) %>%
  mutate(p10_ccc = quantile(ccc, p = 0.10, na.rm = T)) %>%
  filter(ccc <= p10_ccc) %>%
  count(soil_property, model_type, organization)

performance %>%
  group_by(soil_property, model_type) %>%
  mutate(p10_ccc = quantile(ccc, p = 0.10, na.rm = T)) %>%
  filter(ccc <= p10_ccc) %>%
  count(soil_property, model_type, organization) %>%
  group_by(organization) %>%
  summarise(sum = sum(n))

## Lowest performances for SST
performance %>%
  filter(prep_spectra == "SST") %>%
  group_by(soil_property, model_type) %>%
  mutate(p10_ccc = quantile(ccc, p = 0.10, na.rm = T)) %>%
  filter(ccc <= p10_ccc) %>%
  count(soil_property, model_type, organization) %>%
  group_by(organization) %>%
  summarise(sum = sum(n))

performance %>%
  filter(prep_spectra %in% c("SNV", "SST")) %>%
  filter(organization %in% c(4, 5, 9, 12, 13, 17)) %>%
  group_by(organization, soil_property, prep_spectra) %>%
  summarise(ccc = median(ccc), .groups = "drop") %>%
  pivot_wider(names_from = "prep_spectra", values_from = "ccc") %>%
  mutate(difference = (SST/SNV-1)*100)

p.export <- performance %>%
  filter(prep_spectra %in% c("SNV", "SST")) %>%
  filter(organization %in% c(4, 5, 9, 12, 13, 17)) %>%
  group_by(organization, soil_property, prep_spectra) %>%
  summarise(ccc = median(ccc), .groups = "drop") %>%
  pivot_wider(names_from = "prep_spectra", values_from = "ccc") %>%
  mutate(difference = (SST/SNV-1)*100) %>%
  group_by(organization, soil_property)

clipr::write_clip(p.export)

## Dissimilarity

ids.test <- qread("outputs/RT_test_ids.qs")

all.mirspectra.SNV.dissim <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_SNV.csv"))

all.mirspectra.SNV.dissim <- all.mirspectra.SNV.dissim %>%
  filter(sample_id %in% ids.test) %>%
  mutate(prep_spectra = "SNV", .after = 1)

all.mirspectra.SST.dissim <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_SST.csv"))

all.mirspectra.SST.dissim <- all.mirspectra.SST.dissim %>%
  filter(sample_id %in% ids.test) %>%
  select(-ct_subset) %>%
  mutate(prep_spectra = "SST", .after = 1)

all.dissim <- bind_rows(all.mirspectra.SNV.dissim,
                        all.mirspectra.SST.dissim)

all.dissim %>%
  filter(organization %in% c(4, 5, 9, 12, 13, 17)) %>%
  group_by(organization, prep_spectra) %>%
  summarise(dissim = median(distance))
