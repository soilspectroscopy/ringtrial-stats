
## Loading packages
library("tidyverse")
library("lubridate")
library("readr")
library("ggbeeswarm")
library("corrr")
library("multcompView")

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
unique(performance$model_type)

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

# appending metadata

metadata <- read_csv("outputs/instruments_metadata_clean.csv")

clustering <- read_csv("outputs/proportions_clustering.csv")

clustering <- clustering %>%
  select(organization, majority) %>%
  rename(cluster = majority)

clustering %>%
  count(cluster)

performance <- left_join(performance, clustering, by = "organization") %>%
  relocate(cluster, .after = organization)

## Best cluster

performance %>%
  group_by(soil_property, cluster) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(cluster == "C4")

## Worst cluster

performance %>%
  group_by(soil_property, cluster) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(cluster == "C2")

## SNV, Cubist and MBL

performance %>%
  filter(prep_spectra %in% c("SNV")) %>%
  filter(model_type %in% c("MBL", "SST")) %>%
  group_by(soil_property, cluster) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(cluster == "C2")

## SST, Cubist and MBL

performance %>%
  filter(prep_spectra %in% c("SST")) %>%
  filter(model_type %in% c("MBL", "SST")) %>%
  group_by(soil_property, cluster) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(cluster == "C2")

## SST, Cubist and MBL

performance %>%
  filter(prep_spectra %in% c("SST")) %>%
  filter(model_type %in% c("MBL", "SST")) %>%
  group_by(soil_property, cluster) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            perc10th = quantile(ccc, p=0.1, na.rm = T)) %>%
  filter(cluster == "C4")

## Performance gain C2
gain <- performance %>%
  filter(cluster == "C2") %>%
  filter(prep_spectra %in% c("SNV", "SST")) %>%
  group_by(organization, soil_property, prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T)) %>%
  pivot_wider(names_from = "prep_spectra", values_from = "median_ccc") %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(difference = round(((SST/SNV)-1)*100, 0))
View(gain)
