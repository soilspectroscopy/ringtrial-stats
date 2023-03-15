
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
int.cv <- read_csv("outputs/tab_int10CVrep10_PLSR_performance_metrics.csv")

# Overall average Lin's CCC
int.cv %>%
  group_by(soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T))

# Best preprocessing for clay
int.cv %>%
  filter(soil_property == "clay_perc") %>%
  group_by(prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T))

# Overall preprocessing per
int.cv %>%
  group_by(prep_spectra) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T))

# Worst instrument
int.cv %>%
  group_by(organization, soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T)) %>%
  filter(organization %in% c(13, 16))

# Best instrument
int.cv %>%
  group_by(organization, soil_property) %>%
  summarise(median_ccc = median(ccc, na.rm = T),
            iqr_ccc = IQR(ccc, na.rm = T),
            .groups = "drop") %>%
  arrange(soil_property, desc(median_ccc)) %>%
  View()
