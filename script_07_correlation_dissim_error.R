
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

# Distance values

dir.dissimilarity <- paste0(mnt.dir, "dissimilarity/")

dissim.raw <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_raw.csv")) %>%
  mutate(prep_spectra = "raw", .before = 1)

dissim.BOC <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_BOC.csv")) %>%
  mutate(prep_spectra = "BOC", .before = 1)

dissim.SG1stDer <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_SG1stDer.csv")) %>%
  mutate(prep_spectra = "SG1stDer", .before = 1)

dissim.SNV <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_SNV.csv")) %>%
  mutate(prep_spectra = "SNV", .before = 1)

dissim.SNVplusSG1stDer <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_SNVplusSG1stDer.csv")) %>%
  mutate(prep_spectra = "SNVplusSG1stDer", .before = 1)

dissim.wavelet <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_wavelet.csv")) %>%
  mutate(prep_spectra = "wavelet", .before = 1)

dissim.SST <- read_csv(paste0(dir.dissimilarity, "dissim_euclidean_SST.csv")) %>%
  mutate(prep_spectra = "SST", .before = 1)

distance <- bind_rows(dissim.raw, dissim.BOC, dissim.SG1stDer,
                      dissim.SNV, dissim.SNVplusSG1stDer, dissim.wavelet,
                      dissim.SST) %>%
  mutate(prep_spectra = recode(prep_spectra, "SNVplusSG1stDer" = "SNV+SG1stDer"))

distance <- distance %>%
  select(-ct_subset) %>%
  filter(sample_id %in% test.ids)

distance.median <- distance %>%
  group_by(organization, prep_spectra) %>%
  summarise(distance = median(distance), .groups = "drop")

cor.data <- performance %>%
  select(model_type, soil_property, organization, prep_spectra, rmse, ccc) %>%
  left_join(distance.median, by = c("organization", "prep_spectra")) %>%
  filter(!(organization == 16))

## Linear trend

p.lm.rmse <- ggplot(cor.data, aes(x = distance, y = rmse, color = model_type)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = "y ~ x", linewidth = 0.5,
              show.legend = FALSE, se = FALSE) +
  facet_grid(soil_property~prep_spectra, scales = "free") +
  labs(color = "", x = "Median dissimilarity", y = "RMSE") +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.lm.rmse

ggsave("outputs/plot_lm_dissim_rmse.png", p.lm.rmse,
       dpi = 300, units = "in", width = 8, height = 5, scale = 1.25)

p.lm.ccc <- ggplot(cor.data, aes(x = distance, y = ccc, color = model_type)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = "y ~ x", linewidth = 0.5,
              show.legend = FALSE, se = FALSE) +
  facet_grid(soil_property~prep_spectra, scales = "free") +
  labs(color = "", x = "Median dissimilarity", y = "Lin's CCC") +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.lm.ccc

ggsave("outputs/plot_lm_dissim_ccc.png", p.lm.ccc,
       dpi = 300, units = "in", width = 8, height = 5, scale = 1.25)

## Correlation plot

cor.result <- cor.data %>%
  select(-organization) %>%
  group_by(model_type, soil_property, prep_spectra) %>%
  summarise(cor_rmse = cor(distance, rmse),
            cor_ccc = cor(distance, ccc), .groups = "drop") %>%
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

p.cor.rmse <- ggplot(cor.result, aes(x = prep_spectra, y = soil_property, fill = cor_rmse, label = round(cor_rmse, 2))) +
  geom_tile() +
  geom_text(color = "gray20", size = 3) +
  scale_fill_gradient2(midpoint = 0, low = "darkred", mid = "white", high = "darkgreen") +
  facet_wrap(~model_type, ncol = 1) +
  labs(title = "Correlation between RMSE and median dissimilarity", fill = "Correlation", x = NULL, y = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.cor.rmse

ggsave("outputs/plot_cor_dissim_rmse.png", p.cor.rmse,
       dpi = 300, units = "in", width = 6, height = 8, scale = 1)

p.cor.ccc <- ggplot(cor.result, aes(x = prep_spectra, y = soil_property, fill = cor_ccc, label = round(cor_ccc, 2))) +
  geom_tile() +
  geom_text(color = "gray80", size = 3) +
  scale_fill_gradient2(midpoint = 0, low = "darkred", mid = "white", high = "darkgreen") +
  facet_wrap(~model_type, ncol = 1) +
  labs(title = "Correlation between Lin's CCC and median dissimilarity", fill = "Correlation", x = NULL, y = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.cor.ccc

ggsave("outputs/plot_cor_dissim_ccc.png", p.cor.ccc,
       dpi = 300, units = "in", width = 6, height = 8, scale = 1)
