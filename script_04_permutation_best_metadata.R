
## Loading packages
library("tidyverse")
library("lubridate")
library("readr")
library("ggbeeswarm")
library("corrr")
library("rcompanion")

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
                                        "Cubist")))

# best (cubist & SST+SNV) spectra preprocessing

performance <- performance %>%
  filter(prep_spectra %in% c("SST", "SNV")) %>%
  filter(model_type == "Cubist")

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

# plot

f <- function(x) {
  p05 <- quantile(x, probs = 0.10)
  p50 <- quantile(x, probs = 0.50)
  p95 <- quantile(x, probs = 0.90)
  data.frame(ymin = p05, y = p50, ymax = p50)
}

p.dispersion.vert <- performance %>%
  ggplot(aes(x = cluster, y = ccc,
             color = cluster, fill = cluster)) +
  facet_grid(soil_property~prep_spectra) +
  geom_beeswarm(size = 1, cex = 1.2, method = "hex", show.legend = F) +
  stat_summary(fun.data = f, geom = "crossbar", fill = NA,
               linewidth = 0.25, width = 0.75, show.legend = F) +
  scale_y_continuous(limits = c(-0.2, 1.2), breaks = c(0,0.25,0.50,0.75,1.00)) +
  labs(tittle = "",
       y = "Lin's CCC", x = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.dispersion.vert

# final visualization

soil.properties <- pull(performance, soil_property)
prep.spectra <- pull(performance, prep_spectra)

plot.labels.median <- performance %>% # comparison
  count(soil_property, prep_spectra, cluster) %>%
  mutate(ccc = 1.1)

p.cld <- p.dispersion.vert +
  # geom_text(data = plot.labels.median, aes(label = paste0("n=", n)),
  #           color = "gray30", size = 3) +
  labs(title = "Preprocessing comparison for metadata clusters with Cubist",
       caption = paste("Sample sizes: C1=1, C2=6, C3=6, C4=7.")) +
  theme(plot.caption = element_text(size = 8, face = "italic")); p.cld

ggsave("outputs/plot_best_preprocessing_cubist.png", p.cld,
       dpi = 300, units = "in", width = 7, height = 8, scale = 1)
