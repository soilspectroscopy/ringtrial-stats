
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

perf.plsr <- read_csv("outputs/tab_CT-KSSL_PLSR_test_performance.csv")

unique(perf.plsr$prep_spectra)

perf.plsr <- perf.plsr %>%
  mutate(prep_spectra = recode(prep_spectra, "SNVplusSG1stDer" = "SNV+SG1stDer")) %>%
  mutate(prep_spectra = factor(prep_spectra,
                               levels = c("raw",
                                          "BOC",
                                          "SG1stDer",
                                          "SNV",
                                          "SNV+SG1stDer",
                                          "wavelet",
                                          "SST")))

## Example

# plot

f <- function(x) {
  p05 <- quantile(x, probs = 0.10)
  p50 <- quantile(x, probs = 0.50)
  p95 <- quantile(x, probs = 0.90)
  data.frame(ymin = p05, y = p50, ymax = p50)
}

p.plsr.vert <- perf.plsr %>%
  ggplot(aes(x = prep_spectra, y = ccc,
             color = prep_spectra, fill = prep_spectra)) +
  facet_wrap(~soil_property, ncol = 1) +
  geom_beeswarm(size = 1, cex = 1.5, method = "hex", show.legend = F) +
  stat_summary(fun.data = f, geom = "crossbar", fill = NA,
               linewidth = 0.25, width = 0.75, show.legend = F) +
  scale_y_continuous(limits = c(-0.1, 1.2), breaks = c(0,0.25,0.50,0.75,1.00)) +
  labs(tittle = "",
       y = "Lin's CCC", x = NULL) +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.plsr.vert

# permutation test

perf.plsr

perm.plsr <- perf.plsr %>%
  select(soil_property, prep_spectra, ccc)

soil.properties <- pull(distinct(perm.plsr, soil_property), soil_property)
prep.spectra <- pull(distinct(perm.plsr, prep_spectra), prep_spectra)

prep.spectra.pairs <- t(combn(prep.spectra, 2))

prep.spectra.pairs <- tibble("level1" = prep.spectra.pairs[,1],
                             "level2" = prep.spectra.pairs[,2]) %>%
  mutate_all(as.character)

prep.spectra.pairs

# function

permutation.test <- function(sample1, sample2,
                             p.position = 0.50, n.sim = 10000,
                             seed = 1993, hypothesis = "different") {
  
  set.seed(seed)
  
  position.sample1 <- quantile(sample1, p=p.position)
  position.sample2 <- quantile(sample2, p=p.position)
  
  if(hypothesis == "different") {
    
    original.dif <- abs(position.sample1-position.sample2)
  
  }
  
  permuted.dif <- sapply(1:n.sim, function(x) {
  
    permuted.subset <- sample(c(sample1, sample2), size = length(sample1))
    permuted.position <- quantile(permuted.subset, p=p.position)
    
    if(hypothesis == "different") {
      
      abs(position.sample1-permuted.position)
      
    }
    
  })
  
  if(hypothesis == "different") {
    
    p.value <- length(which(permuted.dif > original.dif))/n.sim
    
  }
  
  return(p.value)
  
}

sample1 = rnorm(n = 100, mean = 10, sd = 5)
sample2 = rnorm(n = 100, mean = 12, sd = 5)
permutation.test(sample1 = sample1, sample2 = sample2)

## automated test

permutation.median.list <- list()

i=1
for(i in 1:length(soil.properties)) {
  
  isoil.property <- soil.properties[i]
  
  data <- perm.plsr %>%
    filter(soil_property == isoil.property)
  
  cat(paste0("Started ", isoil.property, " at ", now(), "\n"))
  
  reference.grid <- prep.spectra.pairs
  
  k=1
  for(k in 1:nrow(prep.spectra.pairs)) {
    
    klevel1 = prep.spectra.pairs[[k,1]]
    klevel2 = prep.spectra.pairs[[k,2]]
    
    subset1 <- data %>%
      filter(prep_spectra == klevel1) %>%
      pull(ccc)
    
    subset2 <- data %>%
      filter(prep_spectra == klevel2) %>%
      pull(ccc)
    
    p.value <- permutation.test(sample1 = subset1, sample2 = subset2)
    
    reference.grid[k,"p_value"] <- p.value
    
    cat(paste0("Run for ", klevel1, " & ", klevel2, " at ", now(), "\n"))
    
  }
  
  permutation.median.list[[i]] <- reference.grid %>%
    mutate(soil_property = isoil.property, .before = 1)
  
  cat(paste0("Conclude at ", now(), "\n\n"))
  
}

# long format
permutation.median <- Reduce(bind_rows, permutation.median.list)
permutation.median

# example triangle visualization
permutation.median %>%
  filter(soil_property == "carbon_org_perc") %>%
  corrr::retract(x = level1, y = level2, val = p_value)

# example compact letter display
comparison.data <- permutation.median %>%
  filter(soil_property == "carbon_org_perc") %>%
  mutate(comparison = paste(level1, level2, sep = "-")) %>%
  as.data.frame()

cld.resuls <- cldList(p_value ~ comparison, data = comparison.data,
                      threshold = 0.05,
                      reversed = TRUE)

cld.resuls

# final visualization

comparison.median <- permutation.median %>%
  mutate(comparison = paste(level1, level2, sep = "-")) %>%
  as.data.frame()

plot.labels.median <- comparison.median %>%
  nest_by(soil_property) %>%
  mutate(result = list(cldList(p_value ~ comparison, data = data,
                          threshold = 0.05,
                          reversed = TRUE))) %>%
  select(-data) %>%
  unnest(result) %>%
  select(-MonoLetter) %>%
  rename(prep_spectra = Group, letter = Letter) %>%
  mutate(ccc = 1.1)

p.plsr.vert.cld <- p.plsr.vert +
  geom_text(data = plot.labels.median, aes(label = letter),
            color = "gray30", size = 3) +
  labs(tittle = "Median statistical comparison of preprocessings - PLSR predictions",
       caption = paste("Medians not sharing any letter are significantly different",
                       "by the permutation test at the 5% level of significance.\n",
                       "Box top notch refers to median while the lower notch refers to",
                       "the 10th percentile.")) +
  theme(plot.caption = element_text(size = 8, face = "italic")); p.plsr.vert.cld

ggsave("outputs/plot_cld_example.png", p.plsr.vert.cld,
       dpi = 300, units = "in", width = 7, height = 8, scale = 1)
