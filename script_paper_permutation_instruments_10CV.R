
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

performance <- read_csv("outputs/tab_int10CVrep10_PLSR_performance_metrics.csv")

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
  filter(!(prep_spectra == "wavelet")) %>%
  mutate(organization = as.factor(organization))

## all

# plot

f <- function(x) {
  p05 <- quantile(x, probs = 0.10)
  p50 <- quantile(x, probs = 0.50)
  p95 <- quantile(x, probs = 0.90)
  data.frame(ymin = p05, y = p50, ymax = p50)
}

p.dispersion.vert <- performance %>%
  ggplot(aes(x = organization, y = ccc)) +
  geom_beeswarm(size = 0.35, cex = 0.85, color = "gray30", method = "hex", show.legend = F) +
  stat_summary(fun.data = f, geom = "crossbar", fill = NA, color = "gray40",
               linewidth = 0.25, width = 0.75, show.legend = F) +
  scale_y_continuous(limits = c(-0.2, 1.2), breaks = c(0,0.25,0.50,0.75,1.00)) +
  labs(tittle = "",
       y = "Lin's CCC", x = "Instrument") +
  theme_light() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); p.dispersion.vert

# permutation pairs

performance

perm.models <- performance %>%
  select(organization, ccc)

organizations <- pull(distinct(perm.models, organization), organization)

organization.pairs <- t(combn(organizations, 2))
organization.pairs <- tibble("level1" = organization.pairs[,1],
                             "level2" = organization.pairs[,2]) %>%
  mutate_all(as.character)

organization.pairs

# permutation function

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

# median comparison

data <- perm.models %>%
  filter(!is.na(ccc))

reference.grid <- organization.pairs

k=1
for(k in 1:nrow(organization.pairs)) {
  
  klevel1 = organization.pairs[[k,1]]
  klevel2 = organization.pairs[[k,2]]
  
  subset1 <- data %>%
    filter(organization == klevel1) %>%
    pull(ccc)
  
  subset2 <- data %>%
    filter(organization == klevel2) %>%
    pull(ccc)
  
  p.value <- permutation.test(sample1 = subset1, sample2 = subset2)
  
  reference.grid[k,"p_value"] <- p.value
  
  cat(paste0("Run for ", klevel1, " & ", klevel2, " at ", now(), "\n"))
  
}

permutation.median <- reference.grid
permutation.median

# final visualization

plot.labels.median <- permutation.median %>% # comparison
  mutate(comparison = paste(level1, level2, sep = "-")) %>% # comparison structure A-B
  select(comparison, p_value) %>%
  mutate(soil_property = "all") %>%
  nest_by(soil_property, .key = "significance") %>% # nesting for further analysis
  mutate(significance = list(deframe(significance))) %>% # transf. to named vector
  left_join({performance %>% # median values
      group_by(organization) %>%
      summarise(median = quantile(ccc, p=0.50, na.rm = T), .groups = "drop") %>%
      mutate(soil_property = "all") %>%
      nest_by(soil_property, .key = "median")}, by = "soil_property") %>% # nesting for further analysis
  mutate(median = list(as.data.frame(median))) %>% # mutating to required data.frame format
  mutate(letter = list(multcompLetters3(z = "organization", y = "median", # cld
                                        x = significance, data = median)$monospacedLetters)) %>%
  mutate(letter = list(enframe(letter))) %>% # named vector to table
  mutate(letter = list(rename(letter, organization = name, letter = value))) %>% # renaming
  mutate(join = list(left_join(median, letter, by = "organization"))) %>% # joining the results
  select(-significance, -median, -letter) %>% # cleaning
  unnest(join) %>% # unnesting and original table format
  mutate(ccc = 1.1) %>% # plot y label position
  ungroup()

plot.labels.median <- plot.labels.median %>%
  select(-soil_property) %>%
  mutate(letter = gsub(" ", "_", letter))

plot.labels.median

write_csv({plot.labels.median %>%
    select(-ccc)}, "outputs/tab_statistical_test_instruments_int10CV_ccc.csv")

p.cld <- p.dispersion.vert +
  geom_text(data = plot.labels.median, aes(label = letter),
            color = "gray30", size = 3, angle = 45) +
  # labs(title = "Statistical comparison of internal 10CV of instruments",
  #      caption = paste("Medians not sharing any letter are significantly different",
  #                      "by the permutation test at the 5% level of significance.\n",
  #                      "Box top notch refers to median while the lower notch refers to",
  #                      "the 10th percentile. Preprocessings and soil properties are pooled together.")) +
  theme(plot.caption = element_text(size = 8, face = "italic")); p.cld

ggsave("outputs/plot_cld_instruments_10CV.png", p.cld,
       dpi = 300, units = "in", width = 8, height = 4, scale = 1)
