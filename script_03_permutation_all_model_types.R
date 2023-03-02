
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
                                          "Cubist")))

## all

# plot

f <- function(x) {
  p05 <- quantile(x, probs = 0.10)
  p50 <- quantile(x, probs = 0.50)
  p95 <- quantile(x, probs = 0.90)
  data.frame(ymin = p05, y = p50, ymax = p50)
}

p.dispersion.vert <- performance %>%
  ggplot(aes(x = model_type, y = ccc,
             color = model_type, fill = model_type)) +
  facet_wrap(~soil_property, ncol = 1) +
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

# permutation pairs

performance

perm.models <- performance %>%
  select(soil_property, model_type, ccc)

soil.properties <- pull(distinct(perm.models, soil_property), soil_property)
model.types <- pull(distinct(perm.models, model_type), model_type)

model.types.pairs <- t(combn(model.types, 2))
model.types.pairs <- tibble("level1" = model.types.pairs[,1],
                             "level2" = model.types.pairs[,2]) %>%
  mutate_all(as.character)

model.types.pairs

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

permutation.median.list <- list()

i=1
for(i in 1:length(soil.properties)) {
  
  isoil.property <- soil.properties[i]
  
  data <- perm.models %>%
    filter(soil_property == isoil.property)
  
  cat(paste0("Started ", isoil.property, " at ", now(), "\n"))
  
  reference.grid <- model.types.pairs
  
  k=1
  for(k in 1:nrow(model.types.pairs)) {
    
    klevel1 = model.types.pairs[[k,1]]
    klevel2 = model.types.pairs[[k,2]]
    
    subset1 <- data %>%
      filter(model_type == klevel1) %>%
      pull(ccc)
    
    subset2 <- data %>%
      filter(model_type == klevel2) %>%
      pull(ccc)
    
    p.value <- permutation.test(sample1 = subset1, sample2 = subset2)
    
    reference.grid[k,"p_value"] <- p.value
    
    cat(paste0("Run for ", klevel1, " & ", klevel2, " at ", now(), "\n"))
    
  }
  
  permutation.median.list[[i]] <- reference.grid %>%
    mutate(soil_property = isoil.property, .before = 1)
  
  cat(paste0("Conclude at ", now(), "\n\n"))
  
}

permutation.median <- Reduce(bind_rows, permutation.median.list)
permutation.median

# final visualization

plot.labels.median <- permutation.median %>% # comparison
  mutate(comparison = paste(level1, level2, sep = "-")) %>% # comparison structure A-B
  select(soil_property, comparison, p_value) %>%
  nest_by(soil_property, .key = "significance") %>% # nesting for further analysis
  mutate(significance = list(deframe(significance))) %>% # transf. to named vector
  left_join({performance %>% # median values
      group_by(soil_property, model_type) %>%
      summarise(median = quantile(ccc, p=0.50), .groups = "drop") %>%
      nest_by(soil_property, .key = "median")}, by = "soil_property") %>% # nesting for further analysis
  mutate(median = list(as.data.frame(median))) %>% # mutating to required data.frame format
  mutate(letter = list(multcompLetters3(z = "model_type", y = "median", # cld
                                        x = significance, data = median)$Letters)) %>%
  mutate(letter = list(enframe(letter))) %>% # named vector to table
  mutate(letter = list(rename(letter, model_type = name, letter = value))) %>% # renaming
  mutate(join = list(left_join(median, letter, by = "model_type"))) %>% # joining the results
  select(-significance, -median, -letter) %>% # cleaning
  unnest(join) %>% # unnesting and original table format
  mutate(ccc = 1.1) # plot y label position

write_csv({plot.labels.median %>%
    select(-ccc)}, "outputs/tab_statistical_test_model_types_ccc_all.csv")

p.cld <- p.dispersion.vert +
  geom_text(data = plot.labels.median, aes(label = letter),
            color = "gray30", size = 3) +
  labs(title = "Statistical comparison of model types",
       caption = paste("Medians not sharing any letter are significantly different",
                       "by the permutation test at the 5% level of significance.\n",
                       "Box top notch refers to median while the lower notch refers to",
                       "the 10th percentile. Preprocessings are pooled together.")) +
  theme(plot.caption = element_text(size = 8, face = "italic")); p.cld

ggsave("outputs/plot_cld_model_types_all.png", p.cld,
       dpi = 300, units = "in", width = 7, height = 8, scale = 1)
