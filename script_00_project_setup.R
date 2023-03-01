## Loading packages
library("tidyverse")

## Creating input/output dirs
if(!dir.exists("outputs")){dir.create("outputs")}

## Mounted disk for storing big files
mnt.dir <- "~/projects/mnt-ringtrial/"

## Copying performance results
metadata.dir <- "~/projects/soilspec4gg-mac/ringtrial-metadata/"
modeling.dir <- "~/projects/soilspec4gg-mac/ringtrial-modeling/"

metatada.file <- paste0(metadata.dir, "outputs/instruments_metadata_clean.csv")
file.copy(from = metatada.file, to = paste0("outputs/", basename(metatada.file)), overwrite = T)

list.files(paste0(modeling.dir, "outputs/"))
performance <- list.files(paste0(modeling.dir, "outputs/"), pattern = "performance_metrics.csv|test_performance.csv")
performance

i=1
for(i in 1:length(performance)) {
  ifile <- performance[i]
  file.copy(from = paste0(modeling.dir, "outputs/", ifile),
            to = paste0("outputs/", basename(ifile)), overwrite = T)
}
