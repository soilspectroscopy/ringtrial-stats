
## Loading packages
library("tidyverse")

## Creating input/output dirs
if(!dir.exists("outputs")){dir.create("outputs")}

## Mounted disk for storing big files
mnt.dir <- "~/projects/mnt-ringtrial/"

## Copying other results
metadata.dir <- "~/projects/git/ringtrial-metadata/"
modeling.dir <- "~/projects/git/ringtrial-modeling/"

list.files(paste0(metadata.dir, "outputs/"))

metatada.file <- paste0(metadata.dir, "outputs/instruments_metadata_clean.csv")
file.copy(from = metatada.file, to = paste0("outputs/", basename(metatada.file)), overwrite = T)

clustering.file <- paste0(metadata.dir, "outputs/proportions_clustering.csv")
file.copy(from = clustering.file, to = paste0("outputs/", basename(clustering.file)), overwrite = T)

list.files(paste0(modeling.dir, "outputs/"))

ids.test.file <- paste0(modeling.dir, "outputs/RT_test_ids.qs")
file.copy(from = ids.test.file, to = paste0("outputs/", basename(ids.test.file)), overwrite = T)

performance <- list.files(paste0(modeling.dir, "outputs/"), pattern = "performance_metrics.csv|test_performance.csv")
performance

i=1
for(i in 1:length(performance)) {
  ifile <- performance[i]
  file.copy(from = paste0(modeling.dir, "outputs/", ifile),
            to = paste0("outputs/", basename(ifile)), overwrite = T)
}
