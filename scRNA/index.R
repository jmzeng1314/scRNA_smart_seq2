rm(list = ls())  ## 魔幻操作，一键清空~
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

BiocManager::install('Seurat',ask = F,update = F)




