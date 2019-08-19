# metagenes from ref papers
# for 参考文献ref32：proliferation signature
ref32 <- c("BAG1","ESR1","FOXA1","GPR160","NAT1","NAT1","MAPT","MLPH","PGR")
# for ref31: stroma-related signatures
ref31 <- c('DCN', 'CSPG2', 'CDH11', 'COL3A1', 'FAP', 'PEDF', 'FBN1', 'PDGFRL', 'CTSK', 'HTRA1', 'ASPN', 'SPARC', 'COL5A2', 'LOXL1', 'MMP2', 'SPON1', 'SFRP4', 'ITGBL1', 'CALD1', 'COPZ2', 'MFAP2', 'ANGPTL2', 'PLAU', 'COL1A2', 'LRRC17', 'C1QTNF3', 'SNAI2', 'PCOLCE', 'POSTN', 'ECM2', 'FBLN1', 'ADAM12', 'MMP11', 'AEBP1', 'PDGFRB', 'GAS1', 'COL6A3', 'RARRES2', 'COL6A1', 'TGFB3', 'NDN', 'C1R', 'LRP1', 'COL10A1', 'DPYSL3', 'OLFML2B', 'MMP14', 'DACT1', 'MGC3047', 'THBS2')
# for ref30: stroma-related signatures
# for ref29: endothelial/microvasculature signatures
ref29 <- c('ARAP3','ADCY4','ESAM','ERG','SLC43A3','SOX7','PTPRB','PTPRM','AFAP1L1','MMRN2','TENC1','STARD9','COL4A3','LRRK1','PALD1','NPR3','ROBO4','NOTCH4','TIE1','RASIP1','ACVRL1','RAMP2','FAM110D','EGFL7','SMAD6','FGD5','ENG','CASKIN2','ACKR2','SLC9A3R2','CALCRL','HSPA12B','EPAS1','EHD4','LATS2','ICAM1','HBEGF','PLTP','C1orf54','CTTNBP2NL','MYO1B','SLCO2A1','KIFC1','EPHB4','SOX13','DRAM1','PECAM1','ENTPD1','ICAM2','CLDN5','SDPR','CDH5','GPR116','ELTD1','KDR','HILPDA','NPNT')

# for ref27: endothelial/microvasculature signatures; stroma-related signatures
library(stringr)
ref27_Endothelial_BRCA <- c('ARHGEF15 ARHGEF15 CD34 BCL6B CDH5 CALCRL CXorf36 CD34 ELTD1 CD93 ERG CDH5 ESAM CLEC14A LDB2 CXorf36 MMRN2 ELTD1 MYCT1 ERG TIE1')
ref27_Endothelial_BRCA <- as.character(str_split(ref27_Endothelial_BRCA,' ',simplify = T))

ref27_ECM_BRCA <- c('ADAM12 BNC2 CDH11 COL1A1 COL1A2 COL3A1 COL5A1 COL5A2 COL6A3 DACT1 FAP FBN1 GLT8D2 LUM POSTN SPARC THBS2 VCAN')
ref27_ECM_BRCA <- as.character(str_split(ref27_ECM_BRCA,' ',simplify = T))

# 不过有趣的是作者选取的 ref 27里面的出现在多种癌症的基因，所以是
# ECM是 COL1A1, COL1A2, and COL3A1
# endothelial gene sets 是CDH5, CXorf36, and TIE1



