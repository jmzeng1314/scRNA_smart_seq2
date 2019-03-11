######################################################################################
### 1. Quality control to remove cells ###
######################################################################################
# will read in qc data summary for SS2 data from the single cell platform in Stockholm and print a summary to one file
# and make a list of low quality cells.

# custom cutoffs cand be defined with appropriate flags.
# default cutoffs are defined as 2 standard deviations from the mean
# the exception is number of exonmapping reads, that has default cutoff at 10 000 reads.

# outprefix gives a predefined path/name to each of the output files. (Default is  "qc_filter")

# bulk_file should contain one column with all the sample names of samples that are not single cells, since these will have very different distribution in some of the quality metrics, they are exlcuded in the filtering

# example call:
# Rscript filter_low_qc_samples_SS.R -i SS2_15_0048_qc.txt,SS2_15_0049_qc.txt -o filter_2plates


##################################################
# read input arguments

#save the main path in Path_Main, or set the working directory to where the files are saved (be aware of the folder structure)
install.packages("optparse")
library(optparse)
option_list = list(
  make_option(c("-i","--infiles"), type="character", default=NULL, help="input qc summary files, separated by comma"),
  make_option(c("-o","--outprefix"), type="character", default="qc_filter", help="output prefix, [default qc_filter]"),
  make_option("--bulk_file", type="character", default = NULL, help = "optional, list of samples to exclude from filtering"),
  make_option("--nfail", type="integer", default=2, help = "number of features a cell has to fail to be removed, [default %default]"),
  make_option("--nreads", type="integer", default=10000, help = "cutoff for number of reads mapped to exons, [default %default]"),
  make_option("--uniqmap", type="double", default=NA, help = "cutoff for percent uniquely mapping, [default mean-2*sd]"),
  make_option("--exonic", type="double", default=NA, help = "cutoff for percent exonic mapping, [default mean-2*sd]"),
  make_option("--genes_max", type="integer", default=NA, help = "cutoff for max number of detected genes, [default mean-2*sd]"),
  make_option("--genes_min", type="integer", default=NA, help = "cutoff for min number of detected genes, [default mean+2*sd]"),
  make_option("--correl", type="double", default=NA, help = "cutoff for maximum correlation, [default mean-2*sd]")
)
opt <- parse_args(OptionParser(option_list=option_list))

###############################################
# read each plate
("reading inputs")

infiles = unlist(strsplit(opt$infiles,","))
D<-list()
for (file in infiles){
  D[[file]]<-read.table(file, sep="\t",header=T)
}

D<-Reduce("rbind",D)
samples <- apply(sapply(D[,1:2],as.character),1,paste,collapse="__")

# remove bulk samples if any are present
if (!is.null(opt$bulk_file)){
  bulks<-read.table(bulk.file)
  bulk.idx<-na.omit(match(bulks[m,1],samples))
  D<-D[-bulk.idx,]
  samples<-samples[-bulk.idx]
}

nS<-nrow(D)

# calculate number of exon mapping reads:
D$exonreads = D$reads*D$unique/100*D$exonic/100

##################################################
# plotting and defining cutoffs

plot_types <- c("exonreads","unique","exonic","rpkm1","max.correlation")
cutoffs<-list(exonreads=opt$nreads,unique=opt$uniqmap,exonic=opt$exonic,rpkm1=c(opt$genes_min,opt$genes_max),max.correlation=opt$correl)

# plot overview of all the plates
pdf(paste(opt$outprefix,".plate_overview.pdf",sep=''))
par(mfrow=c(3,2))
for (type in plot_types) {
  data <- split(D[,type],D$experiment)
  boxplot(data,main=type)
}
dev.off()

get_mids<-function(x){
  nb<-length(x)
  rowMeans(cbind(x[-nb],x[-1]))
}

# plot histograms and count number of failed cells per feature
nbin=round(nS/10)
pdf(sprintf("%s.qc_overview.pdf",opt$outprefix))
filtered<-list()
par(mfrow=c(2,3),cex=0.6,oma=c(3,1,2,1))
for (type in plot_types) {
  data <- D[,type]
  m<-mean(na.omit(data))
  s<-sd(na.omit(data))
  if (is.na(cutoffs[[type]][1])) { cutoffs[[type]][1] = m-s*2 }
  filt <- which(data<cutoffs[[type]][1])
  cutname = sprintf("%.2f",cutoffs[[type]][1])
  if (length(cutoffs[[type]])==2) {
    if (is.na(cutoffs[[type]][2])) { cutoffs[[type]][2] = m+s*2 }
    filt <- c(filt,which(data>cutoffs[[type]][2]))
    cutname = sprintf("%.2f & %.2f",cutoffs[[type]][1],cutoffs[[type]][2])
  }
  h<-hist(data,n=nbin,plot=F)
  if (length(filt)>0){
    h1<-hist(data[-filt],breaks=h$breaks,plot=F)
    h2<-hist(data[filt],breaks=h$breaks,plot=F)
    #       b<-barplot(rbind(h1$counts,h2$counts),col=c("red","blue"),las=2,border=NA)
    b<-barplot(rbind(h1$counts,h2$counts),col=c("red","blue"),main=sprintf("%s \n%d cells failed, cutoff %s",type,length(filt),cutname),las=2,border=NA)
    br<-get_mids(h$breaks)
    axis(1,at=b,labels=br)
  }else {
    b<-barplot(h$counts,col="red",main=sprintf("%s \n%d cells failed, cutoff %s",type,length(filt),cutname),las=2,border=NA)
    br<-get_mids(h$breaks)
    axis(1,at=b,labels=br)	  
  }
  filtered[[type]]<-filt
}


# filter out cells with > nfail features too low/high
allfilt<-table(unlist(filtered))
remove <- as.numeric(names(allfilt)[which(allfilt>opt$nfail)])
nFilt<-length(remove)
mtext(sprintf("total %d samples failed QC\nTotal samples: %d",nFilt,nS),side=1,outer=T)

# include a scatterplot with the filtered cells colored and larger
col <- rep("black",nS)
col[remove]<-"blue"
cex <- rep(0.2,nS)
cex[remove]<-1

plot(D[,plot_types],col=col,pch=16,cex=cex)
dev.off()


##########################################
# write a table with failed cells
if (nFilt >0) {
  D2<-D[remove,plot_types]
  reason<-rep('',nFilt)
  for (i in  1:nFilt){
    mF <- unlist(lapply(filtered, function(x) is.element(remove[i],x)))
    reason[i]<-paste(plot_types[mF], collapse=":")
  }
  D2<-cbind(D2,reason)
  rownames(D2)<-samples[remove]
  
  printfile2<-sprintf("%s.filtered_cells.txt",opt$outprefix)
  write.table(D2,file=printfile2,sep="\t")
  (sprintf("%d samples in experiment, %d cells filtered - filtered cells written to %s/%s",nS,nFilt,getwd(),printfile2))
}