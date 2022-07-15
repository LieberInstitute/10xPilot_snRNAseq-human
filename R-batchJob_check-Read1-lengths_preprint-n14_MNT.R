### MNT add 07Apr2021 =======
  # Check Read 1 files for discrepant read lengths, as seen with Br5182-NAc:
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)

FASTQ.dir <- "/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_Tran2021_published/"

### Read in preprint 'samples.manifest.full'
samples.prepr <- read.table("/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-10-15_Tran2021_published/samples.manifest.full",
                            sep="\t", header=F)$V5

# Drop Br5287-DLPFC (poor quality sample; dropped for preprint) and the test sucrose samples
samples.prepr <- samples.prepr[-c(grep("Br5287_DLPFC", samples.prepr),
                                  grep("_suc", samples.prepr))]

R1files <- data.frame(
  sampleName = unlist(sapply(samples.prepr, function(x){
    rep(x,length(list.files(paste0(FASTQ.dir, x), pattern="R1")))}), use.names=F),
  
  R1 = unlist(sapply(samples.prepr,function(x){list.files(paste0(FASTQ.dir,x),
                                                   pattern="R1")}), use.names=F)
)

for(i in 1:nrow(R1files)){
  cat(paste0("Checking R1 length distribution for: ", R1files[i,2], "\n"))
  temp.R1s <- readFastq(paste0(FASTQ.dir, R1files[i,1], "/", R1files[i,2]),
                        withIds=F)
  print(head(sread(temp.R1s), n=4))
  print(table(width(sread(temp.R1s))))
  rm(temp.R1s)
  cat("\n\n")
}
sessionInfo()
