### MNT add 07Apr2021 =======
  # Check Read 1 files in expansion dataset for discrepant read lengths
  #     (as seen with Br5182-NAc in preprint dataset)
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)

FASTQ.dir <- "/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/"

### Read in (2021) 'samples.manifest'
samples.rev <- read.table("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Feb2021/samples.manifest",
                               sep="\t", header=F)$V1

R1files <- data.frame(
  sampleName = c(sapply(samples.rev, function(x){
    rep(x,length(list.files(paste0(FASTQ.dir, x), pattern="R1")))})),
  
  R1 = c(sapply(samples.rev,function(x){list.files(paste0(FASTQ.dir,x),
                                           pattern="R1")}))
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
