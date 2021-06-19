### MNT add 18Jun2021 =======
  # Check Read 1 files for discrepant read lengths, as seen with the original Br5182-NAc:
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)

FASTQ.dir <- "/dcl01/ajaffe/data/lab/singleCell/10x_pilot/FASTQ/Br5182_NAc_reseq/"
files.dirs <- list.files(FASTQ.dir, pattern="^S")

R1files <- data.frame(
  dirName = files.dirs,
  R1 = unlist(sapply(files.dirs,function(x){list.files(paste0(FASTQ.dir,x),
                                                   pattern="R1")}), use.names=F)
)

for(i in 1:nrow(R1files)){
  cat(paste0("Checking R1 length distribution for: ", R1files[i,2], "\n"))
  temp.R1s <- readFastq(paste0(FASTQ.dir, R1files[i,"dirName"], "/", R1files[i,"R1"]),
                        withIds=F)
  print(head(sread(temp.R1s), n=4))
  print(table(width(sread(temp.R1s))))
  rm(temp.R1s)
  cat("\n\n")
}

rm(list=ls())
sessionInfo()