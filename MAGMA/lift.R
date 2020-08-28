# AnJa/MNT segment for lifting coordintes from hg38 > hg19
# 21Aug2020 === === ===

###
library(rtracklayer)
library(GenomicRanges)
BiocManager::install("liftOver")

gene_df = read.delim("GRCh38_Ensembl-93_GENES_all-33538.gene.loc",header=FALSE)
colnames(gene_df)= c("GeneID", "Chr", "Start", "End", "Strand", "Symbol")

# add chr
gene_df$Chr = paste0("chr", gene_df$Chr)
gene_df$Chr[gene_df$Chr == "chrMT"] = "chrM"
gr = makeGRangesFromDataFrame(gene_df,keep=TRUE)
names(gr) = gr$GeneID

## lift
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
lifted_list = range(liftOver(gr, ch))
    table(lengths(lifted_list))
        #  0     1     2     3     4     5     7    11    14    16    21    31
        # 76 33339    99     9     6     2     2     1     1     1     1     1
lifted_list = lifted_list[lengths(lifted_list) == 1]

lifted = unlist(lifted_list)
lifted_df = as.data.frame(lifted)
lifted_df$Symbol = gene_df$Symbol[match(rownames(lifted_df), gene_df$GeneID)]

    # Check:
    lifted_df[lifted_df$Symbol=="GABRQ", ]
        # seqnames     start       end width strand Symbol
        # ENSG00000268089     chrX 151806637 151825999 19363      +  GABRQ
        #   ^ good - in hg19, has a different EnsemblID ("ENSG00000147402")



## Rearrange cols of lifted_df / match format ===
table(rownames(lifted_df) %in% expressedGenes)

table(expressedGenes %in% rownames(lifted_df))
    # FALSE  TRUE
    #   148 30104     ~30,104 still >> than the EnsemblID-matching-approach that retains only 28,462

head(lifted_df)
    #                 seqnames  start    end width strand      Symbol
    # ENSG00000243485     chr1  29554  31109  1556      + MIR1302-2HG
    # ENSG00000237613     chr1  34554  36081  1528      -     FAM138A
    # ENSG00000186092     chr1  65419  71585  6167      +       OR4F5
    # ENSG00000238009     chr1  89295 133723 44429      -  AL627309.1
    # ENSG00000239945     chr1  89551  91105  1555      -  AL627309.3
    # ENSG00000239906     chr1 139790 140339   550      -  AL627309.2

    head(gene_df, n=3)
        #            GeneID  Chr  Start    End Strand      Symbol
        # 1 ENSG00000243485 chr1  29554  31109      + MIR1302-2HG
        # 2 ENSG00000237613 chr1  34554  36081      -     FAM138A
        # 3 ENSG00000186092 chr1  65419  71585      +       OR4F5     <- match this format for MAGMA input

lifted_df$width <- NULL
colnames(lifted_df) <- c("Chr", "Start", "End", "Strand", "Symbol")
    # oh, this was not necessary (won't save colnames)
lifted_df$GeneID <- rownames(lifted_df)
lifted_df <- lifted_df[ ,c(6,1:5)]

# Remove the 'chr' prefix
lifted_df$Chr <- ss(as.character(lifted_df$Chr),"chr",2)

# Finally subset for those 30,104 expressing-and-conserved(?) genes
    #lifted_df <- lifted_df[expressedGenes[expressedGenes %in% rownames(lifted_df)], ]
    # do this way (the above wildly reorders the genes in not-seqlevel-order):
lifted_df <- lifted_df[rownames(lifted_df) %in% expressedGenes, ]


## Write out for MAGMA
write.table(lifted_df, file="./GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc", sep="\t",
            row.names=F, col.names=F, quote=F)









