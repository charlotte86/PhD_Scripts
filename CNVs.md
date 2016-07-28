# PhD_Scripts
PhD Scripts

### this is to be used on the cluster:
/share/apps/R/bin/R
library(ExomeDepth, '~/R/x86_64-unknown-linux-gnu-library/3.1')


## Load databases

data(exons.hg19)
data(genes.hg19)
data(ExomeCount)
data(Conrad.hg19)

print(head(genes.hg19))

## Load bam files and their indexes (bai)

bam.files <- list.files('/cluster/project9/ARVC/bams/prep_7/test/', pattern = ".bam$", full.names = TRUE)

index.files <- list.files('/cluster/project9/ARVC/bams/prep_7/test/', pattern = ".bai$", full.names = TRUE)


## count the reads per exon

my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = bam.files,
                          index.files = index.files,
                          include.chr = FALSE)






ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')



ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space),
                                   pattern = 'chr',
                                   replacement = '')





#### get the annotation datasets to be used later
data(Conrad.hg19)
head(Conrad.hg19.common.CNVs)

ExomeCount.dafr <- AnnotateExtra(x = ExomeCount.dafr,
                           reference.annotation = Conrad.hg19.common.CNVs,
                           min.overlap = 0.5,
                           column.name = 'Conrad.hg19')

exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
                                             IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
                                             names = exons.hg19$name)
                                             
                                             
                                                                                          

### prepare the main matrix of read count data

ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = 'EXOME*')])

nsamples <- ncol(ExomeCount.mat)



### start looping over each sample
for (i in 1:nsamples) {
  #### Create the aggregate reference set for this sample
  my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                     reference.counts = ExomeCount.mat[,-i],
                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                     n.bins.reduced = 10000)
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth',
                    test = ExomeCount.mat[,i],
                    reference = my.reference.selected,
                    formula = 'cbind(test, reference) ~ 1')
                    
  # Now call the CNVs
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = ExomeCount.dafr$space,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = ExomeCount.dafr$names)
  # Now annotate the ExomeDepth object
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = Conrad.hg19.common.CNVs,
                             min.overlap = 0.5,
                             column.name = 'Conrad.hg19')
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.hg19.GRanges,
                             min.overlap = 0.0001,
                             column.name = 'exons.hg19')
  
  
  ##output.file <- paste('Exome_', i, 'csv', sep = '')
  
  names <- names(ExomeCount.dafr[6:ncol(ExomeCount.dafr)])
  
  names <- gsub(names, pattern = ".bam", replacement="")
  
  output.file <- paste(names[i], 'CNVcalls.csv', sep='')
  
  dat <- cbind(names[i], all.exons@CNV.calls)
  
  output.file <- paste(names[i], 'CNVcalls.csv', sep='')
  
  write.table(dat, "/cluster/project9/ARVC/CNV_test.csv", row.names = FALSE, quote=F, sep= ",", append=T, col.names=F)
}


### trying to fix the annotations

exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19$chromosome,
                                             IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
                                             names = exons.hg19$name)
all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = exons.hg19.GRanges,
                           min.overlap = 0.0001,
                           column.name = 'exons.hg19')
