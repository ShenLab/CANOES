# _CANOES_

### Detecting CNVs from exome sequencing data

#### 

<center>[Tutorial](#Tutorial)   |   [Contact](#Contact)   |   [Cite](#Cite)</center>

<a name="Tutorial"></a>

#### Tutorial

To run CANOES, you will need

1.  bedtools, available from [http://bedtools.readthedocs.org/en/latest/](http://bedtools.readthedocs.org/en/latest/)
2.  R, with packages nnls, Hmisc, mgcv and plyr installed. You can install these by executing the following command: `install.packages(c("nnls", "Hmisc", "mgcv", "plyr"))`.
3.  GATK
4.  [CANOES.R](https://github.com/ShenLab/CANOES)

CANOES requires a data frame with the coordinates, GC content and read count per sample for each exome capture region.

If bamlist.txt is a file with the path of a bam on each line and probes.bed is a bed file with the coordinates of the exome capture regions, bedtools can generate read counts with the following command, entered at the LINUX prompt:

    > bedtools multicov -bams `cat bamlist.txt` -bed probes.bed -q 20 > canoes.reads.txt

Make sure the contig names in probes.bed correspond to those in the header of the bam files. The above command will save the read counts in the file canoes.reads.txt.

If FASTA is the path to your reference file and TEMP is a temporary directory, enter the following command at the LINUX prompt to use GATK to calculate the GC content for each exome capture region:

    > java -Xmx2000m -Djava.io.tmpdir=TEMP -jar ./Sting/dist/GenomeAnalysisTK.jar -T GCContentByInterval -L probes.bed -R FASTA -o gc.txt

Again, the contig names in probes.bed should correspond to those in the FASTA.fai file.

We can then use R to stitch together a data frame with all the required information for CANOES. You can follow along with the two sample files [gc.txt](https://github.com/ShenLab/CANOES/blob/master/gc.txt) and [canoes.reads.txt](https://github.com/ShenLab/CANOES/blob/master/canoes.reads.txt), which contain gc content and read count information for probes on chromosome 22 for 26 samples.

First, read in the data:

    > gc <- read.table("gc.txt")$V2
    > canoes.reads <- read.table("canoes.reads.txt")

Now, rename the columns of canoes.reads. Normally, one would use actual sample names instead of S1, S2, …

    > sample.names <- paste("S", seq(1:26), sep="")
    > names(canoes.reads) <- c("chromosome", "start", "end", sample.names)

Create a vector of consecutive target ids:

    > target <- seq(1, nrow(canoes.reads))

Combine the data into one data frame:

    > canoes.reads <- cbind(target, gc, canoes.reads)

Now we can call CNVs in the samples.

    # execute this command from the directory where you have saved CANOES.R
    # first make sure you have the packages nnls, Hmisc, mgcv and plyr installed 
    # you can install these with the command install.packages(c("nnls", "Hmisc", "mgcv", "plyr"))
    > source("CANOES.R")

    # create a vector to hold the results for each sample
    > xcnv.list <- vector('list', length(sample.names))

    # call CNVs in each sample
    > for (i in 1:length(sample.names)){
        xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)
      }

    # combine the results into one data frame
    >  xcnvs <- do.call('rbind', xcnv.list)

    # inspect the first two CNV calls
    > head(xcnvs, 2)
      SAMPLE CNV             INTERVAL     KB CHR   MID_BP    TARGETS NUM_TARG MLCN Q_SOME
    1     S2 DEL 22:25713988-25756059 42.071  22 25735024 1132..1137        6    1     99
    2     S3 DEL 22:24373138-24384231 11.093  22 24378684   936..942        7    0     44

Now we can plot all the CNV calls.

    > pdf("CNVplots.pdf")
    > for (i in 1:nrow(xcnvs)){
        PlotCNV(canoes.reads, xcnvs[i, "SAMPLE"], xcnvs[i, "TARGETS"])
      }
    > dev.off()

Example output: [CNVplots.pdf](https://github.com/ShenLab/CANOES/blob/master/CNVplots.pdf).

Finally, we can genotype a sample to determine the evidence that the sample has a CNV at any interval.

    # genotype all the CNV calls made above in sample S2
    > genotyping.S2 <- GenotypeCNVs(xcnvs, "S2", canoes.reads)
    # inspect the genotype scores for the first two CNV calls
    > head(genotyping.S2, 2)
                                     INTERVAL NQDel SQDel NQDup SQDup
    22:25713988-25756059 22:25713988-25756059     0    99    99     0
    22:24373138-24384231 22:24373138-24384231    86     0    61     0

The SQDel score for the first CNV is high (99) because this CNV was called in sample S2\. The second CNV was called in sample S3, so the SQDel and SQDup quality scores for sample S2 are 0.

* * *

<a name="Contact"></a>

#### Contact

Daniel Backenroth (db2175 at columbia.edu) or Yufeng Shen ( ys2411 at columbia.edu)

* * *

<a name="Cite"></a>

#### How to cite:

Backenroth D, Homsy J, Murillo LR, Glessner J, Lin E, Brueckner M, Lifton R, Goldmuntz E, Chung WK, Shen Y, (2014) CANOES: Detecting rare copy number variants from whole exome sequencing data, Nucleic Acids Research,  doi: 10.1093/nar/gku345 [PMID: 24771342](http://www.ncbi.nlm.nih.gov/pubmed/24771342)

* * *

#### [Shen Lab](http://www.columbia.edu/~ys2411/)


