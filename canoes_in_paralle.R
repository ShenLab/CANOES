#############################################################
# Author: Renjie Tan                                        #
# Date  : 4/18/2022                                         #
# This script is supplied as part of the CANOES pipeline.   #
# It has been modified to take in parameters from the shell # 
# script. It has been made more readable with additional    # 
# comments. It also sources the R script CANOES.R from the  #
# rscripts directory.                                       #
#############################################################

############################################################################
# STEP 1: Fetch the directory locations and file names from the shell script
############################################################################
cat(paste('[',format(Sys.time()),']'," Step1: Fetch the directory locations and file names...\n"))
#install.packages(c("nnls", "Hmisc", "mgcv", "plyr"))
#install.packages(c("doParallel", "foreach", "Rfast" ))

library(doParallel)
library(foreach)
library(Rfast)
library(plyr)

args        <- commandArgs(trailingOnly=TRUE)
scripts_dir <- args[1]
batch_name  <- args[2]
output_dir  <- args[3]
cores       <- as.numeric(args[4])
numrefs     <- as.numeric(args[5])

if(is.na(numrefs)){
    numrefs <- 30
}

if(is.na(cores)){
    cores <- detectCores(logical=F)
    cores <- floor(cores/8)-2
}
print(paste("Script dir:", scripts_dir))
print(paste("Batch name:", batch_name))
print(paste("Output dir:", output_dir))
print(paste("Cores:", cores))
print(paste("Number of reference samples:", numrefs))


setwd(scripts_dir)
source("CANOES.R")
setwd(output_dir)

input_file_name <- paste0(output_dir, '02_merged_rc.csv_TR.bz2') 
output_file_name <- paste0(output_dir,'canoes_calls_',batch_name,'.csv')

raw_cnv_calls_dir <- paste0(output_dir,'/04_raw_cnv_calls/')
log_dir <- paste0(raw_cnv_calls_dir,'/log/')

if(!dir.exists(raw_cnv_calls_dir)){
    cat(paste0('[',format(Sys.time()),']',"\t\tCreat directory: ",raw_cnv_calls_dir,"\n"))
    dir.create(raw_cnv_calls_dir)
}

if(!dir.exists(log_dir)){
    cat(paste0('[',format(Sys.time()),']',"\t\tCreat directory: ",log_dir,"\n"))
	dir.create(log_dir)
}

######################################################
# STEP 2: Import read count file
######################################################
cat(paste('[',format(Sys.time()),']'," Step2: Import read count file ...\n"))

# normalize counts; round so we can use negative binomial
if (file.exists('normalized_counts.RData')){
    cat(paste0('[',format(Sys.time()),']',"\t\tLoading normalized read counts data frame ...\n"))
    load(file=paste('normalized_counts.RData'))
}else{
    # Load/Read the raw read count file
    if (file.exists('canoes.reads.RData')){
        cat(paste0('[',format(Sys.time()),']',"\t\tLoading canoes.reads file\n"))
        load(file=paste('canoes.reads.RData'))
    }else{
        cat(paste0('[',format(Sys.time()),']',"\t\tImporting ",input_file_name,"\n"))
        canoes.reads <- read.table(input_file_name, header=TRUE)

        # Import GC content file
        cat(paste0('[',format(Sys.time()),']',"\t\tLoading gc.txt\n"))
        gc <- read.table("03_gc.txt")$V2

        num_cols <- dim(canoes.reads)[2]
        num_samples <- num_cols - 3
        colnames(canoes.reads)[1]<-"chromosome"
        target <- seq(1, nrow(canoes.reads))
        canoes.reads <- cbind(target, canoes.reads[,1:3], gc, canoes.reads[,-c(1,2,3)])
    }

    ######################################################################################################
    # Note: we moved the following scripts outside of CallCNVs function in CANOES.R.
    if(class(canoes.reads$chromosome) != "numeric"){
        canoes.reads$chromosome <- as.character(canoes.reads$chromosome) #RT modifided
    }

    if (length(setdiff(names(canoes.reads)[1:5], c("target", "chromosome", "start", "end", "gc"))) > 0){
        stop("First five columns of counts matrix must be target, chromosome, start, end, gc")
    }
    if (length(setdiff(unique(canoes.reads$chromosome), seq(1:22))) > 0) {
        # remove sex chromosomes
        cat(paste0('[',format(Sys.time()),']',"\t\tTrying to remove sex chromosomes and 'chr' prefixes\n"))
        canoes.reads <- subset(canoes.reads, !chromosome %in% c("chrX", "chrY", "X", "Y"))
        if (sum(grepl("chr", canoes.reads$chromosome))==length(canoes.reads$chromosome)){
          canoes.reads$chromosome <- gsub("chr", "", canoes.reads$chromosome)
        }
        canoes.reads$chromosome <- as.numeric(canoes.reads$chromosome)
        if (length(setdiff(unique(canoes.reads$chromosome), seq(1:22))) > 0) 
          stop("chromosome must take value in range 1-22 (support for sex chromosomes to come)")
    }

    if(class(canoes.reads$chromosome) != "numeric"){
        canoes.reads$chromosome <- as.numeric(as.character(canoes.reads$chromosome)) #RT modifided
    }
    ######################################################################################################

    canoes.reads <- arrange(canoes.reads, chromosome, start)
    sample.names <- colnames(canoes.reads)[-seq(1,5)]
    # find mean coverage of probes
    cat(paste0('[',format(Sys.time()),']',"\t\tCalculate the overall mean coverage of the matrix (all probes across all the samples)...\n"))
    mean.counts <- mean(apply(canoes.reads[, sample.names], 2, mean))
    cat(paste0('mean.counts: ', mean.counts, "\n"))

    #Read the read count file
    cat(paste0('[',format(Sys.time()),']',"\t\tNormalize counts; round so we can use negative binomial...\n"))

    normalized_counts <- apply(canoes.reads[, sample.names], 2,
        function(x, mean.counts)
                 round(x * mean.counts / mean(x)), mean.counts)
    normalized_counts <- cbind(canoes.reads[,1:5],normalized_counts)
    cat(paste0('[',format(Sys.time()),']',"\t\tSave the canoes.reads counts to file 'canoes.reads.RData' ...\n"))
    save(canoes.reads,file='canoes.reads.RData')
    cat(paste0('[',format(Sys.time()),']',"\t\tSave the normalize counts to file 'normalized_counts.RData' ...\n"))
    save(normalized_counts, file='normalized_counts.RData')
    #rm(canoes.reads)
}

######################################################
# STEP 3: Calculate covariance of read count across samples 
## Note: we removed 'cov' variable to make sure the this variable is shared for all samples.
######################################################
cat(paste('[',format(Sys.time()),']'," Step3: Calculate covariance of read count across samples ...\n"))
if (file.exists(paste('cov.RData'))){
    cat(paste0('[',format(Sys.time()),']',"\t\tLoading covariance file\n"))
    load(file=paste('cov.RData'))
}
if (!exists('covariance')){
    if(!exists('canoes.reads')){
        cat(paste0('[',format(Sys.time()),']',"\t\tLoading file 'canoes.reads.RData' ...\n"))
        load(file=paste0('canoes.reads.RData'))
    }
    if(!exists('sample.names')){
        sample.names <- colnames(canoes.reads)[-seq(1,5)]
    }
    cat(paste0('[',format(Sys.time()),']',"\t\tCalculate covariance of read count across samples by 'cora' (RT modified version)\n"))
    system.time(covariance <- cora(as.matrix(canoes.reads[, sample.names])))

    cat(paste0('[',format(Sys.time()),']',"\t\tSave cov matrix to file ...\n"))
    save(covariance,file=paste('cov.RData'))
}


######################################################
# STEP 4: Run CANOES in parallel 
######################################################
cat(paste('[',format(Sys.time()),']'," Step4: Run CANOES in parallel ...\n"))

clus <- makeCluster(cores)
registerDoParallel(clus, cores=cores)
setwd(output_dir)

sample.names <- colnames(normalized_counts)[-seq(1,5)]
xcnv.list <- vector('list', length(sample.names))

chunk.size <- length(sample.names)/cores
cat(paste0('[',format(Sys.time()),']', " \t\tFor batch ",batch_name, ", start to run CANOES in ", as.character(cores), " cores ...\n"))

foreach(i=1:cores, .combine='rbind') %dopar%{
    log_file <- paste0(log_dir,'/',Sys.info()["nodename"],'_core_',i,'.log')
    write.table(paste0('[',format(Sys.time()),'] cores:', i,' starts ....'), file=log_file, col.names=FALSE, append=FALSE)
    for (x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        if (!file.exists(paste0(raw_cnv_calls_dir,sample.names[x]))){
            write.table(paste0('[',format(Sys.time()),'] [',sample.names[x],'] Touch output file and select reference samples...'), file=log_file, col.names=FALSE, append=TRUE)
            file.create(paste0(raw_cnv_calls_dir,sample.names[x]))
            reference.samples <- setdiff(sample.names, sample.names[x])
            covariances <- covariance[sample.names[x], reference.samples]
            reference.samples <- names(sort(covariances,decreasing=T)[1:min(numrefs, length(covariances))])
            write.table(paste0('[',format(Sys.time()),'] [',sample.names[x],'] Start to call CNVs ...'), file=log_file, col.names=FALSE, append=TRUE)
            xcnv.list[[x]] <- CallCNVs(sample.names[x], reference.samples, normalized_counts, log_file)
            write.table(xcnv.list[[x]], file=paste0(raw_cnv_calls_dir,sample.names[x]), col.names=FALSE, append=FALSE)
            gc()
        }
    }
}
stopImplicitCluster()
stopCluster(clus)

######################################################
# STEP 5: Double check to avoid missing samples 
# Go through all the samples with single thread and pick up the missing ones
######################################################
cat(paste('[',format(Sys.time()),']'," Step5: Double check to avoid missing samples ...\n"))

for (x in 1:length(sample.names)){
    log_file <- paste0(log_dir,'/',Sys.info()["nodename"],'.log')
    if (!file.exists(paste0(raw_cnv_calls_dir,sample.names[x]))){
        file.create(paste0(raw_cnv_calls_dir,sample.names[x]))
        reference.samples <- setdiff(sample.names, sample.names[x])
        covariances <- covariance[sample.names[x], reference.samples]
        reference.samples <- names(sort(covariances,decreasing=T)[1:min(numrefs, length(covariances))])
        print(paste0('[',format(Sys.time()),'] [',sample.names[x],'] Start to call CNVs ...'))
        xcnv.list[[x]] <- CallCNVs(sample.names[x], reference.samples, normalized_counts, log_file)
        write.table(xcnv.list[[x]], file=paste0(raw_cnv_calls_dir,sample.names[x]), col.names=FALSE, append=FALSE)
    }else{
        if(file.size(paste0(raw_cnv_calls_dir,sample.names[x]))>0){
            xcnv.list[[x]] <- read.table(file=paste0(raw_cnv_calls_dir,sample.names[x]), sep=" ")
            xcnv.list[[x]] <- xcnv.list[[x]][-1]
        }
    }

}

################################################################################
# STEP 6: Append output as rows onto a single variable and write the output file
################################################################################
cat(paste('[',format(Sys.time()),']'," Step6: Append output as rows onto a single variable and write the output file ...\n"))

xcnvs <- do.call('rbind', xcnv.list)
colnames(xcnvs) <- c("SAMPLE", "CNV", "INTERVAL", "KB", "CHR", "MID_BP", "TARGETS","NUM_TARG", "MLCN", "SQ", "NQ")
write.csv(xcnvs, file=output_file_name, row.names = FALSE)
