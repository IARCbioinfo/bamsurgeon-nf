#! /usr/bin/Rscript

# this script build a varfile folder based on a bam folder
# we will give each bam file a var-twin
# a varfile name corresponds to associated bam file name
# the resulting var directory can be used in addition to the bam folder to simulate SNV with bamsurgeon

####################################  ARGUMENTS SECTION  ########################################
## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
args <- argsL
samtools="samtools"

## Set an allelic fractions interval if not provided
if(is.null(args$min_true_af)) {args$min_true_af=0.0001} else {args$min_true_af=as.numeric(args$min_true_af)}
if(is.null(args$max_true_af)) {args$max_true_af=1} else {args$max_true_af=as.numeric(args$max_true_af)}
if(is.null(args$n_mut)) {args$n_mut=100} else {args$n_mut=as.numeric(args$n_mut)}
if(is.null(args$hotspot_size)) {args$hotspot_size=25} else {args$hotspot_size=as.numeric(args$hotspot_size)}
if(is.null(args$fasta_ref)) args$fasta_ref = "/appli57/reference/hg19_torrentserver.fasta"
if(is.null(args$reads_thr)) args$reads_thr=5
if(is.null(args$del)) {args$del=FALSE} else {args$del=TRUE} 
if(is.null(args$ins)) {args$ins=FALSE} else {args$ins=TRUE}
if(is.null(args$indel_size)) {args$indel_size=1} else {args$indel_size=as.numeric(args$indel_size)}
if(args$del & args$ins) {cat("ERROR: CANNOT GENERATE INSERTION AND DELETION, NEED TO CHOOSE BEWEEN BOTH\n"); q(save="no")}

## Default setting when no all arguments passed or help needed
if("--help" %in% args | is.null(args$bam_folder) | is.null(args$bed_file)) {
  cat("
      The R Script generate_varfiles.R
      
      Mandatory arguments:
      --bam_folder=somePath          - the explicit path of the bam folder
      --bed_file=somePath            - the explicit path of the target bed file
      --help                         - print this text
      
      Optionnal arguments:
      --fasta_ref=somePath                - the explicit path of the reference fasta file
      --min_true_af=someValue             - the minimum possible value for allelic fraction (default=0.0001)
      --max_true_af=someValue             - the maximum possible value for allelic fraction (default=1)
      --n_mut=someValue                   - the number of mutated bases (default=100)
      --hotspot_size=someValue            - the number of mutated bam for each mutation (default=25)
      --af=someValue                      - the observed allelic fraction (default:10^-(unif(0,4)))
      --reads_thr=someValue               - the value of threshold number of mutated reads to consider the mutation (default=5)
      --hotspot_size=value                - the number of individual being mutated for a same mutation (default=25)
      --del                               - specify that user want to simulate deletions
      --ins                               - specify that user want to simulate insertions
      --del_size=Value                    - size of the deletions
      --ins_size=Value                    - size of the insertions
      
      WARNING : samtools and bedtools have to be in your bin
      the bed must don't have a declaration line and must be sorted and merged
      
      Example:
      ./generate_varfiles_hotspot.R --bam_folder=~/Documents/BAM --var_result_folder=~/Documents/VAR --bed_file myfile.bed \n\n")
  
  q(save="no")
}

#############################  FUNCTIONS  ##############################

bed_to_positions <- function(bed_file, fasta_ref){
  positions=c()
  chrs=c()
  refs=c()
  con <- file(bed_file, open="r")
  while (length(l <- readLines(con, n = 1))) {
    start=as.numeric(strsplit(l,split="\t")[[1]][2])
    end=as.numeric(strsplit(l,split="\t")[[1]][3])
    positions <- c(positions,start:end)
    chrs <- c(chrs,rep(as.character(strsplit(l,split=" ")[[1]][1]),(end-start+1)))
    baseConn = pipe(paste(samtools," faidx ",fasta_ref," ","chr17",":",start,"-",end,sep=""))
    all_bases = readLines(baseConn)[2:length(readLines(baseConn))]
    all_bases2 <- c()
    for(i in 1:length(all_bases)){all_bases2 <- paste(all_bases2, all_bases[i], sep="")}
    refs = paste(refs,all_bases2,sep="")
    close(baseConn)
  }
  close(con)
  refs = unlist(strsplit(refs,""))
  return(list(positions=positions,refs=refs))
}

to_nbReads <- function(indiv_bam, bed){
  baseConn = pipe(paste("cat ",bed," | awk '{print $1\"\t\"$2-1\"\t\"$3}' | ", samtools, " depth -Q 20 -q 20 -d 1000000 -a -b - ",indiv_bam,sep=""))
  cur_base = readLines(baseConn)
  close(baseConn)
  unlist(lapply(cur_base, function(x) as.numeric(strsplit(x,"\t")[[1]][3]))) 
}

#############################   MAIN   #################################
library(snow)

bamnames <- list.files(path=args$bam_folder, pattern=".bam$")
for(bam in bamnames){ assign(paste(bam,".var",sep=""), file(paste(bam,".var",sep=""), open="w"))
                     closeAllConnections() }
bamnames <- bamnames[ !grepl("_2",bamnames) ]; bamnames=gsub(".bam","",bamnames)
bed_info <- bed_to_positions(args$bed_file, args$fasta_ref)
# construction of dataframes of depth (indiv in rows, position in columns)
cl <- makeCluster(10)
clusterExport(cl,ls(), envir=environment()) 
invisible(clusterCall(cl, function() library(plyr)))
DP1_list <- parLapply(cl,paste(args$bam_folder,bamnames,".bam",sep=""),to_nbReads,args$bed_file)
DP1 = data.frame(matrix(unlist(DP1_list), nrow=length(bamnames), byrow=T))
colnames(DP1)=1:length(bed_info$positions);rownames(DP1)=bamnames
DP2_list <- parLapply(cl,paste(args$bam_folder,bamnames,"_2.bam",sep=""),to_nbReads,args$bed_file)
DP2 = data.frame(matrix(unlist(DP2_list), nrow=length(bamnames), byrow=T))
colnames(DP2)=1:length(bed_info$positions);rownames(DP2)=paste(bamnames,"_2")
stopCluster(cl)
closeAllConnections()
# choose the set of indiv to make the hotspot
seen_positions=c(0);position=0
for (mut_id in 1:args$n_mut){
  print(mut_id)
  while(position %in% seen_positions ) {
    if(is.null(args$af)){
      obs_af = 10^(-runif(1,0,4))
      pos_ok=c()
      while(obs_af < args$min_true_af | obs_af > args$max_true_af | length(pos_ok)<1) {
        obs_af = 10^(-runif(1,0,4))
        pos_ok = which(as.numeric(colSums(pmin(DP1,DP2)*obs_af >= args$reads_thr)) >= args$hotspot_size)
      }
    } else {obs_af=args$af}
    position=sample(pos_ok, 1)  
  }
  seen_positions = c(seen_positions, position)
  panel= sample(rownames(DP1)[which((pmin(DP1,DP2)*obs_af)[,position] >= args$reads_thr)],args$hotspot_size)
  if (!args$del & !args$ins) alternative = sample(c("A","T","C","G")[c("A","T","C","G") != bed_info$refs[position]], 1)
  if (args$ins) alternative = sample(c("A","T","C","G"), args$indel_size)
  if (args$del) {
    alternative = ""
    if(is.na(bed_info$positions[position + args$indel_size])) position = position - args$indel_size
    positionS = position
    positionE = position + args$indel_size
  }
  if (!args$del) positionS = positionE = position

  for (indiv in panel) {
    cat("chr17","\t",bed_info$positions[positionS],"\t",bed_info$positions[positionE],"\t",obs_af,"\t", ifelse(args$ins,"INS\t",""), ifelse(args$del,"DEL\t",""), alternative, 
        "\n",sep="",file = paste(indiv,".bam.var",sep=""),append=T)
    cat("chr17","\t",bed_info$positions[positionS],"\t",bed_info$positions[positionE],"\t",obs_af,"\t", ifelse(args$ins,"INS\t",""), ifelse(args$del,"DEL\t",""), alternative, 
        "\n",sep="",file = paste(indiv,"_2.bam.var",sep=""),append=T)
   }
}
