source('functions.R')

library(Biostrings)

processMotifCompare <- function(fileName, pval=0.001) {
    f=read.delim(fileName,header=F,sep='')
    colnames(f)=c("chr","start","stop","strand","score","p_value","motif")
    x=processSingle(f, pval)
    x=stem(x)
    x=subset(x, motif %in% shared$stemmed_motif)
    x=merge(x,shared, by.x='motif', by.y='stemmed_motif')
    x=remove_redundant_TFclass(x)
    out=x[order(x$start),]
    return(out)

}

remove_redundant_TFclass <- function(x) {
    require('GenomicRanges')

    x=subset(x, ! motif=='dl')
    # remove if same motif name and same start position
    x=x[order(x$score),]
    x=x[order(x$start),]
    x=x[!duplicated(x[c('start','class')]),]
    x$id=seq(1:nrow(x))
    # remove duplicate based on overlapping region with and same TF class
    dup=x
    x_gr=with(x, GRanges( chr , IRanges( start , stop) ))
    dup_gr=with(dup, GRanges( chr , IRanges( start , stop) ))
    y=as.data.frame(findOverlaps(x_gr, dup_gr))
    y=subset(y, ! queryHits==subjectHits)
    combine=cbind( x[y$queryHits,], x[y$subjectHits,] )
    colnames(combine)=c('motif','chr', 'start', 'stop', 'strand','score', 'p_value','class','code','id',
    'xmotif','xchr', 'xstart', 'xstop', 'xstrand','xscore', 'xp_value','xclass','xcode','xid')
    combine=subset(combine, class==xclass)
    combine$notkeepid=ifelse(combine$xscore > combine$score, combine$id, combine$xid) 
    x=subset(x, ! id %in% combine$notkeepid)
    return(x)
}



compare_order_rand <- function( seq1, seq2, iter=100) {
    res=c()
    for (i in 1:100){
        samp = paste(sample(as.vector(seq2$code),size=nchar(s2),replace=FALSE),collapse = '')
        o1 = pairwiseAlignment(s1, samp, type="local", scoreOnly = TRUE,gapOpening=5, gapExtension=2)
        o2 = pairwiseAlignment(s1, reverse(samp), type="local", scoreOnly = TRUE,gapOpening=5, gapExtension=2)
        s=max(o1,o2)
        res=c(res, s)
    }
    return(res)

}

#------------------------------#
# Main run
input: motif annotation files
Example compared at TF class
#------------------------------#
iter=100
shared=read.delim('data/eisl_species_sharedmotifs',sep='')
seq1 = processMotifCompare('data/eISLh') 
seq2 = processMotifCompare('data/eISLs') 
s1=paste(as.vector(seq1$code),collapse = '')
s2=paste(as.vector(seq2$code),collapse = '')

#------------------------------#
# Compare forward and reverse
#------------------------------#
fa=pairwiseAlignment(s1, s2, type="local",gapOpening=5, gapExtension=2)
ra=pairwiseAlignment(s1, reverse(s2), type="local",gapOpening=5, gapExtension=2)
max_score=max(fa@score, ra@score)

#------------------------------#
# Scramble motifs
#------------------------------#
out=compare_order_rand(seq1, seq2, iter)
length(out[out>max_score ])/iter 