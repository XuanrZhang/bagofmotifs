library(tm)
library(GenomicRanges)

#' Sorts input frame in ascending order according to start
#' removes duplicates and values in which the p value is > pval  
#' In df$motif, substitutes all :', '-' and '_' chars with 'X'
#' 
#' @param df Input data frame in lasagne output format
#' @param pval maximum pval to keep
#' @return df 
processSingle  <- function(df, pval=0.01){
    df=subset(df, p_value < pval) 
    # order by start
    df = df[!duplicated(df[c("chr","start","motif")]), ]
    df=df[with(df, order(chr, start)), ]
    df$motif = gsub(":", "X", df$motif)
    df$motif = gsub("-", "X", df$motif)
    df$motif = gsub("_", "X", df$motif)
    return(df)
}

stem <- function(raw) {
    raw$motif = gsub("p53", "p53_", raw$motif)
    raw$motif = gsub("NKX6X1", "NKX", raw$motif)
    raw$motif = gsub("NKX3A", "NKX", raw$motif)
    raw$motif = gsub("[[:digit:]]+$", "", raw$motif)
    return(raw)
}

randomize <- function(df){
    
    newdf = data.frame(chr=df$chr, start=df$start, stop=df$stop,
    strand=df$strand, score=df$score, 
        p_value=df$p_value, motif=sample(df$motif, replace = FALSE))
    
    return(newdf)
}

#' Adds width and bin columns to input data frame
#' 
#' @param df Input data frame in lasagne format
#' @param Window window size
#' @param step step size 
#' @return result
slidingWindow <- function( df, window, step){
    seqbegin = df$start[1]
    seqend = df$stop[nrow(df)]
    total = seqend - seqbegin
     if (window >= total ){
        result=df
        result$width = total
        result$bin=sprintf('[%s,%s]', seqbegin, seqend )
        return(result) 
    }
    spots <- seq(from = 0, to = (total - window ), by = step)
    spots <- spots + seqbegin 
    spots_end <- spots + window - 1
    spots_end[length(spots)]=seqend
    bin=sprintf('[%s,%s]', spots,spots_end)
   
    
    df_gr=with( df,  GRanges( rep('a',nrow(df)),IRanges(start, stop),
                              str=strand, score=score,p=p_value, motif=motif, contig=chr ))
    intervals_gr <- GRanges(rep('a',length(spots)),
                            IRanges(start=spots, end=spots_end), names=bin )
    ol=as.data.frame(findOverlaps( df_gr, intervals_gr))
    result = cbind( df[ol$queryHits,],
                    as.data.frame(intervals_gr)[ol$subjectHits,][,c('width','names')])
    names(result)= c("chr", "start", "stop", "strand", "score", "p_value", "motif","width", "bin")
    return(result) 
    }



cutseq <- function(df, binsize=900) {
    seqbegin = df$start[1]
    seqend = df$start[nrow(df)]
    length = seqend - seqbegin
    b = ceiling(length/binsize)
    if (b==1) { df$bin = sprintf('[%s-%s]', seqbegin, seqend) } 
    else{ df$bin =cut(df$start, breaks = b )}# , labels = sprintf("%d", 1:b)) }
    return(df)
}


cossine <- function( m1, m2 ) {
    rxy <- mapply('*', m2, m1[names(m2)])
    dpxy = apply(rxy, 1, sum)
    x = as.matrix(m1)
    dpxx = x %*% t(x)
    y = as.matrix(m2) 
    dpyy = diag( y %*% t(y) )
    denom = sqrt( rep(dpxx, length(dpyy)) * dpyy )
    return( dpxy/denom )

}

weight <- function( df, weightsdf ){
    common = intersect(names(df), names(weightsdf))
    df= df[,common]
    df = mapply('*', df, msig[names(df)])
    return(df)
}


weight_by_bg <- function( df, bg_weights ){
    # downweights common TFBS
    many = bg_weights[ bg_weights < 1 ]
    others = df[, ! names(df) %in% names(many)]
    others = replace(others, ! others==0, 1)[1,]
    common = intersect(names(df), names(many))
    all = c( many[common], as.vector(others) )
    all = unlist(all, use.names=T)
    df = mapply('*', df, all[names(df)])
    return(df)
}


scalar1 <- function(x) {x / sqrt(sum(x^2))}

bg <- function( f, weights=msig, m1 , window=window, step=step, seedval = seedval) {
    
    if ( !missing(seedval)){
        set.seed(42)
    }

    d2_orig=read.delim(f,header=F,sep='')
    colnames(d2_orig)=c("chr","start","stop","strand","score","p_value","motif")
    x=processSingle(d2_orig)
    x=randomize(x)
    #x=cutseq(x, binsize=size)
    x = slidingWindow( x, window, step)
    df2 <- aggregate(motif ~ bin, data = x, paste, collapse = " ")

    corp2 <- Corpus(VectorSource(df2$motif))
    dtm2 <- DocumentTermMatrix(corp2, control=list(tolower=FALSE))
    inspect(dtm2)
    m2 = as.matrix(dtm2) 
    # L2 norm
    m2 = t(apply( m2 , 1, function(x) scalar1(x) ))
    m2 = as.data.frame(m2)

    # weight by msig
    m2<-weight(m2, msig)
    m2=as.data.frame(m2)

    # only use common names
    common = intersect(names(m1), names(m2))
    m1= m1[,common]
    m2= m2[,common]

    # cossine similarity comparing x to all blocks in y
    # m1 is a named vector

    rxy <- mapply('*', m2, m1[names(m2)])
    dpxy = apply(rxy, 1, sum)
    x = as.matrix(m1)
    dpxx = x %*% t(x)
    y = as.matrix(m2) 
    dpyy = diag( y %*% t(y) )
    denom = sqrt( rep(dpxx, length(dpyy)) * dpyy )
    cs = dpxy/denom
    return(cs)

}

assign_p <- function(res, null) {
    pvals=c()
    minpn = nrow(subset(res, score> null$score[1]))
    minp=1/nrow(res)
    pvals = c( pvals, rep(minp, minpn))
    for (i in 1:nrow(res)){
        ind = which(abs(null$score-res$score[i])==min(abs(null$score-res$score[i])))
        if (ind > 1) {
            pvals = c( pvals, ind/nrow(res))
        }
    }
    remain = nrow(res) - length(pvals) 
    if (remain>0){ pvals = c( pvals, rep(1, remain)) }
    return(pvals)
}