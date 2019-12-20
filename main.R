#!/usr/bin/Rscript        

# R version 3.6.2
# Required packages:
# - Bioconducter 3.10
# - tm

# Changelog ~ 
# 17/12/19
# Changed syntax to R 3.6.2 so it talks to current bioconductor lib
# Import GenomocRanges lib for GRanges function
# Changed query and subject paths to absolute to increase portability
# small bug cleanups  
# 18/12/29    
# added optional seed param to bg function 
# Changed script output to the result array 

# TODO 
# take parameters as command line args 
# verify output   

# input params
# bin size
# seed
# p val 
# window size 


source("functions.R")

# get command line args
argv = commandArgs(trailingOnly = TRUE)
if (length(argv) < 2){
  stop("Usage: Main.r <query> <target>")
}

# set query and target
lasagna_out.query = argv[1] # Lasagna motif matching output for query sequence 
lasagna_out.subject = argv[2] # Lasagne motif matching output for target sequence

# look through any remaining command line args for flags
arg_pval = NA
arg_window = NA
arg_seed = NA

i = 3 # first optional arg
while(i <= length(argv)){
  
  # pval, double between 0 and 1
  if( argv[i] == "-pval" || argv[i] == "-p" ){
    #suppressWarnings( (arg_pval = as.double(argv[i+1])) )
    arg_pval = as.double(argv[i+1])
    print(c("argv: ", arg_pval))
    if (is.na(arg_pval) || !is.numeric(arg_pval) || arg_pval < 0 || arg_pval > 1){
      stop("Please check pval, it should be numeric, and 0 <= pval <= 1")
    }
  }
  
  # window, positive integer
  else if( argv[i] == "-window" || argv[i] == "-w"){
    suppressWarnings( (arg_window = as.integer(argv[i+1])) )
    if (is.na(arg_window) || arg_window < 1){
      stop("Please check window, it should be a positive integer")
    }
  }
    
  # seed, any integer 
  else if (argv[i] == "-seed" || argv[i] == "-s"){
    suppressWarnings( (arg_seed = as.integer(argv[i+1])) )
    if (is.na(arg_seed)){
      stop("Please check seed flag, it should be an integer")
    }
  }
  i=i+1 # step  
}

# Absolute path to files
workingdir = getwd()
workingdir = paste(workingdir, "/", sep='')
lasagna_out.query = paste(workingdir, lasagna_out.query, sep='')
lasagna_out.subject = paste(workingdir, lasagna_out.subject, sep = '')

# set parameters

# window
window=900  # set window size larger than length of query
step=900
if (! is.na(arg_window)){ window = arg_window }

# pval
pval_set = 0.001 # default
if (! is.na(arg_pval)){pval_set = arg_pval}

# seed
if (is.na(arg_seed)){seed = NULL}else{seed = arg_seed}

# debug input params
#print(c("Window: ", window))
#print(c("Pval: ", pval_set))
#print(c("Seed", seed))
#stop("Bye")

# creates data frame containing lasagne output of query 
invisible((d1=read.delim(lasagna_out.query,header=F,sep='')))
colnames(d1)=c("chr","start","stop","strand","score","p_value","motif")
x=processSingle(d1, pval=pval_set)
x = slidingWindow( x, window, step)

# creates data frame from query data
# df1[1] is bin 
# df1[2] is motif
df1 <- aggregate(motif ~ bin, data = x, paste, collapse = " ")

# if >1 unique motifs in query take the highest score 
# msig ~ array of scores(double) named with motifs
msig <- aggregate(score ~ motif, data = x, max)
msig <- setNames( msig$score, as.character(msig$motif))
    
# data frame for subject file (lasagne format)
d2=read.delim(lasagna_out.subject,header=F,sep='')
colnames(d2)=c("chr","start","stop","strand","score","p_value","motif")
x=processSingle(d2, pval=0.01) #?
x = slidingWindow( x, window, step)
df2 <- aggregate(motif ~ bin, data = x, paste, collapse = " ")

corp1 <- Corpus(VectorSource(df1$motif))
dtm1 <- DocumentTermMatrix(corp1, control=list(tolower=FALSE))
#inspect(dtm1)
m1 = as.matrix(dtm1) 

# L2 norm
m1 = scalar1(m1)
m1 = as.data.frame(m1)  

# weight by query motif frequencies 
m1<-weight(m1, msig) #names still fine here
m1=as.data.frame(t(m1))

corp2 <- Corpus(VectorSource(df2$motif))
dtm2 <- DocumentTermMatrix(corp2, control=list(tolower=FALSE))
#inspect(dtm2)
m2 = as.matrix(dtm2) 
# L2 norm
m2 = t(apply( m2 , 1, function(x) scalar1(x) ))
m2 = as.data.frame(m2)

# weight by query motif frequencies
m2<-weight(m2, msig)
m2=as.data.frame(m2)

# compare based on common motifs 
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

# bin
res=data.frame( bin=df2$bin, score=cs)
res$chr = d2[1,1]
res=res[order(-res$score),]
res=res[,c(3,1,2)]

# examine top 15
#cs[order(-cs)[1:15]] # prints to stdout

names(cs) = seq(1:length(cs))
#names(cs[order(-cs)[1:15]]) # prints to stdout

  #---------------------# 
# p-val calculations
#---------------------# 

f=lasagna_out.subject
cs_null= bg(f, weights=msig, m1, window=window, step=step, seedval=seed)
  
null=data.frame( bin=df2$bin, score=cs_null)
null$chr = d2[1,1]  
null=null[order(-null$score),] 
null=null[,c(3,1,2)]

    
res$p=assign_p(res,null) 
res$padj = p.adjust(res$p, method='fdr')
      
print(res)

