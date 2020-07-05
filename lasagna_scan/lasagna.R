# This is code to use the Lasagna method to scan for motif matches to the Transfac database given an input FASTA sequence
# Argument file name of sequence .fa file to be scanned
# Note: input sequence must have the .fa suffix and header line must be the genomic position in the format '>sequenceID:start-end'
# Output is a text file containing coordinates, motif names and p-values.
# Requires: python2
# 5th July 2020

filename <- commandArgs(trailingOnly = TRUE)[1]

annot = read.delim('../transfacID2geneName',header=F)
annot$V1 = gsub(".json", "", annot$V1)

modelnames <- Sys.glob("../transfacmatrix/*")


id = strsplit(basename(filename), '.fa')[1]
# lasagna requires the header to be a certain format
file = readLines(filename, -1)
if (grepl( '_', file[1])){
  file[1] = paste('>',strsplit(as.character(id) , split='_')[[1]][2],sep='')
  writeLines(file,filename)
}

dir= dirname(filename)
for (modelname in modelnames) {
  model = strsplit(basename(modelname), '.json')[1]
  motif = subset(annot, V1==model)$V2
  outname = sprintf('%s%s_%s_%s', working, id, model, motif)
  cmd = sprintf('python offline_scanning/scan_UCSC_promoters.py --pvalue 0.01 --compute-pvalue %s %s %s' ,
                modelname, filename, outname )
  system( cmd )
  motifout = read.delim(outname, sep='', header=T)
  if( !nrow(motifout) == 0){
    motifout$motif = motif
    write.table( motifout, sprintf('%s/%s_out', dir,id), quote=F, row.name=F, col.name=F, append=T)
  }
  system( sprintf('rm %s', outname) )
  
}
