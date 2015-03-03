options(stringsAsFactors=FALSE)
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("Biostrings")

# see also: python version of biomart https://pypi.python.org/pypi/biomart/0.4.0

# define biomart to use
newest_mart = biomaRt::useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
# getSequence() is no longer supported for grch37, so not using it
# grch37 = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

buffer = 125 # extra bp on either side of desired position in PCR product


sgrna = read.table('~/d/j/cureffilab/media/2015/03/sgrna_designs.txt',sep='\t',header=TRUE)

# examples for testing
# enst = 'ENST00000379440' # PRNP, + strand
# enst = 'ENST00000389817' # ABCC8, - strand
# sgrna = "TGACCCCAGCCACCACCATG" # against PRNP

tx_genomic_seq = biomaRt::getSequence(id=unique(sgrna$transcript),
                     type='ensembl_transcript_id',
                     seqType='gene_exon_intron',
                     mart=newest_mart)
for (i in 1:dim(sgrna)[1]) {
  spacer = as.character(sgrna$spacer[i])
  genomic_seq = as.character(subset(tx_genomic_seq, ensembl_transcript_id==sgrna$transcript[i], select=gene_exon_intron))
  rel_coordinate = gregexpr(spacer,genomic_seq,fixed=TRUE)[[1]][1]
  if (rel_coordinate == -1) {
    rcomp = as.character(Biostrings::reverseComplement(Biostrings::DNAString(spacer)))
    rel_coordinate = gregexpr(rcomp,genomic_seq,fixed=TRUE)[[1]][1]
  }
  seq_context = substr(genomic_seq,rel_coordinate-buffer,rel_coordinate+buffer)
  sgrna$seq_context[i] = seq_context
}

write.table(sgrna,'~/d/j/cureffilab/media/2015/03/sgrna_designs_with_context.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

# sequence_target = paste(as.character(buffer-nchar(spacer)),",",as.character(2*nchar(spacer)),sep='')
# primer3_input = paste(sep='',
#   "SEQUENCE_ID=","CRISPR_",enst,"_",spacer,"\n",
#   "SEQUENCE_TEMPLATE=",seq_context,"\n",
#   "SEQUENCE_TARGET=",sequence_target,"\n",
#   "PRIMER_TASK=pick_detection_primers","\n",
#   "PRIMER_PICK_LEFT_PRIMER=1","\n",
#   "PRIMER_PICK_INTERNAL_OLIGO=1","\n",
#   "PRIMER_PICK_RIGHT_PRIMER=1","\n",
#   "PRIMER_OPT_SIZE=18","\n",
#   "PRIMER_MIN_SIZE=15","\n",
#   "PRIMER_MAX_SIZE=36","\n",
#   "PRIMER_MAX_NS_ACCEPTED=1","\n",
#   "PRIMER_PRODUCT_SIZE_RANGE=150-200","\n",
#   "P3_FILE_FLAG=1","\n",
#   "SEQUENCE_INTERNAL_EXCLUDED_REGION=",sequence_target,"\n",
#   "PRIMER_EXPLAIN_FLAG=0","\n")
# 
# cat(primer3_input)
