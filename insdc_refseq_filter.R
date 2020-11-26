# The remaining collisions at this point all bear the idenddtical suffix with exception of the numeric value behind the last "." character. We keep the one from INSDC since it seems to be more up to date


args <- commandArgs(trailingOnly = TRUE)
input_filename<-args[1]
output_filename<-args[2]

df<-read.csv(input_filename,sep="\t",header=T)
accs<-df$ACC_NUCCORE

prefices<-sapply(accs,function(x) strsplit(as.character(x),'.',fixed=T)[[1]][1])
potential_collisions<-sapply(prefices,function(x) tail(strsplit(as.character(x),'_',fixed=T)[[1]],n=1))

#1. safety net
checker<-table(potential_collisions)
checker<-names(checker[checker>1])
if (!all(grepl("CP",checker,fixed=T)))
	stop("Remaining collusions of the ACC_NUCCORE column are not trivial.")

#2. safety net
checker<-potential_collisions[duplicated(potential_collisions)]
checker<-names(checker)
if (!all(grepl("NZ_CP",checker,fixed=T)))
	stop("Remaining collusions of the ACC_NUCCORE are not sorted as expected. revisit. insdc_refseq_filter.R")

rowmask<-!duplicated(potential_collisions)
output_df<-df[rowmask,]
write.table(output_df,file=output_filename,quote=F,sep='\t',col.names=T,row.names=F)
