library(dplyr)
#library(rgnparser)
library(taxize)
library(tidytext)
library(stringr)
#Naive improvement script to the table

options(error=traceback)
args = commandArgs(trailingOnly=TRUE)

df<-read.csv(args[1],sep="\t",header=T)

df<-read.csv(args[1],sep="\t",header=F,quote = "")
colnames(df)<-df[1,]
df<-df[-1,]
#Clean up lack of information
missing_information = c("-","na", "n/a", "n.a.", "missing", "none",  "not applicable", "not available", "not collected", "not determined", "not recorded",  "unavailable", "unknown", "unspecified","unidentified", "null", "no", "undocumented", "not defined", "not reported", "not applied")
df[apply(df,2,function(x){ sapply(x,function(y){ tolower(y) %in% missing_information })})]<-""
df$IsolationSource_BIOSAMPLE_processed<-tolower(df$IsolationSource_BIOSAMPLE)
to_correct<-tolower(df$Host_BIOSAMPLE)

#Read in already assigned Input -> Taxa mapping
existing_map<-read.csv(args[2],sep="\t",header=F)
colnames(existing_map)<-c("Input","Assigned")

update_mapping<-function(input_vec=existing_map){
  result<-input_vec$Assigned
  names(result)<-input_vec$Input
  return(result)
}

existing_mapping=update_mapping()




#Use the mapping
apply_existing_map<- function(vector,existing_map=existing_mapping){
  return(sapply(vector,function(x) ifelse(x%in%names(existing_map),existing_map[x],"FALSE")))
}

#Best try at aoutomatically assinging a taxa to a given input
fix_host<-function(host,ask=F,ask_ever=F,bold=T){
      host<-trimws(host)
      if (host %in% names(existing_mapping))
        return(existing_mapping[host])
      result<-get_uid(host,ask=ask)
      if(!is.na(result[1]))
        return(taxize::ncbi_get_taxon_summary(result[1])$name)
      if (bold){
        result<-bold_search(
          gsub(host,pattern = "'",replacement = "")%>% gsub(pattern = '"',replacement = "") 
        )
        if(dim(result)[1]==1&dim(result)[2]>1)
          return(result$taxon)
      }
      if (grepl(host,pattern = 'sp.',fixed=T)){
        return(fix_host(gsub(host,pattern = "sp.",replacement = ""),ask,ask_ever))
      }
      if (grepl(host,pattern = '(',fixed=T)){
        return(fix_host(head(strsplit(host,split = "(",fixed=1)[[1]],n=1),ask,ask_ever))
      }
      if ((!ask)&ask_ever)
        return(fix_host(host,ask = T,ask_ever=T))
      return ("FALSE")
}

separating_characters=paste(strsplit("\"/\\'()[]{};:,.#+~-",split = '')[[1]],collapse = "|\\")

#TODO
get_stems<-function(x){
  trimws(x)
}

#Remove stop words from vector
remove_stopword<-function(x){
  if (x %in% get_stopwords()$word)
    return("")
  return(x)
}

## WITHIN COLUMN CORRECTION

assign_procedure<-function(input_vec){
  #Were they assigned previously?
  to_assign<-unique(input_vec)
  to_assign<-to_assign[apply_existing_map(to_assign)==FALSE]
  #Can we assign them instantly?
  assigned<-sapply(unique(to_assign),fix_host)
  new_assigned=assigned[assigned!="FALSE"]
  to_assign<-assigned[assigned=="FALSE"]%>%names()
  #Can we assign substrings?
  #generate all substrings test and if result is unique: take it
  
  to_assign<-sapply(to_assign,function(x)strsplit(x,split = separating_characters,fixed=F))
  assigned<-lapply(to_assign,FUN=function(x){ 
      result<-sapply(x,FUN= function(y){fix_host(y)})
      result<-result[result!=FALSE]%>%unique()
      if (length(result)==1)
        return (result)
      return ("FALSE")
  })%>%unlist()
  new_assigned=c(new_assigned,assigned[assigned!="FALSE"])
  #TODO check if names is set correct here
  new_assigned<-cbind.data.frame(names(new_assigned),new_assigned)
  if (ncol(new_assigned)==2){
	        colnames(new_assigned)<-c("Input","Assigned")
  }else if (ncol(new_assigned)==0|ncol(new_assigned)==1){
	        new_assigned<-data.frame(Input=character(),
					                Assigned=character())
  }else{  stop("Incorrect dimenesionality in assign procedure")}
  return(new_assigned)
}
print("Assign Biosample")
Host_BIOSAMPLE_assignments<-assign_procedure(to_correct)
existing_map=rbind.data.frame(existing_map,Host_BIOSAMPLE_assignments) %>%unique()
existing_mapping=update_mapping(existing_map)
## CROSS COLUMN CORRECTION
#If no host is specified, try to use the isolation source for host inference
print("Assign IsolationSource")
Iso_source_assignments<-filter(df,Host_BIOSAMPLE=="",df$IsolationSource_BIOSAMPLE_processed!="")$IsolationSource_BIOSAMPLE_processed %>%assign_procedure()
existing_map=rbind.data.frame(existing_map,Iso_source_assignments)%>%unique()
existing_mapping=update_mapping(existing_map)

Host_BIOSAMPLE_corrected<-apply_existing_map(to_correct,existing_mapping)
Host_BIOSAMPLE_inferred<-apply_existing_map(df$IsolationSource_BIOSAMPLE_processed,existing_mapping)
Host_BIOSAMPLE_combined=ifelse((Host_BIOSAMPLE_corrected!="FALSE" & !is.na(Host_BIOSAMPLE_corrected)),Host_BIOSAMPLE_corrected,Host_BIOSAMPLE_inferred)
Host_BIOSAMPLE_final=ifelse(Host_BIOSAMPLE_combined!="FALSE",Host_BIOSAMPLE_combined,"")
df$Host_BIOSAMPLE_processed<-Host_BIOSAMPLE_final
#unchanged<-df$Host_BIOSAMPLE[Host_BIOSAMPLE_final!=df$Host_BIOSAMPLE]
print("Write hosts")
write.table(existing_map,file=args[2],sep="\t",quote=F,row.names=F,col.names=F)
print("Write output")
df%>%select(!IsolationSource_BIOSAMPLE_processed)%>%write.table(file=args[3],sep="\t",quote=F,row.names=F,col.names=T,na = "",qmethod = "escape")

