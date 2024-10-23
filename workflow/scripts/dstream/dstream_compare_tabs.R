#!/usr/bin/Rscript

####################
# LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(data.table))

####################
# ARGS
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--new',  '-n', help='New table with plasmid meta data', required=TRUE)
    parser$add_argument('--old',  '-o', help='Old table with plasmid meta data', required=TRUE)
    parser$add_argument('--new_nonfiltered',  '-f', help='New table before filtering', required=TRUE)
    parser$add_argument('--otab', '-t', help='Output table file')
    parser$add_argument('--olog', '-l', help='Output log file')
    return(parser)
}

####################
# FUNCTIONS
replace_acc_col <- function(df){
    testit::assert('NUCCORE_ACC' %in% colnames(df))
    if( any( !grepl('^.*\\.[0-9]+$', df$NUCCORE_ACC) ) ){
        testit::assert('ACC_FASTA' %in% colnames(df))
        df$NUCCORE_ACC <- df$ACC_FASTA
    }
    return(df)
}

rm_accfasta_col <- function(df){
    if('ACC_FASTA' %in% colnames(df)){
        df <- df[, setdiff(colnames(df), 'ACC_FASTA')]
    }
    return(df)
}

round_coords <- function(df){
    testit::assert('loc_lat' %in% colnames(df))
    testit::assert('loc_lng' %in% colnames(df))
    df$loc_lat <- round(df$loc_lat, 2)
    df$loc_lng <- round(df$loc_lng, 2)
    return(df)
}


####################
# Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

####################
## NUCCORE DATA

newTab <- read.csv(file=args$new, header=TRUE, check.names=FALSE,
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1'))
oldTab <- read.csv(file=args$old, header=TRUE, check.names=FALSE,
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1'))
print(colnames(oldTab))
print(colnames(newTab))

# Unify new_names vs old_names
old_colnames = c("NUCCORE_TaxonID","D1", "D2")
new_colnames = c("TAXONOMY_UID", "UMAP_D1", "UMAP_D2")

# Supressed columns
unique_cols_old_tab <- c(
    "NUCCORE_Genome",
    ## BIOSAMPLE
    "BIOSAMPLE_UID", "BIOSAMPLE_ACC", "BIOSAMPLE_Location", 
    "BIOSAMPLE_Coordinates", "BIOSAMPLE_IsolationSource", "BIOSAMPLE_Host", 
    "BIOSAMPLE_Host_processed", "BIOSAMPLE_Host_processed_source", "BIOSAMPLE_Host_label",
    "BIOSAMPLE_CollectionDate", "BIOSAMPLE_HostDisease", "BIOSAMPLE_SampleType",
    "BIOSAMPLE_HostDisease_processed",
    "loc_lat", "loc_lng", "loc_parsed",
    ## ASSEMBLY
    "ASSEMBLY_UID", "ASSEMBLY_Status", "ASSEMBLY_SeqReleaseDate", 
    "ASSEMBLY_SubmissionDate", "ASSEMBLY_Lastest",
    ## TAXONOMY
    "TAXONOMY_taxon_name", "TAXONOMY_taxon_rank", "TAXONOMY_taxon_lineage", 
    "TAXONOMY_superkingdom", "TAXONOMY_phylum", "TAXONOMY_class",
    "TAXONOMY_order", "TAXONOMY_family", "TAXONOMY_genus", 
    "TAXONOMY_species","TAXONOMY_superkingdom_id", "TAXONOMY_phylum_id", 
    "TAXONOMY_class_id", "TAXONOMY_order_id", "TAXONOMY_family_id",
    "TAXONOMY_genus_id", "TAXONOMY_species_id", 
    ## rMLST
    "rMLST_hits", "rMLST_hitscount",
    ## Others
    "inclusions",
    "has_biosample", "has_assembly", "has_location",
    "Length", "plasmidfinder", "pmlst"
    )

# New columns
unique_cols_new_tab <- c("NUCCORE_has_identical", "VIRAL_category", "VIRAL_putative", "STATUS", "STATUS_note")

OldTab <- setnames(oldTab, 
    old = old_colnames,
    new = new_colnames
)

newAll <- read.csv(file=args$new_nonfiltered, header=TRUE, check.names=FALSE,
                 stringsAsFactors=FALSE, na.strings=c('', 'NA', -1, '-1'))

# Log basic stuff
write(
    x=sprintf(
        'New file: %s (%d x %d)\nOld file: %s (%d x %d)\n\n',
        args$new, nrow(newTab), ncol(newTab),
        args$old, nrow(oldTab), ncol(oldTab)
    ),
    file=args$olog,
    append=FALSE
)
write(
    x=sprintf(
        'Columns:\nOnly in new file: %s\nOnly in old file: %s\n',
        paste(setdiff(colnames(newTab), colnames(oldTab)), collapse=', '),
        paste(setdiff(colnames(oldTab), colnames(newTab)), collapse=', ')
    ),
    file=args$olog,
    append=TRUE
)

write(
    x=sprintf(
        'Old column names: %s \nNew column names: %s \n\n',
        paste(old_colnames, collapse=", "),
        paste(new_colnames, collapse=", ")
    ),
    file=args$olog,
    append=TRUE
)

# Make sure NUCCORE_ACC contains accession w/ version number
# and remove ACC_FASTA column (which was accession + version)
# newTab <- round_coords(rm_accfasta_col(replace_acc_col(newTab)))
oldTab <- round_coords(rm_accfasta_col(replace_acc_col(oldTab)))
#newAll <- round_coords(rm_accfasta_col(replace_acc_col(newAll)))

# Set accessions as row names
rownames(newTab) <- newTab$NUCCORE_ACC
rownames(oldTab) <- oldTab$NUCCORE_ACC
rownames(newAll) <- newAll$NUCCORE_ACC

## COMPARE
# Removed records
rc_rm <- data.frame(
    NUCCORE_ACC=setdiff(oldTab$NUCCORE_ACC, newTab$NUCCORE_ACC),
    Reason=NA,
    stringsAsFactors=FALSE, check.names=FALSE
)
rownames(rc_rm) <- rc_rm$NUCCORE_ACC
for(rc_acc in rownames(rc_rm)){
    is_identical <- !is.na(newTab$Identical) & grepl(sub('\\.', '\\\\.', rc_acc), newTab$Identical)
    is_oldversion <- !is.na(newTab$OldVersion) & grepl(sub('\\.', '\\\\.', rc_acc), newTab$OldVersion)
    if( any( is_identical ) ){
        rc_rm[rc_acc, 'Reason'] <- sprintf(
            'identical to %s',
            paste(newTab$NUCCORE_ACC[is_identical], collapse=';')
        )
    } else if( any( is_oldversion ) ){
        rc_rm[rc_acc, 'Reason'] <- sprintf(
            'old version of %s',
            paste(newTab$NUCCORE_ACC[is_oldversion], collapse=';')
        )
    } else if(rc_acc %in% newAll$NUCCORE_ACC){
        rc_rm[rc_acc, 'Reason'] <- 'removed during filtering'
    } else {
        rc_rm[rc_acc, 'Reason'] <- 'not retrieved'
    }
}

# Added records
rc_add <- setdiff(newTab$NUCCORE_ACC, oldTab$NUCCORE_ACC)

# Changed/updated records
rc_mod <- data.frame(
    NUCCORE_ACC=intersect(newTab$NUCCORE_ACC, oldTab$NUCCORE_ACC),
    Changes="",
    stringsAsFactors=FALSE, check.names=FALSE
)
rownames(rc_mod) <- rc_mod$NUCCORE_ACC
cols <- intersect(colnames(newTab), colnames(oldTab))
cols <- setdiff(cols, c('UMAP_D1', 'UMAP_D2'))

for(rc_acc in rownames(rc_mod)){
    new_v <- as.vector(unlist(newTab[rc_acc, cols]))
    old_v <- as.vector(unlist(oldTab[rc_acc, cols]))
    rc_changes = cols[
        !(
            (is.na(new_v) & is.na(old_v)) |
            (!is.na(new_v) & !is.na(old_v) & new_v == old_v)
        )
    ]
    rc_mod[rc_acc, 'Changes'] <- paste(rc_changes, collapse=';')
}
rc_mod <- rc_mod[rc_mod$Changes != "",]

# Removed/added/changed
all_changes <- data.frame(
    NUCCORE_ACC=c(rownames(rc_rm), rc_add, rownames(rc_mod)),
    Flag=c(
        rep('Removed', nrow(rc_rm)),
        rep('Added', length(rc_add)),
        rep('Changed', nrow(rc_mod))
    ),
    Comment=c(
        rc_rm$Reason,
        rep('', length(rc_add)),
        rc_mod$Changes
    )
)
write.table(x=all_changes, file=args$otab, sep='\t', row.names=FALSE, quote=FALSE)
