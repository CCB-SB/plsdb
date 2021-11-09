#!/usr/bin/Rscript

####################
# LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))

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
    testit::assert('ACC_NUCCORE' %in% colnames(df))
    if( any( !grepl('^.*\\.[0-9]+$', df$ACC_NUCCORE) ) ){
        testit::assert('ACC_FASTA' %in% colnames(df))
        df$ACC_NUCCORE <- df$ACC_FASTA
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
# VARS
# TODO

####################
# Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))
# args <- list(
#     'new'='2018_12_03.tsv',
#     'old'='2018_09_14.tsv',
#     'new_nonfiltered'='../reports/2018_12_03__full2.txt',
#     'otab'='test.tsv',
#     'olog'='test.log'
# )

####################
## DATA
# Read in
newTab <- read.csv(file=args$new, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))
oldTab <- read.csv(file=args$old, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))
newAll <- read.csv(file=args$new_nonfiltered, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))

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

# Make sure ACC_NUCCORE contains accession w/ version number
# and remove ACC_FASTA column (which was accession + version)
newTab <- round_coords(rm_accfasta_col(replace_acc_col(newTab)))
oldTab <- round_coords(rm_accfasta_col(replace_acc_col(oldTab)))
#newAll <- round_coords(rm_accfasta_col(replace_acc_col(newAll)))

# Set accessions as row names
rownames(newTab) <- newTab$ACC_NUCCORE
rownames(oldTab) <- oldTab$ACC_NUCCORE
rownames(newAll) <- newAll$ACC_NUCCORE

## COMPARE
# Removed records
rc_rm <- data.frame(
    ACC_NUCCORE=setdiff(oldTab$ACC_NUCCORE, newTab$ACC_NUCCORE),
    Reason=NA,
    stringsAsFactors=FALSE, check.names=FALSE
)
rownames(rc_rm) <- rc_rm$ACC_NUCCORE
for(rc_acc in rownames(rc_rm)){
    is_identical <- !is.na(newTab$Identical) & grepl(sub('\\.', '\\\\.', rc_acc), newTab$Identical)
    is_oldversion <- !is.na(newTab$OldVersion) & grepl(sub('\\.', '\\\\.', rc_acc), newTab$OldVersion)
    if( any( is_identical ) ){
        rc_rm[rc_acc, 'Reason'] <- sprintf(
            'identical to %s',
            paste(newTab$ACC_NUCCORE[is_identical], collapse=';')
        )
    } else if( any( is_oldversion ) ){
        rc_rm[rc_acc, 'Reason'] <- sprintf(
            'old version of %s',
            paste(newTab$ACC_NUCCORE[is_oldversion], collapse=';')
        )
    } else if(rc_acc %in% newAll$ACC_NUCCORE){
        rc_rm[rc_acc, 'Reason'] <- 'removed during filtering'
    } else {
        rc_rm[rc_acc, 'Reason'] <- 'not retrieved'
    }
}

# Added records
rc_add <- setdiff(newTab$ACC_NUCCORE, oldTab$ACC_NUCCORE)

# Changed/updated records
rc_mod <- data.frame(
    ACC_NUCCORE=intersect(newTab$ACC_NUCCORE, oldTab$ACC_NUCCORE),
    Changes="",
    stringsAsFactors=FALSE, check.names=FALSE
)
rownames(rc_mod) <- rc_mod$ACC_NUCCORE
cols <- intersect(colnames(newTab), colnames(oldTab))
cols <- setdiff(cols, c('D1', 'D2'))
for(rc_acc in rownames(rc_mod)){
    new_v <- as.vector(newTab[rc_acc, cols])
    old_v <- as.vector(oldTab[rc_acc, cols])
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
    ACC_NUCCORE=c(rownames(rc_rm), rc_add, rownames(rc_mod)),
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
