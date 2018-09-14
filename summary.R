#!/usr/bin/Rscript

####################
# LIBS
suppressMessages(library(testit))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(maps))

####################
# ARGS
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--tab', help='Table with plasmid meta data', required=TRUE)
    parser$add_argument('--pdf', help='Output PDF', required=TRUE)
    parser$add_argument('--width', help='Output PDF width', default=6, type='integer')
    parser$add_argument('--height', help='Output PDF height', default=6, type='integer')
    return(parser)
}

####################
# VARS
ranks <- c('species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom')

theme_font <- theme(
    axis.text   =element_text(size=12),
    axis.title  =element_text(size=14, face="bold"),
    legend.text =element_text(size=12),
    legend.title=element_text(size=14, face="bold")
)

####################
# Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

####################
# Data
info <- read.csv(file=args$tab, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))

# proc. data -> year
info$year <- as.numeric(sapply(info$CreateDate_NUCCORE, function(x){ unlist(strsplit(x, '/'))[1] }))

####################
# Data for stat.s and plots

# num. unique taxa per ranks
num_taxa <- data.frame()
for(r in ranks){
    num_taxa <- rbind(num_taxa, data.frame(
        rank=r,
        count=length(unique(na.omit(info[,sprintf('taxon_%s_id', r)])))
    ))
}
num_taxa$rank <- factor(num_taxa$rank, levels=ranks, ordered=TRUE) # for plots

# STs from pMLST
STs <- sapply(info$pmlst, function(x){
    if(is.na(x)){
        return(NA)
    } else {
        return(unlist(strsplit(x, ':'))[1])
    }
})

# Locations: count number of points per locus (same coordinates)
locs_aggr <- aggregate(x=1:nrow(info), by=list(loc_lat=info$loc_lat, loc_lng=info$loc_lng), FUN=function(x){ return(length(x)) })
colnames(locs_aggr)[colnames(locs_aggr)=='x'] <- 'Count'

####################
# Stats

# Number of taxa
print('Number of unique taxa per rank')
print(num_taxa)

# Length and GC summary
print('Length summary')
print(summary(info$Length_NUCCORE))
print('GC summary')
print(summary(info$GC_NUCCORE))

# How many have (mapped) location information
print(sprintf('With location information: %d', sum(!is.na(info$loc_lat))))

# How many have any pMLST hits, how many have an ST
print(sprintf('With an pMLST hit: %d', sum(!is.na(STs))))
print(sprintf('With an pMLST ST: %d',  sum( sapply(STs, function(x){
    if(is.na(x)){
        return(FALSE)
    } else {
        return(!grepl('\\(-\\)$', x))
    }
}) )))

# How many have PlasmidFinder hits
print(sprintf('With an PlasmidFinder hit: %d', sum(!is.na(info$plasmidfinder))))

####################
# Plots

# plot: num. uniq. taxa
plot_utaxa <-   ggplot(data=num_taxa, aes(x=rank, y=count)) +
                geom_col(fill='white', colour='#333333', position="dodge") +
                geom_text(aes(label=count, y=count + 40), position=position_dodge(0.9), angle=0) +
                xlab('') + ylab('Number of unique taxa / rank') +
                theme_bw() +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
                ) +
                theme_font

# plot: dates
plot_year <-    ggplot(data=info, aes(x=year)) +
                geom_histogram(binwidth=1, fill='white', colour='#333333', position="dodge") +
                scale_y_sqrt(breaks=c(100, seq(1000, 10000, 1000))) +
                xlab('Year') + ylab('Number of nuccore records') +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                theme_font

# plot: locations
plot_locs <-    ggplot(data=locs_aggr, aes(x=loc_lng, y=loc_lat, size=Count)) +
                borders("world", colour="#cccccc", fill="#ffffff") +
                geom_point(colour='#333333', shape=4) +
                theme_bw() +
                theme(
                    axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    panel.background=element_blank(),
                    panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank(),
                    legend.position='bottom',
                    legend.direction="horizontal"
                ) +
                theme_font

# plot:
plot_length_gc <-   ggplot(data=info, aes(x=Length_NUCCORE, y=GC_NUCCORE)) +
                    geom_point(shape=20, size=1, colour='#CCCCCC', alpha=0.3) +
                    geom_density_2d(colour='#333333', size=0.5) +
                    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
                    xlab('Sequence length [bp]') + ylab('Sequence GC content [%]') +
                    theme_bw() +
                    theme_font

# save to PDF
pdf(args$pdf, width=args$width, height=args$height)
print(plot_utaxa)
print(plot_year)
print(plot_locs)
print(plot_length_gc)
dev.off()
