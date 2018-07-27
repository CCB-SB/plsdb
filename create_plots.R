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
    parser$add_argument('--path', help='Input files path', required=TRUE)
    parser$add_argument('--pattern', help='Input files pattern', required=TRUE)
    parser$add_argument('--pdf', help='Output PDF', required=TRUE)
    parser$add_argument('--width', help='Output PDF width', default=6, type='integer')
    parser$add_argument('--height', help='Output PDF height', default=6, type='integer')
    return(parser)
}

####################
# VARS
name_map  <- data.frame(
    short=c('insd', 'embl', 'refseq'),
    full=c('INSDC', 'EMBL', 'RefSeq')
)
rownames(name_map) <- name_map$short

theme_font <- theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14, face="bold"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14, face="bold")
)


####################
# Args
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

# Info tables
# files
info_files <- sort(list.files(path=args$path, pattern=sprintf('%s__.*\\.tsv$', args$pattern), full.names=TRUE))
# embedding
info <- NULL
for(info_file in info_files){
    info_name <- sub('\\.tsv', '', unlist(strsplit(basename(info_file), '__'))[2])
    if(info_name %in% rownames(name_map)){
        info_name <- name_map[info_name, 'full']
    }

    tmp <- read.csv(file=info_file, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    tmp$source <- info_name
    tmp$year   <- as.numeric(sapply(tmp$CreateDate_NUCCORE, function(x){ unlist(strsplit(x, '/'))[1] }))

    print(sprintf('Infotable for %s: %d rows', info_name, nrow(tmp)))

    if(is.null(info)){
        info <- tmp
    } else {
        info <- rbind(info, tmp)
    }
}
# plot: dates
plot_year <- ggplot(data=info, aes(x=year, fill=source)) +
             geom_histogram(binwidth=2, colour='white', position = "dodge") +
             scale_y_sqrt() +
             xlab('Year') + ylab('Number of nuccore records') +
             theme_bw() +
             theme_font
# plot: locations
plot_locs <- ggplot(data=info, aes(x=loc_lng, y=loc_lat)) +
             borders("world", colour="white", fill="#CCCCCC") + # map
             geom_point(size=1, alpha=0.5, colour='white', fill='#FF3333', shape=21) +
             facet_wrap(~source, ncol=2) +
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
                plot.background=element_blank()
            ) +
            theme_font
# plot: sequence length distribution
plot_length <-  ggplot(data=info, aes(x=source, y=Length_NUCCORE, fill=source)) +
                geom_violin(scale="width", draw_quantiles=c(0.5), size=1) +
                geom_boxplot(width=.1)
                scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
                guides(fill=FALSE) +
                xlab('') + ylab('Nucleotide sequence length [bp]') +
                theme_bw() +
                theme_font

# save to PDF
pdf(args$pdf, width=args$width, height=args$height)
print(plot_year)
print(plot_locs)
print(plot_length)
dev.off()
