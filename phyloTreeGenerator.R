
if (!require("Hmisc")) install.packages("Hmisc", repos="http://cran.us.r-project.org")
if (!require("corrplot")) install.packages("corrplot", repos="http://cran.us.r-project.org")
if (!require("optparse")) install.packages("optparse", repos="http://cran.us.r-project.org")

library(Hmisc)
library(corrplot)
library(optparse)

option_list = list(
  make_option(c("-i","--input", type="character",
                help="input data file.", metavar="character")),
  make_option(c("-o","--out", type="character",
                help="path to output directory.", metavar="character")),
  make_option(c("-n","--name", type="character",
                help="output results file name.", metavar="character")
));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(length(opt)==0){
  stop("Missing arguements. run --help to see required arguements.")
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("--vcf vcf file input is required to run the program", call.=FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("--out path to output directory is required to run the program", call.=FALSE)
}
if (is.null(opt$name)){
  print_help(opt_parser)
  stop("--name output file name is required to run the program", call.=FALSE)
}

#import numeric genotype data into R
print ("Reading data into R")
in_data <- read.table(opt$input, header = TRUE)

#use rcorr library to get correlation matrix and associated p-values for matching samples
print ("Creating correlation matrix and p-values.")
rcorr_data <- rcorr(as.matrix(in_data), type = c("pearson"))

#function to extract p-values for each matching samples
flattenCorrMatrix <- function(correlation_matrix, p_value_matrix) {
  cor_data_values <- upper.tri(correlation_matrix)
  data.frame(
    row = rownames(correlation_matrix)[row(correlation_matrix)[cor_data_values]],
    column = rownames(correlation_matrix)[col(correlation_matrix)[cor_data_values]],
    cor  =(correlation_matrix)[cor_data_values],
    p = p_value_matrix[cor_data_values]
    )
}

#output correlation values and corresponding p-values for matching samples to a csv file strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][2] == strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', y), ' ')[[1]][2]
correlation_pvalue_data <- flattenCorrMatrix(rcorr_data$r, rcorr_data$P)
print ("Writing correlation and p-values to output file")
write.csv(correlation_pvalue_data[order(-correlation_pvalue_data$cor),], paste(opt$out,"/",opt$name,"_p-values.csv", sep=""))

#create a correlation matrix for phylogenetic tree and heatmap and output into pdf files
correlation_matrix <- cor(in_data, method = "pearson", use = "complete.obs")
print("Creating heatmap...")
correlation_heatmap_filename <- paste(opt$out,"/",opt$name,"_heatmap.pdf", sep="")
#par(mar=c(5,3,10,2)+0.1)
pdf(correlation_heatmap_filename, width=10,height=10)
corrplot(correlation_matrix, title=paste(opt$name,"_heatmap",sep=""), tl.col="black", method = "square", mar=c(0,0,2,0))

#create phylogenetic tree and output to pdf file
print ("Creating phylogenetic tree")
correlation_tree_filename <- paste(opt$out,"/",opt$name,"_tree.pdf", sep="")
correlation_data<-hclust(as.dist(1-cor((in_data), method="pearson")), method="complete")
correlation_data_dendrogram<-as.dendrogram(correlation_data)
pdf(correlation_tree_filename, width=15, height=10)
nodePar <- list(lab.cex = 1.0, pch = c(NA, 19), cex = 1.0, col = "blue")
par(mar = c(20,4,2,2) + 0.1)
plot(correlation_data_dendrogram, cex = 1.0,ylim = c(0,1), main = paste(opt$name,"_Cluster_Dendrogram", sep=""),ylab = "Distance (1-Pearson correlation)",xlab = "", nodePar=nodePar, edgePar = list(col = 2:3, lwd = 2))
abline(h=0.5, col="blue",lty=2)


