# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

# requires package optimx
# http://cran.r-project.org/web/packages/optimx/optimx.pdf
# install.packages('optimx')
library(optimx)

# directory where script is located
args <- commandArgs(trailingOnly = F)
main_dir <- path.expand(dirname(sub("--file=","",args[grep("--file",args)])))

# input data
p_in <- paste(main_dir, '/output/training_preds.map', sep='')
df   <- read.csv(p_in, sep='\t', header=T, stringsAsFactors=F)

# output
p_out <- paste(main_dir, '/output/weights.map', sep='')
out_dat <- as.data.frame(matrix(nrow=0, ncol=2))
row <- 0

# outlier and NA detection
is_outlier <- function(x){
	xy   <- x[!is.na(x)]
	madx <- mad(xy)
	devs <- abs(x - median(xy))
	outs <- ifelse(devs/madx >= 1 | is.na(x), TRUE, FALSE )
	return(outs)
}

# weighted average after outlier removal
weight_preds <- function(preds, weights, fams){
	my_pred  <- 0
	sum_weight <- 0
	outliers <- is_outlier(preds); names(outliers) <- names(preds)
	for (fam in fams){
		if (outliers[fam]) {
			next
		} else if (preds[fam] == 'NA'){
			next
		} else {
			my_pred <- my_pred + ( preds[fam] * weights[fam] )
			sum_weight <- sum_weight + weights[fam]
		}
	}
	return ( my_pred/sum_weight )
}

# median unsigned error
mue <- function(weights){
	unsigned_errors <- c()
	for (rec in ags.list){
		preds           <- rec$preds
		weighted_pred   <- weight_preds(preds, weights, fams)
		unsigned_err    <- abs(rec$truth - weighted_pred)/rec$truth
		unsigned_errors <- c(unsigned_errors, unsigned_err)
	}
	return(median(unsigned_errors))
}

# restructure data object and subset by read length
df_to_list <- function(df, read_length){
	ags.list <- list()
	x        <- df[df$read_length == read_length, ]
	for (genome in x$genome_name){
		 slice <- df[which(x$genome_name == genome),]
		 truth <- slice$true_ags[1]
		 preds <- slice$est_ags
		 names(preds) <- slice$fam
		 rec_data <- list(truth, preds); names(rec_data) <- c('truth', 'preds')
		 rec_name <- paste('x',genome,sep='')
		 ags.list[[rec_name]] <- rec_data
	}
	return(ags.list)
}

# init values
read_lengths  <- sort(unique(df$read_length))

# loop over read lengths
print("Finding optimal weights...")
for (read_length in read_lengths){

	# restructure data and subset by read length
	ags.list <- df_to_list(df, read_length)
	
	# init values
	fams          <- names(ags.list[[1]]$preds)
	lowers        <- rep(0, length(fams)); uppers <- rep(1, length(fams))
	start_weights <- rep(1/length(fams), length(fams)); names(start_weights) <- fams
	
	# find weights to minimize median unsigned error across training libraries (takes ~5 minutes per 150 training libraries)
	optimx_out <- optimx( start_weights, mue, lower=lowers, upper=uppers, method="L-BFGS-B")

	# store output
	for (fam in fams){
		row <- row + 1
		out_dat[row, 1] <- paste(read_length, fam, sep='_')
		out_dat[row, 2] <- as.numeric(optimx_out[fam])
	}		
}
	
# write output
write.table(out_dat, p_out, row.names=F, col.names=F, sep='\t', quote=F)




