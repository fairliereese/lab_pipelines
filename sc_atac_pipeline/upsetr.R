main <- function() { 

	load_packages()
	args <- get_args()

	infile <- args$infile
	df <- read.csv(infile)

	fname <- make_ofile_name(infile)

	png(filename = fname,
        width = 2500, height = 2500, units = "px",
        bg = "white",  res = 300)

	p <- upset(df,
	  mainbar.y.label = 'Chromatin Accessibility Intersections',
	  sets.x.label = 'Cells Per Accessible Gene')
	print(p)

	dev.off()

}

load_packages <- function () {
	suppressPackageStartupMessages(library(UpSetR))
	suppressPackageStartupMessages(library(argparse))
	return
}

get_args <- function() {

	parser <- ArgumentParser()
	parser$add_argument('-f', action='store', dest='infile', help='input file')
 	args <- parser$parse_args()
 	return(args)
}

make_ofile_name <- function(fname) {
	fname <- strsplit(fname, '.csv')[1]
	fname <- paste0(fname, '.png')
	return(fname)
}

main()