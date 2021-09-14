conda activate r_env
Rscript cicero_stuff.R \
	--counts_mat /data/users/freese/mortazavi_lab/data/190628_scATAC/MB1_full_counts_1000_cells_20000_peaks.csv \
	--genome_sizes /data/users/freese/mortazavi_lab/ref/mm10/mm10.chrom.sizes \
	--gtf /data/users/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf \
	--o /data/users/freese/mortazavi_lab/data/190628_scATAC/ \
	--sample MB1
