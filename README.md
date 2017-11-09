# Installation

Install the following libraries are available

	install.packages("gplots", dependencies=TRUE)

	install.packages("ggplot2", dependencies=TRUE)

	install.packages("dynamicTreeCut", dependencies=TRUE)

	source("https://bioconductor.org/biocLite.R")
		
	biocLite('impute')

	biocLite('topGO')

	biocLite('limma')

	biocLite('edgeR')

For pdf combining

	brew install homebrew/tex/pdfjam

## Local installation

You can clone the git repository, or use the archive. The installation location is expected to be `~/git/feminisation_direction`, if you choose a differente location, rename this path in all scripts.

	mkdir ~/git

	cd ~/git

	git clone https://GitHub.com/parisveltsos/feminisation_direction.git

	cd ~/git/feminisation_direction

All output is written in this folder

	mkdir output

# Run analysis scripts

Most analysis is run from the terminal. Occasionally, manually changing or commenting some lines in the scripts is needed.

## Remove rDNA genes from counts

	cd ~/git/feminisation_direction/input

	grep -vFwf rDNA_transcripts.txt count.txt > new_count.txt
	
## Generate the courtship related data and matrices, from the table of full counts

	Rscript ~/git/feminisation_direction/scripts/data_generation.r

## Run the bias analyses (courtship, sex, tissue)

Notice lines 78 and 80 in `bias_analysis.r` allow to exclude sex specific genes. They write the same output files, which can be kept in separate folders in the input folder, for quick access `sex_bias_B`, `sex_bias_B_sex_specific_removed`. The batch script either uses all types of line (E, M, B) or only B (Baseline). Baseline only was used, so that the calls are not influenced by the experimental treatment (ie EMB are for reference only).

	source ~/git/feminisation_direction/scripts/batch_bias.sh

## Run the EM contrast analyses (8 combinations of sex, tissue, courtship)

	source ~/git/feminisation_direction/scripts/batch_constrain.sh

### Compile gene information from each contrast

	cd ~/git/feminisation_direction/output/subsets

	mkdir genes

	for i in $(ls | grep -v pdf | grep -v out | grep -v panels| grep -v genes | grep -v info | grep -v tmp | grep -v chisq); do cp $i/$i\_pseudo_counts.txt genes/; done

	mkdir chisq

	for i in $(ls |grep -v pdf | grep -v out | grep -v panels| grep -v genes | grep -v chisq | grep -v tmp | grep -v info); do cp $i/$i\ chi_tbl.txt chisq/; done

	cd chisq
	
	tail -n +1 * > all_table.txt

### Make panels of interesting graphs, in same order for 8 comparisons

	cd ~/git/feminisation_direction/output/subsets

	mkdir panels

	mkdir tmp

	while read COL_1 COL_2; do echo ">> ";mkdir tmp/$COL_1; for i in $(ls | grep -v pdf | grep EM); do cp $i/$COL_2* tmp/$COL_1; done; pdfjam tmp/$COL_1/*.pdf --landscape --outfile tmp/$COL_1.pdf; pdfjam tmp/$COL_1.pdf '3,7,4,8,1,5,2,6' --nup 2x4 --outfile panels/$COL_1.pdf; done < ~/git/feminisation_direction/scripts/subset_panels.txt

If some heatmaps do not work, replace with the robust version

	pdfjam tmp/Heatmap_DE/*.pdf '1' --landscape --outfile tmp/heatmap_median.pdf; pdfjam tmp/heatmap_median.pdf '3,7,4,8,1,5,2,6' --nup 2x4 --outfile panels/heatmap_median.pdf

	pdfjam tmp/Heatmap_DE/*.pdf '2' --landscape --outfile tmp/heatmap_average.pdf; pdfjam tmp/heatmap_average.pdf '3,7,4,8,1,5,2,6' --nup 2x4 --outfile panels/heatmap_average.pdf

	pdfjam tmp/Heatmap_DE/*.pdf '3' --landscape --outfile tmp/heatmap_wardD.pdf; pdfjam tmp/heatmap_wardD.pdf '3,7,4,8,1,5,2,6' --nup 2x4 --outfile panels/heatmap_wardD.pdf

	rm -rf tmp

	pdfjam MDS_EMB_FB.pdf MDS_EMB_MB.pdf MDS_EMB_FH.pdf MDS_EMB_MH.pdf --landscape --nup 2x2 --outfile MDS_all.pdf

## GO analysis

### For subsets

Go to the folder of the relevant comparisons, and run the analysis for all contrasts

	cd ~/git/feminisation_direction/output/subsets

	for i in $(ls | grep EM | grep -v pdf | grep -v txt); do cd $i/GO_0.1_$i; source ~/git/feminisation_direction/scripts/GO_analysis_10.sh; cd ..; cd ..; done

	mkdir GO_out

	for i in $(ls | grep EM | grep -v pdf | grep -v txt); do cp $i/GO_0.1_$i/Fisher.txt GO_out/Fisher_$i.txt; cp $i/GO_0.1_$i/GO_UP/Fisher.txt GO_out/Fisher_up_$i.txt; cp $i/GO_0.1_$i/GO_DOWN/Fisher.txt GO_out/Fisher_down_$i.txt; done

### For sex bias

Go to the folder of the relevant comparisons, and run the analysis for all contrasts

	cd ~/git/feminisation_direction/output/sex_bias_analysis

	for i in $(ls | grep -v pdf | grep -v txt); do cd $i; for k in $(ls | grep GO_0.05_); do cd $k; source ~/git/feminisation_direction/scripts/GO_analysis_05.sh; cd ..; done; cd ..; done

	mkdir GO_out

	for i in $(ls | grep -v pdf | grep -v GO_out | grep -v txt); do cp $i/GO_0.05_G_M.F_$i\_B/Fisher.txt GO_out/Fisher_$i.txt; cp $i/GO_0.05_G_M.F_$i\_B/GO_UP/Fisher.txt GO_out/Fisher_up_$i.txt; cp $i/GO_0.05_G_M.F_$i\_B/GO_DOWN/Fisher.txt GO_out/Fisher_down_$i.txt; done

### For courtship bias

Go to the folder of the relevant comparisons, and run the analysis for all contrasts

	cd ~/git/feminisation_direction/output/courtship_bias_analysis

	for i in $(ls | grep -v pdf | grep -v txt); do cd $i; for k in $(ls | grep GO_0.1_); do cd $k; source ~/git/feminisation_direction/scripts/GO_analysis_10.sh; cd ..; done; cd ..; done

	mkdir GO_out

	for i in $(ls | grep -v pdf | grep -v GO_out | grep -v txt); do cp $i/GO_0.1_S_V.C_$i\_B/Fisher.txt GO_out/Fisher_$i.txt; cp $i/GO_0.1_S_V.C_$i\_B/GO_UP/Fisher.txt GO_out/Fisher_up_$i.txt; cp $i/GO_0.1_S_V.C_$i\_B/GO_DOWN/Fisher.txt GO_out/Fisher_down_$i.txt; done

## Randomisation

Copy the pseudo_count files from the 8 contrasts into an input folder. Note, more columns are present in those files, than used for randomisation.

The code should select folders `EMFCB, EMFCH, EMFVB, EMFVH, EMMCB, EMMCH, EMMVB, EMMVH`, run the `ls` part only first, to make sure.

	mkdir ~/git/feminisation_direction/input/pseudo_counts

	for i in $(ls ~/git/feminisation_direction/output/subsets/| grep -Ev 'pdf|txt|All|zip|genes|panels|chisq'); do cp ~/git/feminisation_direction/output/subsets/$i/$i\_pseudo_counts.txt ~/git/feminisation_direction/input/pseudo_counts/; done

Run the randomisation scripts. The output files appear in their own `output/results_standard_randomisation` and `output/results_balanced_randomisation` folders.

	Rscript ~/git/feminisation_direction/scripts/StandardisedRandomisation.r

	Rscript ~/git/feminisation_direction/scripts/BalancedRandomisation.r

## Comparison with the McGraw genes

The following makes venn diagrams comparing courted-virgin genes and genes up and down regulated in the McGraw results. There is a chisq test in the end that looks for enrichment.

	source ~/git/feminisation_direction/scripts/batch_mcgraw.sh

## Correlation plot

Plots correlation between homologues from the Veltsos et al and Hollis et al studies.

	source ~/git/feminisation_direction/scripts/correlation_plot.r

The plot is the `/Users/zabameos/git/feminisation_direction/output/correlations/VH_dpse_vs_dmel_EminusM_wo_outliers.pdf` file.

## Chromosomal enrichment of DE genes

Copy DE genes from sex contrasts (sexual selection genes already copied), into their own folder.

	mkdir ~/git/feminisation_direction/output/sex_bias_analysis/genes

	cd ~/git/feminisation_direction/output/sex_bias_analysis

	for i in $(ls | grep -Ev 'pdf|txt|genes|panels|out'); do cp $i/$i\_pseudo_counts.txt genes/; done

Run manually the script `DE_genes.r`

Edit line 23 to appropriate plottype (options appear as comments on the line).
	
Edit line 50-53 to appropriate plottype.

The output genes/*_FDR_chi_tbl.txt make Table 4.

# Downloaded files

`fbgn_fbtr_fbpp_fb_2017_01.tsv` to convert the Hollis microarray transcripts to genes.

The DE genes from Table S3 from Gerrard et al were converted to FBgn IDs using the Flybase ID conversion tool http://flybase.org/static_pages/downloads/IDConv.html 
