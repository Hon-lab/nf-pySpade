#!/usr/bin/env nextflow

/*
* Set the pipeline parameters. Please prepare the following input files:
* transcriptome file from CellRanger output
* sgRNA matrix (boolean), columns are cell ID, rows are sgRNA.
* sgRNA annotation dictionary, indicating what sgRNAs target what region.
* positive control text file to test repression, region and target gene separated by tab.
* Input files details please look at: https://github.com/Hon-lab/pySpade 
* 
* Author: Yihan Wang
* 
*/

//Change this part for different analysis 
params.transcriptome = "/project/GCRB/Hon_lab/s215194/Single_Cell/mini_cm_pilot_02_inducible/_421-424_inducible/_aggr_transcriptome/cm_pilot_mini_02_no_normalization_inducible/outs/count/"
params.sgrna_df = "/project/GCRB/Hon_lab/s215194/Single_Cell/mini_cm_pilot_02_inducible/_421-424_inducible/aggr_dataframe/aggr_combined_df_full.pkl"
//params.sgrna_dict = "/project/GCRB/Hon_lab/s426305/Analysis/IGVF/20240605_WTC11_CMPilot2_PB9/pySpade/sgRNA_dict_hg38.txt"
params.sgrna_dict = "/project/GCRB/Hon_lab/s426305/nextflow_test/subset_sgRNA_dict_hg38.txt"
//params.fc_query = "/project/GCRB/Hon_lab/s426305/Analysis/IGVF/20240605_H9_CMPilot2_inducible/pySpade/promoter_region2.txt"
params.fc_query = '/project/GCRB/Hon_lab/s426305/nextflow_test/promoter_region.txt'

/*
* Default parameters, change if needed.
*/
params.outdir = "./"
params.FDR = 0.1
//FDR = 0.1: 10% false discovery rate.
params.fold_change_cutoff = 0.2
//fold change cutoff = 0.2: fold change must be more than 20%, both up-regulated and down-regulated. 
params.expression_cutoff = 0.05
//expression cutoff = 0.05: genes expressed in more than 5% of cells.

/*
* Test function
*/
//If there are too many perturbation regions (sgrna_dict), the regions can be processed in parallel.
//params.size = 50 

/*
* pySpade functions
*/

process pySpadeprocess {
	executor "slurm"
	queue '256GB,256GBv1,384GB,512GB'
 	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		path transcriptome
		path sgrna_df
		path outdir

	output:
		val 'process_ready'

	script:
		"""
		pySpade process -f $transcriptome/\
						-s $sgrna_df\
						-o $outdir/
		"""
}

process pySpadefc {
    executor "slurm"
    queue '256GB,256GBv1,384GB,512GB'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'process_ready'
        path outdir
		path sgrna_dict
		path fc_query

	output: 
		path "$outdir/pySpade_fc/fold_change.txt"

	script: 
		"""
		mkdir $outdir/pySpade_fc
		pySpade fc -t $outdir/ \
				-d $sgrna_dict \
				-r $fc_query \
				-o $outdir/pySpade_fc/
		"""
}

process randomized_sgrnadf {
	executor "slurm"
    queue 'super'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'process_ready'
		path outdir

	output:
		val 'randomized_sgrna_df_ready'

	script: 
		"""
		mkdir $outdir/FDR
		$outdir/script/randomized_sgrna.py -s $outdir/Singlet_sgRNA_df.h5 \
									 -o $outdir/FDR/Randomized_sgrna_df.h5
		"""
}

process pySpadeDEobs {
	executor "slurm"
    queue '256GB,256GBv1,384GB,512GB'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'process_ready'
		path outdir
		path sgrna_dict
	
	output:
		val 'DEobs_ready'

	script:
		"""
		mkdir $outdir/DEobs/
		pySpade DEobs\
				-t $outdir/Singlet_sub_df.h5\
				-s $outdir/Singlet_sgRNA_df.h5\
				-d $sgrna_dict\
				-n 'cpm'\
				-o $outdir/DEobs/
		"""
}

process pySpadeDEobsFDR {
	executor "slurm"
    queue '256GB,256GBv1,384GB,512GB'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'randomized_sgrna_df_ready'
		path outdir
		path sgrna_dict
	
	output:
		val 'FDR_DEobs_ready'

	script:
		"""
		mkdir $outdir/FDR/DEobs/
		pySpade DEobs\
				-t $outdir/Singlet_sub_df.h5\
				-s $outdir/FDR/Randomized_sgrna_df.h5\
				-d $sgrna_dict\
				-n 'cpm'\
				-o $outdir/FDR/DEobs/
		"""
}

process findDErandRange {
	executor "slurm"
    queue 'super'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'DEobs_ready'
		path outdir
		path sgrna_dict

	output:
		path "$outdir/bin.txt"

	script:
		"""
		mkdir $outdir/DErand/
		$outdir/script/find_DErand_range.py \
				-d $outdir/DEobs/ \
				-s $sgrna_dict \
				-o $outdir/
		"""
}

process pySpadeDErand{
	executor "slurm"
    queue '256GB,256GBv1,384GB,512GB'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		each NUM
		path outdir
		path sgrna_dict

	output:
		val 'DErand_ready', emit: DErand_ready

	script:
		"""
		echo $NUM
		pySpade DErand\
				-t $outdir/Singlet_sub_df.h5\
				-s $outdir/Singlet_sgRNA_df.h5\
				-d $sgrna_dict\
				-n 'cpm'\
				-a 'equal'\
				-o $outdir/DErand/\
				-m $NUM
		"""
}

process pySpadelocal{
	executor "slurm"
    queue 'super'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'DEobs_ready'
		val DErand_ready_list
		path outdir
		path sgrna_dict

	script:
		"""
		pySpade local\
				-f $outdir/ \
				-d $outdir/DEobs/ \
				-s $sgrna_dict \
				-t $outdir/DErand/ \
				-o $outdir/
		"""
}

process pySpadeglobal{
	executor "slurm"
    queue '256GB,256GBv1,384GB,512GB'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'DEobs_ready'
		val DErand_ready_list
		path outdir
		path sgrna_dict

	output:
		val 'global_ready'

	script:
		"""
		pySpade global\
				-f $outdir/ \
				-d $outdir/DEobs/ \
				-s $sgrna_dict \
				-t $outdir/DErand/ \
				-o $outdir/
		"""
}

process pySpadeFDRglobal{
	executor "slurm"
    queue '256GB,256GBv1,384GB,512GB'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'FDR_DEobs_ready'
		val DErand_ready_list
		path outdir
		path sgrna_dict

	output:
		val 'FDR_global_ready'

	script:
		"""
		pySpade global\
				-f $outdir/ \
				-d $outdir/FDR/DEobs/ \
				-s $sgrna_dict \
				-t $outdir/DErand/ \
				-o $outdir/FDR/
		"""
}

process calculateFDR{
	executor "slurm"
    queue 'super'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		val 'global_ready'
		val 'FDR_global_ready'
		path outdir
		val FDR
		val fold_change_cutoff
		val expression_cutoff

	output:
		path "$outdir/Significance_score_cutoff_FDR.txt"

	script:
		"""
		mkdir $outdir/Manhattan_plots
		$outdir/script/calculate_FDR.py \
			-f $outdir/ \
			-g $outdir/unfiltered_global_df.csv \
			-r $outdir/FDR/unfiltered_global_df.csv \
			-c $FDR \
			-cf $fold_change_cutoff \
			-cx $expression_cutoff \
			-o $outdir/
		"""
}

process pySpadeManhattan{
	executor "slurm"
    queue 'super'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		path outdir
		val expression_cutoff
		val fold_change_cutoff
		each significance_score

	output:

	script:
		"""
		echo $significance_score
		pySpade manhattan\
				-f $outdir/ \
				-g $outdir/unfiltered_global_df.csv \
				-cx $expression_cutoff \
				-cf $fold_change_cutoff \
				-cs $significance_score \
				-o $outdir/Manhattan_plots/
		"""
}

process pySpadeFilterLocal{
	executor "slurm"
    queue 'super'
	module 'singularity/3.9.9'
    container './pyspade_v0150.sif'

	input:
		path outdir
		val expression_cutoff
		val fold_change_cutoff
		each significance_score

	output:

	script:
		"""	
		echo $significance_score
		$outdir/script/filtered_local_df.py \
				-f $outdir/ \
				-l $outdir/unfiltered_local_df.csv \
				-cx $expression_cutoff \
				-cf $fold_change_cutoff \
				-cs $significance_score \
				-o $outdir/Manhattan_plots/
		"""
}

/*
* main script workflow
*/

workflow {

log.info """\

	 p y S p a d e ---- N F  P I P E L I N E
	 ========================================
	 transcriptome: ${params.transcriptome}
	 sgrna: ${params.sgrna_df}
	 sgrna_dict: ${params.sgrna_dict}
	 outdir: ${params.outdir}
	 FDR: ${params.FDR}
	 fold change cutoff: ${params.fold_change_cutoff}
	 expression cutoff: ${params.expression_cutoff}
	 """
	 .stripIndent()

//Set parameters  
transcriptome = channel.fromPath(params.transcriptome)
sgrna_df = channel.fromPath(params.sgrna_df)
sgrna_dict = channel.fromPath(params.sgrna_dict)
outdir = channel.fromPath(params.outdir)
fc_query = channel.fromPath(params.fc_query)

FDR = channel.of(params.FDR)
fold_change_cutoff = channel.of(params.fold_change_cutoff)
expression_cutoff = channel.of(params.expression_cutoff)


//Start the process 
process_ch = pySpadeprocess(transcriptome, sgrna_df, outdir)
pySpadefc(process_ch, outdir, sgrna_dict, fc_query)
randomized_ch = randomized_sgrnadf(process_ch, outdir)

//DEobs for real and randomized matrix
DEobs_ch = pySpadeDEobs(process_ch, outdir, sgrna_dict)
FDR_DEobs_ch = pySpadeDEobsFDR(randomized_ch, outdir, sgrna_dict)

//Find the cell number for DErand
bin_text = findDErandRange(DEobs_ch, outdir, sgrna_dict)
process_input = bin_text.splitText().map{it -> it.trim() as Integer}.collect().view()

DErand_ready_channel = pySpadeDErand(process_input, outdir, sgrna_dict)

pySpadelocal(DEobs_ch, DErand_ready_channel.collect(), outdir, sgrna_dict)
Global_ch = pySpadeglobal(DEobs_ch, DErand_ready_channel.collect(), outdir, sgrna_dict)
FDR_global_ch = pySpadeFDRglobal(FDR_DEobs_ch, DErand_ready_channel.collect(), outdir, sgrna_dict)

SS_FDR = calculateFDR(Global_ch, FDR_global_ch, outdir, FDR, fold_change_cutoff, expression_cutoff)
significance_score = SS_FDR.splitText().map{it -> it.trim() as Double}.collect().view()

pySpadeManhattan(outdir, expression_cutoff, fold_change_cutoff, significance_score)
pySpadeFilterLocal(outdir, expression_cutoff, fold_change_cutoff, significance_score)

}
