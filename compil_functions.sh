
# bin/bash

#################################################################################
#################################################################################
   # ____ ____ _  _ ___  _ _        ____ _  _ _  _ ____ ___ _ ____ _  _ ____ #
   # |    |  | |\/| |__] | |        |___ |  | |\ | |     |  | |  | |\ | [__  #
   # |___ |__| |  | |    | |___ ___ |    |__| | \| |___  |  | |__| | \| ___] #   
#################################################################################
#################################################################################

# Script name: compil_functions.sh
# Description: This script contains a collection of functions and for genomic 
#              analysis. 
#
# Usage:
#   Source this script in your shell environment or include it in your scripts
#   to use the dsired functions.
#
# Sections:
#   1. PATHS DEFINITION
#      - Defines paths to required tools, scripts, and data directories.
#        (see config.sh for details and configuration)
#
#   2. BASIC FUNCTIONS
#      - calc: Perform arithmetic calculations.
#      - join_by: Join array elements with a specified delimiter.
#      - multijoin: Join multiple files recursively.
#
#   3. FUNCTIONS
#      - download_SRA: Download SRA files using SRA Toolkit.
#      - mapping_Fastq_bowtie2: Map FASTQ files using Bowtie2.
#      - main_mapping_Fastq: Wrapper for mapping FASTQ files.
#      - peakcalling_MACS3: Perform peak calling using MACS3.
#      - peakcalling_GOPeaks: Perform peak calling using GOPeaks.
#      - main_peakcalling: Wrapper for peak calling with multiple replicates.
#      - replicate_comp: Compare replicates and generate correlation plots.
#      - pairwize_comparison: Perform pairwise comparison of datasets.
#      - analyze_decile: Analyze datasets by dividing into deciles and computing motifs and spacing.
#      - compute_rpkmrip_rpkmril: Compute RPKM for RIP and RIL conditions.
#      - compute_motif: Compute motifs using MEME and TFFM.
#      - compute_DNAShape: Compute DNA shape features.
#      - compute_kmer: Compute k-mer enrichment using GEM.
#      - compute_SeqConv: Compute sequence convolution models.
# 	   - prep_annotation: Prepare genomic annotations for analysis.
#      - compute_NS: Generate negative sets for analysis.
#      - compute_ROCS: Compute ROC curves for model evaluation.
#      - compute_distribution: Compute distribution of binding scores.
#      - compute_space: Analyze spacing and orientation between binding sites.
#      - compute_spaing2TFs: Analyze spacing and orientation between two transcription factors.
#      - compute_spacing_TFBSTSS: Analyze spacing and orientation between binding sites and TSS.
#      - add_coverage: 
#      - add_score: 
#      - cooking_methylation: Analyze methylation impact on binding sites. 
#      - cons_score: Compute conservation scores for binding sites.

IFS=$'\n\t'

################################################################################
########## PATHS CONFIGURATION #################################################
################################################################################
# Please open this file to configure properly the paths to the softwares and scripts
# defined in this file, and the paths to the data directories.
source /home/312.6-Flo_Re/312.6.1-Commun/TF_genomic_analysis/config.sh 

################################################################################
########## BASIC FUNCTIONS #####################################################
################################################################################

#-------------------------------------------------------------------------------
# usage : calc $M1/$M2 or calc $M1/$M2*M3
calc(){
	awk "BEGIN { print "$*" }"; 
}

#-------------------------------------------------------------------------------
# used to create delimited string from array
join_by(){
	local IFS="$1"; shift; echo "$*"; 
}

#-------------------------------------------------------------------------------
# join multiple files recursively
multijoin() {
	out=$1
	shift 1
	cat $1 | awk '{print $1}' > $out
	for f in $*; do join $out $f > tmp; mv tmp $out; done
}

################################################################################
########## FUNCTIONS ###########################################################
################################################################################

# -----------------------------------------------------------------------------
download_SRA(){ 
	# FUNCTION: download_SRA
	# DESCRIPTION:
	#   This function downloads SRA (Sequence Read Archive) files using the SRA Toolkit.
	#   It supports both single-ended and paired-ended reads and compresses the output
	#   into gzipped FASTQ files. The function also allows specifying the number of threads
	#   for parallel processing.
	#
	# USAGE:
	#   download_SRA -f <LIST> -o <PATH> -t [INT] -n <NAMES>
	#
	# ARGUMENTS:
	#   -f <LIST>   : A list of SRA IDs to download (required).
	#   -o <PATH>   : The output directory where the downloaded files will be stored (required).
	#   -t [INT]    : The number of threads to use for parallel processing (optional, default: 1).
	#   -n <NAMES>  : A list of names to rename the downloaded files (optional).
	#   -h, --help  : Display usage information for the function.
	#
	# DEPENDENCIES:
	#   - The SRA Toolkit must be installed and accessible via the $SRAtoolkit environment variable.
	#   - A temporary directory must be set in the $TMPDIR2 environment variable.
	#
	# NOTES:
	#   - Single-ended files are saved as <name>.fastq.gz.
	#   - Paired-ended files are saved as <name>_R1.fastq.gz and <name>_R2.fastq.gz.

	echo "
	___  ____ _ _ _ _  _ _    ____ ____ ___      ____ ____ ____
	|  \ |  | | | | |\ | |    |  | |__| |  \     [__  |__/ |__|
	|__/ |__| |_|_| | \| |___ |__| |  | |__/ ___ ___] |  \ |  |
	"

	local threads=1; 
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-f)
				local list_file=("${!2}")
				echo "-> List of SRA ID set to: ${list_file[@]}";shift 2;;
			-o)
				local dir_out=$2
				echo "-> Name of output directory set to: ${2}";shift 2;;
			-t)
				local threads=$2
				echo "-> Number of threads set to: ${2}";shift 2;;
			-n)
				local name=("${!2}")
				echo "-> Name set to: ${name[@]}";shift 2;;
			-h)
				usage download_SRA ; return;;
			--help)
				usage download_SRA ; return;;
			*)
				echo "Error in arguments"
				echo $1; usage download_SRA ; return;;
		esac
	done

	# Check for required arguments
	local Errors=0
	if [ -z $list_file ]; then echo "ERROR: -f argument needed";Errors+=1;fi
	if [ -z $dir_out ]; then echo "ERROR: -o argument needed";Errors+=1;fi
	if [ $threads -eq 1 ]; then echo "-t argument not used, using 1 thread by default";fi
	if [ $Errors -gt 0 ]; then usage download_SRA ; return 1; fi
	mkdir -p -m 774 $dir_out

	echo "[INFO] - Downloading SRA files to $dir_out"
		for ((i=0;i<${#list_file[@]};i++))
	do
		if [ ! -s $dir_out/${list_file[i]}.fastq.gz ] && [ ! -f $dir_out/${list_file[i]}_R1.fastq.gz ] ; then 
			$SRAtoolkit/prefetch `echo ${list_file[i]}` --max-size 30g --force all -O $dir_out
			$SRAtoolkit/fasterq-dump `echo ${list_file[i]}` -O $dir_out -t $TMPDIR2 --split-3 -p -e $threads

			if [ -s $dir_out/${list_file[i]}.fastq ]; then #Single-ended file

				echo "[INFO] - Downloaded ${list_file[i]} as single-ended"
				gzip -c $dir_out/${list_file[i]}.fastq > $dir_out/${name[i]}.fastq.gz
				rm $dir_out/${list_file[i]}.fastq
				rm -Rf $dir_out/${list_file[i]}
				echo "[INFO] - Download of ${list_file[i]} done, file available at $dir_out/${name[i]}.fastq.gz"
				
			else #Pair-ended files
				echo "[INFO] - Downloading ${list_file[i]} as pair-ended files"
				gzip -c $dir_out/${list_file[i]}_1.fastq > $dir_out/${name[i]}_R1.fastq.gz
				rm $dir_out/${list_file[i]}_1.fastq
				gzip -c $dir_out/${list_file[i]}_2.fastq > $dir_out/${name[i]}_R2.fastq.gz
				rm $dir_out/${list_file[i]}_2.fastq
				rm -Rf $dir_out/${list_file[i]}
				echo "[INFO] - Download of ${list_file[i]} done, files available at $dir_out/${name[i]}_R1.fastq.gz and $dir_out/${name[i]}_R2.fastq.gz"
			fi
		fi
	done
}

#-------------------------------------------------------------------------------
mapping_Fastq_bowtie2 (){ 
	# FUNCTION: mapping_Fastq_bowtie2
	# DESCRIPTION:
	#   This function processes FASTQ files for alignment using Bowtie2. It supports
	#   both single-end and paired-end sequencing data. The function performs the
	#   following steps:
	#     1. Alignment using Bowtie2.
	#     2. Filtering of SAM files to remove low-quality and multimapped reads.
	#     3. Conversion of SAM to BAM format.
	#     4. Sorting and indexing of BAM files.
	#     5. Removal of PCR duplicates.
	#     6. Generation of mapping statistics, including duplication rates.
	#
	# USAGE:
	#   mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s <INT> -pr <INT> -keepDup <"auto"/"all"/INT> -in <PATH> -size <PATH> -debug -multialign
	#
	# ARGUMENTS:
	#   -id       : Input directory containing FASTQ files.
	#   -od       : Output directory for processed files.
	#   -s        : Seed value for Bowtie2 (default: 1254).
	#   -pr       : Number of threads to use (default: 1).
	#   -keepDup  : Duplication level ("auto", "all", or an integer as in MACS3 filterdup).
	#   -in       : Path to Bowtie2 index.
	#   -size     : Path to chromosome size file.
	#   -debug    : Enable debug mode (keeps intermediate files).
	#   -multialign: Allow multimapped reads ("canbe") or filter them out ("none").
	#   -h, --help: Display usage information.
	#
	# DEPENDENCIES:
	#   - Bowtie2
	#   - SAMtools
	#   - MACS2
	#   - FASTQC
	#   - NGmerge
	#   - bedtools
	#
	# NOTES:
	#   - Paired-end files must follow the "_R1" and "_R2" naming convention.

	local proc=1; local seed=1254; local PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"; local mode="PROD"; local sizefile=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize; local multi="none"
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-id)
				local in_dir=$2
				echo "-> Data directory set to: ${2}";shift 2;;
			-od)
				local out_dir=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-s)
				local seed=$2
				echo "-> Seed set to: ${2}";shift 2;;
			-keepDup)
				local keepDup=$2
				echo "-> Level of duplication (auto, all or an integer as in MACS3 filterdup): ${2}";shift 2;;
			-pr)
				local proc=$2
				echo "-> Thread(s) number set to: ${2}";shift 2;;
			-in)
				local PATH_TO_BOWTIE_INDEX=$2
				echo "-> Bowtie index set to: ${2}"; shift 2;;
			-size)
				local sizefile=$2
				echo "-> Sizefile set to: ${2}"; shift 2;;
			-debug)
				local mode="DEBUG"; shift 1;;
			-multialign)
				local multi="canbe"; shift 1;;
			-h)
				usage mapping_Fastq_bowtie2; return;;
			--help)
				usage mapping_Fastq_bowtie2; return;;
			*)
				echo "Error in arguments"
				echo $1; usage mapping_Fastq_bowtie2; return;;
		esac
	done

	# Check for required arguments
	local Errors=0
	if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage mapping_Fastq_bowtie2; return 1; fi

	echo "[INFO] - Processing dataset $in_dir" | tee -a $in_dir/log.txt
	echo $(find $in_dir -name "*.fastq.gz") | tee -a $in_dir/log.txt

	if [ $keepDup == "all" ] | [ $keepDup == "auto" ]; then echo usleep; else keepDup=${keepDup%.*}; fi

	# Process each fastq file
	for fastq in $(find $in_dir -name "*.fastq.gz")
	do 
		if [ -e $fastq ]
		then
			echo "[INFO] - FASTQ found" | tee -a $in_dir/log.txt
			local fastq_file=${fastq##*/}
			local fastq_filename=${fastq_file%.*.gz}
			local filename=${fastq_filename}
			local order=0 #0 for single-end, 1 for R1, 2 for R2
			
			if [[ $fastq_filename == *"_R1"* ]]; then local order=1; local filename=${fastq_filename%_R1*}; fi
			if [[ $fastq_filename == *"_R2"* ]]; then local order=2; local filename=${fastq_filename%_R2*}; fi

			echo "[INFO] - Step 1: alignment using bowtie2" | tee -a $in_dir/log.txt
			local out=$out_dir/${out_dir##*/}  #$fastq_filename
			echo "[INFO] - Out (tmp): " $out | tee -a $in_dir/log.txt

			# ----------- Single End 
			if [ $order -eq 0 ] && [ ! -e "$out.sam" ]
			then
				echo "[INFO] - Single end file detected, mapping..." | tee -a $in_dir/log.txt
				if [ ! -e "$out.sam" ]; 
				then
					$FASTQCmk -t $proc -o $in_dir $fastq
					$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -U $fastq -S $out.sam -p $proc 2>&1 | tee -a $in_dir/log.txt
				fi
				local FORMAT="BAM" 
			else

			# ----------- Paired End 
			echo "[INFO] - Paired end files detected" | tee -a $in_dir/log.txt
			#figure out the order of the pair-ended files
				local FORMAT="BAMPE" 
				local paired_file_order=1
				local fastqIsFirst=false
				if [ "$order" == "$paired_file_order" ]
				then
					local paired_file_order=2
					local fastqIsFirst=true
				fi
			
				local paired_file=${fastq/"R"$order/"R"$paired_file_order} #build up the file name of the pair
				echo "[INFO] - The paired file for $fastq is $paired_file" | tee -a $in_dir/log.txt

				local out=$out_dir/${out_dir##*/} 
				echo "[INFO] - Out: " $out | tee -a $in_dir/log.txt

				local R_NGmerge=$out
				$FASTQCmk -t $proc -o $in_dir $fastq $paired_file 
				echo "[INFO] - Checking paired files..." | tee -a $in_dir/log.txt

				if [ -e ${paired_file} ] && [ ! -e "$out.sam" ]; then
					echo "[INFO] - Paired files found, mapping..." | tee -a $in_dir/log.txt
					if $fastqIsFirst #R1 is first

					then 
						$NGmerge -a -1 $fastq -2 $paired_file -u 41 -g -o $R_NGmerge -n $proc
						$FASTQCmk -t $proc -o $in_dir ${out}_1.fastq.gz ${out}_2.fastq.gz 
						$PATH_TO_BOWTIE2/bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 ${out}_1.fastq.gz -2 ${out}_2.fastq.gz -S $out.sam --dovetail -p $proc 2>&1 | tee -a $in_dir/log.txt
						echo "[INFO] - Processed file ${out}.sam" | tee -a $in_dir/log.txt

					else #R2 is first
						$NGmerge -a -2 $fastq -1 $paired_file -u 41 -g -o $R_NGmerge -n $proc
						$FASTQCmk -t $proc -o $in_dir ${out}_1.fastq.gz ${out}_2.fastq.gz 
						$PATH_TO_BOWTIE2/bowtie2 --seed $seed -x $PATH_TO_BOWTIE_INDEX -1 ${out}_1.fastq.gz -2 ${out}_2.fastq.gz -S $out.sam --dovetail -p $proc 2>&1 | tee -a $in_dir/log.txt
						echo "$PATH_TO_BOWTIE2/bowtie2 -x $PATH_TO_BOWTIE_INDEX -1 $paired_file -2 $fastq -S $fastq.sam"
						echo "[INFO] - Processed file ${out}.sam" | tee -a $in_dir/log.txt
					fi 

				else
					echo "[INFO] - Paired file for $fastq does not exist for data set $(basename $in_dir)" | tee -a $in_dir/log.txt
					continue
				fi

				#cleaning
				rm $fastq
				rm $paired_file
			fi

			echo "[INFO] - Step 2: filtering SAM" | tee -a $in_dir/log.txt
			# filter out reads with more than 2 mismatches, mapping qual <30 and with subalignement
			$PATH_TO_SAMTOOLS/samtools stats $out.sam > $out.sam_stats

			if [ $multi == "none" ]; then # filter out multimapped reads
				# removed the filtering of XS:i reads to keep primary alignements of multimapped reads.
				$PATH_TO_SAMTOOLS/samtools view -F 0x100 -Sh $out.sam | \
			grep -e "^@" -e 'XM:i:[012][^0-9]' | awk '$1~/@/ || $5>30 {print $0}' > $out.filtered.sam 

			elif [ $multi == "canbe" ]; then # do not filter out multimapped reads
				$PATH_TO_SAMTOOLS/samtools view -Sh $out.sam | \
			grep -e "^@" -e 'XM:i:[012][^0-9]' | awk '$1~/@/ || $5>30 {print $0}' > $out.filtered.sam
			else
				echo "[ERROR] - multi argument not recognized, should be none or canbe"
				usage mapping_Fastq_bowtie2; return 1;
			fi

			# cleaning
			rm $out.sam

			echo "[INFO] - Step 3: SAM to BAM conversion" | tee -a $in_dir/log.txt
			$PATH_TO_SAMTOOLS/samtools view -Sh -b $out.filtered.sam \
			> $out.filtered.bam
			rm $out.filtered.sam
			
			# get disribution of duplicates before deduplication
			$PATH_TO_SAMTOOLS/samtools sort -o $out.filtered.sorted.bam $out.filtered.bam
			globsize=$(awk -v sum=0 '{sum+=$2}END{print sum}' $sizefile)
			$PATH_TO_MACS3/macs3 filterdup -i $out.filtered.sorted.bam -f $FORMAT -g $globsize --keep-dup all -o $out.filtered.sorted.bed 
			cat $out.filtered.sorted.bed | uniq -c | sed 's/ //g' | sed 's/chr/\t/g' | cut -f 1 | sort -k1,1rn | uniq -c > $out.dup_distrib_before_deduplication.txt
			
			echo "[INFO] - Step 4: BAM indexing" | tee -a $in_dir/log.txt
			$PATH_TO_SAMTOOLS/samtools index $out.filtered.sorted.bam

			echo "[INFO] - Step 5: printing mapping statistics" | tee -a $in_dir/log.txt
			if [ $order -eq 0 ]; then
				local lines=$(zcat ${out}.fastq.gz | wc -l)
				local reads=$(expr $lines / 4)
				$PATH_TO_SAMTOOLS/samtools stats $out.filtered.sorted.bam > $out.filtered.sorted.stats
				local R1=$reads
				local RMAPPED=$(grep "reads mapped:" $out.filtered.sorted.stats | cut -f 3)
				local R2="NA"
				local MEAN1="NA"
				local SD1="NA"
				local MEAN2="NA"
				local SD2="NA"
			else
				local lines=$(zcat ${out}_1.fastq.gz | wc -l)
				local reads=$(expr $lines / 4)
				local R1=$reads
				local lines=$(zcat ${out}_2.fastq.gz | wc -l)
				local reads=$(expr $lines / 4)
				local R2=$reads
				$PATH_TO_SAMTOOLS/samtools stats $out.filtered.sorted.bam > $out.filtered.sorted.stats
				local RMAPPED=$(grep "reads mapped:" $out.filtered.sorted.stats | cut -f 3)
				echo "[INFO] - Printing insert size (PE only)" | tee -a $in_dir/log.txt
				$PATH_TO_SAMTOOLS/samtools view -f66 ${out}.filtered.sorted.bam  | cut -f9 | awk '{print sqrt($0^2)}' > ${out}.tmpIS
				local MEAN2=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${out}.tmpIS)
				local SD2=$(awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}' ${out}.tmpIS)
				rm ${out}.tmpIS
				local MEAN1=$(grep "insert size average" $out.filtered.sorted.stats | cut -f 3)
				local SD1=$(grep "insert size standard deviation" $out.filtered.sorted.stats | cut -f 3)
			fi
			
			
			echo "[INFO] - Step 6: remove PCR duplicates" | tee -a $in_dir/log.txt
			echo "genome size: $globsize" | tee -a $in_dir/log.txt
			$PATH_TO_MACS3/macs3 filterdup -i $out.filtered.sorted.bam -f $FORMAT -g $globsize --keep-dup $keepDup -o $out.filtered.sorted.dedup.bed | tee -a $in_dir/log.txt
			cat $out.filtered.sorted.dedup.bed | uniq -c | sed 's/ //g' | sed 's/chr/\t/g' | cut -f 1 | sort -k1,1rn | uniq -c > $out.dup_distrib_after_deduplication.txt
			
			# convert the the deduplicated bed to a bam to be used in compute_rpkmrip function
			awk '{print $0"\t"NR}' $out.filtered.sorted.dedup.bed | awk -v OFS="\t" '$2<0{print $1,0,$3,$4,$5,$6,$7; next}{print $1,$2,$3,$4,$5,$6,$7}' > $out.filtered.sorted.dedup.bed4
			bedtools bedtobam -i $out.filtered.sorted.dedup.bed4 -g $sizefile > $out.tmp.bam
			# redo BAM indexing
			$PATH_TO_SAMTOOLS/samtools sort -o $out.filtered.sorted.dedup.bam $out.tmp.bam
			$PATH_TO_SAMTOOLS/samtools index $out.filtered.sorted.dedup.bam
			
			echo '[INFO] - Printing minimal mapping statistics including duplication rate' | tee -a $in_dir/log.txt
			# addition of duplicates stats
			# total numbers of pieces:
			local totFrag=$(cat $out.dup_distrib_before_deduplication.txt | sed 's/^/ /g' | tr ' ' '@' | sed 's/@\+/\t/g' | sed 's/^\t//g' | awk '{ print $1, $2, $1 * $2 }' | awk -F ' ' '{sum+=$3;}END{print sum;}')
			# numbers of unique pieces:
			local uniqFrag=$(cat $out.dup_distrib_before_deduplication.txt | sed 's/^/ /g' | tr ' ' '@' | sed 's/@\+/\t/g' | sed 's/^\t//g' | awk -F '\t' '{sum+=$1;}END{print sum;}')
			# duplication rate:
			local dupRate=$(calc $(calc $totFrag-$uniqFrag)/$totFrag )
			
			echo "R1 R2 RMAPPED MEAN_IS SD_IS MEAN_IS_f66 SD_IS_f66 totFrag uniqFrag dupRate" > $out.minimal.stats
			echo $R1 $R2 $RMAPPED $MEAN1 $SD1 $MEAN2 $SD2 $totFrag $uniqFrag $dupRate >> $out.minimal.stats
			
			# cleaning
			echo "[INFO] - Step 7: removing temporary sam" | tee -a $in_dir/log.txt
			if [ $mode != "DEBUG" ]; then
				rm ${out}.filtered.bam &
				rm ${out}.filtered.sorted.bed &
				rm $out.tmp.bam &
				wait
			fi
			if [ $order -eq 0 ]; then
				rm ${out}.fastq.gz
			else
				rm ${out}_1.fastq.gz
				rm ${out}_2.fastq.gz
			fi
			date | tee -a $in_dir/log.txt
		fi
	done

}

#-------------------------------------------------------------------------------
main_mapping_Fastq(){
	# FUNCTION: main_mapping_Fastq
	# DESCRIPTION:
	#   This function performs the mapping of FASTQ files using Bowtie2. It accepts
	#   various arguments to configure the input/output directories, FASTQ files,
	#   processing parameters, and other options. The function validates the input
	#   arguments, prepares the necessary files, and calls the Bowtie2 mapping
	#   function with the specified parameters.
	#
	# USAGE:
	#   main_mapping_Fastq -fd <PATH> -md <PATH> -f1 <FILE> -f2 [FILE] -n <STRING> 
	#                      -s [INT] -pr [INT] -dedup [INT or "all" or "auto"] 
	#                      -size [FILE] -in [PATH] -debug -multialign
	#
	# ARGUMENTS:
	#   -fd <PATH>       : Path to the directory containing FASTQ files.
	#   -md <PATH>       : Path to the directory where mapping results will be stored.
	#   -f1 <FILE>       : FASTQ file for pair 1 or single-end reads.
	#   -f2 [FILE]       : FASTQ file for pair 2 (optional, for paired-end reads).
	#   -n <STRING>      : Name of the dataset (used for naming output files).
	#   -s [INT]         : Seed for random number generation (default: 1254).
	#   -pr [INT]        : Number of threads to use for processing (default: 1).
	#   -dedup [INT]     : Deduplication setting (1, "all", "auto", or other values).
	#   -size [FILE]     : Path to the size file (e.g., chromosome size file).
	#   -in [PATH]       : Path to the Bowtie2 index files.
	#   -debug           : Enable debug mode for additional logging.
	#   -multialign      : Allow multi-alignment of reads.
	#   -h, --help       : Display usage information for this function.
	# 
	# DEPENDENCIES:
	#   - Bowtie2
	#   - SAMtools
	#   - MACS2
	#   - FASTQC
	#   - NGmerge
	#   - bedtools
	#
	# NOTES:
	#   - The function validates the required arguments and exits with an error
	#     message if any mandatory arguments are missing.
	#   - If the debug mode is enabled, additional logging is performed.
	#   - The function supports both single-end and paired-end FASTQ files.

		echo '
	_  _ ____ _ _  _    _  _ ____ ___  ___  _ _  _ ____
	|\/| |__| | |\ |    |\/| |__| |__] |__] | |\ | | __
	|  | |  | | | \|    |  | |  | |    |    | | \| |__]
	'
	# main_mapping_Fastq -fd <PATH> -md <PATH> -f1 <FILE> -f2 [FILE] -n <STRING> -s [INT] -pr [INT] -dedup [INT or "all" or "auto"] -size [FILE] -in [PATH] -debug -multialign
	local Fastq2="NA"; local deduplication=1; local seed=1254; local proc=1; PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"; local mode="PROD"; sizefile=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize; local multi="none"
	local Argsline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-fd)
				local Fastq_dir=$2
				echo "-> Data directory set to: ${2}";shift 2;;
			-md)
				local Mapping_dir=$2
				echo "-> Mapping directory set to: ${2}";shift 2;;
			-f1)
				local Fastq1=$2
				echo "-> Fastq pair 1 / Fastq non paired end set to: ${2}";shift 2;;
			-f2)
				local Fastq2=$2
				echo "-> Pair 2 set to: ${2}";shift 2;;
			-n)
				local Name=$2
				echo "-> Name of Fastq directory set to: ${2}";shift 2;;
			-s)
				local seed=$2
				echo "-> Seed for random set to: ${2}";shift 2;;
			-pr)
				local proc=$2
				echo "-> Thread(s) number set to: ${2}";shift 2;;
			-dedup)
				local deduplication=$2
				echo "-> Keeping duplicate: ${2}";shift 2;;
			-size)
				local sizefile=$2
				echo "-> Sizefile set to: ${2}"; shift 2;;
			-debug)
				local mode="DEBUG"; shift 1;;
			-multialign)
				local multi="canbe"; shift 1;;
			-in)
				local PATH_TO_BOWTIE_INDEX=$2
				echo "-> Index set to: ${2}";shift 2;;
			-h)
				usage main_mapping_Fastq; return;;
			--help)
				usage main_mapping_Fastq; return;;
			*)
				echo "Error in arguments"
				echo $1; usage main_mapping_Fastq; return;;
		esac
	done
	
	echo "==================="

	# Check for required arguments
	local Errors=0
	if [ -z $proc ]; then echo "-pr argument not used, using 1 processor"; fi
	if [ -z $seed ]; then echo "-s argument not used, using default seed for random (1254)"; fi
	if [ -z $Fastq2 ]; then echo "-f2 argument not used, assuming non paired-end analysis"; fi
	if [ -z $Fastq_dir ]; then echo "ERROR: -fd argument needed"; Errors+=1; fi
	if [ -z $Mapping_dir ]; then echo "ERROR: -md argument needed"; Errors+=1; fi
	if [ -z $Fastq1 ]; then echo "ERROR: -f1 argument needed"; Errors+=1; fi
	if [ -z $Name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage main_mapping_Fastq; return 1; fi

	echo "[INFO] - Searching for FASTQ input files in $Fastq_dir"
	
	# If mapping was already done
	if [[ $Fastq1 == *".done"* ]]; then
		local tmpName=$Fastq2
		if [[ -f $Fastq_dir/${tmpName}_R1.fastq.gz ]]; then
			local Fastq1=$Fastq_dir/${tmpName}_R1.fastq.gz
			if [[ -f $Fastq_dir/${tmpName}_R2.fastq.gz ]]; then
				local Fastq2=$Fastq_dir/${tmpName}_R2.fastq.gz
			fi
		elif [[ -f ${Fastq1%.*}_R1.fastq.gz ]]; then
			if [[ -f ${Fastq1%.*}_R2.fastq.gz ]]; then
				local Fastq2=${Fastq1%.*}_R2.fastq.gz
			fi
			local Fastq1=${Fastq1%.*}_R1.fastq.gz
		fi
		if [[ -f $Fastq_dir/${tmpName}.fastq.gz ]]; then
			local Fastq1=$Fastq_dir/${tmpName}.fastq.gz
		elif [[ -f ${Fastq1%.*}.fastq.gz ]]; then
			local Fastq1=${Fastq1%.*}.fastq.gz
		fi
		if [[ $Fastq1 == *".done"* ]]; then
			echo "ERROR - Fastq not found"
			return
		fi
	fi

	echo $Fastq1
	echo $Fastq2

	echo "==================="

	if [ -d $Mapping_dir/$Name ]; then
	rm -Rf $Mapping_dir/$Name
	fi
	mkdir -m 774 -p $Mapping_dir/$Name

	echo "[INFO] - Starting mapping for dataset $Name" > $Mapping_dir/$Name/log.txt
	echo "main_mapping_Fastq" $Argsline >> $Mapping_dir/$Name/log.txt

	# Copy fastq files to mapping directory
	if [ $(dirname $Fastq1) != $Mapping_dir/$Name ];then
		if [[ $Fastq1 == *"_R1"* ]]; then
			cp $Fastq1 $Mapping_dir/$Name/${Name}_R1.fastq.gz
		else
			cp $Fastq1 $Mapping_dir/$Name/${Name}.fastq.gz
		fi
	fi
	if [ $Fastq2 != "NA" ]; then
		if [ $(dirname $Fastq2) != $Mapping_dir/$Name ];then
			if [[ $Fastq1 == *"_R1"* ]]; then
				cp $Fastq2 $Mapping_dir/$Name/${Name}_R2.fastq.gz
			fi
		fi
	fi

	echo "[INFO] - Calling mapping_Fastq_bowtie2"
	#mapping_Fastq_bowtie2 -id <PATH> -od <PATH> -s <INT> -pr <INT>
	if [ $mode == "DEBUG" ]; then
		if [ $multi == "canbe" ]; then #reads multialignment authorized or not 
			mapping_Fastq_bowtie2 -id $Mapping_dir/$Name -od $Mapping_dir/$Name -s $seed -pr $proc -keepDup $deduplication -in ${PATH_TO_BOWTIE_INDEX} -size $sizefile -debug -multialign
		else
			mapping_Fastq_bowtie2 -id $Mapping_dir/$Name -od $Mapping_dir/$Name -s $seed -pr $proc -keepDup $deduplication -in ${PATH_TO_BOWTIE_INDEX} -size $sizefile -debug
		fi
	else
		if [ $multi == "canbe" ]; then #reads multialignment authorized or not 
			mapping_Fastq_bowtie2 -id $Mapping_dir/$Name -od $Mapping_dir/$Name -s $seed -pr $proc -keepDup $deduplication -in ${PATH_TO_BOWTIE_INDEX} -size $sizefile -multialign
		else
			mapping_Fastq_bowtie2 -id $Mapping_dir/$Name -od $Mapping_dir/$Name -s $seed -pr $proc -keepDup $deduplication -in ${PATH_TO_BOWTIE_INDEX} -size $sizefile
		fi
	fi
}


#-------------------------------------------------------------------------------
peakcalling_MACS3(){
	echo "

	_  _ ____ ____ ____   ____
	|\/| |__| |    [__    ___|
	|  | |  | |___ ___]   ___|

	"
	# peakcalling_MACS3 -id <PATH> -od <PATH> -g <INT> -s <INT> -phs <INT> -t <INT> -kd <STRING> -bl <PATH> -size <PATH> -cs <PATH> -r <BOOL> --llocal [INT] --slocal [INT] --factork [INT] --qval [FLOAT] --minlen [INT] --diag [BOOL] --minfold [INT] --maxfold [INT]
	# defining default parameters
	local Genome_length=120000000; local seedrandom=168159; local phs=200; local threads=1; local keepDup=1; local blacklist="/home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Greenscreen_19012023_merged.bed"; local sizeFile="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize"; local controlstats="NA"; local llocal=5000; local slocal=500; local factork=3; local qval=0.05; local minlen=100; local diag=false; local minfold=5; local maxfold=50
	# parsing arguments
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-id)
				local in_dir=$2
				echo "-> Data directory set to: ${2}";shift 2;;
			-od)
				local out_dir=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-g)
				local Genome_length=$2
				echo "-> Genome length (mappable) set to: ${2}";shift 2;;
			-s)
				local seedrandom=$2
				echo "-> Random seed set to: ${2}"; shift 2;;
			-phs)
				local phs=$2
				echo "-> Peaks half size set to: ${2}"; shift 2;;
			-t)
				local threads=$2
				echo "-> Threads set to: ${2}"; shift 2;;
			-cs)
				local controlstats=$2
				echo "-> Control stats file set to: ${2}"; shift 2;;
			-kd) 
				local keepDup=$2
				echo "-> Level of duplication (auto, all or an integer as in MACS3 filterdup): ${2}"; shift 2;;
			-bl)
				local blacklist=$2
				echo "-> Blacklist set to: ${2}"; shift 2;;
			-size)
				local sizeFile=$2
				echo "-> SizeFile set to: ${2}"; shift 2;;
			-r)
				local redo_analysis=$2
				echo "-> redo_analysis set to (true/false): ${2}"; shift 2;;
			--llocal)
				local llocal=$2
				echo "-> llocal set to: $2"; shift 2;;
			--slocal)
				local slocal=$2
				echo "-> slocal set to: $2"; shift 2;;
			--factork) 
				local factork=$2
				echo "-> factork set to: $2"; shift 2;;
			--qval) 
				local qval=$2
				echo "-> qval set to: $2"; shift 2;;
			--minlen) 
				local minlen=$2
				echo "-> minlen set to: $2"; shift 2;;
			--minfold)
				local minfold=$2
				echo "-> minfold set to: $2"; shift 2;;
			--maxfold)
				local maxfold=$2
				echo "-> maxfold set to: $2"; shift 2;;
			-diag)
				local diag=$2
				echo "-> Diagnostic set to (true/false): $2"; shift 2;;
			-h)
				usage peakcalling_MACS3; return;;
			--help)
				usage peakcalling_MACS3; return;;
			*)
				echo "[ERROR] - in argument $1 $2";
				usage peakcalling_MACS3; return;;
		esac
	done
	# checking arguments for errors
	local Errors=0
	if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage peakcalling_MACS3; return 1; fi
	echo "[INFO] - Processing folder: ${in_dir##*/}" # treating replicate folder
	local out_dir=${out_dir}/${in_dir##*/}
	
	local folderdone=true
	# checking if replicate folder already has essential results files (Skipping the treatment of this folder if this is true
	if [ "$redo_analysis" == "true" ]; then 
		local folderdone=false
		echo "[INFO] - Redoing replicate analysis as requested"
	elif [ ! -f ${out_dir}/${in_dir##*/}_peaks.bed ] || [ ! -f ${out_dir}/${in_dir##*/}_cov.bdg ]; then
		local folderdone=false
		echo "[WARNING] - Missing output files, (re)doing replicate analysis"
	fi

	echo "[INFO] - Is folder done ? ${folderdone}"
	mkdir -m 774 -p $out_dir

	local log=$out_dir/${in_dir##*/}.log;
	echo "" > $log

	# Saving Arguments used for this analysis
	echo $@ >> $log
	local controls=()
	local replicates=()
	# creating list of bam files for multiple controls and multiple replicates at once (usually only one replicate here as multiple replicate are dealt with mspc)
	for bam in $(find $in_dir -name "*.filtered.sorted.bam")
	do
		if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" == *"control"* ]]; then
			controls+=("$bam")
		else
			replicates+=("$bam")
		fi
	done
	local all_ok=true;

	if [ "$folderdone" == "false" ]; then
		mkdir -p -m 774 $out_dir/temp

		local R2_value="NA"
		local R2_value_control="NA"	
		for bam in $(find $in_dir -name "*.filtered.sorted.bam")
		do
			if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
				R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${bam%.filtered.sorted.bam}.minimal.stats)
			fi
		done
		if [ $controlstats == "NA" ]; then
			local R2_value_control=$R2_value
		else
			R2_value_control=$(awk -v FS=" " 'NR>1{print $2}' $controlstats)
		fi
		if [ $R2_value != "NA" ]; then
			local formatBed="BEDPE"
			local formatBAM="BAMPE"
		else
			local formatBed="BED"
			local formatBAM="BAM"
		fi
		echo "[INFO] - Preparing BAM files from sample and control..." | tee -a $log
		$PATH_TO_MACS3/macs3 randsample -i ${replicates[@]} -f $formatBAM -p 100 -o $out_dir/temp/replicate.bed >> $log 2>&1;
		if [ $R2_value_control != "NA" ]; then
			local formatBAM="BAMPE"
		else
			local formatBAM="BAM"
		fi

		if [ ${#controls[@]} -eq 0 ]; then
			all_ok=false
			controls=("NA")
			echo "[WARNING] - No control found for data set ${in_dir##*/}, running MACS3 without control" | tee -a $log
		else
			$PATH_TO_MACS3/macs3 randsample -i ${controls[@]} -f $formatBAM -p 100 -o $out_dir/temp/control.bed >> $log 2>&1;
		fi
		
		

		if [ $all_ok == true ]; then
			echo "[INFO] - Started MACS with controls" | tee -a $log
			# call summit parameters is used by default as the summit of peaks are very usefull data for MSPC analysis afterwards
			$PATH_TO_MACS3/macs3 callpeak -t $out_dir/temp/replicate.bed -c $out_dir/temp/control.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --keep-dup $keepDup --outdir $out_dir --llocal $llocal --slocal $slocal --min-length $minlen -q $qval --seed $seedrandom --mfold $minfold $maxfold --call-summits >> $log 2>&1;

			if [ "$diag" == "true" ]; then
				# controlling slocal and llocal used to confirm visually (ex on IGB) that none is overload the biais in the earlier callpeak.
				#first fragment length is computed
				echo "[INFO] - Diagnostic mode active, manual analysis in progress for llocal and slocal evaluation" | tee -a $log
				echo "[STEP 1] - Computing fragment length from control" | tee -a $log
				$PATH_TO_MACS3/macs3 predictd -f BEDPE -i  $out_dir/temp/control.bed > $out_dir/temp/tmplog.txt 2>&1;
				local fragment_length=$(cat $out_dir/temp/tmplog.txt | grep "Average insertion length of all pairs " | awk -v OFS="\t" '{print $17}')
				local fragment_length_halved=$(cat $out_dir/temp/tmplog.txt | grep "Average insertion length of all pairs " | awk -v OFS="\t" '{print int($17/2+0.5)}')
				# d background
				echo "[STEP 2] - Computing d background with fragment length of $fragment_length" | tee -a $log
				$PATH_TO_MACS3/macs3 pileup -i  $out_dir/temp/control.bed -B --extsize ${fragment_length_halved} -o $out_dir/d_bg.bdg -f BEDPE
				
				# small local background (slocal)
				echo "[STEP 3] - Computing slocal background with slocal of $slocal" | tee -a $log
				$PATH_TO_MACS3/macs3 pileup -i  $out_dir/temp/control.bed -B --extsize $slocal -o $out_dir/temp/slocal_bg.bdg -f BEDPE
				# normalised to be comparable to d background
				echo "[STEP 3b] - Normalising slocal background to be comparable to d background" | tee -a $log
				local normslocal=$(calc $fragment_length/$slocal)
				$PATH_TO_MACS3/macs3 bdgopt -i $out_dir/temp/slocal_bg.bdg -m multiply -p $normslocal -o $out_dir/slocal_bg_norm.bdg

				# large local background (llocal)
				echo "[STEP 4] - Computing llocal background with llocal of $llocal" | tee -a $log
				$PATH_TO_MACS3/macs3 pileup -i  $out_dir/temp/control.bed -B --extsize $llocal -o $out_dir/temp/llocal_bg.bdg -f BEDPE
				# normalised to be comparable to d background
				echo "[STEP 4b] - Normalising llocal background to be comparable to d background" | tee -a $log
				local normllocal=$(calc $fragment_length/$llocal)
				$PATH_TO_MACS3/macs3 bdgopt -i $out_dir/temp/llocal_bg.bdg -m multiply -p $normllocal -o $out_dir/llocal_bg_norm.bdg

				# merging all backgrounds by setteing the final background value as the maximum of the 3 in each genomic position
				echo "[STEP 5] - Merging d background, slocal and llocal normalised backgrounds" | tee -a $log
				$PATH_TO_MACS3/macs3 bdgcmp -m max -t $out_dir/slocal_bg_norm.bdg -c $out_dir/llocal_bg_norm.bdg -o $out_dir/sllocal_bg_norm.bdg
				$PATH_TO_MACS3/macs3 bdgcmp -m max -t $out_dir/sllocal_bg_norm.bdg -c $out_dir/d_bg.bdg -o $out_dir/d_sllocal_bg_norm.bdg
				$bdg2bwig $out_dir/d_sllocal_bg_norm.bdg $sizeFile $out_dir/d_sllocal_bg_norm.bw
			fi

		elif [ ${controls[0]} == "NA" ]; then
			echo "[INFO] - Started MACS without controls" | tee -a $log
			$PATH_TO_MACS3/macs3 callpeak -t $out_dir/temp/replicate.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --keep-dup $keepDup --outdir $out_dir --seed $seedrandom --call-summits --llocal $llocal --slocal $slocal --min-length $minlen -q $qval >> $log 2>&1;
		else
		   echo "[ERROR] - Not all needed input present." | tee -a $log
		   return 1
		fi
		if [ $formatBed == "BEDPE" ]; then
			fragment_length=$(cat $log | grep "fragment size = " | awk -v OFS="\t" '{print $13}')
		else
			fragment_length=$(cat $log | grep "predicted fragment length is" | awk -v OFS="\t" '{print $14}') 
		fi
		echo "[INFO] - Recomputing coverages along the genome" | tee -a $log
		# Using size of fragment and bam file to recompute coverages along the genome.
		# We do not use the bdegraph generated by macs3 as they are normalize by size of smallest bam
		for bam in $(find $in_dir -name "*.filtered.sorted.dedup.bam")
		do
			if [[ ! -L $bam && -f $bam ]]; then # checking if file is not a link and exists
				if [ ${R2_value} == "NA" ]; then # checking for PE, or single ended reads
					echo "[INFO] - Processing SE reads for coverage: $bam" | tee -a $log
					if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
						bedtools bamtobed -i $bam > $in_dir/${in_dir##*/}.bamtobed.bed
						awk -v fraglen=${fragment_length} -v OFS="\t" '$6=="+"{print $1,$2,$3+fraglen,$4,$5,$6}$6=="-"{print $1,$2-fraglen,$3,$4,$5,$6}' $in_dir/${in_dir##*/}.bamtobed.bed | awk -v OFS="\t" '$2<=0{print $1,1,$3,$4,$5,$6;next}{print $0}'  > $in_dir/${in_dir##*/}.ext.bed
						if [ -f ${bam%.filtered*}.minimal.stats ]; then
							echo "[INFO] - Using ${bam%.filtered*}.minimal.stats to get library size"
							local libsize=$(awk 'NR>1{print $9}' ${bam%.filtered*}.minimal.stats)
						else
							local libsize=$($PATH_TO_SAMTOOLS/samtools view -f 0 -c $bam)
						fi
						local scale=$(calc 1000000/$libsize)
						
						tmpsave=$LC_COLLATE;LC_COLLATE=C # to avoid problems with different locale using sort
						bedtools genomecov -bga -scale $scale -i $in_dir/${in_dir##*/}.ext.bed -g $sizeFile | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_cov.bdg
						LC_COLLATE=$tmpsave
						if [ -f $in_dir/${in_dir##*/}.bamtobed.bed ]; then rm $in_dir/${in_dir##*/}.bamtobed.bed; fi
						if [ -f $in_dir/${in_dir##*/}.ext.bed ]; then rm $in_dir/${in_dir##*/}.ext.bed; fi
					fi
				else
					echo "[INFO] - Processing PE reads for coverage: $bam" | tee -a $log
					if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
						if [ -f ${bam%.filtered*}.minimal.stats ]; then
							local libsize=$(awk 'NR>1{print $9}' ${bam%.filtered*}.minimal.stats)
						else
							local libsize=$($PATH_TO_SAMTOOLS/samtools view -f 0 -c $bam)
						fi
						local scale=$(calc 1000000/$libsize)
						tmpsave=$LC_COLLATE;LC_COLLATE=C
						bedtools genomecov -bga -scale $scale -ibam $bam | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_cov.bdg
						LC_COLLATE=$tmpsave
					fi
				fi
				if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
					echo "[INFO] - Checking for chr/Chr problem in $bam" | tee -a $log
					if [ $(head -1 $out_dir/${in_dir##*/}_cov.bdg | grep "Chr" | wc -l) -gt 0 ] && [ $(head -1 $sizeFile | grep "chr" | wc -l) -gt 0 ]; then
						sed -i 's/Chr/chr/g' $out_dir/${in_dir##*/}_cov.bdg
					elif [ $(head -1 $out_dir/${in_dir##*/}_cov.bdg | grep "chr" | wc -l) -gt 0 ] && [ $(head -1 $sizeFile | grep "Chr" | wc -l) -gt 0 ]; then
						sed -i 's/chr/Chr/g' $out_dir/${in_dir##*/}_cov.bdg
					fi
					# converting to bigwig
					$bdg2bwig $out_dir/${in_dir##*/}_cov.bdg $sizeFile $out_dir/${in_dir##*/}_cpm.bw
				fi
			fi
		done
		local nbpeaks=$(wc -l $out_dir/${in_dir##*/}_peaks.narrowPeak | awk '{print $1}')
		# checking if peaks are detected, as next computation doesn't like no peaks cases
		if [ $nbpeaks -eq 0 ]; then
			echo "[ERROR] - No peaks in ${in_dir##*/}_peaks.narrowPeak created, MACS3 must have FAIL to compute peaks" | tee -a $log
		elif [ ${#controls[@]} -eq 0 ]; then
			if [ -f $blacklist ]; then
				echo "[INFO] - No controls but a blacklist has been submitted, filtering using these" | tee -a $log
				bedtools intersect -v -a $out_dir/${in_dir##*/}_peaks.narrowPeak -b $blacklist -wa > $out_dir/${in_dir##*/}_filtered.narrowPeak
			else
				echo "[WARNING] - No controls or blacklist submitted... filtering step omitted" | tee -a $log
				cp $out_dir/${in_dir##*/}_peaks.narrowPeak $out_dir/${in_dir##*/}_filtered.narrowPeak
			fi
			local basename=${in_dir##*/}
		else
			echo "[INFO] - Filtering out peaks detected in input using $factork on 7th column of narrowpeaks from MACS3" | tee -a $log
			# filtering out peaks that are not k times higher in control than in sample.
			sort -k1,1 -k2,2n $out_dir/${in_dir##*/}_peaks.narrowPeak | awk -v k=$factork -v OFS='\t' '$7>k{print $0}' > $out_dir/temp/${in_dir##*/}_filteredtmp1.narrowPeak
			# removing blacklit
			if [ -f $blacklist ]; then
				echo "[INFO] - Filtering out blacklist regions using $blacklist" | tee -a $log
				cat $blacklist | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_forbidden_regions.bed
				bedtools intersect -v -a $out_dir/temp/${in_dir##*/}_filteredtmp1.narrowPeak -b $out_dir/${in_dir##*/}_forbidden_regions.bed -wa | awk -v OFS="\t" '{print $0}' > $out_dir/temp/${in_dir##*/}_filteredtmp2.narrowPeak
			else
				echo "[INFO] - No blacklist submitted - skipping this step" | tee -a $log
				cat $out_dir/temp/${in_dir##*/}_filteredtmp1.narrowPeak > $out_dir/temp/${in_dir##*/}_filteredtmp2.narrowPeak
			fi
			# computing control coverages on the peaks
			controlbdg=()
			for ctrlbam in $(find $in_dir -name "control_*.filtered.sorted.bam")
			do
				# this whole section remove peaks that have k times mean genome coverage of input
				bamname=${ctrlbam##*/}
				## Scaling factor for single-end data, counting every mapped read (bitwise flag = 0)
				ScaleTotalMappedReads=$(bc <<< "scale=6;1000000/$($PATH_TO_SAMTOOLS/samtools view -f 0 -c $ctrlbam)")
				tmpsave=$LC_COLLATE;LC_COLLATE=C
				bedtools genomecov -bga -ibam $ctrlbam -scale $ScaleTotalMappedReads | sort -k1,1 -k2,2n > $out_dir/temp/${bamname%.filtered.sorted.bam}.bdg
				LC_COLLATE=$tmpsave
				controlbdg+=("$out_dir/temp/${bamname%.filtered.sorted.bam}.bdg")
			done
			
			if [ ${#controlbdg[@]} -lt 2 ]; then
				toreformat=${controlbdg[0]}
			else
				$PATH_TO_MACS3/macs3 cmbreps -i ${controlbdg[@]} -m max -o $out_dir/temp/control_mean.bdg
				toreformat=$out_dir/temp/control_mean.bdg
			fi
			$bdg2bwig $toreformat $sizeFile $out_dir/control.bw
			
			bedtools map -a $out_dir/temp/${in_dir##*/}_filteredtmp2.narrowPeak -b $toreformat -c 4,4 -o max,mean > $out_dir/temp/${in_dir##*/}_filtered_ctrlcov.bdg
			bedtools map -a $out_dir/temp/${in_dir##*/}_filteredtmp2.narrowPeak -b $out_dir/${in_dir##*/}_cov.bdg -c 4,4 -o max,mean |awk -v OFS="\t" '{print $11,$12}' > $out_dir/temp/${in_dir##*/}_filtered_cov.bdg

			paste $out_dir/temp/${in_dir##*/}_filtered_ctrlcov.bdg $out_dir/temp/${in_dir##*/}_filtered_cov.bdg | awk -v OFS="\t" -v k=$factork '{if(($11*k)<$13&&($12*k)<$14){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > $out_dir/${in_dir##*/}_filtered.narrowPeak
			echo "~==================================================~" | tee -a $log
			echo "[REPORT] - Initial number of peaks detected by MACS3: $(wc -l $out_dir/${in_dir##*/}_peaks.narrowPeak | awk '{print $1}') " | tee -a $log
			echo "[REPORT] - Peaks after filtering on MACS3 SignalValue using th=$factork $(wc -l $out_dir/temp/${in_dir##*/}_filteredtmp1.narrowPeak | awk '{print $1}')" | tee -a $log
			echo "[REPORT] - Peaks after filtering out forbidden region using provided blacklist $(wc -l $out_dir/temp/${in_dir##*/}_filteredtmp2.narrowPeak | awk '{print $1}')" | tee -a $log
			echo "[REPORT] - Peaks after final filtering using custom computation of control coverages $(wc -l $out_dir/${in_dir##*/}_filtered.narrowPeak | awk '{print $1}')" | tee -a $log
			echo "~==================================================~" | tee -a $log
			## peaks to be used in case of multiple replicates per sample
			awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' $out_dir/${in_dir##*/}_filtered.narrowPeak > $out_dir/${in_dir##*/}_peaks.bed

			## peaks used in case there is no replicates
			awk -v OFS="\t" '{print $1, $2+$10-"'$phs'", $2+$10+"'$phs'"}' $out_dir/${in_dir##*/}_filtered.narrowPeak > $out_dir/${in_dir##*/}_narrow.bed
		fi
		rm -Rf $out_dir/temp
	fi
}

#-------------------------------------------------------------------------------
peakcalling_GOPeaks(){
	echo "

	____ ____ ___  ____ ____ _  _ ____
	| __ |  | |__] |___ |__| |_/  [__ 
	|__] |__| |    |___ |  | | \_ ___]

	"
	# peakcalling_GOPeaks -id <PATH> -od <PATH> -s <INT> -phs <INT> -t <INT> -size <PATH>
	# defining default parameters
	local seedrandom=6318; local phs=200; local threads=2; local sizeFile="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize"
	# parsing arguments
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-id)
				local in_dir=$2
				echo "-> data directory set to: ${2}";shift 2;;
			-od)
				local out_dir=$2
				echo "-> output directory set to: ${2}";shift 2;;
			-s)
				local seedrandom=$2
				echo "-> random seed set to: ${2}"; shift 2;;
			-phs)
				local phs=$2
				echo "-> peaks half size set to: ${2}"; shift 2;;
			-t)
				local threads=$2
				echo "-> threads set to: ${2}"; shift 2;;
			-size)
				local sizeFile=$2
				echo "-> sizeFile set to: ${2}"; shift 2;;
			-h)
				usage peakcalling_GOPeaks; return;;
			--help)
				usage peakcalling_GOPeaks; return;;
			*)
				echo "[ERROR] - in argument: $1 $2";
				usage peakcalling_GOPeaks; return;;
		esac
	done
	# checking arguments for errors
	local Errors=0
	if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage peakcalling_MACS3; return 1; fi


	echo "[INFO] - Processing folder: ${in_dir##*/}"  # treating replicate folder
	local out_dir=${out_dir}/${in_dir##*/}

	# Environment à tester AJ
	bamPEFragmentSize=/home/312.6-Flo_Re/312.6.1-Commun/Conda_forge/Flore/miniforge3/envs_dirs/TFgenomics_DL/bin/bamPEFragmentSize

	mkdir -m 774 -p $out_dir
	local log=$out_dir/${in_dir##*/}.log;
	echo "" > $log
	# Saving Arguments used for this analysis
	echo $@ >> $log
	local controls=()
	local replicates=()
	# creating list of bam files for multiple controls and multiple replicates at once (usually only one replicate here as multiple replicate are dealt with mspc)
	echo "[INFO] - Bam files found:" | tee -a $log 
	for bam in $(find $in_dir -name "*.filtered.sorted.bam")
	do
		echo "$bam" | tee -a $log
		if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" == *"control"* ]]; then
			controls+=("$bam")
		else
			replicates+=("$bam")
		fi
	done
	
    mkdir -p -m 774 $out_dir/controls
    mkdir -p -m 774 $out_dir/replicates
    
    cd $out_dir/controls
    # checking if control files are present, not blocking the analysis if not, just changing macs2 launching parameters
    if [ ${#controls[@]} -eq 0 ]; then
        echo "[WARNING] - No control files present" | tee -a $log
    elif [ ${#controls[@]} -eq 1 ]; then
        echo "[INFO] - Only one bam in control, no merge needed" | tee -a $log
        cp ${controls[0]} control.bam
    else
        # pooling controls bam into a unique bam
        printf "%s\n" "${controls[@]}" > listbams.txt
        echo "[INFO] - Concatenating controls: ${controls[@]}" | tee -a $log
        if [ ! -f $out_dir/controls/control.bam ]; then
            $PATH_TO_SAMTOOLS/samtools merge --threads $(calc ${threads}-1) -b listbams.txt control.bam
        fi
    fi
    if [ ! -f $out_dir/controls/control.bam.bai ]; then
		echo "[INFO] - Indexing control bam" | tee -a $log
        $PATH_TO_SAMTOOLS/samtools index $out_dir/controls/control.bam 
    fi
    cd $out_dir/replicates
    # same check for replicates, but merging should not be usefull as we use MSPC to handle multiple replicates ( WARNING not tested with multiple replicates)
    if [ ${#replicates[@]} -eq 0 ]; then
        echo "[ERROR] - No replicates present" | tee -a $log
		return 1
    elif [ ${#replicates[@]} -eq 1 ]; then
        echo "[INFO] - Only one bam as replicate, no merge needed" | tee -a $log
        if [ ! -f $out_dir/replicates/replicate.bam ]; then
            cp ${replicates[0]} replicate.bam
        fi
    else
        echo "[INFO] - Concatenating replicates" | tee -a $log
        printf "%s\n" "${replicates[@]}" > listbams.txt
        $PATH_TO_SAMTOOLS/samtools merge --threads $(calc ${threads}-1) -b listbams.txt replicate.bam
    fi
    if [ ! -f $out_dir/replicates/replicate.bam.bai ]; then
        $PATH_TO_SAMTOOLS/samtools index $out_dir/replicates/replicate.bam
    fi
    R2_value="NA"
    # using minimal stats file of replicate to decide if analysis is paired ended or single ended as macs2 needs this info
    for bam in $(find $in_dir -name "*.filtered.sorted.bam")
    do
        if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
        	R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${bam%.filtered.sorted.bam}.minimal.stats)
        fi
    done
    cd $out_dir
	echo "[INFO] - Starting GOPeaks" | tee -a $log
    conda run -n order_66 gopeaks -p 1e-100 -m 20 -w $phs -b $out_dir/replicates/replicate.bam -c $out_dir/controls/control.bam -o $out_dir/${in_dir##*/} -s $sizeFile

    if [ ${R2_value} != "NA" ]; then
        # PE dataset
        $bamPEFragmentSize -b $out_dir/replicates/replicate.bam --table $out_dir/metrics.txt
        fragment_length=$(cat $out_dir/metrics.txt | awk 'NR==2{print $6}')
    else
        # SE dataset
        macs2 predictd -i $out_dir/replicates/replicate.bam -f BAM --out_dir $outdir >> $log 2>&1
        fragment_length=$(cat $log | grep "predicted fragment length is" | awk -v OFS="\t" '{print $14}')
    fi

    echo "[REPORT] - Fragment length: $fragment_length" | tee -a $log
    echo "[INFO] - Recomputing coverages along the genome" | tee -a $log
    # Using size of fragment and bam file to recompute coverages along the genome.
    # We do not use the bedgraph generated by macs2 as they are normalize by size of smallest bam
    for bam in $(find $in_dir -name "*.filtered.sorted.bam")
    do
        echo $bam
        if [[ ! -L $bam && -f $bam ]]; then
            if [ ${R2_value} == "NA" ]; then # checking for PE, or single ended reads
                if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
                    bedtools bamtobed -i $bam > $in_dir/${in_dir##*/}.bamtobed.bed
                    awk -v fraglen=${fragment_length} -v OFS="\t" '$6=="+"{print $1,$2,$3+fraglen,$4,$5,$6}$6=="-"{print $1,$2-fraglen,$3,$4,$5,$6}' $in_dir/${in_dir##*/}.bamtobed.bed | awk -v OFS="\t" '$2<=0{print $1,1,$3,$4,$5,$6;next}{print $0}'  > $in_dir/${in_dir##*/}.ext.bed
                    local libsize=$(samtools view -f 0 -c $bam) # counting every mapped read (bitwise flag = 0)
                    local scale=$(calc 1000000/$libsize) # scaling factor for single-end data to get CPM normalization
                    
                    tmpsave=$LC_COLLATE;LC_COLLATE=C
                    bedtools genomecov -bga -scale $scale -i $in_dir/${in_dir##*/}.ext.bed -g $sizeFile | sed 's/Chr/chr/g' | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_cov.bdg
                    LC_COLLATE=$tmpsave
                    $bdg2bwig $out_dir/${in_dir##*/}_cov.bdg $sizeFile $out_dir/${in_dir##*/}_cpm.bw
                fi
            else
                if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
                    local libsize=$(samtools view -f 0 -c $bam) 
                    local scale=$(calc 1000000/$libsize) # scaling factor for paired-end data to get CPM normalization
                    tmpsave=$LC_COLLATE;LC_COLLATE=C
                    bedtools genomecov -bga -scale $scale -ibam $bam | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_cov.bdg
                    LC_COLLATE=$tmpsave
                    $bdg2bwig $out_dir/${in_dir##*/}_cov.bdg $sizeFile $out_dir/${in_dir##*/}_cpm.bw
                fi
            fi
        fi
    done
    
	bedtools genomecov -bg -ibam $out_dir/replicates/replicate.bam -g $sizeFile | sortBed -i > $out_dir/${in_dir##*/}_covtmp.bed
    bedtools intersect -b $out_dir/${in_dir##*/}_covtmp.bed -a $out_dir/${in_dir##*/}_peaks.bed -wo > $out_dir/peaks_tmp.bed      
    awk -v OFS="\t" -v max=0 -v ch=0 -v pss=0 -v pts=0 -v id_old=0 '(id=$1$2$3) {
        if(id!=id_old){
            if(id_old==0){
                pos=int(($5+$6)/2)
                max=$7
                ch=$1
                pss=$2
                pts=$3
            } else {
                print ch,pss,pts,max,pos
                pos=int(($5+$6)/2)
                max=$7
                ch=$1
                pss=$2
                pts=$3
            }
        }else{
            if($7>max){
                pos=int(($5+$6)/2)
                max=$7
                ch=$1
                pss=$2
                pts=$3
            }
        } id_old=id
    } END{print ch,pss,pts,max,pos}' peaks_tmp.bed  > $out_dir/${in_dir##*/}_max_peaks.bed
    # cat $out_dir/${in_dir##*/}_max_peaks.bed >> all_pics_tmp.bed
	# removing temporary folders
    rm $out_dir/peaks_tmp.bed
	rm -R $out_dir/controls;
    rm -R $out_dir/replicates;
}

#-------------------------------------------------------------------------------
main_peakcalling(){
	# FUNCTION: main_peakcalling
	# DESCRIPTION:
	#   This function performs peak calling for genomic data analysis. It supports
	#   multiple peak callers (e.g., MACS3, GOPeaks) and handles various input 
	#   parameters to customize the analysis. The function processes input data, 
	#   applies filters, and generates consensus peaks and coverage files. 
	#   Consensus peaks from replicate are dtermined with MSPC tool and resized.
	# 
	#      Example for peak resizing:
	#      </>: extremities for the peaks
	#      | : maximum of peaks, determined by MACS3
    # 
	# 	       replicate 1:         <-------------|------------------->
	#        replicate 2  <---------|----------------------->
	#      replicate 3:                   <------------------------|--------->
    # 
	#      consensus MSPC: <------------------------------------------------->
	#      peaks resized:  <-----------|-----------> <-----------|----------->
	# 
	#
	# USAGE:
	#   main_peakcalling -id <PATH> -cd <PATH> -od <PATH> -nc <STR> -g <INT> -top <INT> 
	#                    -ps <INT> -size <PATH> -bl <PATH> -mspc <INT%> -pm <STR> 
	#                    -s <INT> -t <INT> --keep-dup <INT> -r --diag --llocal <INT> 
	#                    --factork <FLOAT> --qval <FLOAT> --minfold <INT> --maxfold <INT> 
	#                    --minlen <INT> --weak <FLOAT> --strong <FLOAT>
	#
	# PARAMETERS:
	#   -id       : List of path(s) to the input data directory (required, mapping directory).
	#   -cd       : List of path to the control data directory (optional).
	#   -od       : Path to the output directory (required).
	#   -nc       : Name for the consensus directory (required).
	#   -g        : Genome length (mappable) (default: 120000000).
	#   -top      : Maximum number of peaks to retain (default: 0, no limit).
	#   -ps       : Peak size (default: 200).
	#   -size     : Path to the genome size file.
	#   -bl       : Path to the blacklist file.
	#   -mspc     : Percentage for MSPC (default: 100%).
	#   -pm       : Peak caller to use (e.g., MACS3, GOPeaks) (default: MACS3).
	#   -s        : Random seed for reproducibility (default: 168159).
	#   -t        : Number of threads to use (default: 8).
	#   --keep-dup: Keep duplicate reads (default: 1).
	#   -r        : Redo analysis (default: false).
	#   --diag    : Enable diagnostic mode (default: false).
	#   --llocal  : Large local background size to discrimine from 
	#               noise signal(default: 5000).
	#   --factork : Factor k for custom control filtering, coverage 
	#               difference between control and peaks (default: 3).
	#   --qval    : q-value threshold for peak detection (default: 0.05).
	#   --minfold : Minimum fold change for peak detection (default: 5).
	#   --maxfold : Maximum fold change for peak detection (default: 50).
	#   --minlen  : Minimum length for peak detection (default: 100).
	#   --weak    : Weak p-value threshold for MSPC (default: 1e-5).
	#   --strong  : Strong p-value threshold for MSPC (default: 1e-9).
	#
	# DEPENDENCIES:
	#   - bedtools
	#   - samtools
	#   - MACS3
	#   - GOPeaks
	#
	# NOTES:
	#   - The output directory will be overwritten if it already exists.

	echo "	

	_  _ ____ _ _  _     ___  ____ ____ _  _ ____ ____ _    _    _ _  _ ____
	|\/| |__| | |\ |     |__] |___ |__| |_/  |    |__| |    |    | |\ | | __
	|  | |  | | | \| ___ |    |___ |  | | \_ |___ |  | |___ |___ | | \| |__]	

		 "

	local top=0; local mspcpc=100; local phs=200; local Genome_length=120000000; local seedrandom=168159; local keepDup=1;  local threads=8; local redo_analysis="false"; local input_dir=("NA"); local sizeFile="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize"; local blacklist="/home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Greenscreen_19012023_merged.bed"; local peakcaller="MACS3"; local llocal=5000; local factork=3; local qval=0.05; local minlen=100; local diag="false"; local weak=1e-5; local strong=1e-9; local minfold=5; local maxfold=50; local SDfilter=0;
	local Argsline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-id)
				local in_dir=("${!2}")
				echo "-> Data directory set to: ${in_dir[@]}";shift 2;;
			-cd)
				local input_dir=("${!2}")
				echo "-> Control directory set to: ${input_dir[@]}";shift 2;;
			-od)
				local out_dir=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-nc)
				local name_cons=$2
				echo "-> Name for consensus directory set to: ${2}";shift 2;;
			-g)
				local Genome_length=${2%.*}
				echo "-> Genome length (mappable) set to: ${2}";shift 2;;
			-top)
				local top=$2
				echo "-> Maximum number of peaks set to: ${2}";shift 2;;
			-ps)
				local phs=$(calc $2/2.0) # $2 is the peak size! but we take only half, hence peak half size...
				echo "-> Peak size: ${2}";shift 2;;
			-size)
				local sizeFile=$2; 
				echo "-> Size file for genome is: $2"; shift 2;;
			-bl)
				local blacklist=$2 ; 
				echo "-> Blacklist set to: $2"; shift 2 ;;
			-mspc)
				local mspcpc=$2
				echo "-> Percentage for mspc set to: $2"; shift 2;;
			-pm)
				local peakcaller=${2^^}
				echo "-> Selected peakcaller: ${2^^}"; shift 2;;
			-s)
				local seedrandom=$2
				echo "-> Random seed set to: ${2}"; shift 2;;
			-t)
				local threads=$2
				echo "-> Threads set to: ${2}"; shift 2;;
			--keep-dup)
				local keepDup=$2
				echo "-> Keep-dup set to: ${2}"; shift 2;;
			--SD-filtering)
				local SDfilter=1
				echo "-> filtering peaks with SD of coverages"; shift 1;;
			-r)
				local redo_analysis="true"
				echo "-> Redo analysis activated"; shift 1;;
			-diag)
				local diag="true"
				echo "-> Diagnostic mode activated"; shift 1;;
			--llocal)
				local llocal=${2%.*}
				echo "-> Large local background set to: ${2}"; shift 2;;
			--factork) 
				local factork=$2
				echo "-> Factor k for custom control filtering set to: ${2}"; shift 2;;
			--qval) 
				local qval=$2
				echo "-> q-value for peak detection set to: ${2}"; shift 2;;
			--minfold)
				local minfold=$2
				echo "-> Minimum fold change for peak detection set to: ${2}"; shift 2;;
			--maxfold)
				local maxfold=$2
				echo "-> Maximum fold change for peak detection set to: ${2}"; shift 2;;
			--minlen) 
				local minlen=${2%.*}
				echo "-> Minimum length for peak detection set to: ${2}"; shift 2;;
			--strong) 
				local strong=$2
				echo "-> Strong p-value for MSPC set to: ${2}"; shift 2;;
			--weak) 
				local weak=$2
				echo "-> Weak p-value for MSPC set to: ${2}"; shift 2;;
			-h)
				usage main_peakcalling; return;;
			--help)
				usage main_peakcalling; return;;
			*)
				echo "Error in arguments"
				echo $1"\t"$2; usage main_peakcalling; return;;
		esac
	done

	# Checking arguments for errors
	local Errors=0
	if [ -z $in_dir ]; then echo "ERROR: -id argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ -z $name_cons ]; then echo "ERROR: -nc argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage main_peakcalling; return 1; fi
	local list_peaks_mspc=()
	local list_bdg=()
	local list_rm=()

	if [[ $mspcpc != *"%"* ]]; then local mspcpc=$mspcpc"%" ; fi
	if [[ -d $out_dir/$name_cons ]]; then rm -Rf $out_dir/$name_cons; fi
	
	# logs
	mkdir -p -m 774 $out_dir/$name_cons/temp
	# Saving Arguments used for this analysis
	echo "main_peakcalling" $Argsline > $out_dir/$name_cons/log.txt
	tmpsave=$LC_COLLATE;LC_COLLATE=C

	# ------- PARSE INPUT DIRECTORIES
	for dir in ${in_dir[@]};
	do
		# link bamfile of the controls into the directory of the sample
		bam_ech=$(find $dir -name "*.filtered.sorted.bam")
		local i=1
		controlstats="NA"
		if [ ${input_dir[0]} != "NA" ]; then
			for ctrldir in ${input_dir[@]};
			do
				for bam in $(find $ctrldir -name "*.filtered.sorted.bam")
				do
					controlstats=${bam%.filtered.sorted*}.minimal.stats
					if [ ! -f $dir/control_$i.filtered.sorted.bam ] || [[ $bam -nt $dir/control_$i.filtered.sorted.bam ]]; then
						ln -s $bam $dir/control_$i.filtered.sorted.bam
					fi
					local i=$(($i+1))
				done
			done
		fi

		local name=${dir##*/}
		
		# -------GOPPeaks
		if [ $peakcaller == "GOPEAKS" ]; then
			echo "[INFO] - Using GOPeaks for peakcalling"
			echo "[INFO] - Launching GOPeaks peakcalling..." | tee -a $out_dir/$name_cons/log.txt
			# expected outputs of function peakcalling_GOPeaks per replicate
			list_peaks_mspc+=("$out_dir/$name/${name}_max_peaks.bed") 
			list_bdg+=("$out_dir/$name/${name}_cov.bdg")
			peakcalling_GOPeaks -id $dir -od $out_dir -phs $phs -t $threads -size $sizeFile -s $seedrandom
			cat $out_dir/$name/${name}_max_peaks.bed | tee -a $out_dir/$name_cons/temp/all_pics_tmp.bed

		# ------- MACS3
		elif [ $peakcaller == "MACS3" ]; then
			echo "[INFO] - Launching MACS3 peakcalling..." | tee -a $out_dir/$name_cons/log.txt
			# expected outputs of function peakcalling_MACS2 per replicate
			list_peaks_mspc+=("$out_dir/$name/${name}_peaks.bed") 
			list_bdg+=("$out_dir/$name/${name}_cov.bdg")
			if [ -f $controlstats ] || [ $controlstats != "NA" ]; then
				peakcalling_MACS3 -id $dir -od $out_dir -g $Genome_length -s $seedrandom -phs $phs -kd $keepDup -t $threads -size $sizeFile -r $redo_analysis -bl $blacklist -cs $controlstats --llocal $llocal --factork $factork --qval $qval --minlen $minlen -diag $diag --minfold $minfold --maxfold $maxfold
			else
				peakcalling_MACS3 -id $dir -od $out_dir -g $Genome_length -s $seedrandom -phs $phs -kd $keepDup -t $threads -size $sizeFile -r $redo_analysis -bl $blacklist --llocal $llocal --factork $factork --qval $qval --minlen $minlen -diag $diag --minfold $minfold --maxfold $maxfold
			fi

			echo "[INFO] - Computing FRiP and other stats" | tee -a $out_dir/$name_cons/log.txt
			for elt in ${bam_ech[@]};
			do
				if [ ! -f $out_dir/$name/${name}_stats.txt ] || [[ $out_dir/$name/${name}_peaks.bed -nt $out_dir/$name/${name}_stats.txt ]]; then # checking if this step is already done to not do it again
					if [[ $elt != *"control"* ]]; then # filter out all bam that contain control in his name.
						# computing Fresquency of Reads in Peaks (FRiP)
						total=$(samtools view -c $elt)
						inpeak=$(bedtools sort -i $out_dir/$name/${name}_peaks.bed | bedtools merge -i stdin | bedtools intersect -u -a $elt -b stdin -ubam | samtools view -c)
						FRIP=$(calc $inpeak/$total*100)
						echo "$name $FRIP% $inpeak / $total" | tee -a $out_dir/$name_cons/log.txt
						echo "$name $FRIP% $inpeak / $total" > $out_dir/$name/${name}_FreqReadInPeak.txt
						R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${elt%.filtered.sorted.bam}.minimal.stats)
						if [ $R2_value == "NA" ]; then # number of tags/fragments, way of retrieving this info is PE/SE dependant
							totTags=$(grep "total tags in treatment" $out_dir/$name/${name}.log | awk '{print $NF}')
							if [ $keepDup == "all" ]; then
								filtTags=$totTags
							else
								filtTags=$(grep "tags after filtering in treatment" $out_dir/$name/${name}.log | awk '{print $NF}')
							fi
						fi
						if [ $R2_value != "NA" ]; then
							totTags=$(grep "total fragments in treatment" $out_dir/$name/${name}.log | awk 'NR==1{print $NF}' )
							if [ $keepDup == "all" ]; then
								filtTags=$totTags
							else
								filtTags=$(grep "fragments after filtering in treatment" $out_dir/$name/${name}.log | awk 'NR==1{print $NF}' )
							fi
						fi
						filtPeaks=$(wc -l $out_dir/$name/${name}_filtered.narrowPeak | cut -d " " -f 1)
						echo "Sample totalTags filtTags filtPeaks FRIP" > $out_dir/$name/${name}_stats.txt
						echo $name $totTags $filtTags $filtPeaks $FRIP >> $out_dir/$name/${name}_stats.txt
					fi
				fi
			done
		else
			echo "[ERROR] - undefined peakcaller selected" | tee -a $out_dir/$name_cons/log.txt
			echo "->Correct values are: MACS2 GOPEAKS MACS3" | tee -a $out_dir/$name_cons/log.txt
		fi
	done

	# ------- MSPC TO DETERMINE CONSENSUS PEAKS AMONG REPLICATES
	if [ $peakcaller == "MACS2" ] || [ $peakcaller == "MACS3" ]; then
		if [ "${#in_dir[@]}" -gt 1 ]; then 
			# merging all coverages into one file, to be used for MSPC
			
			# p_val=$(awk 'BEGIN {print 10**-30}') #TODO implement this for MSPC (pass lower values than 1e-09 to MACS3) ## JL: maybe keep this as an idea for future improvement. 
			${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w $weak -s $strong -c ${mspcpc} -o $out_dir/$name_cons/MSPC -d 5 | tee -a $out_dir/$name_cons/log.txt
		
			# filtering max number of peaks or not, depending on user choice. Filter is done on score
			if [ $top -eq 0 ]; then  
				sed '1d' $out_dir/$name_cons/MSPC/ConsensusPeaks.bed | awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}.bed
			else
				sed '1d' $out_dir/$name_cons/MSPC/ConsensusPeaks.bed | sort -k5,5nr | awk -v OFS="\t" '{print $1,$2,$3}' > $out_dir/$name_cons/${name_cons}_comp.bed
				head -${top} $out_dir/$name_cons/${name_cons}_comp.bed | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}.bed
			fi
			local i=0
			local files=()
			local reppeaks=()

			# prep files for next big filter, here we compute positions for the maximum of each peaks, for each replicate
			for dir in ${in_dir[@]};
			do
				local name=${dir##*/}
				if [ $i -eq 0 ]; then
					# $5 start=filtered narropeak & $13 position of maximum 
					bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13,$5,$6,"1"}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4" "$5" "$6" "$7}END{print save}' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed
				else
					bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v rep=$i -v OFS="\t" '{print $1,$2,$3,$5+$13,$5,$6,rep+1}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4" "$5" "$6" "$7}END{print save}' | awk -v OFS=" " '{$1=$2=$3="";print $0}' | sed 's/   //' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed
				fi
				reppeaks+=("$out_dir/$name/${name}_filtered.narrowPeak")
				files+=("$out_dir/$name_cons/temp/tmp_peaks_$i.bed")
				local i=$(($i+1))
			done

			cat ${reppeaks[@]} | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,"customMade","boundRegion",$2,$3,$7,".",".","ID="$4; print $1,"customMade","Peak",$2,$3,$7,".",".","ID="$4"."NR";Parent="$4; print $1,"customMade","Max",$2+$10,$2+$10+1,$7,".",".","ID="$4"."NR"."NR";Parent="$4"."NR;}' > $out_dir/$name_cons/${name_cons}_replicates_peaks.gff
			rm $out_dir/$name_cons/$name_cons.bed # cleaning unused misleading file

			# ------- PEAK RESIZING BASED ON MSPC CONSENSUS AND REPLICATE MAXIMUM POSITIONS
			# The maximum of each peaks is compared among each replicates. if maximum are close enough (peak 2/5 of peak size or 300 bp (lowest value)) they are considered as same maximum, otherwize multiple peaks are created to insure that 1 maximum => 1 peaks rule is respected.
			echo ${#files[@]} ${mspcpc//%} | tee -a $out_dir/$name_cons/log.txt

			paste ${files[@]} | awk -v OFS="\t" '{
				nb=0
				for(i=4;i<=NF;i+=4){
					if($i in subpeaks){
						subpeaks[$i]=subpeaks[$i]"\t"$(i)"\t"$(i+1)"\t"$(i+2)"\t"$(i+3)
					}else{
						subpeaks[$i]=$(i)"\t"$(i+1)"\t"$(i+2)"\t"$(i+3)
						ordersub[nb++]=$i
					}
				}
				n=asort(ordersub,sortedindex)
				save=""
				for (i = 1; i <= n; i++) {
					if(save!=""){
						save=save"\t"subpeaks[sortedindex[i]]
					}else{
						save=subpeaks[sortedindex[i]]
					}
				}
				print $1,$2,$3,save; delete subpeaks; delete ordersub
			}' > $out_dir/$name_cons/temp/${name_cons}_narrow.bed
			
			thpc=$(calc ${#files[@]}*${mspcpc//%}/100 )
			echo "${thpc%.*},$phs" | tee -a $out_dir/$name_cons/log.txt
			
			# awk of the death, blame JL 
			awk -v SIZE=$phs -v th=${thpc%.*} -v FS="[ \t]" -v OFS="\t" 'function abs(v) {return v < 0 ? -v : v}
			{
				for(i=4;i<=NF;i+=4){ # creating an array of array with subpeaks from all replicates
					if($i!=-1){
						arraylen=0
						if($(i+3) in listsubpeaks){
							arraylen=length(listsubpeaks[$(i+3)]["max"])
						}
						arraylen++
						listsubpeaks[$(i+3)]["max"][arraylen]=$i
						listsubpeaks[$(i+3)]["start"][arraylen]=$(i+1)
						listsubpeaks[$(i+3)]["stop"][arraylen]=$(i+2)
					}
				}
				nbmax=0; 
				if( (($3-$2)*2/5)>300){
					compareValue=300
				}else{
					compareValue=($3-$2)*2/5
				}
				for (replicate in listsubpeaks){
					# print "==",replicate,"==" #DEBUG
					if (nbmax==0){ # no local maximum defined yet
						for (subpeakmax in listsubpeaks[replicate]["max"]){ 
							# for each subpeaks in this replicate; generate a local maximum
							listlocalmaxi[nbmax+1]["max"]=listsubpeaks[replicate]["max"][subpeakmax];
							listlocalmaxi[nbmax+1]["moy"]=listsubpeaks[replicate]["max"][subpeakmax];
							listlocalmaxi[nbmax+1]["nb"]=1;
							listlocalmaxi[nbmax+1]["start"]=listsubpeaks[replicate]["start"][subpeakmax];
							listlocalmaxi[nbmax+1]["stop"]=listsubpeaks[replicate]["stop"][subpeakmax];
							listlocalmaxi[nbmax+1]["rep"]=replicate;
							nbmax=nbmax+1
						}
					}else{ # at leats one local maximum defined
						candidate=0;
						for (subpeakmax in listsubpeaks[replicate]["max"]){
							# maximum taken from subpeaks
							ValueSubpeakmax=listsubpeaks[replicate]["max"][subpeakmax] 
							for (indexlocalmax in listlocalmaxi){
								repUsed=0
								# maximum taken from registered local maximum
								Valuelocalmax=listlocalmaxi[indexlocalmax]["max"]
								# print indexlocalmax,"Valuelocalmax",Valuelocalmax,"ValueSubpeakmax",ValueSubpeakmax,"-" #DEBUG
								if ( abs(Valuelocalmax - ValueSubpeakmax) <= compareValue ){
									# maximum from subpeak is close enough from registered maximum
									# print " "," ","close enough" #DEBUG
									for(oldsubpeakmax in listOfSubpeakmax){
										# print " "," "," ",listOfSubpeakmax[oldsubpeakmax],indexlocalmax #DEBUG
										if(listOfSubpeakmax[oldsubpeakmax]==indexlocalmax){
											repUsed=1; 
											# print " "," ","index used",candidate #DEBUG

											if (candidate!=0){
												#test if new subpeak is  closer than old one
												oldcandidateValuelocalmax=listlocalmaxi[candidate]["max"]
												oldValueSubpeakmax=listsubpeaks[replicate]["max"][oldsubpeakmax]
												# print abs(Valuelocalmax - ValueSubpeakmax), abs(Valuelocalmax - oldValueSubpeakmax)
												if( abs(Valuelocalmax - ValueSubpeakmax) < abs(Valuelocalmax - oldValueSubpeakmax) ){
													#new subpeak is better
													# print " "," "," ","passed"
													if(abs(Valuelocalmax - ValueSubpeakmax) < abs(oldcandidateValuelocalmax - ValueSubpeakmax) ){
														# new candidate is better
														# print " "," ","taking over" #DEBUG
														candidate=indexlocalmax
														delete listOfSubpeakmax[oldsubpeakmax]
													}
												}
											}else{
												oldValueSubpeakmax=listsubpeaks[replicate]["max"][oldsubpeakmax]
												if( abs(Valuelocalmax - ValueSubpeakmax) < abs(Valuelocalmax - oldValueSubpeakmax) ){
													candidate=indexlocalmax
													delete listOfSubpeakmax[oldsubpeakmax]
												}
											}
										}
									}
									if (repUsed==0){
										if (candidate!=0){
											# subpeak is already candidate to participate in a localmax
											oldcandidateValuelocalmax=listlocalmaxi[candidate]["max"]
											if ( abs(Valuelocalmax - ValueSubpeakmax) < abs(oldcandidateValuelocalmax - ValueSubpeakmax) ){
												# print " "," ","replacing max" #DEBUG
												# checking if new localmax is better than the older one
												candidate=indexlocalmax
											}
										}else{
											# print " "," ","adding to this max" #DEBUG
											# subpeak is candidate for a local max for the first time
											candidate=indexlocalmax
										}
									}
								}
							}
							# for a replicate, saving all relation subpeak - candidate
							listOfSubpeakmax[subpeakmax]=candidate; candidate=0;
						}
						for (subpeakmax in listOfSubpeakmax){
							candidate=listOfSubpeakmax[subpeakmax]
							if (candidate==0){
								listlocalmaxi[nbmax+1]["max"]=listsubpeaks[replicate]["max"][subpeakmax];
								listlocalmaxi[nbmax+1]["moy"]=listsubpeaks[replicate]["max"][subpeakmax];
								listlocalmaxi[nbmax+1]["nb"]=1
								listlocalmaxi[nbmax+1]["start"]=listsubpeaks[replicate]["start"][subpeakmax];
								listlocalmaxi[nbmax+1]["stop"]=listsubpeaks[replicate]["stop"][subpeakmax];
								listlocalmaxi[nbmax+1]["rep"]=replicate
								nbmax=nbmax+1
							}else{
								listlocalmaxi[candidate]["nb"]+=1
								listlocalmaxi[candidate]["moy"]+=listsubpeaks[replicate]["max"][subpeakmax];
								listlocalmaxi[candidate]["rep"]=listlocalmaxi[candidate]["rep"]","replicate;
								# if (listlocalmaxi[candidate]["start"] > listsubpeaks[replicate]["start"][subpeakmax]){
								# 	listlocalmaxi[candidate]["start"]=listsubpeaks[replicate]["start"][subpeakmax]
								# }
								# if (listlocalmaxi[candidate]["stop"] < listsubpeaks[replicate]["stop"][subpeakmax]){
								# 	listlocalmaxi[candidate]["stop"]=listsubpeaks[replicate]["stop"][subpeakmax]
								# }
								listlocalmaxi[candidate]["start"]+=listsubpeaks[replicate]["start"][subpeakmax]
								listlocalmaxi[candidate]["stop"]+=listsubpeaks[replicate]["stop"][subpeakmax]
							}
						} 
						delete listOfSubpeakmax;
					}
				}
				for (indexlocalmax in listlocalmaxi){
					if(listlocalmaxi[indexlocalmax]["nb"] >= th){
						moy=int(listlocalmaxi[indexlocalmax]["moy"]/listlocalmaxi[indexlocalmax]["nb"])
						print $1,moy-SIZE,moy+SIZE, int(listlocalmaxi[indexlocalmax]["start"]/listlocalmaxi[indexlocalmax]["nb"]), int(listlocalmaxi[indexlocalmax]["stop"]/listlocalmaxi[indexlocalmax]["nb"]), listlocalmaxi[indexlocalmax]["rep"]
						listlocalmaxi[indexlocalmax]["max"]=0
					}
				};
				# saving localmaxi close to each other with too few replicates to support them
				for (indexlocalmax1 in listlocalmaxi){
					for(indexlocalmax2 in listlocalmaxi){
						if((indexlocalmax2!=indexlocalmax1)&&(listlocalmaxi[indexlocalmax1]["max"]!=0)&&(listlocalmaxi[indexlocalmax2]["max"]!=0)){
							moy1=int(listlocalmaxi[indexlocalmax1]["moy"]/listlocalmaxi[indexlocalmax1]["nb"])
							moy2=int(listlocalmaxi[indexlocalmax2]["moy"]/listlocalmaxi[indexlocalmax2]["nb"])
							nbtot=listlocalmaxi[indexlocalmax1]["nb"]+listlocalmaxi[indexlocalmax2]["nb"]
							allreps=listlocalmaxi[indexlocalmax1]["rep"]","listlocalmaxi[indexlocalmax]["rep"]
							split(allreps,listallreps,",");
							ind=0;
							for(rep in listallreps){
								if(ind==0){
									listconfirmedRep[listallreps[rep]]=listallreps[rep]
									# confirmedRep=listallreps[rep]
									ind+=1
								}else{
									if(listallreps[rep] in listconfirmedRep);else{
										listconfirmedRep[listallreps[rep]]=listallreps[rep]
										# confirmedRep=confirmedRep","listallreps[rep]
										ind+=1
									}
								}
							}
							asort(listconfirmedRep);ind2=0
							for(elt in listconfirmedRep){
								if(ind2==0){
									confirmedRep=listconfirmedRep[elt]
									ind2+=1
								}else{
									confirmedRep=confirmedRep","listconfirmedRep[elt]
								}
							}
							delete listconfirmedRep; delete listallreps
							if((abs(moy1-moy2)<SIZE)&&(ind>=th)){
								moy=int((listlocalmaxi[indexlocalmax1]["moy"] + listlocalmaxi[indexlocalmax2]["moy"]) / (listlocalmaxi[indexlocalmax1]["nb"] + listlocalmaxi[indexlocalmax2]["nb"]))
								start=int((listlocalmaxi[indexlocalmax1]["start"] + listlocalmaxi[indexlocalmax2]["start"]) / (listlocalmaxi[indexlocalmax1]["nb"] + listlocalmaxi[indexlocalmax2]["nb"]))
								stop=int((listlocalmaxi[indexlocalmax1]["stop"] + listlocalmaxi[indexlocalmax2]["stop"]) / (listlocalmaxi[indexlocalmax1]["nb"] + listlocalmaxi[indexlocalmax2]["nb"]))
								print $1,moy-SIZE, moy+SIZE, start, stop,confirmedRep
								listlocalmaxi[indexlocalmax1]["max"]=0
								listlocalmaxi[indexlocalmax2]["max"]=0
							}
						}
					}
				}
				delete listlocalmaxi; delete listsubpeaks;
			}' $out_dir/$name_cons/temp/${name_cons}_narrow.bed | sed 's/\t\+/\t/g;s/^\t//' | awk -v OFS="\t" '{print $0}' > $out_dir/$name_cons/${name_cons}_wextendedpos.bed
			awk -v OFS="\t" '{print $1,$2,$3}' $out_dir/$name_cons/${name_cons}_wextendedpos.bed > $out_dir/$name_cons/${name_cons}_narrow.bed
			wc -l $out_dir/$name_cons/${name_cons}_narrow.bed | tee -a $out_dir/$name_cons/log.txt

			# finally the mean coverage values of replicates is used as consensus coverage values
			echo "[INFO] - Computing mean coverage for consensus coverage. This might take a while...." | tee -a $out_dir/$name_cons/log.txt
			$PATH_TO_MACS3/macs3 cmbreps -i ${list_bdg[@]} -m mean -o $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg
		else
			# ------- NO REPLICATES, SIMPLY COPYING FILES
			# no step of replicates "merging" is necessary when you have only one replicate
			local name=${in_dir[0]##*/}
			cp $out_dir/$name/${name}_narrow.bed $out_dir/$name_cons/${name_cons}_narrow.bed
			cp ${list_bdg[0]} $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg

		fi

		#  ------- FINALIZING PEAKS COVERAGE AND PEAKS FILES
		awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$4}' $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}_cov.bdg
		rm $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg
		$bdg2bwig $out_dir/$name_cons/${name_cons}_cov.bdg $sizeFile $out_dir/$name_cons/${name_cons}_cov.bw
		local i=0
		local files=()
		for dir in ${in_dir[@]}; # retrieving coverage values at maximum position to use as filter if necessary.
		do
			local name=${dir##*/}
			if [ $i -eq 0 ]; then
				bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cov.bdg -wao | awk -v OFS="\t" '{print $1":"$2"|"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 > $out_dir/$name_cons/temp/tmp_max_${i}.txt
			else
				bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cov.bdg -wao | awk -v OFS="\t" '{print $1":"$2"|"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 | awk '{print $2}' > $out_dir/$name_cons/temp/tmp_max_${i}.txt
			fi
			files+=("$out_dir/$name_cons/temp/tmp_max_${i}.txt")
			local i=$(($i+1))
		done
		paste ${files[@]} | sed "s/[|:]/\t/g" > $out_dir/$name_cons/${name_cons}_max.bed
		awk -v OFS="\t" '{mean=0;for(i=4;i<=NF;i++) {mean+=$i};print $1,$2,$3,mean/(i-3)}' $out_dir/$name_cons/${name_cons}_max.bed | sort -k4,4nr > $out_dir/$name_cons/${name_cons}_maxMean.bed
		if [ $SDfilter -eq 1 ]; then
			echo "Applying standard deviation filter on maximum coverage values" | tee -a $out_dir/$name_cons/log.txt
			# filtering peaks with a maximum value too different from the mean of maximum values (more than 2*SD away from mean)
			# awk '{print $0, sqrt((1/3)*((($4^2)+($5^2)+($6^2)) - ((($4+$5+$6)^2)/3))),($4+$5+$6)/3}' $out_dir/$name_cons/${name_cons}_max.bed | awk '{print $0,$8-($7*2),$8+($7*2)}' | awk '{keep=0;if($4<$9||$4>$10){keep=keep+1};if($5<$9||$5>$10){keep=keep+1};if($6<$9||$6>$10){keep=keep+1};if(keep<2){print $0}}' > $out_dir/$name_cons/${name_cons}_maxMean_SDfiltered.bed
			# modified on 2025-09-11 to be functional with any number of replicates, not only 3. #TOTEST 
			awk -v n=${#in_dir[@]} '{sum=0; sumsq=0; for(i=4;i<=NF;i++) {sum+=$i; sumsq+=($i)^2}; mean=sum/(i-4); sd=sqrt((sumsq/(i-4)) - (mean^2)); print $0, sd, mean}' $out_dir/$name_cons/${name_cons}_max.bed | awk '{print $0,$8-($7*2),$8+($7*2)}' | awk '{keep=0; for(i=4;i<=NF-3;i++){ if($i<$9||$i>$10){keep=keep+1} }; if(keep<('"${#in_dir[@]}"'/2)){print $0}}' > $out_dir/$name_cons/${name_cons}_maxMean_SDfiltered.bed
		else
			cp $out_dir/$name_cons/${name_cons}_maxMean.bed $out_dir/$name_cons/${name_cons}_maxMean_SDfiltered.bed
		fi
		rm -Rf $out_dir/$name_cons/*/

	# ------- MSPC FOR GOPEAKS 
	elif [ $peakcaller == "GOPEAKS" ]; then
		cat $out_dir/$name_cons/temp/all_pics_tmp.bed | sortBed -i > $out_dir/$name_cons/temp/all_peaks_sorted_tmp.bed

		# ------- Cluster peaks to get consensus positions
		bedtools cluster -i $out_dir/$name_cons/temp/all_peaks_sorted_tmp.bed | awk -v OFS='\t' -v id=0 -v save="" '{
			if(id==$6){
				save=save"\t"$2"\t"$3"\t"$5
			}else{
				if(save==""){
					id=$6;
					save=$1"\t"$2"\t"$3"\t"$5;
				}else{
					print save;
					save=$1"\t"$2"\t"$3"\t"$5;
					id=$6
				}
			}
		}END{print save}' > $out_dir/$name_cons/temp/ligned_peaks.bed

		#  ------- Get all replicats peaks 
		echo "[INFO] - Get all replicats peaks" | tee -a $out_dir/$name_cons/log.txt
		awk -v OFS='\t' '{if(NF==10)print $0}' $out_dir/$name_cons/temp/ligned_peaks.bed > $out_dir/$name_cons/temp/ligned_peaks_3R.bed

		while read -r line; do 
			echo "$line" | awk '{chr=$1; for (i=2; i<=NF; i+=3) printf "%s\t%s\t%s\t%s\n", chr, $(i), $(i+1), $(i+2)}' > $out_dir/$name_cons/temp/line.bed
			bedtools merge -i $out_dir/$name_cons/temp/line.bed -c 4 -o mean  | awk -v OFS='\t' '{print $1,$2,$3,$4}' >> $out_dir/$name_cons/temp/consensus_pics_3MA.bed
		done < $out_dir/$name_cons/temp/ligned_peaks_3R.bed  

		#  ------- Get complex peaks (potential subpeaks)
		echo "[INFO] - Get complexe peaks"	 | tee -a $out_dir/$name_cons/log.txt
		awk -v OFS='\t' '{if(NF>10) print $0 }' $out_dir/$name_cons/temp/ligned_peaks.bed > $out_dir/$name_cons/temp/ligned_peaks_4AMR.bed

		# ------- Process complexe peaks to get consensus
		while read -r line; do 
			echo "$line" | awk '{chr=$1; for (i=2; i<=NF; i+=3) printf "%s\t%s\t%s\t%s\n", chr, $(i), $(i+1), $(i+2)}' > $out_dir/$name_cons/temp/line.bed
			awk -v OFS="\t" '{print $4}' $out_dir/$name_cons/temp/line.bed | awk -v OFS="\t" '{printf "%s\t", $1} END {print ""}' > $out_dir/$name_cons/temp/tmp_max.bed
			bedtools merge -i $out_dir/$name_cons/temp/line.bed  | awk -v OFS='\t' '{print $1,$2,$3}' > $out_dir/$name_cons/temp/tmp_coordinate.bed
			paste $out_dir/$name_cons/temp/tmp_coordinate.bed $out_dir/$name_cons/temp/tmp_max.bed > $out_dir/$name_cons/temp/peak.bed
			awk -v SIZE=$phs -v th=3 -v FS="[ \t]" -v OFS="\t" 'function abs(v) {return v < 0 ? -v : v} {nbmax=0; 
				for(i=4;i<=NF;i++) { 
					if($i!= -1){ 
						if(nbmax!=0){
							ok=1;
							for(lmax in max){
								localmax=max[lmax];
								if(abs(localmax-$i)<=SIZE){ 
									moy[localmax]+=$i; nb[localmax]+=1; ok=0; break 
								}
							} 
							if(ok==1){ 
								max[i-3]=$i; moy[$i]=$i;nb[$i]=1
							} 
						} 
						if(nbmax==0){
							max[i-3]=$i; moy[$i]=$i;nb[$i]=1;nbmax=1
						} 
					}
				}; 
				for(kmax in max){
					localmax=max[kmax]; 
					print($1, int(moy[localmax]/nb[localmax])-SIZE, int(moy[localmax]/nb[localmax])+SIZE)
					delete max[kmax];delete nb[localmax]; delete moy[localmax]
				}
				}' $out_dir/$name_cons/temp/peak.bed | sed 's/\t\+/\t/g;s/^\t//' | awk -v OFS="\t" '{if($2>-SIZE){print $0}}' | sortBed -i >> $out_dir/$name_cons/temp/consensus_pics_4MA.bed
		done < $out_dir/$name_cons/temp/ligned_peaks_4AMR.bed 
			
		cat $out_dir/$name_cons/temp/consensus_pics_3MA.bed $out_dir/$name_cons/temp/consensus_pics_4MA.bed | awk -v OFS="\t" '{print $1,$2,$3}'| sortBed -i > $out_dir/$name_cons/temp/LFY_CT_consensus_tmp.bed   
		
		# ------- Filter peaks via blacklist
		echo "Filtre via la blacklist" | tee -a $out_dir/$name_cons/log.txt
		bedtools intersect -v -a $out_dir/$name_cons/temp/LFY_CT_consensus_tmp.bed  -b $blacklist -wa > $out_dir/$name_cons/temp/LFY_CT_consensus_filtered.bed

		# ------- Compute normalized coverage per million reads in peaks and generate consensus files
		RIP_list=()
		BDG_list=()
		for dir in ${in_dir[@]};
			do
			local sample=${dir##*/}
			echo "scaling coverages for: $sample" | tee -a $out_dir/$name_cons/log.txt
			# ------- Compute RIP in peaks
			samplebam=$dir/${sample}.filtered.sorted.bam
			conspeaks=$out_dir/$name_cons/temp/LFY_CT_consensus_filtered.bed
			inpeaks=$(bedtools sort -i $conspeaks | bedtools merge -i stdin | bedtools intersect -u -a $samplebam -b stdin -ubam | samtools view -c)
			r_factor=$(calc 1000000/$inpeaks) 
			bedtools coverage -a $conspeaks -b $samplebam -sorted -F 1 | awk -v OFS="\t" -v a=$r_factor '{print $1,$2,$3,$4,$5=$5*a}' >  $out_dir/${sample}/${sample}_RPMinpeaks.bed

			# ------- Compute scaled bedgraph
			local libsize=$(samtools view -f 0 -c $bam)
			local scale=$(calc 1000000/$libsize)
			bedtools genomecov -ibam $samplebam -bga -scale $scale | sed 's/Chr/chr/g' > $out_dir/${sample}/${sample}_CPM.bdg
			RIP_list+=("$out_dir/${sample}/${sample}_RPMinpeaks.bed")
			BDG_list+=("$out_dir/${sample}/${sample}_CPM.bdg")
		done

		# ------- Verify peak consistency accross replicates and compute final consensus values (2.5 factor RiP)
		paste "${RIP_list[@]}" | awk -v OFS="\t" '{print $1,$2,$3,$5,$10,$15}' | awk -v OFS="\t" '{if($4/$5 > 0.4 && $4/$5 < 2.5 && $4/$6 > 0.4 && $4/$6 < 2.5){print $1,$2,$3,$4=($5+$10+$15)/3}}' > $out_dir/$name_cons/temp/LFY_CT_consensus_filtered_verified.bed

		# ------- Generate final consensus bedgraph and bigwig files
		echo "Merging bedgraphs for consensus" | tee -a $out_dir/$name_cons/log.txt
		time $PATH_TO_MACS/macs2 cmbreps -i ${BDG_list[@]} -m mean -o $out_dir/$name_cons/temp/tmp.bdg
		cat $out_dir/$name_cons/temp/tmp.bdg | sortBed -i > $out_dir/$name_cons/${name_cons}_cpm.bdg
		$bgtbw_path/bedGraphToBigWig $out_dir/$name_cons/${name_cons}_cpm.bdg /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt $out_dir/$name_cons/${name_cons}_cpm.bw

		# ------ Center peaks on maximum coverage position
		echo "peaks are been centralized on maximum" | tee -a $out_dir/$name_cons/log.txt
		bedtools intersect -b $out_dir/$name_cons/${name_cons}_cpm.bdg -a $out_dir/$name_cons/temp/LFY_CT_consensus_filtered_verified.bed -wo > $out_dir/$name_cons/temp/peaks_tmp.bed
		awk -v OFS="\t" -v max=0 -v ch=0 -v pss=0 -v pts=0 -v id_old=0 '(id=$1$2$3) {
			if(id!=id_old){
				if(id_old==0){
					max=$8
					ch=$1
					pss=$6
					pts=$7
				} else {
					print ch,pss,pts,max
					max=$8
					ch=$1
					pss=$6
					pts=$7
				}
			}else{
				if($8>max){
					max=$8
					ch=$1
					pss=$6
					pts=$7
				}
			} id_old=id
		} END{print ch,pss,pts,max}' $out_dir/$name_cons/temp/peaks_tmp.bed | awk -v SIZE=$phs -v OFS='\t' '(nt=int(($2+$3)/2)) {print $1,nt-SIZE,nt+SIZE,$4}' | sort -k4,4nr > $out_dir/$name_cons/${name_cons}_maxMean.bed

		awk -v OFS="\t" '{print $1,$2,$3}' $out_dir/$name_cons/${name_cons}_maxMean.bed | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}_narrow.bed
	fi
	LC_COLLATE=$tmpsave

	# # ------- SUPPLEMNTAL FILE WITH NUMBER OF PEAKS
	peaksNb=$(wc -l < "$out_dir/$name_cons/${name_cons}_narrow.bed")
	echo "$peaksNb" > "$out_dir/$name_cons/peaksNb.txt"
	
}

#-------------------------------------------------------------------------------
replicates_comparisons (){
	# FUNCTION: replicates_comparisons
	# DESCRIPTION:
	#   This function compares replicates for peak calling and coverage computation. 
	#   It supports two normalization methods: "Not Input Normalized" and "Input Re-Normalized". 
	#
	# USAGE:
	#   replicates_comparisons -p <PATH> -nb <INT> -od <PATH>
	#                          -gd <PATH> -n <STRING> -c [STRING]
	#
	# ARGUMENTS:
	#   -p          : Path to the consensus peaks file. Required.
	#   -nb         : Number of replicates in the sample. Required.
	#   -od         : Output directory for results. Required.
	#   -gd         : General directory containing `Peakcalling` and `Mapping` subdirectories. Required.
	#   -n          : Name of the sample. Required.
	#   -c          : Color for dots in plots (hex code, e.g., `#FF0000` or `FF0000`). Optional. Default: `#000000`.
	#   -h, --help  : Display usage information for this function.
	# 
	# DEPENDENCIES:
	#   - compute_rpkmrip_rpkmril
	#   - Pairwise_comps (R script)

	echo "
	____ ____ ___      ____ ____ _  _ ___  ____ ____ _ ____ ____ _  _
	|__/ |___ |__]     |    |  | |\/| |__] |__| |__/ | [__  |  | |\ |
	|  \ |___ |    ___ |___ |__| |  | |    |  | |  \ | ___] |__| | \|
	"
	#TODO: add option to show consensus peaks and/or all peaks from all replicates. (JL 09/2025; from comments of FP) ## JL: I'll do a branch for this soon
	local color="#000000";
	local Argsline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Consensus peaks file set to: ${2}";shift 2;;
			-nb)
				local nbrep=$2
				echo "-> Number of replicates in the sample set to: ${2}";shift 2;;
			-od)
				local result=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-gd)
				local general_dir=$2
				echo "-> General result directory set to: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> Name of sample set to: ${2}";shift 2;;
			-c)
				local color=$2;
				echo "-> Color set for dots: $color"; shift 2;;
			-h)
				usage replicates_comparisons ; return;;
			--help)
				usage replicates_comparisons ; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage replicates_comparisons ; return;;
		esac
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument is needed"; Errors+=1; fi
	if [ -z $nbrep ]; then echo "ERROR: -nb argument is needed"; Errors+=1; fi
	if [ -z $result ]; then echo "ERROR: -od argument is needed"; Errors+=1; fi
	if [ -z $name ]; then echo "ERROR: -n argument is needed"; Errors+=1; fi
	if [ -z $general_dir ]; then echo "ERROR: -gd argument is needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage replicates_comparisons ; return 1; fi

	[[ ${color} =~ ^#.* ]] || color="#${color}"

	# initialize lists
	local rep_list=()
	local bdg_list=()
	local remove_list=()
	mkdir -p -m 774 $result/$name/NotInputNormalized $result/$name/InputReNormalized
	local log=$result/$name/log.txt

	echo "" > $log
	echo $Argsline >> $log

	# prepare lists of replicates files
	for rep in $(seq $nbrep)
	do
		rep_list+=("${name}rep${rep}")
		bdg_list+=("$general_dir/Peakcalling/${name}rep${rep}/${name}rep${rep}_cpm.bdg")
		rep_peaks_list+=("$general_dir/Peakcalling/${name}rep${rep}/${name}rep${rep}_narrow.bed")
		remove_list+=("$result/$name/NotInputNormalized/${name}rep${rep}_filt.cov.bed" "$result/$name/InputReNormalized/${name}rep${rep}_filt.cov.bed" "$result/$name/NotInputNormalized/${name}rep${rep}_cpmrip.bdg" "$result/$name/InputReNormalized/${name}rep${rep}_cpmrip.bdg")
	done

	#concatenate peaks from all replicates, keeping track of their origin
	list_data=("$general_dir/Peakcalling")
	list_bamdir=("$general_dir/Mapping")

	# Not Input Normalized method
	echo "[INFO] - Computing coverages with RPKM-RiP method (not Input normalized)" | tee -a $log
	compute_rpkmrip_rpkmril -p $peaks -bd list_bamdir[@] -pd list_data[@] -sn rep_list[@] -o $result/$name/NotInputNormalized/ -m "inPeaks"
	$PATHRscript $Pairwise_comps $result/$name/NotInputNormalized/ $result/$name/NotInputNormalized/peaks_perSample_rpkminPeaks.txt $color

	# Input ReNormalized method
	echo "[INFO] - Computing coverages with RPKM-RiL method (Input Re-Normalized)" | tee -a $log
	compute_rpkmrip_rpkmril -p $peaks -bd list_bamdir[@] -pd list_data[@] -sn rep_list[@] -o $result/$name/InputReNormalized/ -m "inLibs"
	$PATHRscript $Pairwise_comps $result/$name/InputReNormalized/ $result/$name/InputReNormalized/peaks_perSample_rpkminLibs.txt  $color

	# Cleaning 
	echo "[INFO] - Cleaning temporary files" | tee -a $log
	for f in ${remove_list[@]}; do
		if [ -f $f ]; then
			rm $f
		else
			echo "No $f to remove !" | tee -a $log
		fi
	done
	# rm ${remove_list[@]}
	rm $result/$name/*/tmp_peaks
	rm $result/$name/NotInputNormalized/tmpFiltTags.txt || echo "No $result/$name/NotInputNormalized/tmpFiltTags.txt to remove !"
	rm $result/$name/InputReNormalized/tmpFiltTags.txt || echo "No $result/$name/InputReNormalized/tmpFiltTags.txt to remove !"
}

#-------------------------------------------------------------------------------
pairwize_comparison(){
	# FUNCTION: pairwize_comparison
	# DESCRIPTION:
	#   This function performs pairwise comparison of two datasets (e.g., ChIP-seq replicates or conditions).
	#   It merges peaks, computes coverage (RPKM) in peaks and libraries, and generates comparative plots.
	#   The function supports filtering by coverage/height, extended peak mode, and custom colors for visualization.
	#   It also handles replicates and control samples for normalization.
	#
	# USAGE:
	# pairwize_comparison -n1 <STRING> -n2 <STRING> -od <PATH> -id <PATH> 
	#        -f [INT] -he [FLOAT] -bd [PATH]  -rep1 [LIST of STRING] 
	#        -rep2 [LIST of STRING] -p1 [PATH] -p2 [PATH] -gcov [STRING] -ctrl1 [PATH]
	#        -ctrl2 [PATH] -ext [BOOL] -c1 [STRING] -c2 [STRING] -cn [STRING] 
	#        -he [INT] -id [PATH] -id2 [PATH] -bd2 [PATH] -th [FLOAT] -h
	#
	# ARGUMENTS:
	#   -n1         : Name of the first dataset. Required.
	#   -n2         : Name of the second dataset. Required.
	#   -od         : Directory where results will be saved. Required.
	#   -id         : General data directory (peak calling directory) for the first dataset. Required.
	#   -id2        : General data directory for the second dataset. Optional. Default: same as -id.
	#   -f          : Minimum coverage required in both samples to consider a peak. Optional. Default: 0.
	#   -he         : Minimum height required in both samples to consider a peak. Optional. Default: 0.
	#   -bd         : General BAM directory for the first dataset. Optional.
	#   -p1         : Path to the peaks file for the first dataset. Optional. Default: "NA".
	#   -p2         : Path to the peaks file for the second dataset. Optional. Default: "NA".
	#   -bd2        : General BAM directory for the second dataset. Optional. Default: same as -bd.
	#   -gcov       : Whether to compute reads count from BAM at each peak ("yes" or "no"). Optional. Default: "yes".
	#   -rep1       : List of replicate names for the first dataset. Optional.
	#   -rep2       : List of replicate names for the second dataset. Optional.
	#   -ctrl1      : List of control BAM files for the first dataset. Optional. Default: ("null").
	#   -ctrl2      : List of control BAM files for the second dataset. Optional. Default: ("null").
	#   -ext        : Activate extended peaks mode. Optional.
	#   -c1         : Color for the first dataset. Optional. Default: "#2A9D8F".
	#   -c2         : Color for the second dataset. Optional. Default: "#B8475C".
	#   -cn         : Neutral color for common peaks. Optional. Default: "#BEBAB7".
	#   -h, --help  : Display usage information for this function.
	#   -th		    : Threshold for fold change in differential analysis. Optional. Default: 1.
	#
	# DEPENDENCIES:
	#   - compute_rpkmrip_rpkmril
	#   - Pairwise_comps (R script)
	#
	# NOTES:
	#   - The function requires the following tools: `bedtools`
	#   - If replicates are not provided, the function attempts to find them in the BAM directory.
	#   - If extended mode is activated, the function uses extended peak positions.
	#   - The function supports filtering peaks by coverage or height.
	#   - If both datasets have replicates, edgeR analysis is performed for differential analysis.

	echo "
	___  ____ _ ____ _ _ _ _ ___  ____     ____ ____ _  _ ___  ____ ____ _ ____ ____ _  _ 
	|__] |__| | |__/ | | | |   /  |___     |    |  | |\/| |__] |__| |__/ | [__  |  | |\ | 
	|    |  | | |  \ |_|_| |  /__ |___ ___ |___ |__| |  | |    |  | |  \ | ___] |__| | \|                                                      
	"
	local get_bedtools_cov="yes"; local filterCov=0; local filterHeight=0; local bamdir2="NA"; local data2="NA"; local peaksext1="NA"; local peaksext2="NA"; local extendedmode=false; local list_ctrl1=("null"); local list_ctrl1=("null"); local color1="2A9D8F"; local color2="B8475C"; local colorN="BEBAB7"; local thresholdFC=1;
	local Argsline=$@; local keep_original="no"
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-n1)
				local name1=$2
				echo "-> Name of directory for dataset 1 is: ${2}";shift 2;;
			-n2)
				local name2=$2
				echo "-> Name of directory for dataset 2 is: ${2}";shift 2;;
			-od)
				local result=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-id)
				local data=$2
				echo "-> General data directory (i.e. peakcalling directory) set to: ${2}";shift 2;;
			-id2)
				local data2=$2
				echo "-> General data directory (i.e. peakcalling directory) for second dataset set to: ${2}";shift 2;;
			-f)
				local filterCov=$2
				echo "-> Coverage must be at least >${2} in both sample to be considered as peaks";shift 2;;
			-he)
				local filterHeight=$2
				echo "-> Height must be at least >${2} in both sample to be considered as peaks";shift 2;;
			-bd)
				local bamdir=$2
				echo "-> General bam directory set to: ${2}";shift 2;;
			-p1)
				local peaksext1=$2
				echo "-> Peaks to use for dataset 1 set to: ${2}";shift 2;;
			-p2)
				local peaksext2=$2
				echo "-> Peaks to use for dataset 2 set to: ${2}";shift 2;;
			-bd2)
				local bamdir2=$2
				echo "-> General bam directory for second dataset set to: ${2}";shift 2;;
			-gcov)
				local get_bedtools_cov=$2
				echo "-> 'yes' or 'no' to compute reads count from bam at each peak: ${2}";shift 2;;
			-rep1)
				local list_rep1=("${!2}")
				echo "-> List of replicates names for dataset 1 in bam directory set to: ${list_rep1[@]}";shift 2;;
			-rep2)
				local list_rep2=("${!2}")
				echo "-> List of replicates names for dataset 2 in bam directory set to: ${list_rep2[@]}";shift 2;;
			-ctrl1)
				local list_ctrl1=("${!2}"); 
				echo "-> List of control .bam files for dataset 1: ${list_ctrl1[@]}"; shift 2;;
			-ctrl2)
				local list_ctrl2=("${!2}"); 
				echo "-> List of control .bam files for dataset 2: ${list_ctrl2[@]}"; shift 2;;
			-ext)
				local extendedmode=true; 
				echo "-> Extended peaks mode: true "; shift 1;;
			-c1)
				local color1=$2; 
				echo "-> Color for set1: "$color1; shift 2;;
			-c2)
				local color2=$2;
				echo "-> Color for set2: "$color2; shift 2;;
			-cn)
				local colorN=$2;
				echo "-> Neutral color for common peaks set to: "$colorN; shift 2;;
			-keep)
				local keep_original=$2;
				echo "-> Keeping original peaks to compute additional RiP before merging ("yes" or "no"): "$keep_original; shift 2;;
			-th)
				local thresholdFC=$2;
				echo "-> Threshold for fold change to consider a peak as specific to one of the two datasets (e.g., 2 for 2-fold change): "$thresholdFC; shift 2;;
			-h)
				usage pairwize_comparison; return;;
			--help)
				usage pairwize_comparison; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage pairwize_comparison; return;;
		esac
	done
	local Errors=0
	if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
	if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
	if [ -z $result ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ -z $data ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage pairwize_comparison; return 1; fi
	if [ $data2 == "NA" ]; then
		echo "[WARNING] - id2 not used, assuming same Peakcalling directory for both datasets"
		local data2=$data
	fi
	if [ $bamdir2 == "NA" ]; then
		echo "[WARNING] - bd2 not used, assuming same bam directory for both datasets"
		local bamdir2=$bamdir
	fi

	# Validating color codes
	[[ ${color1} =~ ^#.* ]] || color1="#${color1}"
	[[ ${color2} =~ ^#.* ]] || color2="#${color2}"
	[[ ${colorN} =~ ^#.* ]] || colorN="#${colorN}"

	local out_dir=$result/${name1}_${name2}
	mkdir -p -m 774 $out_dir
	# touch $out_dir/log.txt
	local log=$out_dir/log.txt
	echo $Argsline > $log

	# ---------- If more than 1 replicate, get the wextendedpos.bed instead of the narrow.bed file
	if (( ${#list_rep1[@]} > 1 && ${#list_rep2[@]} > 1 )); then #checking if more than 1 replicat in both dataset 
		echo "[INFO] - More than one replicate detected: switching to extended peaks positions for first sample" | tee -a $log
		local changed=0
		if [ -f "${peaksext1%/*}/${name1}_wextendedpos.bed" ]; then #check the existence of an extended bedfile (dataset 1)
			peaksext1="${peaksext1%/*}/${name1}_wextendedpos.bed"
			changed=$changed+1 #switch to extended mode (dataset 1)
		fi
		echo "[INFO] - Peaks file path used for ${name1}: ${peaksext1}" | tee -a $log
		echo "[INFO] - More than one replicate detected: switching to extended peaks positions for second sample" | tee -a $log
		if [ -f "${peaksext2%/*}/${name2}_wextendedpos.bed" ]; then #check the existence of an extended bedfile (dataset 1)
			peaksext2="${peaksext2%/*}/${name2}_wextendedpos.bed"
			changed=$changed+1 #switch to extended mode (dataset 2)
		fi
		echo "[INFO] - Peaks file path used for ${name2}: ${peaksext2}" | tee -a $log

		# switch to extended mode or not 
		if [ $changed -eq 2 ]; then
			extendedmode=true
			echo "[INFO] - Extended mode activated" | tee -a $log
		else
			extendedmode=false
			echo "[INFO] - Extended mode deactivated (at least one of the extended peaks file not found)" | tee -a $log
		fi
	else
		extendedmode=false
		echo "[INFO] - Extended mode deactivated" | tee -a $log
	fi

	# Preparing peaks files
	local cov1=$data/$name1/${name1}_cov.bdg
	local cov2=$data2/$name2/${name2}_cov.bdg

	# Selecting peaks on their height:
	if [ $peaksext1 == "NA" ]; then
		if [ $filterHeight -eq 0 ];then
			local peaks1=$data/$name1/${name1}_narrow.bed # narrow bed as default if no filter applied on peaks height
		else #filter based on chosen height
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
			local peaks1=$data/$name1/${name1}_narrow.heightFiltered.bed
		fi
	else
		local peaks1=$peaksext1
	fi
	if [ $peaksext2 == "NA" ]; then
		if [ $filterHeight -eq 0 ];then
			local peaks2=$data2/$name2/${name2}_narrow.bed
		else
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data2/$name2/${name2}_max.bed > $data2/$name2/${name2}_narrow.heightFiltered.bed
			local peaks2=$data2/$name2/${name2}_narrow.heightFiltered.bed
		fi
	else
		local peaks2=$peaksext2
	fi

	echo "[INFO] - Preparing peaks files for comparison" | tee -a $log
	if [ $extendedmode == true ]; then
		# extended peaks mode
		awk -v OFS="\t" -v name=$name1 'NF==5{print $1,$4,$5,name,$2,$3;next}{print $1,$2,$3,name,$2,$3}' $peaks1  > $out_dir/${name1}_peaks.bed
		awk -v OFS="\t" -v name=$name2 'NF==5{print $1,$4,$5,name,$2,$3;next}{print $1,$2,$3,name,$2,$3}' $peaks2  > $out_dir/${name2}_peaks.bed
	else
		# normal peaks mode
		awk -v OFS="\t" -v name=$name1 '{print $1,$2,$3,name}' $peaks1  > $out_dir/${name1}_peaks.bed
		awk -v OFS="\t" -v name=$name2 '{print $1,$2,$3,name}' $peaks2  > $out_dir/${name2}_peaks.bed
	fi

	# ---------- Merging peaks
	cat $out_dir/${name1}_peaks.bed $out_dir/${name2}_peaks.bed | sort -k1,1 -k2,2n > $out_dir/${name1}_${name2}_peaks.bed
	echo "[INFO] - Merging peaks files" | tee -a $log
	local peak_file=$out_dir/${name1}_${name2}_peaks_processed.bed
	
	if [ $extendedmode == true ]; then
		# extended peaks mode
		python ${PATH_TO_COMPIL}/bin/merge_all_peaks_vextended.py -f1 $name1 -f2 $name2 -o $result
		cat $out_dir/${name1}_${name2}_peaks_merged.bed <(bedtools intersect -a $out_dir/${name1}_${name2}_peaks_uniques.bed -b $out_dir/${name1}_${name2}_peaks_merged.bed -v -f 0.90 -wa) | sort -k1,1 -k2,2n | uniq > $out_dir/${name1}_${name2}_peaks_processedext.bed
		awk -v OFS="\t" '{print $1,$5-200,$5+200,$4}' $out_dir/${name1}_${name2}_peaks_processedext.bed | sort -k1,1 -k2,2n > $peak_file
	else
		# normal peaks mode
		python $merge_peaks -f1 $name1 -f2 $name2 -o $result
		cat $out_dir/${name1}_${name2}_peaks_merged.bed <(bedtools intersect -a $out_dir/${name1}_${name2}_peaks_uniques.bed -b $out_dir/${name1}_${name2}_peaks_merged.bed -v -f 0.90 -wa) | sort -k1,1 -k2,2n | uniq > $out_dir/${name1}_${name2}_peaks_processed.bed
	fi
	
	# ---------- Computing coverages in peaks (RPKM-RiP) for both samples
	local peak_file=$out_dir/${name1}_${name2}_peaks_processed.bed
	paste $peak_file | sed "1ichr\tbegin\tend\tname" >  $out_dir/table_${name1}_${name2}.csv #init table file
	
	if [ ! -z $bamdir ] && [ -z $list_rep1 ]; then # Bamdir defined but not list_rep
		echo "[INFO] - Searching replicates in bam directory if not provided" | tee -a $log
		local list_rep1=()
		local list_rep2=()
		for rep in {1..10}; do 
			local foundrep1=0
			local foundrep2=0
			if [[ -d $bamdir/${name1}rep${rep} ]]; then
				echo "[INFO] - found $bamdir/${name1}rep${rep} as replicate for $name1" | tee -a $log
				list_rep1+=("${name1}rep${rep}")
				foundrep1=1
			fi
			if [[ -d $bamdir/${name1}_rep${rep} ]]; then
				echo "[INFO] - found $bamdir/${name1}_rep${rep} as replicate for $name1" | tee -a $log
				list_rep1+=("${name1}_rep${rep}")
				foundrep1=1
			fi

			if [[ -d $bamdir2/${name2}rep${rep} ]]; then
				echo "[INFO] - found $bamdir2/${name2}rep${rep} as replicate for $name2" | tee -a $log
				list_rep2+=("${name2}rep${rep}")
				foundrep2=1
			fi
			if [[ -d $bamdir2/${name2}_rep${rep} ]]; then
				echo "[INFO] - found $bamdir2/${name2}_rep${rep} as replicate for $name2" | tee -a $log
				list_rep2+=("${name2}_rep${rep}")
				foundrep2=1
			fi
			if [[ $foundrep1 -eq 0 ]] && [[ $foundrep2 -eq 0 ]]; then
				break # exit the loop if no more replicates found for both samples
			fi
		done
	
	fi
	if [ -z $bamdir ] && [ -z $list_rep1 ]; then # Bamdir & list_rep not defined
		$PATHRscript $do_plot_quant $result $name1 $name2 $out_dir/table_${name1}_${name2}.csv RiL 
	else
		echo "[INFO] - Saved replicates for ${name1}: ${list_rep1[@]}" | tee -a $log
		echo "[INFO] - Saved replicates for ${name2}: ${list_rep2[@]}" | tee -a $log

		if [ ${list_ctrl1[0]} != "null" ]; then
			echo "[INFO] - Saved control bam files for ${name1}: ${list_ctrl1[@]}" | tee -a $log
			local list_all=("${list_rep1[@]}" "${list_rep2[@]}" "${list_ctrl1[@]}" "${list_ctrl2[@]}")
		else
			echo "[INFO] - No control bam files provided for ${name1}" | tee -a $log
			local list_all=("${list_rep1[@]}" "${list_rep2[@]}")
		fi
		local list_bamdir=("$bamdir" "$bamdir2")
		local list_data=("$data" "$data2")

		# ---------- Compute RPKM in libraries
		echo "[INFO] - Computing coverages with RPKM-RiL method" | tee -a $log
		compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd list_bamdir[@] -pd list_data[@] -sn list_all[@] -m "inLibs" -o $out_dir -g $get_bedtools_cov

		# Computing RPKM in libraries with control re-normalization if controls provided
		if [ ${list_ctrl1[0]} != "null" ]; then 
			last_col_n1="$((${#list_rep1[@]}-1+4))"
			first_col_n2="$(($last_col_n1+1))"
			first_ctrl_1="$(($first_col_n2+${#list_rep2[@]}))"
			first_ctrl_2="$(($first_ctrl_1+${#list_ctrl1[@]}))"
			# echo $last_col_n1 $first_col_n2 $first_ctrl_1 $first_ctrl_2 # DEBUG
			awk -v OFS="\t" -v c=$last_col_n1 -v fc1=$first_ctrl_1 -v lc1=$(($first_ctrl_2)) '(NR==1 && $2!="begin" && $2!="start") || NR>1{moy=0;uncor=0;nb=0;moyc=0;nbc=0;for(i=fc1;i<lc1;i++){nbc++;moyc+=$i};for(i=4;i<=c;i++){nb++;moy+=($i-moyc/nbc);uncor+=$i};val=((moy>0)?moy:(0.01*nb))/nb;print $1,$2,$3,val,uncor/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name1}_RPKMril.txt
		
			awk -v OFS="\t" -v c=$first_col_n2 -v fc1=$(($first_ctrl_1 - 1)) -v lc1=$(($first_ctrl_2)) '(NR==1 && $2!="begin" && $2!="start") || NR>1 {moy=0;uncor=0;nb=0;moyc=0;nbc=0;for(i=lc1;i<=NF;i++){nbc++;moyc+=$i};for(i=c;i<=fc1;i++){nb++;moy+=($i-moyc/nbc);uncor+=$i};val=((moy>0)?moy:(0.01*nb))/nb;print $1,$2,$3,val,uncor/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name2}_RPKMril.txt
		else
			last_col_n1="$((${#list_rep1[@]}-1+4))"
			first_col_n2="$(($last_col_n1+1))"
			awk -v OFS="\t" -v c=$last_col_n1 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name1}_RPKMril.txt
			awk -v OFS="\t" -v c=$first_col_n2 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name2}_RPKMril.txt
			
		fi

		local list_all=("${list_rep1[@]}" "${list_rep2[@]}")

		# ---------- Compute RPKM in peaks
		echo "[INFO] - Computing coverages with RPKM-RiP method" | tee -a $log
		compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd list_bamdir[@] -pd list_data[@] -sn list_all[@] -m "inPeaks" -o $out_dir -g "no"
		awk -v OFS="\t" -v c=$last_col_n1 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminPeaks.txt > $out_dir/${name1}_RPKMrip.txt
		awk -v OFS="\t" -v c=$first_col_n2 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminPeaks.txt > $out_dir/${name2}_RPKMrip.txt

		# ---------- For those who want to compute original and not merged peaks RiP
	if [[ $keep_original == "yes" ]]; then 

		mkdir -p -m 774 $out_dir/${name1}_peaks/${name1}_rpkm/
		mkdir -p -m 774 $out_dir/${name2}_peaks/${name2}_rpkm/
		mkdir -p -m 774 $out_dir/${name1}_peaks/${name2}_rpkm/
		mkdir -p -m 774 $out_dir/${name2}_peaks/${name1}_rpkm/

		list_bamdir_1=("$bamdir")
		list_data_1=("$data")
		list_rep1=("${list_rep1[@]}")
		list_bamdir_2=("$bamdir2")
		list_data_2=("$data2")
		list_rep2=("${list_rep2[@]}")

		#  1) Peaks dataset1 + coverage dataset1 
		compute_rpkmrip_rpkmril -p $out_dir/${name1}_peaks.bed -bd list_bamdir_1[@] -pd list_data_1[@] -sn list_rep1[@] -m "inPeaks" -o $out_dir/${name1}_peaks/${name1}_rpkm/ -g "yes"
		awk -v OFS="\t" 'NR==1 {print $1,$2,$3,"mean_RPKM"; next}{s=0;n=0; for(i=4;i<=NF;i++){if($i!="" && $i!="NA"){s+=$i;n++}} m=(n>0?s/n:0); print $1,$2,$3,m }' $out_dir/${name1}_peaks/${name1}_rpkm/peaks_perSample_rpkminPeaks.txt > $out_dir/${name1}_peaks/${name1}_rpkm/${name1}_RPKMrip.txt
		#  2) Peaks dataset1 + coverage dataset2 
		compute_rpkmrip_rpkmril -p $out_dir/${name1}_peaks.bed -bd list_bamdir_2[@] -pd list_data_2[@] -sn list_rep2[@] -m "inPeaks" -o $out_dir/${name1}_peaks/${name2}_rpkm/ -g "yes"
		awk -v OFS="\t" 'NR==1 {print $1,$2,$3,"mean_RPKM"; next}{s=0;n=0; for(i=4;i<=NF;i++){if($i!="" && $i!="NA"){s+=$i;n++}} m=(n>0?s/n:0); print $1,$2,$3,m }' $out_dir/${name1}_peaks/${name2}_rpkm/peaks_perSample_rpkminPeaks.txt > $out_dir/${name1}_peaks/${name2}_rpkm/${name2}_RPKMrip_in_${name1}_peaks.txt
		#  3) Peaks dataset2 + coverage dataset2 
		list_bamdir_2=("$bamdir2")
		list_data_2=("$data2")
		list_rep2=("${list_rep2[@]}")
		compute_rpkmrip_rpkmril -p $out_dir/${name2}_peaks.bed -bd list_bamdir_2[@] -pd list_data_2[@] -sn list_rep2[@] -m "inPeaks" -o $out_dir/${name2}_peaks/${name2}_rpkm/ -g "yes"
		awk -v OFS="\t" 'NR==1 {print $1,$2,$3,"mean_RPKM"; next}{s=0;n=0; for(i=4;i<=NF;i++){if($i!="" && $i!="NA"){s+=$i;n++}} m=(n>0?s/n:0); print $1,$2,$3,m }' $out_dir/${name2}_peaks/${name2}_rpkm/peaks_perSample_rpkminPeaks.txt > $out_dir/${name2}_peaks/${name2}_rpkm/${name2}_RPKMrip.txt
		#  4) Peaks dataset2 + coverage dataset1 
		compute_rpkmrip_rpkmril -p $out_dir/${name2}_peaks.bed -bd list_bamdir_1[@] -pd list_data_1[@] -sn list_rep1[@] -m "inPeaks" -o $out_dir/${name2}_peaks/${name1}_rpkm/ -g "yes"
		awk -v OFS="\t" 'NR==1 {print $1,$2,$3,"mean_RPKM"; next}{s=0;n=0; for(i=4;i<=NF;i++){if($i!="" && $i!="NA"){s+=$i;n++}} m=(n>0?s/n:0); print $1,$2,$3,m }' $out_dir/${name2}_peaks/${name1}_rpkm/peaks_perSample_rpkminPeaks.txt > $out_dir/${name2}_peaks/${name1}_rpkm/${name1}_RPKMrip_in_${name2}_peaks.txt

	fi

	# ---------- Generating comparison plots
	echo '~###################################'
		if [ -f $out_dir/${name1}_RPKMrip.txt ] && [ -f $out_dir/${name2}_RPKMrip.txt ] && [ -f $out_dir/${name1}_RPKMril.txt ] && [ -f $out_dir/${name2}_RPKMril.txt ]; then
			# Both RiP and RiL available
			if [ ${#list_rep1[@]} -gt 1 ] && [ ${#list_rep2[@]} -gt 1 ]; then 
				echo "[INFO] - Both datasets have replicates, launching edgeR analysis" | tee -a $log
				$PATHRscript $StatEdgeR $out_dir/peaks_perSample_rpkminPeaks.txt $out_dir ${#list_rep1[@]} ${#list_rep2[@]} 
				paste <(awk 'NR==1 && $2!="begin"{print $0}NR>1{print $0}' $out_dir/table_${name1}_${name2}.csv) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMrip.txt) <(awk -v OFS="\t" 'NR==1 && $2!="start"{print $4,$5}NR>1{print $4,$5}' $out_dir/${name1}_RPKMril.txt) <(awk -v OFS="\t" 'NR==1 && $2!="start"{print $4,$5}NR>1{print $4,$5}' $out_dir/${name2}_RPKMril.txt) <(awk -v OFS="\t" 'NR!=1 {print $2,$6}' $out_dir/results_edgeR.tsv) | sed "1ichr\tbegin\tend\tname\t${name1}_RiP\t${name2}_RiP\t${name1}_RiL\t${name1}_RiLuc\t${name2}_RiL\t${name2}_RiLuc\tlogFC\tFDR" > $out_dir/table_${name1}_${name2}_RiL_RiP.tsv
				paste <(awk -v OFS="\t" '{print $1,$2,$3,$4}' $out_dir/table_${name1}_${name2}.csv) <(awk '{s = ""; for (i = 4; i <= NF; i++){if(s==""){s=$i}else{s = s"\t"$i}}; print s}' $out_dir/peaks_perSample_rpkminPeaks.txt) <(awk '{s = ""; for (i = 4; i <= NF; i++){if(s==""){s=$i}else{s = s"\t"$i}}; print s}' $out_dir/peaks_perSample_rpkminLibs.txt) <(awk -v OFS="\t" '{print $2,$6}' $out_dir/results_edgeR.tsv) > $out_dir/table_${name1}_${name2}_RiL_RiP2.tsv
			# Only one or none have replicates
			else
				paste <(awk 'NR==1 && $2!="begin"{print $0}NR>1{print $0}' $out_dir/table_${name1}_${name2}.csv) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMril.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMril.txt) | sed "1ichr\tbegin\tend\tname\t${name1}_RiP\t${name2}_RiP\t${name1}_RiL\t${name2}_RiL" > $out_dir/table_${name1}_${name2}_RiL_RiP.tsv
			fi
		fi
		# Generating plots
		echo '~###################################'
		# TODO compaisonplotmaker Rscript is not defineda t all ? 
		$PATHRscript $comparisonPlotMaker -od $result/${name1}_${name2} -n1 $name1 -n2 $name2 -t $out_dir/table_${name1}_${name2}_RiL_RiP.tsv -c1 $color1 -c2 $color2 -cn $colorN -th ${thresholdFC}
	fi
	if [ ! -z $bamdir ] || [ ! -z $list_rep1 ]; then
		rm $out_dir/*.bdg
	fi

}

#-------------------------------------------------------------------------------
analyze_decile(){
	# FUNCTION: analyze_decile
	# DESCRIPTION:
	#   This function performs a decile-based analysis of peak datasets, focusing on motif discovery
	#   and spacing analysis between peaks. It divides the peaks into deciles based on their
	#   normalized coverage ratio (RiP or RiL), computes motifs for each decile, and analyzes
	#   the spacing between peaks using provided matrices. The results are visualized using heatmaps.
	#
	# USAGE:
	# analyze_decile -n1 <STRING> -n2 <STRING> -o <PATH> -g <FILE> 
	#        -m <LIST of FILE> -t STRING [-th <FLOAT> -thf <PATH>] -ol [INT] -or [INT]
	#        -mins [INT] -maxs [INT] [-sym] -d [INT] -mim [INT] -mam [INT] -l [INT] 
	#        -nm [INT] -s [INT]
	#
	# ARGUMENTS:
	#   -n1          : Name of the first dataset. Required.
	#   -n2          : Name of the second dataset. Required.
	#   -o           : Directory where results will be saved. Required.
	#   -t           : Type of normalization to use ("RiP" or "RiL"). Optional. Default: "RiP".
	#   -g           : Path to the genome FASTA file. Optional.
	#   -m           : List of matrix files for spacing analysis. Optional.
	#   -th          : List of thresholds for spacing analysis. Optional. Required if -thf is not provided.
	#   -thf         : Directory containing threshold files. Optional. Required if -th is not provided.
	#   -ol          : Left offset for spacing analysis. Optional. Default: 0.
	#   -or          : Right offset for spacing analysis. Optional. Default: 0.
	#   -mins        : Minimum spacing between peaks. Optional. Default: 0.
	#   -maxs        : Maximum spacing between peaks. Optional. Default: 30.
	#   -sym         : Use symmetric matrix mode for spacing analysis. Optional.
	#   -d           : Number of deciles to divide the peaks into. Optional. Default: 10.
	#   -mim         : Minimum motif length for motif discovery. Optional. Default: 5.
	#   -mam         : Maximum motif length for motif discovery. Optional. Default: 30.
	#   -l           : Size of the learning set for motif analysis. Optional. Default: 600.
	#   -nm          : Number of motifs to discover. Optional. Default: 5.
	#   -s           : Random seed for motif analysis. Optional. Default: 42.
	#   -h, --help   : Display usage information for this function.
	# 
	# DEPENDENCIES:
	#   - compute_motif
	#   - compute_space
	#
	# NOTES:
	#   - The function requires the input file `table_${name1}_${name2}_RiL_RiP.tsv`, which should be generated
	#     by the `pairwize_comparison` function beforehand.
	#   - If thresholds are not provided, the function attempts to retrieve them from the specified directory.
	#   - The function uses `compute_motif` for motif discovery and `compute_space` for spacing analysis.
	#   - The function supports both symmetric and asymmetric matrix modes for spacing analysis.
	#   - Heatmaps are generated using R scripts (`Heatmap.R` or `HeatmapAR.R`).

	echo "
	____ _  _ ____ _    _   _ ___  ____     ___  ____ ____ _ _    ____
	|__| |\ | |__| |     \_/    /  |___     |  \ |___ |    | |    |___
	|  | | \| |  | |___   |    /__ |___ ___ |__/ |___ |___ | |___ |___
	"
	local typeofnorm="RiP"; local offset_left=0; local offset_right=0; local min_spacing=0; local max_spacing=30; local thresholds_dir="null"; local thre=(0); local divide_in=10; local matrix_type="ASYMMETRIC"; local minmotiflength=5; local maxmotiflength=30; local learningsetSize=600; local numberofmotifs=5; local seed=42
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-n1)
				local name1=$2
				echo "-> Name of directory for dataset 1 is: ${2}";shift 2;;
			-n2)
				local name2=$2
				echo "-> Name of directory for dataset 2 is: ${2}";shift 2;;
			-o)
				local dir_comp=$2
				echo "-> Ouput directory set to: ${2}"; shift 2;;
			-t)
				local typeofnorm=$2
				echo "-> Type of normalization to use set to: ${2}"; shift 2;;
			-g)
				local genome=$2
				echo "-> Genome file set to: ${2}"; shift 2;;
			-m)
				local list_matrices=("${!2}")
				echo "-> Matrix file: ${list_matrices[@]}"; shift 2;;
			-th)
				local thre=("${!2}")
				echo "-> Thresholds set to: ${thre[@]}";shift 2;;
			-thf)
				local thresholds_dir=$2
				echo "-> Thresholds directory set to: ${2}"; shift 2;;
			-ol)
				local offset_left=$2
				echo "-> Offset left set to: ${2}"; shift 2;;
			-or)
				local offset_right=$2
				echo "-> Offset right set to"; shift 2;;
			-mins)
				local min_spacing=$2;
				echo "-> Minimum length of spacing: ${2}"; shift 2;;
			-maxs)
				local max_spacing=$2;
				echo "-> Maximum length of spacing: ${2}"; shift 2;;
			-sym)
				local matrix_type="SYMMETRIC"
				echo "-> Symmetric matrix mode activated"; shift 1;;
			-d)
				local divide_in=$2
				echo "-> Dividing cloud of peaks in ${2} groups"; shift 2;;
			-mim)
				local minmotiflength=${2};
				echo "-> Minimum motif length: ${2}"; shift 2;;
			-mam)
				local maxmotiflength=${2};
				echo "-> Maximum motif length: ${2}"; shift 2;;
			-l)
				local learningsetSize=${2};
				echo "-> Size of learning set for motif analysis: ${2}"; shift 2;;
			-nm)
				local numberofmotifs=${2};
				echo "-> Number of motifs to find: ${2}"; shift 2;;
			-s)
				local seed=${2};
				echo "-> Random seed for motif analysis: ${2}"; shift 2;;
			-h)
				usage analyze_decile; return;;
			--help)
				usage analyze_decile; return;;
			*)
				echo "Error in arguments"
				echo $1; usage analyze_decile; return;;
		esac
	done
	local Errors=0
	if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
	if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
	if [ -z $dir_comp ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
	if (( $(echo "${thre[0]} == 0.0" |bc -l) )) && [ ${thresholds_dir} == "null" ]; then echo "ERROR either -th or -thd arguments needed; -th arguments needs to be a list"; Errors+=1; fi
	if [ ${thresholds_dir} != "null" ] && [ ! -d ${thresholds_dir} ]; then echo "ERROR directory for thresholds acquisition does not exist, please check your argument -thd"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage analyze_decile; return 1; fi

	# Get RiL_RiP table
	mkdir -p -m 774 ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/
	if [ ! -f "${dir_comp}/table_${name1}_${name2}_RiL_RiP.tsv" ]; then
		echo "[ERROR] - File ${dir_comp}/table_${name1}_${name2}_RiL_RiP.tsv not found, please run pairwize_comparison first"
		exit 1
	else
		echo "[INFO] - Copying file ${dir_comp}/table_${name1}_${name2}_RiL_RiP.tsv to ${dir_comp}/${name1}_${name2}/"
		cp "${dir_comp}/table_${name1}_${name2}_RiL_RiP.tsv" "${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv"
	fi

	local log=${dir_comp}/${name1}_${name2}/${name1}vs${name2}/log.txt
	echo $Argsline > $log

	# Recover thresholds if not provided
	if (( $(echo "${thre[0]} == 0.0" |bc -l) )); then
		echo "[WARNING] - Thresholds not provided, trying to get them from directory ${thresholds_dir}" | tee -a $log
		local name=${name1}
		thre=() 
		for i in 50 70 90; do 
			file="scores_99_${i}_info.txt"
			if [ -f ${thresholds_dir}/${file} ]; then 
				thre+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then 
				thre+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
			if [ -f ${thresholds_dir}/${name}/${file} ]; then 
				thre+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then 
				thre+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
		done
		if (( $(echo "${thre[0]} == 0.0" |bc -l) )); then
			echo "[ERROR] - Thresholds recovery from ${thresholds_dir} unsuccessful" | tee -a $log; exit 1
		fi
	fi
	
	# Dividing peaks in deciles
	local total=$(wc -l ${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv | awk '{print $1-1}')
	local total=$(wc -l ${dir_comp}/table_${name1}_${name2}_RiL_RiP.tsv | awk '{print $1-1}')
	local decile=$(calc $total/$divide_in | awk '{print int($1+0.5)}')
	echo "$total $decile" | tee -a $log

	if [[ $typeofnorm == *"RiP"* ]] || [[ $typeofnorm == *"RIP"* ]] || [[ $typeofnorm == *"rip"* ]]; then
		# CFR from RiP normalisation
		awk 'NR!=1&&$7!=0&&$8!=0{print $1,$2,$3,$7,$8,$7/$8,$4;next}NR!=1&&$7!=0{print $1,$2,$3,$7,$8,$7/0.001,$4;next}NR!=1&&$8!=0{print $1,$2,$3,$7,$8,0.001,$4/$8}' ${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv | sort -k6,6n | awk -v decile=$decile -v OFS="\t" '{printf "%s\t%d\t%d\t%f\t%f\t%f\t%d\t%s\n", $1,$2,$3,$4,$5,$6,((NR-1)/decile)+1,$7}' | awk -v maxidiv=$divide_in -v OFS="\t" '$7==maxidiv+1{print $1,$2,$3,$4,$5,$6,maxidiv,$8; next}{print $0}' > ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv
	else
		# CFR from RiL normalisation
		awk 'NR!=1&&$9!=0&&$10!=0{print $1,$2,$3,$9,$10,$9/$10,$4;next}NR!=1&&$9!=0{print $1,$2,$3,$9,$10,$9/0.001,$4;next}NR!=1&&$10!=0{print $1,$2,$3,$9,$10,0.001/$10,$4}' ${dir_comp}/${name1}_${name2}/table_${name1}_${name2}_RiL_RiP.tsv | sort -k6,6n | awk -v decile=$decile -v OFS="\t" '{printf "%s\t%d\t%d\t%f\t%f\t%f\t%d\t%s\n", $1,$2,$3,$4,$5,$6,((NR-1)/decile)+1,$7}' | awk -v maxidiv=$divide_in -v OFS="\t" '$7==maxidiv+1{print $1,$2,$3,$4,$5,$6,maxidiv,$8; next}{print $0}' | sort -k1,1 -k2,2n | uniq > ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv
		
	fi
	if [ $matrix_type == "SYMMETRIC" ]; then
		echo -e "Spacing\tDecile\tAR" > ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv
	else
		echo -e "Spacing\tDecile\tER\tDR\tIR" > ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv
	fi

	# Analyzing each decile
	for ((l=1;l<=$divide_in;l++)); do
		echo "[INFO] - Analyzing decile $l" | tee -a $log
		awk -v OFS="\t" -v th=$l '$7==th{print $0}' ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv > ${dir_comp}/${name1}_${name2}/tmp_peaks.tsv
		echo "[INFO] - Computing motifs for decile $l" | tee -a $log
		compute_motif -p ${dir_comp}/${name1}_${name2}/tmp_peaks.tsv -n decile${l} -g $genome -od ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Motif -ls $learningsetSize -s $seed -nm $numberofmotifs -mim $minmotiflength -mam $maxmotiflength -c 4
		
		list_peaks=("${dir_comp}/${name1}_${name2}/tmp_peaks.tsv")
		list_name=("decile${l}")
		echo "[INFO] - Computing spacing for decile $l" | tee -a $log
		if [ $matrix_type == "SYMMETRIC" ]; then
			compute_space -p list_peaks[@] -m list_matrices[@] -n list_name[@] -od ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Spacing -maxs $max_spacing -mins $min_spacing -ol ${offset_left} -or ${offset_right} -g ${genome} -th thre[@] -sym
			awk -v OFS="\t" -v decile=$l 'NR!=1{print $1,decile,$3}' ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Spacing/decile${l}/Zscore_stats_F2.tsv >> ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv
		else
			compute_space -p list_peaks[@] -m list_matrices[@] -n list_name[@] -od ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Spacing -maxs $max_spacing -mins $min_spacing -ol ${offset_left} -or ${offset_right} -g ${genome} -th thre[@]
			awk -v OFS="\t" -v decile=$l 'NR!=1{print $1,decile,$3,$11,$7}' ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Spacing/decile${l}/Zscore_stats_F2.tsv >> ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv
		fi
		
	done
	echo "[INFO] - Making heatmap" | tee -a $log
	if [ $matrix_type == "SYMMETRIC" ]; then
		$PATHRscript $HeatmapSymmetric ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv ${dir_comp}/${name1}_${name2}/${name1}vs${name2} $name1 $name2 ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv CFR $divide_in
	else
		$PATHRscript $HeamapClassic ${dir_comp}/${name1}_${name2}/${name1}vs${name2}/Recap_F2.tsv ${dir_comp}/${name1}_${name2}/${name1}vs${name2} $name1 $name2 ${dir_comp}/${name1}_${name2}/tmp_table_${name1}_${name2}.tsv CFR $divide_in
	fi
}

#-------------------------------------------------------------------------------
compute_rpkmrip_rpkmril(){ 
	# FUNCTION: compute_rpkmrip_rpkmril
	# DESCRIPTION:
	#   This function computes the coverage of reads per peak using bedtools and normalizes the results
	#   according to the specified mode ("inPeaks" or "inLibs"). It processes BAM files for each sample,
	#   calculates the Fraction of Reads In Peaks (FRIP), and generates output files for further analysis.
	#
	# USAGE:
	#   compute_rpkmrip_rpkmril -p <PATH> -bd <LIST of PATH> -pd <LIST of PATH>
	#                            -sn <LIST of STRING> -m <STRING> -o <PATH> -g [STRING]
	#
	# ARGUMENTS:
	#   -p         : Path to the file containing consensus peaks. Required.
	#   -bd        : List of directories containing BAM files for each sample. Required.
	#   -pd        : List of directories containing peak calling results for each sample. Required.
	#   -sn        : List of sample names to process. Required.
	#   -m         : Normalization mode: "inPeaks" or "inLibs". Required.
	#   -o         : Directory where output files will be saved. Required.
	#   -g         : Whether to compute coverage ("yes" or "no"). Optional. Default: "yes".
	#   -h, --help : Display usage information for this function.
	#
	# NOTES:
	#   - The function assumes the existence of the following tools: `bedtools`, `samtools`
	
	echo "
	____ ___  _  _ _  _ ____ _ ___      ____ ___  _  _ _  _ ____ _ _   
	|__/ |__] |_/  |\/| |__/ | |__]     |__/ |__] |_/  |\/| |__/ | |   
	|  \ |    | \_ |  | |  \ | |    ___ |  \ |    | \_ |  | |  \ | |___
		"
	local getCov="yes";
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> File containing peaks set to ${2}";shift 2;;
			-bd)
				local bam_dir=("${!2}")
				echo "-> Path to mapping directory set to ${bam_dir[@]}";shift 2;;
			-pd)
				local peaks_dir=("${!2}")
				echo "-> Path to peaks calling directory set to ${peaks_dir[@]}";shift 2;;
			-sn)
				local samples_names=("${!2}")
				echo "-> List of sample's name set to: ${samples_names[@]}";shift 2;;
			-m)
				local mode=$2
				echo "-> Normalization mode set to 'inPeaks' or 'inLibs' (number of reads retained by MACS3) ${2}";shift 2;;
			-o)
				local out_dir=$2
				echo "-> Out directory set to : ${2}";shift 2;;
			-g)
				local getCov=$2
				echo "-> Computation of coverage set to (yes/no): ${2}";shift 2;;
			-h)
				usage compute_rpkmrip_rpkmril; return;;
			--help)
				usage compute_rpkmrip_rpkmril; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				echo $1; usage compute_rpkmrip_rpkmril; return;;
		esac	
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
	if [ -z $bam_dir ]; then echo "ERROR: -bd argument needed"; Errors+=1; fi
	if [ -z $peaks_dir ]; then echo "ERROR: -pd argument needed"; Errors+=1; fi
	if [ -z $samples_names ]; then echo "ERROR: -sn argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
	if [ -z $mode ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
	if [ -z $getCov ]; then echo "reads count will be computed using bedtools"; fi
	if [ $Errors -gt 0 ]; then usage compute_rpkmrip; return 1; fi	

	mkdir -p -m 774 $out_dir
	local log=${out_dir}/log.txt
	echo $Argline > $log
	touch $out_dir/tmpTotalTags.txt

	# Format peaks file for bedtools
	awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3}NR>1{print $1,$2,$3}' $peaks | sort -k1,1 -k2,2n > $out_dir/tmp_peaks
	printf "pr -mts <(cut -f -3 $out_dir/tmp_peaks) " > $out_dir/tmp_run.sh

	# Compute coverage per peak for each sample
	if [ getCov=="yes" ]; then 
		## if tmpFiltTags file is already present, remove it
		if [[ -f $out_dir/tmpTotalTags.txt ]]; then
			rm $out_dir/tmpTotalTags.txt
			touch $out_dir/tmpTotalTags.txt
		fi
		if [[ -f $out_dir/tmpFiltTags.txt ]]; then
			rm $out_dir/tmpFiltTags.txt
			touch $out_dir/tmpFiltTags.txt
		fi
		# same thing for FRIP file
		if [[ -f ${out_dir}/RIP.txt ]]; then
			rm ${out_dir}/RIP.txt
		fi
	fi

	# Compute coverage per peak for each sample
	printf "\n[INFO] - Computing reads coverage per peak using bedtools coverage\n" | tee -a $log
	for SAMP in  ${samples_names[@]}
	do
		echo "[INFO] - Treating sample: $SAMP" | tee -a $log
		printf "<(cut -f 4 ${out_dir}/${SAMP}_filt.cov.bed) " >> $out_dir/tmp_run.sh
		
		bam_file="NA"

		# Find bam file in provided directories
		for bams in ${bam_dir[@]}; do
			if [[ -f "${bams}/$SAMP/$SAMP.filtered.sorted.dedup.bam" ]]; then
				local bam_file=${bams}/$SAMP/$SAMP.filtered.sorted.dedup.bam
				break
			elif [[ -f "${bams}/$SAMP/$SAMP.filtered.sorted.bam" ]]; then
				local bam_file=${bams}/$SAMP/$SAMP.filtered.sorted.bam
				break
			fi
			if [[ -f "$bam_file" ]]; then
				break
			fi
			
		done
		if [ "$bam_file" == "NA" ] || [ ! -f "$bam_file" ]; then
			echo "[ERROR] - BAM file for sample $SAMP not found, please check your -bd argument"; return 1
		else
			echo "[INFO] - BAM file for sample $SAMP found: $bam_file" | tee -a $log
		fi

		if [ "$getCov" == "yes" ]; then
			echo "[INFO] - Computing coverage for sample $SAMP" | tee -a $log
			bedtools coverage -a $out_dir/tmp_peaks -b $bam_file > ${out_dir}/${SAMP}_filt.cov.bed #compute coverage
			
			# Compute FRIP
			inpeaks=$(bedtools sort -i $out_dir/tmp_peaks | bedtools merge -i stdin | bedtools intersect -u -a $bam_file -b stdin -ubam | $PATH_TO_SAMTOOLS/samtools view -c)
			total=$($PATH_TO_SAMTOOLS/samtools view -c $bam_file)
			FRIP=$(calc $inpeaks/$total*100)
			echo "[REPORT] - ${FRIP}% of reads in peaks (${inpeaks}/${total})" | tee -a $log
			echo "${inpeaks}" >> ${out_dir}/RIP.txt
			
			#format filtered tags (for in libs normalization)
			local found=0 #init
			# Look for the stats file in the provided peaks calling directories
			for Peaksdir in ${peaks_dir[@]}
			do
				if [ -f $peaks_dir/$SAMP/${SAMP}_stats.txt ]; then
					# Extract totalTags and filtTags from stats file
					cut -d " " -f 2 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v totalTags >> $out_dir/tmpTotalTags.txt
					cut -d " " -f 3 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v filtTags >> $out_dir/tmpFiltTags.txt
					local found=1
					break
				fi
			done
			
			if [ $found -eq 0 ] && [[ $bam_file == *"control"* ]]; then
				echo $(($total)) >> $out_dir/tmpTotalTags.txt
				echo $(($total)) >>  $out_dir/tmpFiltTags.txt
			fi
		else
			printf "\n[INFO] - Skipping bedtools coverage because reads coverage has already been computed\n" | tee -a $log
		fi
	done
	
	# Finalize command to get all reps reads count
	printf "> ${out_dir}/allReps_RC.txt" >> $out_dir/tmp_run.sh
	bash $out_dir/tmp_run.sh
	rm $out_dir/tmp_run.sh
	#cp ${out_dir}/allReps_RC.txt
	echo "[INFO] - R script to normalize reads count in peaks or in library/mapped started" | tee -a $log
	echo '~###################################' | tee -a $log
	printf "\n"
	$PATHRscript $rpkmrip_rpkmril ${out_dir}/allReps_RC.txt $mode $out_dir ${samples_names[@]}
	printf "\n[INFO] - End of rpkmrip_rpkmril\n"  | tee -a $log
}

#-------------------------------------------------------------------------------
compute_motif(){
	# FUNCTION: compute_motif
	# DESCRIPTION:
	#   This function performs motif discovery and analysis using MEME and TFFM tools.
	#   It prepares sequence datasets, runs MEME for motif discovery, and generates
	#   Position Frequency Matrices (PFM) and TFFM models. It supports both palindromic
	#   and non-palindromic motifs, and allows for customization of motif length,
	#   learning set size, and other parameters. 
	#
	# USAGE:
	#   compute_motif -p <PATH> -n <STRING> -od <PATH>
	#                 -c [INT] -g [PATH] -ls [INT] -nm [INT] [-t] -mim [INT] -mam [INT]
	#                  [-pal] -stffm [INT] -s [INT] -top [INT]
	#                  [-neg <neg_file>] [-cutoff <fracCutoff>]
	#
	# ARGUMENTS:
	#   -p         : Path to the BED file containing peaks. Required.
	#   -n         : Name for the sub-directory and output files. Required.
	#   -od        : Directory where results will be saved. Required.
	#   -c         : Column number in the input BED file that gives sequence coverage (RPKM). Optional. Default: 0.
	#   -g         : Path to the genome FASTA file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas".
	#   -ls        : Size of the learning set. Optional. Default: 600.
	#   -nm        : Number of PWM motifs to generate. Optional. Default: 1.
	#   -t         : Transpose the PWM motif if required. Optional.
	#   -mim       : Minimum length for motifs. Optional. Default: 6.
	#   -mam       : Maximum length for motifs. Optional. Default: 15.
	#   -pal       : Activate palindromic mode for motif discovery. Optional.
	#   -stffm     : PFM index used for TFFM computation. Optional. Default: 0.
	#   -s         : Seed for random number generation. Optional. Default: 1254.
	#   -top       : Maximum number of peaks to consider for the testing set. Optional. Default: 0.
	#   -neg       : Path to the control file (negative set). Optional. Default: "NA".
	#   -cutoff    : Fraction of sequences having at least one motif to consider a motif as valid. Optional. Default: 0.333.
	#   -h, --help : Display usage information for this function.
    # 
	# DEPENDENCIES:
	#   - bedtools
	#   - MEME
	#   - TFFM 
	#
	# NOTES:
	#   - If the learning set size is too large compared to the input BED file, the entire file is used for both learning and testing.
	#   - The function supports both palindromic and non-palindromic motifs.

	echo "
	____ ____ _  _ ___  _  _ ___ ____     _  _ ____ ___ _ ____
	|    |  | |\/| |__] |  |  |  |___     |\/| |  |  |  | |___
	|___ |__| |  | |    |__|  |  |___ ___ |  | |__|  |  | |   
	
	"
	local pal=false; local motifstffm=0; local learning_size=600; local min_motif=6; local max_motif=15; local seed=1254; local nb_motifs=1; local col_coverage=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local neg_file="NA"; local selection_tffm=0; local transpose=false; local topPeaks=0; local fracCutoff=0.333;
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Peak file set to: ${2}"; shift 2;;
			-c)
				local col_coverage=$2
				echo "-> Column number, in the input bed file (e.g. output by initial_comp), that gives sequences coverage (RPKM): ${2}"; shift 2;;
			-n)
				local name=$2
				echo "-> Name for sub-directory set to: ${2}"; shift 2;;
			-g)
				local genome=$2
				echo "-> Genome to use set to: ${2}"; shift 2;;
			-od)
				local results_motif=$2
				echo "-> Result directory set to: ${2}"; shift 2;;
			-ls)
				local learning_size=${2%.*}
				echo "-> Size for learning set: ${2}"; shift 2;;
			-nm)
				local nb_motifs=${2%.*}
				echo "-> Number of PWM motif generated: ${2}"; shift 2;;
			-t)
				local transpose=true
				echo "-> PWM motif used will be transposed as required"; shift 1;;
			-mim)
				local min_motif=$2
				echo "-> Minimum length for motif set to: ${2}"; shift 2;;
			-mam)
				local max_motif=$2
				echo "-> Maximum length for motif set to: ${2}"; shift 2;;
			-pal)
				local pal=true
				echo "-> Palindromic mode activated";shift 1;;
			-stffm)
				local selection_tffm=$2
				echo "-> PFM used for tffm computation set to: ${2}"; shift 2;;
			-s)
				local seed=$2
				echo "-> Seed for random set to: ${2}"; shift 2;;
			-top)
				local topPeaks=${2%.*}
				echo "-> Maximum number of peaks to consider: $2"; shift 2;;
			-neg)
				local neg_file=$2
				echo "-> Control file set to: ${2}"; shift 2;;
			-cutoff)
				local fracCutoff=$2;
				echo "-> Fraction of sequences having at least one motif to consider a motif as valid: ${2}"; shift 2;;
			-h)
				usage compute_motif; return;;
			--help)
				usage compute_motif; return;;
			*)
				echo "Error in arguments"
				echo $1; usage compute_motif; return;;
		esac
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; 
	elif [ ! -f $peaks ]; then echo "ERROR: peak file $peaks not found"; Errors+=1;	fi
	if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
	if [ -z $results_motif ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage compute_motif; return 1; fi

	# Check motif sizes
	if [ $max_motif -lt $min_motif ]; then
		echo "[WARNING] - Things got mixed up, max size of motif can't be lower than min size."
		echo -e "\t[INFO] - Using min size as max size "
		local max_motif=$(calc $max_motif+$min_motif) # max motif = total
		local min_motif=$(calc $max_motif-$min_motif) # new min_motif = old max_motif (total - old min = old max)
		local max_motif=$(calc $max_motif-$min_motif) # new max_motif = old min_motif (total - new min = old min)
	fi 

	local lowlimit=$(calc $learning_size/2)
	mkdir -p -m 774 $results_motif/$name/sets $results_motif/$name/meme $results_motif/$name/tffm
	local log=$results_motif/$name/log.txt
	echo $Argline > $log

	# Prepare learning and testing sets
	if [ $(wc -l $peaks | awk '{print $1}') -le $(calc $learning_size+$lowlimit) ]; then
		echo "[WARNING] - Not enough sequences for proper training/testing, using whole set for testing!" | tee -a $log
		awk -v OFS='\t' '{print $1, $2, $3}' $peaks > $results_motif/$name/sets/${name}_testingset.bed
		awk -v OFS='\t' '{print $1, $2, $3}' $peaks | head -n $learning_size > $results_motif/$name/sets/${name}.bed
		if [ $neg_file != "NA" ]; then
			cat $neg_file > $results_motif/$name/sets/tmpneg${name}.bed
		fi
	else

		if [ $col_coverage != 0 ]; then
			echo "[INFO] - BED file will be sorted according to column $col_coverage" | tee -a $log
			# Check if there are any non-empty values in the $col_coverage column
			# awk -v col="$col_coverage" 'NR>1 && $col != "" && $col != 0' "$peaks" | head
			# echo "="
			# nbcolFile=$(head -2 $peaks | tail -1 | tr -cd '\t' | wc -c | awk '{print $1+1}')
			# echo "="
			# if [[ $nbcolFile -ge $col_coverage ]] ; then
			# 	# head $peaks
			# 	# awk 
			# 	echo "[ERROR] - You asked to sort peaks according to coverage in column $col_coverage but no values are found in this column. Exiting."
			# 	exit 1
			# fi

			# Select top learning_size peaks for learning set
			cat $peaks | cut -f 1,2,3,$col_coverage \
			| sort -k4,4rn > $results_motif/$name/sets/tmp${name}.bed
			head -n $learning_size $results_motif/$name/sets/tmp${name}.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
			cat $peaks | cut -f 1,2,3,$col_coverage \
			| sort -k4,4rn | sed "1,${learning_size}d" > $results_motif/$name/sets/tmp${name}.bed
			
			# Select top peaks for testing set
			if [ $topPeaks != 0 ]; then
				echo "[INFO] - using only top $topPeaks peaks for testing set" | tee -a $log
				head -n $topPeaks $results_motif/$name/sets/tmp${name}.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			else
				# Use all remaining peaks for testing set
				cat $results_motif/$name/sets/tmp${name}.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			fi

			# Prepare negative set if provided
			if [ $neg_file != "NA" ]; then
				cat $neg_file | head -n $learning_size > $results_motif/$name/sets/tmpneg${name}.bed
			fi
			# clean tmp file
			rm $results_motif/$name/sets/tmp${name}.bed


		# If no coverage column is specified, randomly shuffle the BED file
		else
			echo "[INFO] - BED file will be randomly shuffled to select learning set" | tee -a $log
			shuf $peaks | head -n $learning_size | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}.bed
			if [ $topPeaks != 0 ]; then
				sed "1,${learning_size}d" $peaks | head -n $topPeaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			else
				sed "1,${learning_size}d" $peaks | awk -v OFS='\t' '{print $1, $2, $3}' > $results_motif/$name/sets/${name}_testingset.bed
			fi
			
		fi
	fi

	# Extract sequences for learning set
	awk -v OFS="\t" '$2>=1{print $1,$2,$3;next}{print $1,"1",$3}' $results_motif/$name/sets/${name}.bed > $results_motif/$name/sets/${name}_resized.bed
	# get sequences
	bedtools getfasta -fi $genome -fo $results_motif/$name/sets/${name}.fas -bed $results_motif/$name/sets/${name}_resized.bed

	# Run MEME for motif discovery
	echo "[INFO] - Starting MEME analysis" | tee -a $log
	if [ $neg_file == "NA" ]; then
		if $pal ; then
			$meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-pal -meme-maxsize $(calc $learning_size*2000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed -fimo-skip -ccut 0 -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR/JASPAR2024_CORE_plants_non-redundant_pfms_meme.meme 
			local type_pal="sym"
		else
			$meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-maxsize $(calc $learning_size*2000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed -fimo-skip -ccut 0 -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR/JASPAR2024_CORE_plants_non-redundant_pfms_meme.meme 
			local type_pal="asym"
		fi
	else
		bedtools getfasta -fo $results_motif/$name/sets/tmpneg${name}.fa -fi $genome -bed $results_motif/$name/sets/tmpneg${name}.bed
		if $pal ; then
			$meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-pal -meme-maxsize $(calc $learning_size*2000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed -fimo-skip -ccut 0 -neg $results_motif/$name/sets/tmpneg${name}.fa -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR/JASPAR2024_CORE_plants_non-redundant_pfms_meme.meme
			local type_pal="sym"
		else
			$meme_prog -oc $results_motif/$name/meme -nmeme $learning_size -meme-maxsize $(calc $learning_size*2000) -meme-minw $min_motif -meme-maxw $max_motif -meme-nmotifs $nb_motifs -dreme-m 0 -noecho $results_motif/$name/sets/${name}.fas -seed $seed -fimo-skip -ccut 0 -neg $results_motif/$name/sets/tmpneg${name}.fa -db /home/312.6-Flo_Re/312.6.1-Commun/data/meme_db/motif_databases/JASPAR/JASPAR2024_CORE_plants_non-redundant_pfms_meme.meme
			local type_pal="asym"
		fi
	fi

	# Get MEME mini format 
	$meme2meme $results_motif/$name/meme/meme_out/meme.txt > $results_motif/$name/meme/meme_out/meme_mini.txt
	
	local path_to_meme_mini=$results_motif/$name/meme/meme_out/meme_mini.txt
	
	#Filter out motif (MEME-1 motif found in less than one third of the learning set)
	local nsitesM1=$(cat $path_to_meme_mini | grep "MEME-1" -A 2 | tail -n 1 | awk -F 'nsites= ' ' {print $2}' | awk '{print $1}')
	local fracM1=$(calc $nsitesM1/$learning_size)
	if awk "BEGIN {exit !($fracM1 < $fracCutoff)}"; then
		echo "[WARNING] - MEME-1 motif is found in less than one third of the learning set!" | tee -a $log
	fi
	
	# Convert to PFM
	bash $meme2pfm $path_to_meme_mini $name ${type_pal} #> $results_motif/$name/${name}.pfm

	# Generate TFFM model: preprataion 
	if $pal ; then
		echo "[INFO] - generate TFFM learning set" | tee -a $log
		python3 $prepMEMEforPalTFFM -m $results_motif/$name/meme/meme_out/meme.txt -f $results_motif/$name/sets/${name}_tffm_learningset.fas -p ${selection_tffm} # THIS SCRIPT NEEDS TO BE CORRECTED TODO ##JL: still doesn't know what is wrong with it ?
		local learning_set_tffm=$results_motif/$name/sets/${name}_tffm_learningset.fas
	# 		local learning_set_tffm=$results_motif/$name/sets/${name}.fas
	# 		cat $learning_set_tffm | grep ">" -v | tr "A" "@" | tr 'T' '@' | tr 'C' '@' | tr 'G' '@'  | sed 's/@//g' | tr '\n' 'U'
	else
		local learning_set_tffm=$results_motif/$name/sets/${name}.fas
	fi

	# Generate TFFM model: run TFFM
	echo "[INFO] - generate TFFM model from MEME output" | tee -a $log
	$Python_TFFM $pfmTOtffm -r $results_motif/$name/tffm/ -f $learning_set_tffm -m $results_motif/$name/meme/meme_out/meme.txt -p ${selection_tffm}
	# Copy TFFM model to final output
	local pfmtokeep=$(calc $selection_tffm+1) # keeping the pfm used for tffm for spacing and ROCs computation
	cp $results_motif/$name/tffm/tffm_first_order.xml $results_motif/$name/${name}_tffm.xml

	# Generate logo and reverse complement of motif
	echo "[INFO] - generate logo and reverse complement of motif" | tee -a $log
	Headermotif=$(head -1 $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm | awk '{print $4}' )
	head -1 $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm > $results_motif/$name/${name}_rc.pfm
	echo -e 'A\tC\tG\tT' >> $results_motif/$name/${name}_rc.pfm
	tail -n +3 $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm | tac | awk -v OFS="\t" '{print $4,$3,$2,$1}' >> $results_motif/$name/${name}_rc.pfm
	if awk "BEGIN {exit !($fracM1 > $fracCutoff)}"; then cp $results_motif/$name/meme/meme_out/logo_rc${pfmtokeep}.png $results_motif/$name/${name}_logorc.png; fi

	if awk "BEGIN {exit !($fracM1 > $fracCutoff)}"; then cp $results_motif/$name/meme/meme_out/logo${pfmtokeep}.png $results_motif/$name/${name}_logo.png; fi
	cp $results_motif/$name/meme/meme_out/Motif_MEME_seperateFiles/Motif_${pfmtokeep}.pfm $results_motif/$name/${name}.pfm
}

#-------------------------------------------------------------------------------
prep_annotation(){
	# FUNCTION: prep_annotation
	# DESCRIPTION:
	#   This function prepares genome annotation files by extracting and formatting
	#   promoter, intronic, intergenic, and exonic regions from a GFF3 file.
	#   It generates a comprehensive annotation file for downstream analysis,
	#   such as peak annotation or genomic feature enrichment.
	#   The function also supports generating a FAI (Fasta Index) file if not provided.
	#
	# USAGE:
	# 	prep_annotation -n <FILE> -p <INT> -g <FILE> -s <FILE> -o <PATH>
	#
	# ARGMENTS:
	#   -n          : Path to the GFF3 annotation file of the genome. Required.
	#   -o          : Path to the output annotation file. Required.
	#   -p          : Length of promoter regions to extract (in bp). Optional. Default: 1000.
	#   -g          : Path to the genome FASTA file. Optional. Required if -s is not provided.
	#   -s          : Path to the genome size file (FAI format). Optional. Required if -g is not provided.
	#   -h, --help  : Display usage information for this function.
	#
	# NOTES:
	#   - The function requires `bedtools` and `samtools` to be installed and available in the system PATH.
	#   - If `-s` is not provided but `-g` is, the function will generate a FAI file from the genome FASTA.
	#   - The function filters the GFF3 file to retain only primary RNA features (e.g., genes, mRNA).
	#   - Promoter regions are extracted upstream of genes (adjustable with `-p`).
	#   - Intronic and intergenic regions are computed using `bedtools complement`.
	#   - The output file is sorted and formatted for downstream analysis.

	local promoterLength=1000; local Size="NA"; local Genome="NA"
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-n)
				local gff=$2
				echo "-> GFF file of the genome set to: ${2}";shift 2;;
			-p)
				local promoterLength=$2
				echo "-> Size for the promoters set to: ${2}";shift 2;;
			-g)
				local Genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-s)
				local Size=$2
				echo "-> Size file of the genome (fai) set to: ${2}";shift 2;;
			-o)
				local out_file=$2
				echo "-> Output file & directory set to: ${2}";shift 2;;
			-h)
				usage prep_annotation; return;;
			--help)
				usage prep_annotation; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage prep_annotation; return;;
		esac
	done

	local Errors=0
	if [ -z $gff ]; then echo "ERROR: -n argument needed";Errors+=1;fi
	if [ -z $out_file ]; then echo "ERROR: -o argument needed";Errors+=1;fi
	if [ $Genome == "NA" ] && [ $Size == "NA" ]; then 
		echo "ERROR: -g or -s argument needed";Errors+=1;
	else
		if [ $Size == "NA" ]; then echo "[WARNING] - No size file (fai) given (-s parameter); if no fai file in the very \n[WARNING] - same directory as the used genome file; a fai file will be generated."; local Size=${Genome}.fai ;fi
	fi
	if [ $Errors -gt 0 ]; then usage prep_annotation; return 1; fi

	local out_dir=$(dirname $out_file)
	mkdir -p -m 774 $out_dir
	local log=$out_dir/log.txt
	echo $Argline > $log

	# this awk should do the work to adapt to 'most' gff3, tested with A.thaliana
	awk -v OFS="\t" '{split($9,a,";"); split(a[2],b,"."); if(b[2]~"1"){print $0}else{if($3=="gene"){print $0}}}' $gff | uniq > $out_dir/tmpPrimaryRNA.gff3

	# Generate chrom sizes file
	if [ $Size == ${Genome}.fai ]; then
		echo "[INFO] - Generating fai file from genome fasta" | tee -a $log
		$PATH_TO_SAMTOOLS/samtools faidx $Genome
		awk 'OFS="\t" {print $1, $2}' ${Genome}.fai | sort -k1,1 -k2,2n > $out_dir/chromSizes.bed
	else
		awk 'OFS="\t" {print $1, $2}' $Size | sort -k1,1 -k2,2n > $out_dir/chromSizes.bed
	fi

	echo "[INFO] - Generating promoter regions" | tee -a $log
	# Promoter extraction
	if [ $promoterLength -gt 0 ]; then
		awk -v OFS="\t" -vFS="[=\t]" -v lengthprom=$promoterLength '$1 ~ /^#/ {next} $7=="+" && $3=="gene" {print $1,$2,$3,$4,$5,$6,$7,$8,$9"="$10"="$11; print $1,".","promoter",$4-lengthprom-1,$4-1,".","+",".","ID="$11".1.TAIR10.PROMOTER;Parent="$11".TAIR10" } $7=="-" && $3=="gene" {print $1,$2,$3,$4,$5,$6,$7,$8,$9"="$10"="$11; print $1,".","promoter",$5+1,$5+1+lengthprom,".","-",".","ID="$11".1.TAIR10.PROMOTER;Parent="$11".TAIR10" }' $out_dir/tmpPrimaryRNA.gff3 | awk -v OFS="\t" '$3 != "gene" {print $0;next}' | awk -v OFS="\t" '$4<0{print $1,".","promoter","1",$5,$6,$7,$8,$9;next}{print $0}' | awk -v OFS="\t" '$5>$4{print $0}' > $out_dir/promoters.gff3
		cat $out_dir/tmpPrimaryRNA.gff3 $out_dir/promoters.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > $out_dir/in_sorted.gff
		bedtools complement -i <(cat $out_dir/tmpPrimaryRNA.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}') -g $out_dir/chromSizes.bed > $out_dir/intergenicTMP.bed
	else
		# Option exist in case you don't want to use promoters in NSG. 
		cat $out_dir/tmpPrimaryRNA.gff3 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > $out_dir/in_sorted.gff
	fi

	echo "[INFO] - Generating intergenic regions" | tee -a $log
	bedtools complement -i $out_dir/in_sorted.gff -g $out_dir/chromSizes.bed > $out_dir/intergenic.bed

	echo "[INFO] - Generating intronic regions" | tee -a $log
	# Exon extraction
	awk -v OFS="\t" '$1 ~ /^#/ {next} $3!= "gene" && $3!="mRNA" && $3!="promoter" {print $1, $4-1, $5+1;next}' $out_dir/in_sorted.gff > $out_dir/exon.bed
	if [ $promoterLength -gt 0 ]; then
		bedtools complement -i <(cat $out_dir/exon.bed $out_dir/intergenicTMP.bed | sort -k1,1 -k2,2n) -g $out_dir/chromSizes.bed > $out_dir/intron.bed
	else
		bedtools complement -i <(cat $out_dir/exon.bed $out_dir/intergenic.bed | sort -k1,1 -k2,2n) -g $out_dir/chromSizes.bed > $out_dir/intron.bed
	fi

	echo "[INFO] - Generating final annotation file" | tee -a $log
	# Final annotation file
	cat <(awk -v OFS="\t" '{print $0,"intron"}' $out_dir/intron.bed) <(awk -v OFS="\t" '{print $0,"intergenic"}' $out_dir/intergenic.bed) <(awk -v OFS="\t" '$1 ~ /^#/ {next} $3!= "gene" && $3!="mRNA" && $3!="rRNA"{print $1, $4, $5,$3;next}' $out_dir/in_sorted.gff) | awk -v OFS="\t" '$2<0{print $1,"1",$3,$4;next}{print $0}' | awk -v OFS="\t" '$4=="CDS" {print $1,$2,$3,"exon";next}{print $0}' | awk -v OFS="\t" '$2!="0" && $4=="intergenic"{print $1,$2+1,$3,$4;next}{print $0}' | uniq | sort -k1,1 -k2,2n | awk -v OFS="\t" '$3>$2{print $0}' > $out_file
	
	echo "[INFO] - Cleaning temporary files" | tee -a $log
	[ -f $out_dir/tmpPrimaryRNA.gff3 ] && rm $out_dir/tmpPrimaryRNA.gff3
	[ -f $out_dir/promoters.gff3 ] && rm $out_dir/promoters.gff3
	[ -f $out_dir/in_sorted.gff ] && rm $out_dir/in_sorted.gff
	[ -f $out_dir/chromSizes.bed ] && rm $out_dir/chromSizes.bed
	[ -f $out_dir/intergenicTMP.bed ] && rm $out_dir/intergenicTMP.bed
	[ -f $out_dir/intergenic.bed ] && rm $out_dir/intergenic.bed
	[ -f $out_dir/exon.bed ] && rm $out_dir/exon.bed
	[ -f $out_dir/intron.bed ] && rm $out_dir/intron.bed
}

#-------------------------------------------------------------------------------
compute_NS(){
	# FUNCTION: compute_NS
	# DESCRIPTION:
	#   This function generates negative sets of genomic regions for comparative analysis.
	#   It uses a Python script (`negative_set_script`) to create negative sets that match
	#   the GC content and other characteristics of the input peak regions.
	#   The function supports both a standard mode and a simplified mode (`-st` or `--timeout`),
	#   which is faster but may be less precise.
	#
	# USAGE:
	# compute_NS -p <FILE> -n <STRING> -g <FILE> -od <PATH> -anf <FILE> 
	#        -nb [INT] -ws [INT] -gc [INT] -lt [INT] -s [INT] -mp [INT] -c [INT] 
	#        [--timeout] [-st]
	#
	# ARGUMENTS:
	#   -p          : Path to the BED file containing peak regions. Required.
	#   -n          : Prefix for the output files. Required.
	#   -g          : Path to the genome FASTA file. Required.
	#   -od         : Directory where results will be saved. Required.
	#   -anf        : Path to the annotation file (BED format). Required.
	#   -nb         : Number of negative sets to generate. Optional. Default: 1.
	#   -ws         : Window size for negative set generation. Optional. Default: 250.
	#   -gc         : Maximum allowed difference in GC content between positive and negative sets. Optional. Default: 0.03.
	#   -lt         : Number of regions to generate per annotation type. Optional. Default: 1000.
	#   -mp         : Maximum number of peaks to retain for negative set generation. Optional. Default: 15000.
	#   -c          : Column number in the input BED file that gives sequence coverage (RPKM). Optional. Default: 0.
	#   --timeout   : Activate timeout mode. If the standard mode takes too long, it switches to the simplified mode. Optional.
	#   -st         : Activate simplified mode for faster execution. Optional.
	#   -s          : Random seed for reproducibility. Optional. Default: 1254.
	#   -h, --help  : Display usage information for this function.
	#
	# NOTES:
	#   - The function requires the Python script `negative_set_script` to be available and executable.
	#   - If the `-c` option is provided, the input BED file is sorted by the specified coverage column.
	#   - The function supports two modes:
	#     - **Standard mode**: More precise but potentially slower.
	#     - **Simplified mode (`-st` or `--timeout`)**: Faster but may be less precise.
	#   - If the standard mode times out (after 15 minutes per negative set), the function automatically switches to the simplified mode.
	#   - The function checks that the number of positive and negative sequences is the same in the output files.

	echo "

	____ ____ _  _ ___  _  _ ___ ____     _  _ ____
	|    |  | |\/| |__] |  |  |  |___     |\ | [__ 
	|___ |__| |  | |    |__|  |  |___ ___ | \| ___]

	"
	local window_size=250; local number_of_NS=1; local seed=1254; local deltaGC=0.03; local limit_type=1000; local maxPeaks=15000;local col_coverage=0;
	local timeout=false; local simple_type=false; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas";
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Peaks file set to: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> Prefix for results set to: ${2}";shift 2;;
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-od)
				local results=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-anf)
				local annotation_file=$2
				echo "-> Annotation file set to: ${2}";shift 2;;
			-nb)
				local number_of_NS=$2
				echo "-> Number of negative sets to create set to: ${2}";shift 2;;
			-ws)
				local window_size=$2
				echo "-> Window size set to: ${2}";shift 2;;
			-gc)
				local deltaGC=$2
				echo "-> Delta GC allowed set to: ${2}";shift 2;;
			-lt)
				local limit_type=$2
				echo "-> Number of region by origin set to: ${2}";shift 2;;
			-s)
				local seed=$2
				echo "-> Seed for random set to: ${2}";shift 2;;
			-mp)
				local maxPeaks=$2
				echo "-> Maximum number of peaks retained for NS generation: ${2}";shift 2;;
			-c)
				local col_coverage=$2
				echo "-> Column number, in the input bed file (e.g. output by initial_comp), that gives sequences coverage (RPKM): ${2}";shift 2;;
			--timeout)
				local timeout=true
				echo "-> Timeout option active: NSG will start as default, but if it takes too long (>15 min per negative set demanded) it will be launched again with simple_type option";shift 1;;
			-st)
				local simple_type=true
				echo "-> Simple_type option active: compute_NS will start with -st option";shift 1;;
			-h)
				usage compute_NS; return;;
			--help)
				usage compute_NS; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage compute_NS; return;;
		esac
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
	if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
	if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ -z $annotation_file ]; then echo "ERROR: -anf argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage compute_NS; return 1; fi
	if $timeout ; then echo "[INFO] - timeout mode activated"; fi
	if $simple_type ; then echo "[INFO] - simple_type mode activated"; fi


	mkdir -p -m 774 ${results}
	local log=$results/log.txt
	echo $Argline > $log

	# If coverage column specified, sort peaks by coverage and keep top maxPeaks peaks
	if [ $col_coverage != 0 ]; then
		echo "[INFO] - BED file will be sorted according to column $col_coverage" | tee -a $log
		cut -f 1,2,3,$col_coverage $peaks | sort -k4,4rn | awk -v top=$maxPeaks 'NR <= top {print}' > $results/tmp_peaks.bed 
	else
		
		# Romain 3/12/25 Exit if the number of peaks exceeds $maxPeaks because we did not sorted peaks so some very good peaks may end up in negative set.
		line_count=$(wc -l < "$peaks")
		if [ "$line_count" -gt "$maxPeaks" ]; then
			echo "Error: The number of lines in $peaks ($line_count) exceeds the maximum allowed ($maxPeaks). Because peaks are not sorted by coverage or q values, sequences corresponding to good peaks may end up in your negative set." | tee -a $log
			exit 1
		fi
		# If no coverage column specified, randomly select maxPeaks peaks
		cat $peaks | awk -v top=$maxPeaks 'NR <= top {print}' > $results/tmp_peaks.bed
		
	fi
	peaks=$results/tmp_peaks.bed
	
	# Run negative set generator
	if $timeout ; then # if --timeout option specified in compute_NS
		echo "[INFO] - Running NSG in default mode, with timeout in case of failure, as requested" | tee -a $log
		local exit_status=0
		local time="$(calc 15*60*$number_of_NS)s" # timeout at 15 minutes per set
		echo $time 
		timeout $time python3 $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size || local exit_status=124
		# Check that output files were created and contain same number of positive and negative sequences
		if [[ -f ${results}/${name}_pos.bed ]] && [[ -f ${results}/${name}_1_neg.bed ]]; then
			lenpos=$( wc -l ${results}/${name}_pos.bed | awk '{print $1}' )
			lenneg=$( wc -l ${results}/${name}_1_neg.bed | awk '{print $1}' )
		else
			echo "[WARNING] - NSG did not produce output files" | tee -a $log
			local lenpos=0
			local lenneg=1
		fi
		# If command timed out (exit status 124) or if output files not created or if different number of positive and negative sequences, run again with -st option
		if [[ $exit_status -eq 124 ]] || [[ $lenpos -ne $lenneg ]] ; then
			echo -e "[INFO] - Running NSG in default mode took too long/failed; running with --simpletype faster option instead" | tee -a $log
			python3 $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size -st
		else
			echo "[REPORT] - NSG successfully ran with standard method" | tee -a $log
		fi
	# Otherwise, run according to simple_type option
	elif $simple_type ; then # if -st option specified in compute_NS
		echo "[INFO] - Running NSG with simple_type as requested" | tee -a $log
		python3 $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size -st
	else # default mode
		echo "[INFO] - Running NSG in default mode" | tee -a $log
		python3 $negative_set_script -pos $peaks -of ${name} -od ${results} -fas  $genome -bed $annotation_file -n $number_of_NS -r $seed -GC $deltaGC -l $limit_type -bs $window_size
	fi

	# Final check that output files were created and contain same number of positive and negative sequences
	if [[ -f ${results}/${name}_pos.bed ]] && [[ -f ${results}/${name}_1_neg.bed ]]; then
			lenpos=$( wc -l ${results}/${name}_pos.bed | awk '{print $1}' )
			lenneg=$( wc -l ${results}/${name}_1_neg.bed | awk '{print $1}' )
			if [[ $lenpos -ne $lenneg ]]; then
				echo "[ERROR] - NSG did not produce same number of positive and negative sequences" | tee -a $log
				return 1
			else
				echo "[REPORT] - NSG successfully produced ${lenpos} positive and negative sequences" | tee -a $log
			fi
	else
		echo "[ERROR] - NSG did not produce output files" | tee -a $log
		return 1
	fi
	rm $results/tmp_peaks.bed
}

#-------------------------------------------------------------------------------
compute_DNAshape(){
	# FUNCTION: compute_DNAshape
	# DESCRIPTION:
	#   This function prepares training and testing datasets for DNA shape analysis.
	#   It splits input peak regions into training and testing sets, generates negative sets,
	#   and extracts FASTA sequences for both foreground and background regions.
	#   It then trains a DNA shape classifier using either a TFFM (Transcription Factor Flexible Model)
	#   or a PWM (Position Weight Matrix) and generates heatmaps of DNA shape features.
	#
	# USAGE:
	# compute_DNAshape -f <FILE> -m <FILE> -o <PATH> -ls <INT> -n <STRING> 
	#        -g <FILE> -anf <FILE> -e [INT] -c [INT] -ws [INT] -gc [FLOAT] -lt [INT]
	#        -top [INT] -s [INT] 
	#
	# Arguments:
	#   -f           : Path to the BED file containing foreground peak regions. Required.
	#   -m           : Path to the matrix file (PFM or TFFM format). Required.
	#   -o           : Directory where results will be saved. Required.
	#   -e           : Number of base pairs to extend peaks. Optional. Default: 0.
	#   -ls          : Number of sequences for the learning set. Optional. Default: 250.
	#   -n           : Prefix for output files. Optional.
	#   -c           : Column number in the input BED file that gives sequence coverage (RPKM). Optional. Default: 0.
	#   -g           : Path to the genome FASTA file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas".
	#   -anf         : Path to the annotation file (BED format). Optional.
	#   -ws          : Window size for negative set generation. Optional. Default: 250.
	#   -gc          : Maximum allowed difference in GC content between positive and negative sets. Optional. Default: 0.03.
	#   -lt          : Number of regions to generate per annotation type. Optional. Default: 1000.
	#   -top         : Maximum number of peaks to consider for the testing set. Optional. Default: 0.
	#   -s           : Random seed for reproducibility. Optional. Default: 1254.
	#   -h, --help   : Display usage information for this function.
	#
	# NOTES:
	#   - The function requires `bedtools`, `compute_NS`, and Python scripts (`$Python_TFFM`) for DNA shape analysis.
	#   - If the `-c` option is provided, the input BED file is sorted by the specified coverage column.
	#   - The function splits the input peaks into training and testing sets, either randomly or by coverage.
	#   - Negative sets are generated using the `compute_NS` function.
	#   - The function supports both TFFM and PWM matrices for training DNA shape classifiers.

	echo "
	___  _  _ ____     ____ _  _ ____ ___  ____
	|  \ |\ | |__|     [__  |__| |__| |__] |___
	|__/ | \| |  | ___ ___] |  | |  | |    |___
	"

	local window_size=250; local seed=1254; local deltaGC=0.03; local limit_type=1000;
	local timeout=false; local extend=0; local col_coverage=0; local topPeaks=0
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-f)
				local foreground=${2};
				echo "-> BED file set to: ${2}";shift 2;;
			-m)
				local matrix=${2};
				echo "-> Motif matrix (pfm or tffm) set to: ${2}";shift 2;;
			-e)
				local extend=${2};
				echo "-> Extend peaks by: ${2} bp";shift 2;;
			-o)
				local output=${2};
				echo "-> Output directory set to: ${2}";shift 2;;
			-ls)
				local learning_size=${2%.*};
				echo "-> Number of sequences for learning set: ${2}";shift 2;;
			-n)
				local name=${2};
				echo "-> Prefix for results set to: ${2}";shift 2;;
			-c)
				local col_coverage=${2};
				echo "-> Column number, in the input bed file (e.g. output by initial_comp), that gives sequences coverage: ${2}";shift 2;;
			-g)
				local genome=$2;
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-anf)
				local annotation_file=$2
				echo "-> Annotation file set to: ${2}";shift 2;;
			-ws)
				local window_size=$2;
				echo "-> Window size set to: ${2}";shift 2;;
			-gc)
				local deltaGC=$2;
				echo "-> Delta GC allowed set to: ${2}";shift 2;;
			-lt)
				local limit_type=$2
				echo "-> Number of region by origin set to: ${2}";shift 2;;
			-top)
				local topPeaks=${2%.*}
				echo "-> Maximum number of peaks to consider: $2";shift 2;;
			-s)
				local seed=$2
				echo "-> Seed for random set to: ${2}";shift 2;;
			--help)
				usage compute_DNAshape; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage compute_DNAshape; return;;
		esac
	done
	
	local Errors=0
	if [ -z $foreground ]; then echo "ERROR: no foreground/peak bed file specified, use -f"; Errors+=1; fi
	if [ -z $matrix ]; then echo "ERROR: no matrix given, use -m (only pfm or TFFM accepted)"; Errors+=1; fi
	if [ -z $output ]; then echo "ERROR: no output prefix (PATH+prefix) given, use -o"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage compute_DNAshape; return 1; fi
	
	local lowlimit=$(calc $learning_size/2)
	mkdir -p -m 774 $output/$name/testing $output/$name/training
	local log=$output/$name/log.txt
	echo $Argline > $log

	if [ $(wc -l $foreground | awk '{print $1}') -le $(calc $learning_size+$lowlimit) ]; then
		echo "[WARNING] - Not enough sequences for proper training/testing (need at least learning_size + learning_size/2 sequences)" | tee -a $log
		awk -v OFS='\t' '{print $1, $2, $3}' $foreground > $output/$name/testing/foreground_testing.bed		
		awk -v OFS='\t' '{print $1, $2, $3}' $foreground > $output/$name/training/foreground_training.bed
	else
		# Coverage column is specified
		if [ $col_coverage != 0 ]; then
			echo "[INFO] - Sorting peaks by coverage before splitting in training and testing sets" | tee -a $log
			if [ $topPeaks != 0 ]; then
				# Filter top peaks
				cat $foreground | cut -f 1,2,3,$col_coverage \
				| sort -k4,4rn > $output/$name/training/tmpOrdered.bed
			
				# Create training and testing sets
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				
				sed "1,${learning_size}d" $output/$name/training/tmpOrdered.bed > $output/$name/training/tmpOrdered2.bed
				head -n $topPeaks $output/$name/training/tmpOrdered2.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
				rm $output/$name/training/tmpOrdered2.bed
			else
				# No filtering of top peaks
				cat $foreground | cut -f 1,2,3,$col_coverage \
				| sort -k4,4rn > $output/$name/training/tmpOrdered.bed

				# Create training and testing sets
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				cat $foreground | cut -f 1,2,3,$col_coverage \
				| sort -k4,4rn | sed "1,${learning_size}d" | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
			fi
			rm $output/$name/training/tmpOrdered.bed
		
		# No coverage column specified
		else
			echo "[INFO] - Randomly splitting peaks in training and testing sets" | tee -a $log
			shuf $foreground > $output/$name/training/tmpOrdered.bed
			if [ $topPeaks != 0 ]; then
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				sed "1,${learning_size}d" $output/$name/training/tmpOrdered.bed | head -n $topPeaks | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
			else
				head -n $learning_size $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/training/foreground_training.bed
				sed "1,${learning_size}d" $output/$name/training/tmpOrdered.bed | awk -v OFS='\t' '{print $1, $2, $3}' > $output/$name/testing/foreground_testing.bed
			fi
			rm $output/$name/training/tmpOrdered.bed
		fi
	fi

	# Generate negative sets and sequences
	echo "[INFO] - Generating negative sets and sequences for training set" | tee -a $log
	compute_NS -p $output/$name/training/foreground_training.bed -n $name -g $genome -od $output/$name/training -anf $annotation_file -nb 1 -ws $window_size -gc $deltaGC -lt $limit_type -s $seed --timeout
	
	# Reformat bed files to add name column
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/training/${name}_pos.bed > $output/$name/training/foreground_training.bed
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/training/${name}_1_neg.bed > $output/$name/training/background_training.bed
	
	echo "[INFO] - Generating negative sets and sequences for testing set" | tee -a $log
	compute_NS -p $output/$name/testing/foreground_testing.bed -n $name -g $genome -od $output/$name/testing -anf $annotation_file -nb 1 -ws $window_size -gc $deltaGC -lt $limit_type -s $seed --timeout
	
	# Reformat bed files to add name column
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/testing/${name}_pos.bed > $output/$name/testing/foreground_testing.bed
	awk -v OFS='\t' '{print $1,$2,$3,$1":"$2"-"$3}' $output/$name/testing/${name}_1_neg.bed > $output/$name/testing/background_testing.bed
	
	# Extract FASTA sequences
	bedtools getfasta -fo $output/$name/testing/background_testing.fas -fi $genome -bed $output/$name/testing/background_testing.bed
	bedtools getfasta -fo $output/$name/training/background_training.fas -fi $genome -bed $output/$name/training/background_training.bed
	bedtools getfasta -fo $output/$name/testing/foreground_testing.fas -fi $genome -bed $output/$name/testing/foreground_testing.bed
	bedtools getfasta -fo $output/$name/training/foreground_training.fas -fi $genome -bed $output/$name/training/foreground_training.bed
	
	# -------- Specified matrix is TFFM
	if [[ $matrix == *".xml"* ]]; then # TFFM scores computation
		echo "[INFO] - Training a first order TFFM + DNA shape classifier.";
		$Python_TFFM $ComputeDNAshaped trainTFFM -T $matrix \
			-i $output/$name/training/foreground_training.fas -I $output/$name/training/foreground_training.bed \
			-b $output/$name/training/background_training.fas -B $output/$name/training/background_training.bed \
			-o $output/$name/${name}_fo_classifier -t first_order \
			-1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;
		echo "[INFO] - Generating DNA shape heatmap from TFFM classifier.";
		$Python_TFFM $HeatmapDNAshape -c $output/$name/${name}_fo_classifier.pkl -2
		inkscape -o $output/$name/${name}_ShapeheatmapTFFM.png -w 1250 -h 1000 $output/$name/${name}_fo_classifier.pkl.svg
		cp $matrix  $output/$name/usedTFFM_fo.xml
	fi

	# -------- Specified matrix is PFM
	if [[ $matrix == *".pfm"* ]]; then  # PWM scores
		echo "[INFO] - Matrix in pfm format detected" | tee -a $log
		echo "[INFO] - Converting pfm to JASPAR format" | tee -a $log
		awk 'NR==1{print ">MA0000.1 "$4}' $matrix > $output/$name/JASPAR.pfm
		awk 'NR!=1{print $0}NR==2{print "[","[","[","["}END{print "]","]","]","]"}' $matrix | awk '
		{ 
			for (i=1; i<=NF; i++)  {
				if($i ~ /^[0-9]+$/){
					if ($i ~ /^[0-9]{1}$/){a[NR,i] = "  "$i}
					if ($i ~ /^[0-9]{2}$/){a[NR,i] = " "$i}
					if ($i ~ /^[0-9]{3}$/){a[NR,i] = $i}
					if ($i ~ /^[0-9]{4}$/){a[NR,i] = $i} 
				} else{ if ($i == "[") {a[NR,i] = " "$i} else {a[NR,i] = $i}}
			}
		}
		NF>p { p = NF }
		END {
			for(j=1; j<=p; j++) {
				str=a[1,j]
				for(i=2; i<=NR; i++){
					str=str" "a[i,j];
				}
				print str
			}
		}' >> $output/$name/JASPAR.pfm
		echo "[INFO] - Training a PSSM + DNA shape classifier." | tee -a $log
		$Python_TFFM $ComputeDNAshaped trainPSSM -f $output/$name/JASPAR.pfm \
			-i $output/$name/training/foreground_training.fas -I $output/$name/training/foreground_training.bed \
			-b $output/$name/training/background_training.fas -B $output/$name/training/background_training.bed \
			-o $output/$name/${name}_PSSM_classifier \
			-1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n;
		echo "[INFO] - Generating DNA shape heatmap from PSSM classifier." | tee -a $log
		$Python_TFFM $HeatmapDNAshape -c $output/$name/${name}_PSSM_classifier.pkl -2
		inkscape -o $output/$name/${name}_ShapeheatmapPSSM.png -w 1250 -h 1000 $output/$name/${name}_PSSM_classifier.pkl.svg
	fi
}

#-------------------------------------------------------------------------------
compute_SeqConv(){
	# FUNCTION: compute_SeqConv
	# DESCRIPTION:
	#   This function prepares a dataset of peak sequences, trains a SeqConv model for motif discovery,
	#   and generates a motif logo from the model's predictions.
	#   It uses a Python-based pipeline to prepare the data, train the model, and predict motifs.
	#
	# USAGE:
	# compute_SeqConv -p <FILE> -n <STRING> -od <PATH> -ps [INT] -g [FILE] 
	#        -sf [FILE] -s [INT]
	#
	# ARGUMENTS:
	#   -p          : Path to the BED file containing peak coordinates. Required.
	#   -ps         : Size of peaks in base pairs. Required. Default: 200.
	#   -ns         : Path to the BED file containing negative set coordinates. Required.
	#   -n          : Prefix for output directory and files. Required.
	#   -g          : Path to the genome FASTA file. Required.
	#   -od         : Directory where results will be saved. Optional.
	#   -sf         : Path to the chromosome size file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize".
	#   -s          : Random seed for reproducibility. Optional. Default: 1254.
	#   -split      : Percentage of data to use for training set. Optional. Default: 80.
	#   -h, --help  : Display usage information for this function.
	# 
	# DEPENDENCIES
	#  - SeqConv_train.py 
	#  - SeqConv_motif.py
	#  - bedtools
	#  - plot_logo.R

	echo "
	____ ____ ____ ____ ____ _  _ _  _
	[__  |___ |  | |    |  | |\ | |  |
	___] |___ |_\| |___ |__| | \|  \/ 
	
	"
	local peak_size=200; local seed=1254; local genome=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas ; local sizefile=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize; local split=80; local window_size=250;
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Name of bed files with peaks coordinates: ${2}";shift 2;;
			-ps)
				local peak_size=$2
				echo "-> Peaks size (200bp by default): ${2}";shift 2;;
			-ns)
				local negative_set=$2
				echo "-> Name of bed file to use as negative input for training and testing: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> A name used as prefix for output directory and files: ${2}";shift 2;;
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-od)
				local outdir=$2
				echo "-> Name of out directory: ${2}";shift 2;;
			-sf)
				local sizefile=$2
				echo "-> Chromosome size file: ${2}";shift 2;;
			-split)
				local split=$2
				echo "-> Splitting data percentage for training set: ${2}";shift 2;;
			-s)
				seed=$2
				echo "-> Seed for random set to: ${seed}"; shift 2;;
			--help)
				usage compute_SeqConv; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage compute_SeqConv; return;;
		esac
	done

	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed";Errors+=1;fi
    if [ -z $peak_size ]; then echo "ERROR: -ps argument needed";Errors+=1;fi
	if [ -z $name ]; then echo "ERROR: -n argument needed";Errors+=1;fi
	if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
	if [ $Errors -gt 0 ]; then usage compute_SeqConv ; return 1; fi

	mkdir -p -m 774 $outdir/$name/train $outdir/$name/model $outdir/$name/loss $outdir/$name/roc $outdir/$name/motif $outdir/$name/predict
	local log=$outdir/$name/log.txt
	echo $Argline > $log

	rm -f "$outdir/${name}/train/${name}_train.txt"
	touch "$outdir/${name}/train/${name}_train.txt"

	local train_file=$outdir/${name}/train/${name}_train.txt

	# ----- prepare data for SeqConv
	echo "[INFO] - Splitting positive and negative peaks into train and test sets" | tee -a $log
	
	awk -v seed="$seed" 'BEGIN{srand(seed)} {print rand(), $0}' "$peaks" | sort -k1,1n | cut -d' ' -f2- > "$outdir/${name}/tmp_pos_shuf.bed"
	awk -v seed="$seed" 'BEGIN{srand(seed)} {print rand(), $0}' "$negative_set" | sort -k1,1n | cut -d' ' -f2- > "$outdir/${name}/tmp_neg_shuf.bed"

	total_pos=$(wc -l < $outdir/${name}/tmp_pos_shuf.bed)
	total_neg=$(wc -l < $outdir/${name}/tmp_neg_shuf.bed)
	if [[ $total_pos -ne $total_neg ]]; then
		echo "[ERROR] - Different number of positive and negative sequences detected: ${total_pos} positive vs ${total_neg} negative"
		# exit 0
	fi
	train_pos=$((total_pos * split / 100))
	train_neg=$((total_neg * split / 100))
	# train splitting
	head -n $train_pos $outdir/${name}/tmp_pos_shuf.bed > $outdir/${name}/train/pos_train.bed
	head -n $train_neg $outdir/${name}/tmp_neg_shuf.bed > $outdir/${name}/train/neg_train.bed
	# test splitting 
	tail -n +$((train_pos + 1)) $outdir/${name}/tmp_pos_shuf.bed > $outdir/${name}/predict/pos_test.bed
	tail -n +$((train_neg + 1)) $outdir/${name}/tmp_neg_shuf.bed > $outdir/${name}/predict/neg_test.bed

	echo "[INFO] - Extracting sequences and formatting train.txt" | tee -a $log
	# - positive input 
	seqtk subseq -l 0 "$genome" "$outdir/${name}/train/pos_train.bed" | grep -v ">" > "$outdir/${name}/train/pos_seq.txt"
	paste "$outdir/${name}/train/pos_seq.txt" "$outdir/${name}/train/pos_train.bed" | awk '{print "1\t"$1"\t"$2"\t"$3"\t"$4}' > "$outdir/${name}/train/${name}_train.txt"
	# - negative input 
	seqtk subseq -l 0 "$genome" "$outdir/${name}/train/neg_train.bed" | grep -v ">" > "$outdir/${name}/train/neg_seq.txt"
	paste "$outdir/${name}/train/neg_seq.txt" "$outdir/${name}/train/neg_train.bed" | awk '{print "0\t"$1"\t"$2"\t"$3"\t"$4}' >> "$outdir/${name}/train/${name}_train.txt"

	# Format needed for training (positive input paste negative input):
	# 1	sequence	chrom	start	end
	# ...
	# 0	sequence	chrom	start	end

	# exit 0

	#----- train model
	echo "[INFO] - Training SeqConv model" | tee -a $log
	python $SeqConv_train ${outdir}/${name}/train/${name}_train.txt $peak_size $name $outdir $seed
	echo "-> Done! Model file availible at: ${outdir}/${name}/model/${name}_model.keras"

	#----- motif prediction
	echo "[INFO] - For motif prediction, extracting 500 first sequences" | tee -a $log
	head -n 500 ${outdir}/${name}/train/${name}_train.txt > ${outdir}/${name}/train/${name}_train_500.txt
	python $SeqConv_motif ${outdir}/${name}/train/${name}_train_500.txt $peak_size $name $outdir ${outdir}/${name}/model/${name}_model.keras $seed
	
	#----- plotting motif 
	echo "[INFO] - Plotting motif" | tee -a $log
	$PATHRscript $plot_logo ${outdir}/${name}/motif/${name}_seq.meme ${outdir}/${name}/motif/${name}_logo.png

	# clean 
	rm $outdir/${name}/tmp_neg_shuf.bed $outdir/${name}/tmp_pos_shuf.bed

	echo "[ALL DONE]"
}

#-------------------------------------------------------------------------------
generate_weighted_fasta(){
	# FUNCTION: generate_weighted_fasta
	# DESCRIPTION:
	#   This function generates a weighted FASTA file from a BED file containing genomic regions.
	#   It extracts sequences from a reference genome and optionally weights them based on a coverage column.
	#   The resulting FASTA file can be used for downstream analyses such as motif discovery or k-mer analysis.
	#
	# USAGE:
	# generate_weighted_fasta -f <FILE> -c [INT] -s <INT> -o <PATH> 
	#        -p <STRING> - g <FILE>
	#
	# ARGUMENTS:
	#   -f          : Path to the BED file containing genomic regions. Required.
	#   -s          : Number of sequences to include in the weighted FASTA file. Required.
	#   -o          : Directory where results will be saved. Required.
	#   -p          : Prefix for output files. Required.
	#   -c          : Column number in the BED file that contains coverage information. Optional. Default: 0.
	#   -g          : Path to the reference genome FASTA file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas".
	#   -h, --help  : Display usage information for this function.
	#
	# NOTES:
	#   - The function requires `bedtools` to extract sequences from the reference genome.
	#   - If the `-c` option is provided, the BED file is sorted by the specified coverage column before selecting the top sequences.
	#   - If the coverage column is used, the FASTA headers are modified to include coverage information.
	echo "
	_ _ _ ____ _ ____ _  _ ___ _ _  _ ____     ____ ____ ____ ___ ____
	| | | |___ | | __ |__|  |  | |\ | | __     |___ |__| [__   |  |__|
	|_|_| |___ | |__] |  |  |  | | \| |__] ___ |    |  | ___]  |  |  |

	" 
	local colcov=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; 
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-f)
				local file=$2
				echo "-> Analyzed bed file is: ${2}";shift 2;;
			-c)
				local colcov=$2
				echo "-> Column with corresponding coverage is: ${2}";shift 2;;
			-s)
				local size=$2
				echo "-> Number of seqs (that will be used by KMAC) for the training set is $2"; shift 2;;
			-o)
				local out_dir=$2
				echo "-> Output directory is: ${2}";shift 2;;
			-p)
				local prefix=$2
				echo "-> Prefix is: ${2}";shift 2;;
			-g)
				local genome=$2
				echo "-> Reference genome to extract seq from is:$2"; shift 2;;
			-h)
				usage generate_weighted_fasta; return;;
			--help)
				usage generate_weighted_fasta; return;;
			*)
				echo "Error in arguments"
				echo $1; usage; return;;
		esac
	done

	local Errors=0
	if [ -z $file ]; then echo "ERROR: -f argument needed"; Errors+=1; fi
	if [ -z $size ]; then echo "ERROR: -s argument needed"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
	if [ -z $prefix ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage generate_weighted_fasta; return 1; fi

	mkdir -p -m 774 $out_dir/$prefix
	local log=$out_dir/$prefix/log.txt
	echo $Argline > $log
	echo -e "[INFO] - Generating weighted fasta" | tee -a $log
	
	# Generate weighted bed and fasta files
	if [ $colcov == 0 ]; then
		cat $file | awk -v top=$size 'NR <= top {print}' > $out_dir/$prefix/${prefix}.weighted.bed
		bedtools getfasta -fi $genome -bed $out_dir/$prefix/${prefix}.weighted.bed -fo $out_dir/$prefix/tmp.fas -fo $out_dir/$prefix/${prefix}.weighted.fas
	else
		echo "[INFO] - Using coverage column $colcov to weight sequences" | tee -a $log
		echo "[WARNING] - THIS WEIGHTING OF SEQUENCE HAS TO BE TRIPLE CHECKED" #TODO 17/09/2025 JL to RBM: what does this mean exactly?
		# AJ do we need to sort the file by coverage to be sure here ? TO 
		cat $file | cut -f 1,2,3,$colcov \
			| sort -k4,4rn | awk -v top=$size 'NR <= top {print}' > $out_dir/$prefix/${prefix}.weighted.bed
		bedtools getfasta -fi $genome -bed $out_dir/$prefix/${prefix}.weighted.bed -fo $out_dir/$prefix/tmp.fas
		cat $out_dir/$prefix/tmp.fas | awk -F ":" '{print ">"$3,":"$4,$1}' \
			| sed 's/^> : //g' | sed 's/ >/ /g' \
			| sed 's/ :/:/g' > $out_dir/$prefix/${prefix}.weighted.fas
	fi
}

#-------------------------------------------------------------------------------

compute_kmer(){
	# FUNCTION: compute_kmer
	# DESCRIPTION:
	#   This function performs k-mer analysis using the KMAC program
	#   (part of the GEM suite) to identify enriched k-mer sets in a
	#   positive sequence set, optionally compared against a negative
	#   set. It splits the positive (and, if provided, negative) BED
	#   file into training and testing sets, either by a training
	#   percentage (-split) or by taking the top N sequences
	#   (-ls, using the original, non-shuffled order). It then
	#   extracts FASTA sequences for the training sets with
	#   bedtools, runs KMAC via the GEM jar, and moves/renames the
	#   KMAC output directory into the function's output directory.
	#
	# USAGE:
	#   compute_kmer -b1 <FILE> -n1 <STRING> -kwi <INT> -kmi <INT>
	#                -kma <INT> -o <PATH> (-split <INT> | -ls <INT>)
	#                -g [FILE] -b2 [FILE] -ktop [INT] -r [STRING]
	#
	# ARGUMENTS:
	#   -b1       : BED file for positive sequences (required).
	#   -n1       : Name of the first condition, usually the positive
	#               set (required).
	#   -kwi      : K-mer window size (required).
	#   -kmi      : Minimum k-mer length (required).
	#   -kma      : Maximum k-mer length (required).
	#   -o        : Output directory (required).
	#   -split    : Percentage of sequences to use for the training
	#               set. Mutually exclusive with -ls; exactly one of
	#               -split or -ls must be provided.
	#   -ls       : Number of top sequences (original, non-shuffled
	#               order) to use for the training set. Mutually
	#               exclusive with -split; exactly one of -split or
	#               -ls must be provided.
	#   -g        : Path to the genome FASTA file (optional, has an
	#               internal default).
	#   -b2       : BED file for negative sequences (optional). If
	#               omitted, KMAC is run without a negative set.
	#   -ktop     : Number of top k-mers to report (default: 5).
	#   -r        : Run ID, used to name the KMAC outputs.
	#   -h, --help: Display usage information for this function.
	#
	# DEPENDENCIES:
	#   - bedtools
	#   - java, and the KMAC program via the $GEM environment
	#     variable (must point to the GEM jar file).
	#
	# NOTES:
	#   - Exactly one of -split or -ls must be provided; providing
	#     both, or neither, aborts the function.
	#   - In SPLIT mode, the positive (and negative, if any) BED file
	#     is shuffled before splitting into train/test.
	#   - In LEARNING_SIZE mode (-ls), the ORIGINAL (non-shuffled)
	#     order of the BED file is used to pick the top N lines as
	#     the training set.
	#   - If a negative set is provided and -ls exceeds the number of
	#     negative peaks, the function aborts.
	#   - The KMAC output directory is expected to be named
	#     "<run_id>_outputs" in the current working directory; it is
	#     moved to "<out_dir>/<mode>_<background>/KMAC_outputs" and
	#     any absolute paths in its ksm_list.txt file are rewritten
	#     accordingly.
	#   - All actions are logged to a log.txt file in the final
	#     output directory.

    echo "
        ____ ____ _  _ ___  _  _ ___ ____     _  _ _  _ ____ ____
        |    |  | |\/| |__] |  |  |  |___     |_/  |\/| |___ |__/
        |___ |__| |  | |    |__|  |  |___ ___ | \_ |  | |___ |  |

	"
	local ktop=5; local annotation="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.bed"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local bed2=""; local split=""; local learning_size="";
	local Argline=$@
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -g)     local genome="$2"
					echo "-> Reference genome to extract seq from is:$2"; shift 2;;
            -b1)    local bed1="$2"
					echo "-> BED file for positive set: ${2}"; shift 2;;
            -b2)    local bed2="$2"
					echo "-> BED file for negative set: ${2}"; shift 2;;
            -n1)    local name1="$2"
					echo "-> Name for positive set: ${2}"; shift 2;;
            -kwi)   local kmerWindow="$2"
					echo "-> K-mer window size: ${2}"; shift 2;;
            -kmi)   local kmerMin="$2"
					echo "-> Minimum k-mer size: ${2}"; shift 2;;
            -kma)   local kmerMax="$2"
					echo "-> Maximum k-mer size: ${2}"; shift 2;;
            -ktop)  local ktop="$2"
					echo "-> Number of top k-mers to report: ${2}"; shift 2;;
            -o)     local out_dir="$2"
					echo "-> Output directory: ${2}"; shift 2;;
            -split) local split="$2"
					echo "-> Using SPLIT mode with ${2}% of data for training"; shift 2;;
            -ls)    local learning_size="${2%.*}"
					echo "-> Using LEARNING_SIZE mode with top ${2} sequences for training"; shift 2;;
			-r)
				local RUN=$2
				echo "-> Run ID for worflows: ${2}";shift 2;;
			-h)
				usage compute_kmer; return;;
			--help)
				usage compute_kmer; return;;
			*)
				echo "Error in arguments"
				echo $1; usage; return;;
        esac
    done

    # -----------------------------
    # Validate arguments
    # -----------------------------
    Errors=0
    [[ -z "$bed1" ]]        && echo "[ERROR] -b1 missing"       && Errors=1
    [[ -z "$name1" ]]       && echo "[ERROR] -n1 missing"       && Errors=1
    [[ -z "$kmerWindow" ]]  && echo "[ERROR] -kwi missing"      && Errors=1
    [[ -z "$kmerMin" ]]     && echo "[ERROR] -kmi missing"      && Errors=1
    [[ -z "$kmerMax" ]]     && echo "[ERROR] -kma missing"      && Errors=1
    [[ -z "$out_dir" ]]     && echo "[ERROR] -o missing"        && Errors=1

	# log for agruments 
	local log=$out_dir/log.txt
	echo $Argline > $log

    # Split vs learning_size mutual exclusivity
    if [[ -n "$split" && -n "$learning_size" ]]; then
        echo "[ERROR] - You cannot specify both --split and --learning_size" | tee -a $log
        Errors=1
    fi
    if [[ -z "$split" && -z "$learning_size" ]]; then
        echo "[ERROR] - One of --split or --learning_size must be provided" | tee -a $log
        Errors=1
    fi

    if (( Errors > 0 )); then
        echo "[ERROR] - Aborting compute_kmer" | tee -a $log
        return 1
    fi

    echo "[INFO] - Preparing positive BED" | tee -a $log

    # -----------------------------
    # Define mode and background
    # -----------------------------
    if [[ -n "$split" ]]; then
        mode="split_${split}"
    else
        mode="ls_${learning_size}"
    fi

    if [[ -n "$bed2" ]]; then
        background="providedNS"
    else
        background="dinuclNS"
    fi

    # -----------------------------
    # Prepare directories
    # -----------------------------
    out_dir="${out_dir}/${mode}_${background}"
	
	# log for chosen mode
	local log=$out_dir/log.txt
	echo $Argline > $log

    echo "[INFO] - Output directory set to: $out_dir" | tee -a $log

    mkdir -p "${out_dir}/Infiles"

    # -----------------------------
    # Prepare positive set
    # -----------------------------

    # Shuffle once (needed for split mode; harmless for ls mode)
    shuf "$bed1" > "${out_dir}/Infiles/tmp_pos_shuf.bed"
    total_pos=$(wc -l < "${out_dir}/Infiles/tmp_pos_shuf.bed")

    # MODE A: SPLIT
    if [[ -n "$split" ]]; then
        echo "[INFO] - Using SPLIT mode: ${split}% train" | tee -a $log
        echo "[INFO] - Split mode shuffles peaks before splitting" | tee -a $log

        train_pos=$(( total_pos * split / 100 ))
        if (( train_pos <= 0 )); then
            echo "[WARNING] - Computed train_pos <= 0; forcing train_pos=1" | tee -a $log
            train_pos=1
        fi

        head -n "$train_pos" "${out_dir}/Infiles/tmp_pos_shuf.bed" > "${out_dir}/Infiles/pos_train.bed"
        tail -n +$((train_pos + 1)) "${out_dir}/Infiles/tmp_pos_shuf.bed" > "${out_dir}/Infiles/pos_test.bed"

    # MODE B: LEARNING_SIZE
    elif [[ -n "$learning_size" ]]; then
        echo "[INFO] - Using LEARNING_SIZE mode: top ${learning_size} sequences" | tee -a $log

        if (( learning_size > total_pos )); then
            echo "[WARNING] - learning_size (${learning_size}) exceeds number of peaks (${total_pos})" | tee -a $log
        fi

        # Use ORIGINAL order of $bed1 (not shuffled) for LS
        head -n "$learning_size" "$bed1" > "${out_dir}/Infiles/pos_train.bed"
        tail -n +$((learning_size + 1)) "$bed1" > "${out_dir}/Infiles/pos_test.bed"
    fi

    # Create fasta from training set
    fastaPOS="${out_dir}/Infiles/pos_train.fa"
    bedtools getfasta -fi "$genome" -bed "${out_dir}/Infiles/pos_train.bed" -fo "$fastaPOS"

    rm -f "${out_dir}/Infiles/tmp_pos_shuf.bed"

    # -----------------------------
    # Prepare negative set (optional)
    # -----------------------------
    if [[ -n "$bed2" ]]; then

        echo "[INFO] Preparing negative BED using same mode: $mode" | tee -a $log

        # Shuffle once (used for split; harmless for ls)
        shuf "$bed2" > "${out_dir}/Infiles/tmp_neg_shuf.bed"
        total_neg=$(wc -l < "${out_dir}/Infiles/tmp_neg_shuf.bed")

        # Warn if pos/neg size mismatch (not fatal)
        if (( total_neg != total_pos )); then
            echo "[WARNING] Positive and negative sets have different sizes:" | tee -a $log
            echo "          POS = ${total_pos}, NEG = ${total_neg}" | tee -a $log
        fi

        # MODE A: SPLIT
        if [[ -n "$split" ]]; then
            echo "[INFO] Negative set: SPLIT mode (${split}% train)" | tee -a $log
 
            train_neg=$(( total_neg * split / 100 ))
            if (( train_neg <= 0 )); then
                echo "[WARNING] Computed train_neg <= 0; forcing train_neg=1" | tee -a $log
                train_neg=1
            fi

            head -n "$train_neg" "${out_dir}/Infiles/tmp_neg_shuf.bed" > "${out_dir}/Infiles/neg_train.bed"
            tail -n +$((train_neg + 1)) "${out_dir}/Infiles/tmp_neg_shuf.bed" > "${out_dir}/Infiles/neg_test.bed"

        # MODE B: LEARNING_SIZE
        elif [[ -n "$learning_size" ]]; then
            echo "[INFO] Negative set: LEARNING_SIZE mode (top ${learning_size})" | tee -a $log

            if (( learning_size > total_neg )); then
                echo "[ERROR] learning_size (${learning_size}) exceeds number of negative peaks (${total_neg})"
                rm -f "${out_dir}/Infiles/tmp_neg_shuf.bed"
                return 1
            fi

            head -n "$learning_size" "$bed2" > "${out_dir}/Infiles/neg_train.bed"
            tail -n +$((learning_size + 1)) "$bed2" > "${out_dir}/Infiles/neg_test.bed"
        fi

        # Build FASTA for the negative train set
        fastaNEG="${out_dir}/Infiles/neg_train.fa"
        bedtools getfasta -fi "$genome" -bed "${out_dir}/Infiles/neg_train.bed" -fo "$fastaNEG"

        rm -f "${out_dir}/Infiles/tmp_neg_shuf.bed"
    fi

    # -----------------------------
    # Run KMAC
    # -----------------------------
    echo "[INFO] Running KMAC..." | tee -a $log

    if [[ -z "$GEM" ]]; then
        echo "[ERROR] GEM variable not set. export GEM=/path/to/gem.jar" | tee -a $log
        return 1
    fi

    local out_name="${mode}_${background}"

	#below in the command I could have done --out_name \"$out_name\" \
	CMD="java -Xmx8G -jar \"$GEM\" KMAC \
		--pos_seq \"$fastaPOS\" \
		--k_win \"$kmerWindow\" \
		--k_min \"$kmerMin\" \
		--k_max \"$kmerMax\" \
		--k_top \"$ktop\" \
		--out_name $RUN \
		--print_aligned_seqs ON \
		--gc -1"

	if [[ -n "$bed2" ]]; then
		CMD="$CMD --neg_seq \"$fastaNEG\""
	fi

    # Run the command
    eval $CMD

    # -----------------------------
    # Move results to final location
    # -----------------------------
    # Paths
	local kmac_out="${RUN}_outputs"
    local kmac_dir="$out_dir/KMAC_outputs"

    # current path
    local tmpPWD=$(pwd)
    echo "[INFO] Current working directory: $tmpPWD" | tee -a $log

    # Update output file in current directory to outdir path (in case KMAC outputs absolute paths in its result files)
    if [[ -f "$kmac_out/${RUN}.ksm_list.txt" ]]; then
        echo "[INFO] Updating paths in $kmac_out/${RUN}.ksm_list.txt..." | tee -a $log
        sed -i "s~$tmpPWD/$kmac_out~$kmac_dir~g" "$kmac_out/${RUN}.ksm_list.txt" 
        echo "[INFO] Updated paths:" | tee -a $log
    else
        echo "[WARNING] File $kmac_out/KMAC.ksm_list.txt not found. Skipping path update."
    fi

	# Move output directory from current location to final output directory
	if [[ -d "$kmac_out" ]]; then
		if [[ -d "$kmac_dir" ]]; then
			echo "[WARNING] Output directory $kmac_dir already exists. It will be overwritten." | tee -a $log
			rm -rf "$kmac_dir"
		fi
		mv "$kmac_out" "$kmac_dir"
		echo "[INFO] Moved KMAC output to: $kmac_dir" | tee -a $log
	else
		echo "[ERROR] Expected KMAC output directory '$kmac_out' not found. Please check if KMAC ran successfully." | tee -a $log
		return 1
	fi

    printf "\nKMAC has finished\n\n"
    echo "[REPORT] - KMAC results are stored in: $kmac_dir" | tee -a $log
}

#-------------------------------------------------------------------------------
compute_ROCS(){
	# FUNCTION: compute_ROCS
	# DESCRIPTION:
	#   This function computes and plots Receiver Operating Characteristic (ROC) curves
	#   for various types of sequence analysis, including TFFM, K-mer, DNAshape, PWM, and SeqConv models.
	#   It evaluates the performance of these models in distinguishing between positive and negative sequence sets.
	#   The function supports multiple matrices, peak files, and negative sets, and allows customization
	#   of various parameters such as offsets, peak size, and scoring methods.
	#
	# USAGE:
	#   compute_ROCS -p <LIST of FILE> -ns <LIST of FILE> -m <LIST of FILE>
	#        -n <LIST of STRING> -od <PATH> -g [FILE] [-pc] -color [LIST of STRING]
	#        -nKSM [INT] -sKSM [STRING]
	#
	# ARGUMENTS:
	#   -p          : List of BED or FASTA files containing peak sequences. Required.
	#   -ns         : List of BED or FASTA files containing negative sequences. Required.
	#   -m          : List of matrix files (TFFM, K-mer, PWM or SeqConv). Required.
	#   -n          : List of names associated with each peak/negative set/matrix combination. Required.
	#   -g          : Path to the genome FASTA file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas".
	#   -od         : Directory where results will be saved. Required.
	#   -pc         : Activate Pocc mode for PWM analysis. Optional.
	#   -d          : Dimer conformation for PWM analysis. Optional. Default: "Mono".
	#   -ol         : Left offset for peak extension. Optional. Default: 0.
	#   -or         : Right offset for peak extension. Optional. Default: 0.
	#   -nKSM       : The k-mer set (m0 or m1 or m2 ect) to use in KSM search. Optional. Default: 0.
	#   -sKSM       : K-mer score method (SUM, BEST, or MEAN). Optional. Default: "BEST".
	#   -ps         : Size of peaks in base pairs. Optional. Default: 200.
	#   -color      : List of colors for ROC curves. Optional. Default: A predefined list of colors.
	#   -fpr_lim    : Limit for False Positive Rate in ROC curves. Optional. Default: 1 (full ROC).
	#   -h, --help  : Display usage information for this function.
	#
	# DEPENDENCIES:
		# - ../bin/parse_KSM_scores.py
		# - ../bin/compute_POcc.py
		# - ../bin/get_best_score_tffm.py
		# - ../bin/scores.py
		# - ../bin/plots_ROCS_multiple.py

	echo "
	____ ____ _  _ ___  _  _ ___ ____     ____ ____ ____ ____
	|    |  | |\/| |__] |  |  |  |___     |__/ |  | |    [__ 
	|___ |__| |  | |    |__|  |  |___ ___ |  \ |__| |___ ___]
	
	"
	local pocc=false; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544' '#244c63' '#566324' '#4db352' '#4d99b3' '#7c4db3' '#b34da8' '#4db396' '#b34d75' '#FF4000' '#240B3B' '#61380B' '#7A2513' '#7DAAAD' '#F8DB4B' '#B4F84B' '#687574' '#9FE1DB'); local number_of_ksm=0; local ksm_score_meth="BEST"; local dimer="Mono"; local offset_left=0; local offset_right=0; local peaksize=200; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local seed=1234; local fpr_lim=1;
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=("${!2}")
				echo "-> Peaks files set to: ${peaks[@]}";shift 2;;
			-ns)
				local negative_sets=("${!2}")
				echo "-> Negative files set to: ${negative_sets[@]}";shift 2;;
			-m)
				local matrices=("${!2}")
				echo "-> Matrix files set to: ${matrices[@]}";shift 2;;
			-n)
				local names=("${!2}")
				echo "-> Names associated set to: ${names[@]}";shift 2;;
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-od)
				local results=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-pc)
				local pocc=true
				echo "-> Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
			-d)
				local dimer=$2
				echo "-> PFM will be tested as a dimer too: conformations: ${2}"; shift 2;;
			-ol)
				local offset_left=${2%.*}
				echo "-> Offset on the left set to: ${2%.*}";shift 2;;
			-or)
				local offset_right=${2%.*}
				echo "-> Offset on the right set to: ${2%.*}";shift 2;;
			-nKSM)
				local number_of_ksm=${2%.*}
				echo "-> The k-mer set (m0 or m1 or m2 ect) to use in the ksm search: ${2}";shift 2;;
			-sKSM)
				local ksm_score_meth=$2
				echo "-> Kmer score method (SUM, BEST or MEAN): ${2}";shift 2;;
			-ps)
				local peaksize=$2
				echo "-> Peaks size: ${2}";shift 2;;
			-color)
				local colors=("${!2}")
				echo "-> Colors set to: ${colors[@]}";shift 2;;
			-s)
				local seed=$2
				echo "-> Seed for shuffling peaks and negative sets: ${2}";shift 2;;
			-fpr_lim) 
				local fpr_lim=$2
				echo "-> You asked early ROCS version with a FPR limit set to: ${2} (if lower than 1, computing partial early ROCs)";shift 2;;
			-h)
				usage compute_ROCS; return;;
			--help)
				usage compute_ROCS; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage compute_ROCS; return;;
		esac
	done

	# Create log file and print arguments
	log=$results/log.txt
	echo $Argline > $log
	
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed or needs to be a LIST"; Errors+=1; fi
	if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed or needs to be a LIST"; Errors+=1; fi
	if [ -z $matrices ]; then echo "ERROR: -m argument needed or needs to be a LIST"; Errors+=1; fi
	if [ -z $names ]; then echo "ERROR: -n argument needed or needs to be a LIST"; Errors+=1; fi
	if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	
	# checking that lists are of same length
	if [ ${#names[@]} -ne ${#peaks[@]} ] || [ ${#names[@]} -ne ${#matrices[@]} ] || [ ${#names[@]} -ne ${#negative_sets[@]} ] || [ ${#names[@]} -gt ${#colors[@]} ]; then
		if [ ${#peaks[@]} -eq ${#negative_sets[@]} ] && [ ${#peaks[@]} -ne ${#matrices[@]} ] && [ ${#names[@]} -eq ${#peaks[@]} ] && [ ${#colors[@]} -gt ${#matrices[@]} ] ; then
			for ((i=1;i<${#matrices[@]};i++)); do
				peaks+=("${peaks[0]}")
				negative_sets+=("${negative_sets[0]}")
				names+=("${names[0]}")
			done
		else
			echo ${#peaks[@]}
			echo ${#negative_sets[@]}
			echo ${#matrices[@]}
			echo ${#names[@]}

			if [ ${#peaks[@]} -eq ${#negative_sets[@]} ] && [ ${#peaks[@]} -ne ${#matrices[@]} ] && [ ${#names[@]} -eq ${#matrices[@]} ] && [ ${#colors[@]} -gt ${#matrices[@]} ] ; then
				for ((i=1;i<${#matrices[@]};i++)); do
				peaks+=("${peaks[0]}")
				negative_sets+=("${negative_sets[0]}")
			done
			else
				echo "[ERROR] - -n, -m & -ns have to be lists of same length";
			fi
		fi
	fi
	if [ $Errors -gt 0 ]; then usage compute_ROCS; return 1; fi
	
	# making sure colors are in hex format
	for ((i=0;i<${#colors[@]};i++)); do
		[[ ${colors[$i]} =~ ^#.* ]] || colors[$i]="#${colors[$i]}"
	done

	
	local scores=()
	local endnames=()
	i=0
	mkdir -p -m 774 ${results}/scores
	mkdir -p -m 774 ${results}/sequences

	local log=$results/log.txt
	echo $Argline > $log

	# processing each peak/negative set/matrix combination
	for ((i=0;i<${#peaks[@]};i++))
	do
		local negative_set=${negative_sets[i]}
		local name=${names[i]// /_}
		local matrice=${matrices[i]}
		local peak=${peaks[i]}
		echo -e "\tpeaks $peak" | tee -a $log
		echo -e "\tNegset $negative_set" | tee -a $log
		echo -e "\tname $name" | tee -a $log
		echo -e "\tmatrice $matrice" | tee -a $log
		echo "[INFO] - Preparing data" | tee -a $log

		if [[ $peak != *".fa"* ]]; then # if bed file, convert to fasta
			if [[ $matrice == *".xml"* ]]; then # TFFM
				echo "[INFO] - Preparing data for TFFM" | tee -a $log
				awk -v OFS="\t" '{print $1,$2-2,$3}' $peak > $results/sequences/pos_dnashape_tffm_reworked.bed
				awk -v OFS="\t" '{print $1,$2-2,$3}' $negative_set > $results/sequences/neg_dnashape_tffm_reworked.bed
				local peak=$results/sequences/pos_dnashape_tffm_reworked.bed
				local negative_set=$results/sequences/neg_dnashape_tffm_reworked.bed
			fi
			if [[ $matrice == *".pkl"* ]]; then # DNAshape
				echo "[INFO] - Preparing data for DNAshape" | tee -a $log
				awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $peak > $results/sequences/pos_dnashape_tffm_reworked.bed
				awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $negative_set > $results/sequences/neg_dnashape_tffm_reworked.bed
				local peak=$results/sequences/pos_dnashape_tffm_reworked.bed
				local negative_set=$results/sequences/neg_dnashape_tffm_reworked.bed
			fi
			if [[ $matrice == *".keras"* ]]; then # SeqConv
				echo "[INFO] - Preparing data for SeqConv" | tee -a $log
				bedtools getfasta -fi $genome -bed $peak | awk -v OFS="\t" '$1 !~ /^>/ {print 1,$1}' > $results/seq_pos.txt
				bedtools getfasta -fi $genome -bed $negative_set | awk -v OFS="\t" '$1 !~ /^>/ {print 0,$1}' > $results/seq_neg.txt
				# concatenate 
				paste $results/seq_pos.txt $peak > $results/seq_train_pos.txt
				paste $results/seq_neg.txt $negative_set > $results/seq_train_neg.txt
				cat $results/seq_train_pos.txt $results/seq_train_neg.txt > $results/sequences/seqconv_pos_neg.txt
				rm $results/seq*.txt # clean
				seqconv_dataset=$results/sequences/seqconv_pos_neg.txt #to feed the model
			fi

			# Get sequences
			echo "[INFO] - Extracting sequences from bed files" | tee -a $log
			bedtools getfasta -fi $genome -bed $peak -name -fo $results/sequences/PWM_pos_set.fa
			sed -i "/^>/! {s/n/-/g; s/\(.*\)/\U\1/g}" $results/sequences/PWM_pos_set.fa
			local peak=$results/sequences/PWM_pos_set.fa
			bedtools getfasta -fi $genome -bed $negative_set -name -fo $results/sequences/PWM_neg_set.fa
			sed -i "/^>/! {s/n/-/g; s/\(.*\)/\U\1/g}" $results/sequences/PWM_neg_set.fa
			local negative_set=$results/sequences/PWM_neg_set.fa
		fi

		# Computing scrore for each matrix type
		echo "[INFO] - Computing scores" | tee -a $log
		if [[ $matrice == *".xml"* ]]; then # TFFM scores computation
			echo "===== TFFM =====" | tee -a $log
			echo "[INFO] - Computing TFFM scores" | tee -a $log
			$Python_TFFM $tffmscores -o ${results}/scores/TFFM_scores_pos.tsv -pos $peak -t ${matrice}
			echo "[INFO] - Computing TFFM scores on negative set" | tee -a $log
			$Python_TFFM $tffmscores -o ${results}/scores/TFFM_scores_neg.tsv -pos $negative_set -t ${matrice}
			paste <(awk '$1!="None"{print $8;next}{print "0.0"}' ${results}/scores/TFFM_scores_pos.tsv) <(awk '$1!="None"{print $8;next}{print "0.0"}' ${results}/scores/TFFM_scores_neg.tsv) >"${results}/scores/TFFM_${name}_scores_roc.tsv"
			# adding score files to the list
			scores+=("${results}/scores/TFFM_${name}_scores_roc.tsv")
			endnames+=("${name}_TFFM")
		fi

		if [[ $matrice == *".txt"* ]]; then # K-mer scores computation
			echo "===== KSM =====" | tee -a $log

			if [ "$number_of_ksm" -eq 0 ]; then
				head -n 1 "$matrice" > "${results}/scores/tmp_ksm_list"
			else
				sed -n "$((number_of_ksm + 1))p" "$matrice" > "${results}/scores/tmp_ksm_list"
			fi		
			
			cat $peak | tr ':' '_' > ${results}/scores/tmpPos.fsa # coz gem does not like ":" ...
			echo "[INFO] - Scanning positive sequences with KSM" | tee -a $log
			java -Xmx8G -jar $GEM KSM --fasta ${results}/scores/tmpPos.fsa --ksm ${results}/scores/tmp_ksm_list --out ${results}/scores/ksm.scan.pos
			cat $negative_set | tr ':' '_' > ${results}/scores/tmpNeg.fsa
			echo "[INFO] - Scanning negative sequences with KSM" | tee -a $log
			java -Xmx8G -jar $GEM KSM --fasta ${results}/scores/tmpNeg.fsa --ksm ${results}/scores/tmp_ksm_list --out ${results}/scores/ksm.scan.neg
			rm ${results}/scores/tmpNeg.fsa ${results}/scores/tmpPos.fsa 
			#rm ${results}/scores/tmp_ksm_list
			python3 $parse_KSM_scores -p ${results}/scores/ksm.scan.pos.motifInstances.txt -n ${results}/scores/ksm.scan.neg.motifInstances.txt -o ${results}/scores/KSM_${name}_${ksm_score_meth}_scores_roc.tsv -m $ksm_score_meth
			mv ${results}/scores/ksm.scan.pos.motifInstances.txt ${results}/scores/KSM_${ksm_score_meth}_scores_pos.txt
			mv ${results}/scores/ksm.scan.neg.motifInstances.txt ${results}/scores/KSM_${ksm_score_meth}_scores_neg.txt
			
			# adding score files to the list
			scores+=("${results}/scores/KSM_${name}_${ksm_score_meth}_scores_roc.tsv")
			endnames+=("${name}_KSM")
		fi
	
		if [[ $matrice == *".pkl"* ]]; then # DNAshape
			echo "===== DNAshape =====" | tee -a $log
			if [[ $matrice == *"_fo_classifier"* ]]; then # DNAshape done with first order TFFM
				echo "[INFO] - Applying the trained detailed TFFM + DNA shape classifier on foreground sequences." | tee -a $log
				$Python_TFFM $ComputeDNAshaped applyTFFM -T $(dirname $matrice)/usedTFFM_fo.xml -i $peak -I $results/sequences/pos_dnashape_tffm_reworked.bed -c $matrice -o ${results}/scores/DNAshape_scores_pos.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
				echo "[INFO] - Applying the trained detailed TFFM + DNA shape classifier on background sequences." | tee -a $log
				$Python_TFFM $ComputeDNAshaped applyTFFM -T $(dirname $matrice)/usedTFFM_fo.xml -i $negative_set -I $results/sequences/neg_dnashape_tffm_reworked.bed -c $matrice -o ${results}/scores/DNAshape_scores_neg.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
				paste <(awk 'NR!=1{print $6}' ${results}/scores/DNAshape_scores_pos.txt | sort -nr ) <(awk 'NR!=1{print $6}' ${results}/scores/DNAshape_scores_neg.txt | sort -nr ) >"${results}/scores/DNAshape_${name}_scores_roc.tsv"
				# adding score files to the list
				scores+=("${results}/scores/DNAshape_${name}_scores_roc.tsv")
				endnames+=("${name}_DNAshape")
			fi
			if [[ $matrice == *"_PSSM_classifier"* ]]; then # DNAshape done with PSSM/PFM
				echo "[INFO] - Applying the trained detailed PSSM + DNA shape classifier on foreground sequences." | tee -a $log
				$Python_TFFM $ComputeDNAshaped applyPSSM -f $(dirname $matrice)/JASPAR.pfm -i $peak -I $results/sequences/pos_dnashape_tffm_reworked.bed -c $matrice -o ${results}/scores/DNAshapePSSM_pred_pos.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
				echo "[INFO] - Applying the trained detailed PSSM + DNA shape classifier on background sequences." | tee -a $log
				$Python_TFFM $ComputeDNAshaped applyPSSM -f $(dirname $matrice)/JASPAR.pfm -i $negative_set -I $results/sequences/neg_dnashape_tffm_reworked.bed -c $matrice -o ${results}/scores/DNAshapePSSM_pred_neg.txt -1 $helt $prot $mgw $roll -2 $helt2 $prot2 $mgw2 $roll2 -n -v 0.000001
				paste <(awk 'NR!=1{print $6}' ${results}/scores/DNAshapePSSM_pred_pos.txt | sort -nr ) <(awk 'NR!=1{print $6}' ${results}/scores/DNAshapePSSM_pred_neg.txt | sort -nr ) >"${results}/scores/tab_${name}_PSSMshape.tsv"
				# adding score files to the list
				scores+=("${results}/scores/tab_${name}_PSSMshape.tsv")
				endnames+=("${name}_PSSMShape")
			fi
		fi

		if [[ $matrice == *".jaspar"* ]]; then # PWM from JASPAR to PFM conversion
			matrixname=$(basename $matrice .jaspar)
			echo "[INFO] - Converting JASPAR to PFM for ${name}" | tee -a $log
			echo "MATRIX COUNT ASYMMETRIC ${matrixname} SIMPLE" > ${results}/scores/${matrixname}.pfm
			awk -F"\t" -v OFS="\t" '!/^>/' $m \
				| awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $i}}
				NF>p { p = NF }
				END {for(j=1; j<=p; j++) {str=a[1,j];for(i=2; i<=NR; i++){str=str" "a[i,j];}print str}}'\
				| awk '!/^[\[\]]/' \
				| sed 's/ /\t/g' >> ${results}/scores/${matrixname}.pfm
			local matrice=${results}/scores/${matrixname}.pfm
		fi
	
		if [[ $matrice == *".pfm"* ]]; then  # PWM scores
			echo "===== PWM =====" | tee -a $log
			echo "[INFO] - Computing PWM scores for ${name}"
			echo $peak
			python3 $scores_prog -m ${matrice} -f $peak -o ${results}/scores/
			echo "[INFO] - Computing PWM scores for ${name} negative set"
			python3 $scores_prog -m ${matrice} -f $negative_set -o ${results}/scores/
			paste <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $peak).scores | sort -u -k1,1 | awk '{print $8}'  | sort -nr )  <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $negative_set).scores | sort -u -k1,1 | awk '{print $8}' | sort -nr ) | awk '{print $0}' >"${results}/scores/PWM_${name}_scores_roc.tsv"
			# adding score files to the list
			scores+=("${results}/scores/PWM_${name}_scores_roc.tsv")
			endnames+=("${name}_PWM")
			
			if $pocc ; then # PWM Pocc
				echo "===== Pocc =====" | tee -a $log
				echo "[INFO] - Computing Pocc scores for ${name}" 
				python3 $pocc_pfm -s ${results}/scores/$(basename $peak).scores -o ${results}/scores/Pocc_pos.pocc
				echo "[INFO] - Computing Pocc scores for ${name} negative set" 
				python3 $pocc_pfm -s ${results}/scores/$(basename $negative_set).scores -o ${results}/scores/Pocc_neg.pocc
				paste <( sort -nr ${results}/scores/Pocc_pos.pocc )  <( sort -nr ${results}/scores/Pocc_neg.pocc ) > "${results}/scores/${name}_Pocc_scores.tsv"
				# adding score files to the list
				scores+=("${results}/scores/${name}_Pocc_scores.tsv")
				endnames+=("${name}_Pocc")
			fi
	
			if [ ${dimer} != "Mono" ]; then # Dimer PFM scores
				local dist=${dimer:2}
				local conf=${dimer:0:2}
				local dimer_name="${name}_$dimer"
				echo "[INFO] - Computing dimer matrix for ${name} with conformation ${conf} and distance ${dist}" | tee -a $log
				bash ${dimer_builder} ${matrice} ${dimer_name} $dist $conf $offset_left $offset_right > ${results}/scores/$(basename $matrice)_${dimer}.pfm
				echo "[INFO] - Computing dimer PWM scores for ${name} with conformation ${conf} and distance ${dist}" | tee -a $log
				python3 $scores_prog -m ${results}/scores/$(basename $matrice)_${dimer}.pfm -f $peak -o ${results}/scores/
				echo "[INFO] - Computing dimer PWM scores for ${name} negative set with conformation ${conf} and distance ${dist}" | tee -a $log
				python3 $scores_prog -m ${results}/scores/$(basename $matrice)_${dimer}.pfm -f $negative_set -o ${results}/scores/
				paste <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $peak).scores | sort -u -k1,1 | awk '{print $8}'  | sort -nr )  <(sort -u -k1,1 -k8,8nr ${results}/scores/$(basename $negative_set).scores | sort -u -k1,1 | awk '{print $8}' | sort -nr ) | awk '{print $0}' >"${results}/scores/tab_${name}_${dimer}.tsv"
				# adding score files to the list
				scores+=("${results}/scores/tab_${name}_${dimer}.tsv")
				endnames+=("${name}_${dimer}")
			fi
		fi

		if [[ $matrice == *".keras"* ]]; then # SeqConv scores
			echo "===== SeqConv =====" | tee -a $log
			seqconv_dataset=$results/sequences/seqconv_pos_neg.txt
			awk 'NF >= 5' $results/sequences/seqconv_pos_neg.txt > $results/sequences/seqconv_pos_neg.txt.tmp && mv $results/sequences/seqconv_pos_neg.txt.tmp $results/sequences/seqconv_pos_neg.txt

			mkdir -p -m 774 $results/roc
			echo "[INFO] - Computing SeqConv scores for ${name}" | tee -a $log
			python $SeqConv_predict $seqconv_dataset $peaksize $name $results $matrice $seed
			cp $results/roc/${name}_scores.txt $results/scores

			# Formatting scores for ROC plotting
			awk -v OFS='\t' '$1 == 1.0 {print $6}' $results/scores/${name}_scores.txt | sort -g -r > $results/scores/tmp_pos_${name}.txt
			awk -v OFS='\t' '$1 == 0.0 {print $6}' $results/scores/${name}_scores.txt | sort -g -r > $results/scores/tmp_neg_${name}.txt
			paste $results/scores/tmp_pos_${name}.txt $results/scores/tmp_neg_${name}.txt > $results/scores/SeqConv_${name}_scores_roc.tsv

			rm $results/scores/tmp_pos_${name}.txt $results/scores/tmp_neg_${name}.txt #clean
			rm -r ${results}/roc #seqconv secondary outputs

			scores+=("$results/scores/SeqConv_${name}_scores_roc.tsv")
			endnames+=("${name}_SeqConv")
		fi
	done

	echo ${colors[@]}
	echo "[INFO] - Plotting ROC curves" | tee -a $log

	# Plotting ROC curves
	python3 $plot_ROCS_prog -s ${scores[@]} -n ${endnames[@]} -o ${results} -of ROC.svg -c ${colors[@]} -lim 1
	inkscape --export-type=png -w 1000 -h 1000 ${results}/ROC.svg 2>/dev/null
	rm $results/nKSM_AUC.txt

	echo "[INFO] - FPR limits set to ${fpr_lim}, plotting additional curves" | tee -a $log

	# Subresults directory if FPR_limit is below 1 (partial ROC)
	if (( $(echo "$fpr_lim < 1" | bc -l) )); then
		echo "[INFO] - FPR limit is below 1, creating subdirectory for partial ROC curves" | tee -a $log
		results=${results}/FPR_${fpr_lim}
		mkdir -p -m 774 $results
		python3 $plot_ROCS_prog -s ${scores[@]} -n ${endnames[@]} -o ${results} -of ROC.svg -c ${colors[@]} -lim $fpr_lim
		inkscape --export-type=png -w 1000 -h 1000 ${results}/ROC.svg 2>/dev/null
		rm $results/nKSM_AUC.txt
	fi
}

#-------------------------------------------------------------------------------

compute_distribution(){
	# FUNCTION: compute_distribution
	# DESCRIPTION:
	#   This function calculates and visualizes the distribution of motif scores (TFFM or PWM)
	#   across a genome. It processes each chromosome in parallel to generate scores,
	#   samples a subset of these scores, and creates an interactive distribution plot.
	#   This is particularly useful for analyzing the global occurrence of motifs in a genome 
	#   and determining threshold for downstream analysis.
	#
	# USAGE:
	#   compute_distribution -m <file> -n <STRING> -od <PATH> -g [FILE] -r [INT]
	#         -c [LIST of STRING]
	#
	# ARGUMENTS:
	#   -m          : Path to the matrix file (TFFM in `.xml` format or PWM in `.pfm` format). Required.
	#   -od         : Directory where results will be saved. Required.
	#   -n          : Prefix for output files and directories. Required.
	#   -g          : Path to the genome FASTA file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas".
	#   -r          : Identifier for the interactive plot request. Optional. Default: 0.
	#   -c          : List of chromosomes to analyze. Optional. Default: ("chr1" "chr2" "chr3" "chr4" "chr5").
	#   --all       : Flag to use all sequences for distribution plot without sampling. Optional.
	#   -h, --help  : Display this help message.
	#
	# DEPENDENCIES 
		# - ../bin/scores.py 
		# - get_all_score_tffm.py
		# - Environment with pyhon version 2.7 for TFFM
		# - ../bin/InteractiveDistrib.r
	# NOTES:
	#   - For TFFM matrices, the function computes scores for the entire genome at once.
	#   - For PWM matrices, the function processes each chromosome in parallel to speed up computation.
	#   - The function samples 10% of the scores from each chromosome to create a manageable dataset.
	echo "
	____ ____ _  _ ___  _  _ ___ ____     ___  _ ____ ___ ____ _ ___ 
	|    |  | |\/| |__] |  |  |  |___     |  \ | [__   |  |__/ | |__]
	|___ |__| |  | |    |__|  |  |___ ___ |__/ | ___]  |  |  \ | |__]

	"
	local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; 
	local chrset=("chr1" "chr2" "chr3" "chr4" "chr5"); 
	local request=0
	local allseq=false
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-m)
				local matrice=$2
				echo "-> Matrix files set to: ${matrice}";shift 2;;
			-od)
				local results=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> Names associated set to: ${name}";shift 2;;
			-r)
				local request=$2;shift 2;;
			-c)
				local chrset=("${!2}"); echo ${!2}
				shift 2;;
			--all)
				local allseq=true
				echo "-> All sequences will be used for the distribution plot, no sampling will be done"; shift 1
				;;
			-h)
				usage compute_distribution; return;;
			--help)
				usage compute_distribution; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage compute_distribution; return;;
		esac
	done

	# Log 
	log=$results/$name/log.txt
	echo $Argline > $log

	local Errors=0
	if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; fi

	if [ $Errors -gt 0 ]; then usage compute_space; return 1; fi
	mkdir -p -m 774 $results/${name} $results/${name}/scores

	fainame="${genome}.fai"
	if [ ! -f $fainame ]; then
		#not tested
		samtools faidx $genome
	fi

	if [ $allseq == true ]; then
		local fastatoUse=$genome
	else
		bedtools random -l 500 -n 10000 -g $fainame | sort -k1,1 -k2,2n | bedtools merge -i - > $results/${name}/subseq_genome.bed
		bedtools getfasta -bed $results/${name}/subseq_genome.bed -fi $genome -fo $results/${name}/subseq_genome.fa
		local fastatoUse=$results/${name}/subseq_genome.fa
	fi

	# length of each sequence in the fasta file. then do the sum
	local Sumoflength=$(awk '/^>/{if (seqlen){print seqlen}; seqlen=0; next; } {seqlen+=length($0)} END{print seqlen}' $fastatoUse | awk '{sum+=$1} END{print sum}') 
	echo "Total length of sequences used for score distribution: $Sumoflength"
	if [[ $matrice == *".xml"* ]]; then
		$Python_TFFM $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $genome -t ${matrice}
		local scores=${results}/${name}/scores/${name}_tffm_scores_pos.tsv
	elif [[ $matrice == *".pfm"* ]]; then
		# local scores=${results}/${name}/scores/subseq_genome.fa.scores
		local scores=${results}/${name}/scores/$(basename $fastatoUse).scores
		if [[ ! -f $scores ]] || [[ ${matrice} -nt $scores ]] ; then
			python3 ${scores_prog} -m ${matrice} -f $fastatoUse -o ${results}/${name}/scores/
		fi
	fi
	$PATHRscript $interactivedistrib $scores $results/${name}/ $Sumoflength
	#relpace 120000000 par le nombre de tfbs du fichier fasta.
}	



#-------------------------------------------------------------------------------
compute_space(){
	# FUNCTION: compute_space
	# 
	# DESCRIPTION:
	#   This function computes the spacing between motif occurrences (TFFM or PWM) 
	#   in a set of peaks and visualizes the distribution of enriched conformations
	#   (based on a Z-score calculation).
	# 
	# USAGE:
	#  compute_space -p [LIST of FILE] -m [LIST of FILE] -n [LIST of STRING] -od <PATH> 
	#                -th [LIST of FLOAT] -thf <PATH> -maxy [FLOAT] -miny [FLOAT] -maxs [INT] 
	#                -mins [INT] -ol [INT] -or [INT] -g [FILE] -sym -color [LIST of STRING]
	# 
	# ARGUMENTS:
	#   -p           : List of peak files (BED or FASTA format). Required.
	#   -m           : List of matrix files (TFFM in `.xml` format or PWM in `.pfm` format). Required.
	#   -n           : List of names associated with each peak/matrix pair. Required.
	#   -od          : Directory where results will be saved. Required.
	#   -th          : List of thresholds for motif score to consider an occurrence. Optional. Default: (0).
	#   -thf         : Directory containing threshold files for each matrix. Optional. Default: "null".
	#   -maxy        : Maximum enrichment to display on the plot. Optional. Default: 0 (auto).
	#   -miny        : Minimum enrichment to display on the plot. Optional. Default: 0 (auto).
	#   -maxs        : Maximum spacing to compute. Optional. Default: 50.
	#   -mins        : Minimum spacing to compute. Optional. Default: 0.
	#   -ol          : Offset on the left motif for spacing computation. Optional. Default: 0.
	#   -or          : Offset on the right motif for spacing computation. Optional. Default: 0.
	#   -g           : Path to the genome FASTA file. Optional.
	#   -sym         : Flag to indicate if the matrix is symmetric (i.e., same motif on both sides). 
	#                  Optional. Default: false.
	#   -color        : List of colors for the plot in hex format. Optional.
	#   -h, --help   : Display this help message.
	#
	# DEPENDENCIES:
	#   - ../bin/scores.py or ../bin/get_all_score_tffm.py depending on the matrix type
	#   - ../bin/get_interdistances.py
	#   - ../bin/Zacing.R
	# 
	# NOTES:
	#  - The function checks for the presence of either -th or -thf arguments to determine the thresholds 
	#    for motif occurrences (they can be computed with scores_disctribution function).

    echo '
    ____ ____ _  _ ___  _  _ ___ ____     ____ ___  ____ ____ ____
    |    |  | |\/| |__] |  |  |  |___     [__  |__] |__| |    |___
    |___ |__| |  | |    |__|  |  |___ ___ ___] |    |  | |___ |___
    '
	local Errors=0; local matrix_type="ASYMMETRIC"; local thresholds_dir="null"; local thresholds=(0); local maxy=0; local miny=0; local maxSpace=50; local minSpace=0; local offset_left=0; local offset_right=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544' '#244c63' '#566324' '#4db352' '#4d99b3' '#7c4db3' '#b34da8' '#4db396' '#b34d75' '#FF4000' '#240B3B' '#61380B' '#7A2513' '#7DAAAD' '#F8DB4B' '#B4F84B' '#687574' '#9FE1DB')
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=("${!2}")
				echo "-> Peaks files set to: ${peaks[@]}";shift 2;;
			-m)
				local matrices=("${!2}")
				echo "-> Matrix files set to: ${matrices[@]}";shift 2;;
			-n)
				local names=("${!2}")
				echo "-> Names associated set to: ${names[@]}";shift 2;;
			-th)
				local thresholds=("${!2}")
				echo "-> Thresholds set to: ${thresholds[@]}";shift 2;;
			-thf)
				local thresholds_dir=$2
				echo "-> Thresholds directory set to: ${2}"; shift 2;;
			-od)
				local results=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-maxy)
				local maxy=$2
				echo "-> Maximum enrichment to display set to: ${2}";shift 2;;
			-miny)
				local miny=$2
				echo "-> Minimum enrichment to display set to: ${2}";shift 2;;
			-maxs)
				local maxSpace=${2%.*}
				echo "-> Maximum spacing to compute set to: ${2%.*}";shift 2;;
			-mins)
				local minSpace=${2%.*}
				echo "-> Minimum spacing to compute set to: ${2%.*}";shift 2;;
			-ol)
				local offset_left=${2%.*}
				echo "-> Offset on the left set to: ${2%.*}";shift 2;;
			-or)
				local offset_right=${2%.*}
				echo "-> Offset on the right set to: ${2%.*}";shift 2;;
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-sym)
				local matrix_type="SYMMETRIC"; shift 1;;
			-color)
				local colors=("${!2}")
				echo "-> Colors set to: ${colors[@]}";shift 2;;
			-h)
				usage compute_space ;Errors+=1; shift 1;;
			--help)
				usage compute_space ;Errors+=1; shift 1;;
			*)
				echo "Error in arguments"
				echo $1; usage compute_space; return;;
		esac
	done

    # ----- Checking arguments -----
    if [ -z $peaks ]; then echo "[ERROR] - -p argument needed or may not be a LIST"; Errors+=1; fi
    if [ -z $matrices ]; then echo "[ERROR] - -m argument needed  or may not be a LIST"; Errors+=1; fi
    if [ -z $names ]; then echo "[ERROR] - -n argument needed  or may not be a LIST"; Errors+=1; fi
    if [ -z $results ]; then echo "[ERROR] - -od argument needed"; Errors+=1; fi
    if (( $(echo "${thresholds[0]} == 0.0" | bc -l) )) && [ ${thresholds_dir} == "null" ]; then
        echo "[ERROR] - either -th or -thd arguments needed; -th arguments needs to be a list"
        Errors+=1
    fi
    if [ ${thresholds_dir} != "null" ] && [ ! -d ${thresholds_dir} ]; then
        echo "[ERROR] - directory for thresholds acquisition does not exist, please check your argument -thd"
        Errors+=1
    fi
    if [[ ${#peaks[@]} -ne ${#matrices[@]} ]] || [[ ${#peaks[@]} -ne ${#names[@]} ]]; then
        echo "[ERROR] - LISTs for -p, -n and -m arguments should have the same length"
        Errors+=1
    fi
    if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; fi
    if [ -z $miny ]; then echo "-miny argument not used, no limits"; fi
    if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; fi
    if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; fi
    if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; fi
    if [ -z $offset_right ]; then echo "-or argument not used, no offset"; fi
    if [ -z $genome ]; then
        echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"
    fi
    if [ $Errors -gt 0 ]; then usage compute_space; return 1; fi

	# log
	log=$results/$name/log.txt
	echo $Argline > $log

    # ----- Before launching -----
    for ((i=0;i<${#colors[@]};i++)); do
        [[ ${colors[$i]} =~ ^#.* ]] || colors[$i]="#${colors[$i]}"
    done

    if [ $maxSpace -lt $minSpace ]; then
        echo "[WARNING] - minimum and maximum spacing limits seems mixed up. Correcting..."
        local tmp=$minSpace
        local minSpace=$maxSpace
        local maxSpace=$tmp
    fi
    
    if [ ${matrix_type} == "SYMMETRIC" ]; then
        echo "- Palindromic mode enabled -"
    fi

    echo "[INFO] - Everything looks good, starting computation of spacing between TFBSs"

	# ----- Main loop -----
    for ((i=0;i<${#peaks[@]};i++)); do
        local name=${names[i]// /_}
        local matrice=${matrices[i]}
        local peak=${peaks[i]}
        local length_mat=0

        # ----- Scores computation -----
        echo "[INFO] - Compute TFBSs scores for ${name} under given peaks"
        mkdir -p -m 774 $results/${name} $results/${name}/scores

		if [[ $peak == *"score"* ]];then
			local pos_file=$peak # no preparation needed
		elif [[ $peak != *".fa"* ]]; then
			echo "[INFO] - Transforming bedfile into fasta"
			bedtools getfasta -fi $genome -bed $peak -fo $results/${name}/${name}_pos_set.fa # need to tranform in fasta
			local peak=$results/${name}/${name}_pos_set.fa
		fi
		if [[ $peak == *".fa"* ]]; then #for FASTA sequences only
			if [[ $matrice == *".xml"* ]]; then #if TFFM 
				if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_pos.tsv ];then
					$Python_TFFM $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $peak -t ${matrice}
				fi
                if (( $(echo "${thresholds[0]} == 0.0" |bc -l) )); then
					# Get thresholds from distribution if not set by user
					echo "searching for thresholds in directory provided"
                    thresholds_TF=() 
                    for i in 50 70 90; do 
                        file="scores_99_${i}_info.txt"
                        if [ -f ${thresholds_dir}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
                        if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
                        if [ -f ${thresholds_dir}/${name}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
                        if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
                    done
				else
					thresholds_TF=(${thresholds[@]})
                fi
				# result name (minus .fa)
				tmpresultName=${peak%.fa}

				local pos_file=${results}/${name}/scores/${name}_tffm_scores_pos.tsv

			elif [[ $matrice == *".pfm"* ]]; then #if PWM 
				local length_mat=$(awk -v OFS="[\t ]" '{if(NR==1){if($5=="SIMPLE"){typM="simple"} else {typM="dependency"};count=-1;next};if(typM=="simple"){count++};if(typM=="dependency"){if($1=="DEPENDENCY"){count-=2;exit};count++}}END{print count}' $matrice)
				if [ ! -f ${results}/${name}/scores/$(basename $peak).scores ] || [ ${results}/${name}/scores/$(basename $peak).scores -ot $peak ];then
					python3 ${scores_prog} -m ${matrice} -f $peak -o ${results}/${name}/scores/
				fi

					# searching distribution scores if th are set to 0 
                if (( $(echo "${thresholds[0]} == 0.0" |bc -l) )); then
				    echo "User set thresholds not detected, searching for thresholds in compute_distribution folder"
                    thresholds_TF=() 
					# Get thresholds from distribution if not set by user
                    for i in 50 70 90; do 
                        file="scores_99_${i}_info.txt"
                        if [ -f ${thresholds_dir}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
                        if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
                        if [ -f ${thresholds_dir}/${name}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
                        if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
						echo ${thresholds_TF[@]}
					done
				else
					thresholds_TF=(${thresholds[@]})
                fi 
				local pos_file=${results}/${name}/scores/$(basename $peak).scores
			fi
		fi
		
		# ----- compute spacing conformations in peaks -----
		echo $spacing_mk
		python3 $spacing_mk -o $results/${name}/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th ${thresholds_TF[@]} -ol $offset_left -or $offset_right -lm $length_mat
		
		# ----- compute Z-score and plots -----
		$PATHRscript $Zscore_spacing -f ${results}/${name}/${name}_spacing.tsv -m $matrix_type -th ${thresholds_TF[@]} -od ${results}/${name} --maxy $maxy --miny $miny -c ${colors[@]}
	done

}

#-------------------------------------------------------------------------------
spacing_2TFs () {
	# FUNCTION: spacing_2TFs
	# 
	# DESCRIPTION:
	#   This function computes the spacing between occurrences of two different motifs (TFFM or PWM)
	#   in a set of peaks and visualizes the distribution of enriched conformations (based on a Z-score calculation).
	#
	# USAGE:
	#  spacing_2TFs -p [FILE] -ma [FILE] -mb [FILE} -n [STRING] -od <PATH>
	#                -th [LIST of FLOAT] -thf [PATH] -maxy [FLOAT] -miny [FLOAT] -maxs [INT]
	#                -mins [INT] -ol [INT] -or [INT] -ol2 [INT] -or2 [INT] -g [FILE] -sym -color [LIST of STRING]
	#
	# ARGUMENTS:
	#   -p           : Peak file (BED or FASTA format). Required.
	#   -ma          : The first matrix files (TFFM in `.xml` format or PWM in `.pfm` format). Required.
	#   -mb          : The second matrix files (TFFM in `.xml` format or PWM in `.pfm` format). Required.
	#   -n           : Name associated with the analysis, used as prefix for output files. Required.
	#   -od          : Directory where results will be saved. Required.
	#   -tha         : List of thresholds for motif scores to consider an occurrence for the first TF. Optional. Default: (0) for both TFs.
	#   -thb         : List of thresholds for motif scores to consider an occurrence for the second TF. Optional. Default: (0) for both TFs.
	#   -thfa        : Directory containing threshold files for the first matrix. Optional. Default: "null".
	#   -thfb        : Directory containing threshold files for the second matrix. Optional. Default: "null".
	#   -maxy        : Maximum enrichment to display on the plot. Optional. Default: 0 (auto).
	#   -miny        : Minimum enrichment to display on the plot. Optional. Default: 0 (auto).
	#   -maxs        : Maximum spacing to compute. Optional. Default: 50.
	#   -mins        : Minimum spacing to compute. Optional. Default: 0.
	#   -ol          : Offset on the left motif for spacing computation for the first TF. Optional. Default: 0.
	#   -or          : Offset on the right motif for spacing computation for the first TF. Optional. Default: 0.
	#   -ol2         : Offset on the left motif for spacing computation for the second TF. Optional. Default: 0.
	#   -or2         : Offset on the right motif for spacing computation for the second TF. Optional. Default: 0.
	#   -g           : Path to the genome FASTA file. Optional.
	#   -sym         : Flag to indicate if the matrices are symmetric (i.e., same motif on both sides). Optional. Default: false.
	#   -color        : List of colors for the plot in hex format. Optional.
	#   -h, --help   : Display this help message.
	#
	# DEPENDENCIES:
	#   - ../bin/scores.py or ../bin/get_all_score_tffm.py depending on the matrix type
	#   - ../bin/get_interdistances.py
	#   - ../bin/Zacing_2TFs.R
	#
	# NOTES:
	#  - The function checks for the presence of either -tha/-thb or -thfa/-thfb arguments to determine the thresholds
	#    for motif occurrences (they can be computed with scores_disctribution function).

    echo '
____ ___  ____ ____ _ _  _ ____     ____ _____ ____  ___
[__  |__] |__| |    | |\ | | __     ___|   |   |__  |___ 
___] |    |  | |___ | | \| |__] ___ |___   |   |    ___|
    '

	local Errors=0; local matrix_type="ASYMMETRIC"; local thresholds_dir="null"; local thresholds_a=(0); local thresholds_b=(0); local maxy=0; local maxSpace=50; local minSpace=0; local offset_left=0; local offset_right=0; local offset_left2=0; local offset_right2=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Name of bed files with peaks coordinates: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> Name used as prefix for output directory and files: ${2}";shift 2;;
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-ma)
				local PFMa=$2
				echo "-> PFM file for first TF: ${2}";shift 2;;
			-mb)
				local PFMb=$2
				echo "-> PFM file for second TF: ${2}";shift 2;;
			-od)
				local outdir=$2
				echo "-> Name of out directory: ${2}";shift 2;;
			-tha)
				local thresholds_a=("${!2}")
				echo "-> Array of three PWM score thresholds for first TF: ${thresholds_a[@]}";shift 2;;
			-thb)
				local thresholds_b=("${!2}")
				echo "-> Array of three PWM score thresholds for second TF: ${thresholds_b[@]}";shift 2;;
			-thfa)
				local thresholds_dira=$2
				echo "-> Thresholds directory set to: ${2}"; shift 2;;
			-thfb)
				local thresholds_dirb=$2
				echo "-> Thresholds directory set to: ${2}"; shift 2;;
			-maxs)
				local maxSpace=$2
				echo "-> Maximum spacing to compute set to: ${2}";shift 2;;
			-mins)
				local minSpace=$2
				echo "-> Minimum spacing to compute set to: ${2}";shift 2;;
			-ol)
				local offset_left=$2
				echo "-> 1st matrix offset on the left set to: ${2}";shift 2;;
			-or)
				local offset_right=$2
				echo "-> 1st matrix offset on the right set to: ${2}";shift 2;;
			-ol2)
				local offset_left2=$2
				echo "-> 2nd matrix offset on the left set to: ${2}";shift 2;;
			-or2)
				local offset_right2="$2"
				echo "-> 2nd matrix offset on the right set to: ${2}";shift 2;;
			-sym)
				local matrix_type="SYMMETRIC"; shift 1;;
			-h)
				usage spacing_2TFs ; return;;
			--help)
				usage spacing_2TFs ; return;;
			*)
				echo "Error in arguments"
				echo $1; usage spacing_2TFs; return;;
		esac
	done
	# checking arguments
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed";Errors+=1;fi
	if [ -z $name ]; then echo "ERROR: -n argument needed";Errors+=1;fi
	if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
	if [ -z $PFMa ]; then echo "ERROR: -ma argument needed";Errors+=1;fi
	if [ -z $PFMb ]; then echo "ERROR: -mb argument needed";Errors+=1;fi
	if [ -z $outdir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
	if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; fi
	if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; fi
	if [ $Errors -gt 0 ]; then usage spacing_2TFs ; return 1; fi

	# log
	log=$outdir/log.txt
	echo $Argline > $log
	
	# ----- Thresholds acquisition -----
	echo "[INFO] - Checking thresholds for both TFs"

	if (( $(echo "${thresholds_a[0]} == 0.0" |bc -l) )); then
		echo "[INFO] - User set thresholds TFa not detected, searching for thresholds in Scores_Distribution folder"
		thresholds_TFa=()
		local thresholds_dir=$thresholds_dira
		for i in 50 70 90; do 
			file="scores_99_${i}_info.txt"
			if [ -f ${thresholds_dir}/${file} ]; then thresholds_TFa+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then thresholds_TFa+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
			if [ -f ${thresholds_dir}/${name}/${file} ]; then thresholds_TFa+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then thresholds_TFa+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
		done
	else
		thresholds_TFa=(${thresholds_a[@]})
	fi 

	if (( $(echo "${thresholds_b[0]} == 0.0" |bc -l) )); then
		echo "[INFO] - User set thresholds TFb not detected, searching for thresholds in Scores_Distribution folder"
		thresholds_TFb=()
		local thresholds_dir=$thresholds_dirb
		for i in 50 70 90; do 
			file="scores_99_${i}_info.txt"
			if [ -f ${thresholds_dir}/${file} ]; then thresholds_TFb+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then thresholds_TFb+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
			if [ -f ${thresholds_dir}/${name}/${file} ]; then thresholds_TFb+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then thresholds_TFb+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
		done
	else
		thresholds_TFb=(${thresholds_b[@]})
	fi 

	# ----- Reporting thresholds ------
	echo "[REPORT] - Thresholds listed for TFa:"
	echo ${thresholds_TFa[@]}
	length_TFa=${#thresholds_TFa[@]}
	echo "[REPORT]  -Thresholds listed for TFb:"
    echo ${thresholds_TFb[@]}
	length_TFb=${#thresholds_TFb[@]}

	# Checking thresholds sizes
	if [[ "$length_TFa" -ne "$length_TFb" ]]; then echo "[ERROR] - TFa and TFb thresholds list need to be the same size !"; exit 1; fi
	
	# ----- Preparing scores files ------
	mkdir -p $outdir/scores_TFa $outdir/scores_TFb
	bedtools getfasta -fi $genome -bed  $peaks > $outdir/${name}.fasta #get peaks fasta for scores 
	peaks=$outdir/${name}.fasta

	# ------ compute scores for both TFs ------
	echo "[INFO] - Computing scores for both TFs on given peaks"
	if [ ! -f $outdir/scores_TFa/${name}.fasta.scores ]; then 
		if [[ $PFMa == *".xml"* ]]; then
			$Python_TFFM $tffm_all_scores -o $outdir/scores_TFa/${name}.fasta.scores -pos $peaks -t ${PFMa}
		elif [[ $PFMa == *".pfm"* ]]; then
			$Python_TFFM $scores_prog -m $PFMa -f $peaks -o $outdir/scores_TFa; 
		fi
	fi

	if [ ! -f $outdir/scores_TFb/${name}.fasta.scores ]; then 
		if [[ $PFMb == *".xml"* ]]; then
			$Python_TFFM $tffm_all_scores -o $outdir/scores_TFb/${name}.fasta.scores -pos $peaks -t ${PFMb}
		elif [[ $PFMb == *".pfm"* ]]; then
			$Python_TFFM $scores_prog -m $PFMb -f $peaks -o $outdir/scores_TFb; 
		fi
	fi

	# clean 
	rm $peaks

	# ----- compute spacings -----
	mkdir -p -m 774 $outdir/${name}

	echo "[INFO] - Checking motif lengths"
	lm_a=$(cat $outdir/scores_TFa/${name}.fasta.scores | cut -f 7 | head -1)
	printf "\nLength of motif $lm_a\n"
	lm_b=$(cat $outdir/scores_TFb/${name}.fasta.scores | cut -f 7 | head -1)
	printf "\nLength of motif $lm_b\n"

	echo "[INFO] - Get interdistance for spacing 2 TFs"
	python3 $get_interdistances_2TF -pos $outdir/scores_TFa/${name}.fasta.scores -ol=$offset_left -or=$offset_right -ol2=$offset_left2 --offset_right2TF=$offset_right2 -neg $outdir/scores_TFa/${name}.fasta.scores -pos2 $outdir/scores_TFb/${name}.fasta.scores -neg2 $outdir/scores_TFb/${name}.fasta.scores -th ${thresholds_TFa[@]} -th2 ${thresholds_TFb[@]} -smin $minSpace -smax $maxSpace -lm $lm_a -lm2 $lm_b -wi -o $outdir/${name}/out

	echo "[INFO] - Preparing spacing_pos.tsv"
	sed "1d" $outdir/${name}/out_spacing_pos21.tsv | awk -v OFS="\t" '{print $1,$2,$4,$3,$6,$5,$8,$7}' > $outdir/${name}/out_spacing_pos21_converted.tsv
	cat $outdir/${name}/out_spacing_pos12.tsv $outdir/${name}/out_spacing_pos21_converted.tsv |sed "1d" | sed "s/-/:/" | sed "s/:/\t/g" | sed "s/_/\t/" | awk -v OFS="\t" '{print $0,$3-$2+1}' | sort -k1,1 -k2,2n -k8,8n | uniq | sed "1ichr\tStart\tEnd\tConf\tSpace\tScore1\tScore2\tmatricePosition1\tmatricePosition2\tcorrectedPosition1\tcorrectedPosition2\tSize" > $outdir/${name}/out_spacing_pos.tsv
	
	# ------ compute Z-score and plots ------
	echo "[INFO] - Launching Zacing_2TFs.R to compute statistics and plots"
	$PATHRscript $Zacing_2TFs $outdir/${name}/out_spacing_pos.tsv $matrix_type ${thresholds_TFa[@]} ${thresholds_TFb[@]} $outdir/${name} $length_TFa $length_TFb
	awk -v OFS="\t" 'NR!=1{
	if($11>$10){print $1,$2+$10,$2+$11,$4$5"|"$6"|"$7,$6+$7,"+"}else{print $1,$2+$11,$2+$10,$4$5"|"$6"|"$7,$6+$7,"-"
	}
	}' $outdir/${name}/out_spacing_pos.tsv > $outdir/${name}/out_spacing_pos.bed
	bedtools getfasta -name+ -fi $genome -bed $outdir/${name}/out_spacing_pos.bed -fo $outdir/${name}/out_spacing_pos.fas

	# ----- compute PFM for enriched conformations -----
	echo "[INFO] - Computing PFM for enriched conformations"
	nbCandidatConfo=$(wc -l $outdir/${name}/candidates.txt | awk '{print $1}') #extract enriched conformations from Zacing_2TFs output
	if [ $nbCandidatConfo -gt 0 ]; then
		while read f; do
			echo $f
			trueconfo=$(echo $f | grep -o "[a-Z]*")
			dist=$(echo $f | grep -o "[0-9]*")

			echo "[INFO] - Computing PFM for conformation ${trueconfo} with distance ${dist}"
			bash $builder_2matrix $PFMa $offset_left $offset_right $PFMb $offset_left2 $offset_right2 $f $trueconfo ${dist} > $outdir/${name}/$f.pssm;

			tail -n +3 $outdir/${name}/$f.pssm > ${outdir}/${name}/pfm2transpose.tsv
			if [ -f ${outdir}/${name}/${f}.pfm ]; then rm ${outdir}/${name}/${f}.pfm; fi
			col="$(head -1 ${outdir}/${name}/pfm2transpose.tsv | wc -w)"
			for i in $(seq 1 $col); do
				awk '{ print $'$i' }' ${outdir}/${name}/pfm2transpose.tsv | paste -s -d "\t" >> ${outdir}/${name}/${f}.pfm
			done
			rm ${outdir}/${name}/pfm2transpose.tsv
			
		done <$outdir/${name}/candidates.txt

		echo "[INFO] - Converting PFM to MEME format and generating logos for enriched conformations"
		/home/prog/meme/meme_4.12.0/bin/jaspar2meme -pfm ${outdir}/${name} > ${outdir}/${name}/memePFM.txt
		rm $outdir/${name}/*.pfm
		/home/prog/meme/meme_4.12.0/bin/meme2images -png -rc ${outdir}/${name}/memePFM.txt ${outdir}/${name}
		rm ${outdir}/${name}/memePFM.txt
		while read f; do
			mv ${outdir}/${name}/$f.pssm ${outdir}/${name}/$f.pfm
		done <$outdir/${name}/candidates.txt
	else
		echo "[REPORT]- No enrichment found"
	fi
}

#-------------------------------------------------------------------------------
spacing_TFBSTSS () { 
	# FUNCTION: spacing_TFBSTSS
	#
	# DESCRIPTION:
	#   This function computes the spacing between occurrences of a motif (TFFM or PWM) and TSSs 
	#   in a set of peaks and visualizes the distribution of enriched conformations (based on a Z-score calculation). 
	#
	# USAGE:
	#  spacing_TFBSTSS -p [FILE] -m [FILE] -n [STRING] -od <PATH>
	#                -thtf [LIST of FLOAT] -thf [PATH] -thtss [LIST of FLOAT] -maxy [FLOAT] -miny [FLOAT] -maxs [INT]
	#                -mins [INT] -ol [INT] -or [INT] -g [FILE] -sym -color [LIST of STRING]
	#                -w [INT] -f [INT]
	#
	# ARGUMENTS:
	#   -p           : Peak file (BED or FASTA format). Required.
	#   -ma          : Matrix file (TFFM in `.xml` format or PWM in `.pfm` format). Required.
	#   -n           : Name associated with the analysis, used as prefix for output files. Required.
	#   -od          : Directory where results will be saved. Required.
	#   -thtf        : List of thresholds for motif scores to consider an occurrence for TFBS. Optional. Default: (0).
	#   -thf         : Directory containing threshold files for the TFBS. Optional. Default: "null".
	#   -thtss       : List of thresholds for motif scores to consider an occurrence for TSS. Optional. Default: (-2) for all three thresholds.
	#   -maxy        : Maximum enrichment to display on the plot. Optional. Default: 0 (auto).
	#   -miny        : Minimum enrichment to display on the plot. Optional. Default: 0 (auto).
	#   -maxs        : Maximum spacing to compute. Optional. Default: 1000.
	#   -mins        : Minimum spacing to compute. Optional. Default: 0.
	#   -ol          : Offset on the motif for spacing computation. Optional. Default: 0.
	#   -or          : Offset on the motif for spacing computation. Optional. Default: 0.
	#   -g           : Path to the genome FASTA file. Optional.
	#   -ga          : Path to the annotated genome GFF file. Optional.
	#   -w           : Window size to count orientations. Optional. Default: 10.
	#   -f           : Filter to consider only peaks at proximity from TSSs in bp range. Optional. Default: 1000.
	#   -sym		 : Flag to indicate if the matrix is symmetric. Optional. Default: false.
	# 
	# DEPENDENCIES:
	#   - ../bin/scores.py or ../bin/get_all_score_tffm.py depending on the matrix type
	#   - ../bin/get_interdistances.py
	#   - ../bin/Zacing_TFBSTSS.R
	#
	# NOTES:
	#  - The function checks for the presence of either -thtf or -thf arguments to determine the thresholds
	#    for motif occurrences (they can be computed with scores_disctribution function).
	#  - TSS positions are determined based on the annotated genome GFF file, and the closest TSS to each peak is 
    #    considered for spacing computation.
	#  - Thresholds for TSS (thtss) ar false thresholds only used for distance s=computations and not 
	#    for filtering TSSs, as all TSSs are considered in the analysis.

	echo '
	____ ___  ____ ____ _ _  _ ____     ___ ____ ___  ____    ___ ____ ____
	[__  |__] |__| |    | |\ | | __      |  |___ |__] [__  __  |  [__  [__ 
	___] |    |  | |___ | | \| |__] ___  |  |    |__] ___]     |  ___] ___]
	'

	local Errors=0; local matrix_type="ASYMMETRIC"; local filter=1000; local windows=10; local thresholds_dir="null"; local thresholds=(0); local thresholds_TSS=("-2" "-2" "-2");local maxy=0; local maxSpace=1000; local minSpace=0; local offset_left=0; local offset_right=0; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome_a="/home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Phytozome/PhytozomeV12/Athaliana/annotation/Athaliana_167_TAIR10.gene.gff3"
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Name of BED files with peaks coordinates: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> A name used as prefix for output directory and files: ${2}";shift 2;;
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-ga)
				local genome_a=$2
				echo "-> Annotated GFF genome set to: ${2}";shift 2;;
			-ma)
				local PFM=$2
				echo "-> PFM file for TF: ${2}";shift 2;;
			-f) 
				local filter=$2
				echo "-> Filter only peaks at proximity from TSSs in bp range of:";shift 2;;
			-w) 
				local windows=$2
				echo "-> Windows size to counts spacings set to:";shift 2;;
			-od)
				local outdir=$2
				echo "-> Name of out directory: ${2}";shift 2;;
			-thtf)
				local thresholds=("${!2}")
				echo "-> An array of three TF PWM score thresholds: ${thresholds_TF[@]}";shift 2;;
			-thtss)
				local thresholds_TSS=("${!2}")
				echo "-> An array of three false PWM score thresholds for TSS: ${thresholds_TSS[@]}";shift 2;;
			-thf)
				local thresholds_dir=$2
				echo "-> Thresholds directory set to: ${2}"; shift 2;;
			-maxs)
				local maxSpace=${2%.*}
				echo "-> Maximum spacing to compute set to: ${2}";shift 2;;
			-mins)
				local minSpace=${2%.*}
				echo "-> Minimum spacing to compute set to: ${2}";shift 2;;
			-ol)
				local offset_left=${2%.*}
				echo "-> Matrix offset on the left set to: ${2}";shift 2;;
			-or)
				local offset_right=${2%.*}
				echo "-> Matrix offset on the right set to: ${2}";shift 2;;
			-sym)
				local matrix_type="SYMMETRIC"; shift 1;;
			-h)
				usage  ; return;;
			--help)
				usage  ; return;;
			*)
				echo "Error in arguments"
				echo $1; usage spacing_TFBSTSS ; return;;
		esac
	done

	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed";Errors+=1;fi
	if [ -z $name ]; then echo "ERROR: -n argument needed";Errors+=1;fi
	if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
	if [ -z $genome_a ]; then echo "ERROR: -ga argument needed";Errors+=1;fi
	if [ -z $PFM ]; then echo "ERROR: -ma argument needed";Errors+=1;fi
	if [ -z $outdir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
	if [ -z $thresholds ]; then echo "ERROR: -thtf argument needed or needs to be a LIST";Errors+=1;fi
	if [ -z $thresholds_TSS ]; then echo "ERROR: -thtss argument needed or needs to be a LIST";Errors+=1;fi
	if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; fi
	if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; fi
	if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; fi
	if [ -z $offset_right ]; then echo "-or argument not used, no offset"; fi
	if [ $Errors -gt 0 ]; then usage spacing_TFBS2TSS ; return 1; fi

	# log
	log=$outdir/log.txt
	echo $Argline > $log

	# ----- Preparing scores files -----
	# Get thresholfds for TFBS from distribution if not set by user
	echo "[INFO] - User TFBS set thresholds not detected, searching for thresholds in Scores_Distribution folder"
	if (( $(echo "${thresholds[0]} == 0.0" |bc -l) )); then
		thresholds_TF=() 
		for i in 50 70 90; do 
			file="scores_99_${i}_info.txt"
			if [ -f ${thresholds_dir}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
			if [ -f ${thresholds_dir}/${name}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then thresholds_TF+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
		done
	else thresholds_TF=(${thresholds[@]})
	fi 

	length_thtf=${#thresholds_TF[@]}
	length_thtss=${#thresholds_TSS[@]}

	echo "[REPORT] - Thresholds listed for TFBS:"
	echo ${thresholds_TF[@]}
	echo "[REPORT] - Thresholds listed for TSS:"
	echo ${thresholds_TSS[@]}
	
	echo "[INFO] - Filtering peaks at proximity of TSSs"

	# ----- Filtering peaks -----
	field_nb=$(head -n1 $peaks | awk '{print NF}') #counting number of columns and get only 3 first fields 
	if [[ $field_nb > 3 ]]; then awk -v OFS="\t" '{print $1, $2, $3}' $peaks > temp_peaks && mv temp_peaks $peaks ; fi 
	# Extarting TSS + converting them into peaks
	awk -v OFS="\t" '$3 == "gene" && ($7 == "+" || $7 == "-") {print $1, ($7 == "+" ? $4 : $5), $7}' $genome_a | awk -v OFS="\t" '$1 != "ChrM" && $1 != "ChrC" {print $1, $2, $3}' | awk -v OFS="\t" '{print tolower($1),$2" "$2,$3}' | awk -v s=$filter -v OFS="\t" '{print $1, $2-s, $3+s, $4}' > $outdir/${name}/TSS_into_peaks.bed  
	# Peaks and TSS peaks intersect 
	bedtools intersect -a $peaks -b $outdir/${name}/TSS_into_peaks.bed -f 1 -wo > $outdir/${name}/peaks_filtered.bed 	

	# ----- Preparing scores files for TSSs -----
	mkdir -p $outdir/${name}/scores_TF
    mkdir -p $outdir/${name}/scores_TSS

	echo "[INFO] - Preparing scores files for TSSs"

	#determine TSS position within peaks (substracting filter) + fake score header
	awk -v s=$filter -v OFS="\t" '{print $1":"$2"-"$3,($6-s)-$2,($6-s)-$2, $7}' $outdir/${name}/peaks_filtered.bed > $outdir/${name}/TSS_score_1.bed 		
	#searching for TSS first base in with fasta file with tmp fasta TSS coordonees 
	awk -v s=$filter -v OFS="\t" '{print $1,$6-s,($6-s)+1}' $outdir/${name}/peaks_filtered.bed > $outdir/${name}/tmp.bed 
	bedtools getfasta -fi $genome -bed $outdir/${name}/tmp.bed > $outdir/${name}/tmp2.bed 
	#extracting from FASTA nucleotides and add last fields of fake scores file in a TSS_score_2.bed 
	awk -v OFS="\t" '$1~ /^[^>]/ {print $1}' $outdir/${name}/tmp2.bed | awk -v OFS="\t" '{print $1, $2="PWM",$3=1,$4=-1}'> $outdir/${name}/TSS_score_2.bed 
	# assembling
	paste $outdir/${name}/TSS_score_1.bed $outdir/${name}/TSS_score_2.bed > $outdir/${name}/scores_TSS/${name}.fasta.scores 
	# Example : 
	# chr1:37969-38169	-98	-98	-	A	PWM	1	-1
	rm $outdir/${name}/tmp.bed $outdir/${name}/tmp2.bed $outdir/${name}/TSS_score_1.bed $outdir/${name}/TSS_score_2.bed $outdir/${name}/peaks_filtered.bed #cleaning

	# ----- prepare peaks FASTA files for scoring -----
	bedtools getfasta -fi $genome -bed  $peaks > $outdir/${name}/${name}.fasta
	peaks=$outdir/${name}/${name}.fasta
	
	# ----- Preparing PFM files for spacings computations -----
	# If SYMMETRIC matrix = faking a asymmetric one (easier with Z-scores)
	# => modifying .pfm file header from "SYMMETRIC" to "ASYMMETRIC" and use this file for scoring 
	echo "[INFO] - Computing scores for peaks for spacing computation"

	if [[ "$matrix_type" == "SYMMETRIC" ]]; then
		cat $PFM | sed 's/SYMMETRIC/ASYMMETRIC/' > $outdir/${name}/${name}_asymmetric.pfm
		echo "[INFO] - Your matrix is SYMMETRIC, faking ASYMMETRIC"
		$Python_TFFM $scores_prog -m $outdir/${name}/${name}_asymmetric.pfm -f $peaks -o $outdir/${name}/scores_TF
	else
		$Python_TFFM $scores_prog -m $PFM -f $peaks -o $outdir/${name}/scores_TF
	fi	
	
	mkdir -p -m 774 $outdir/${name}

	# ----- Compute spacings -----
	# get motif length 
	lm_TF=$(cat $outdir/${name}/scores_TF/${name}.fasta.scores | head -1 | cut -f 7 || true)
	printf "\n[INFO] - Length of TF motif $lm_TF\n"
	lm_TSS=$(cat $outdir/${name}/scores_TSS/${name}.fasta.scores | head -1 | cut -f 7  || true)
	printf "[INFO] - Length of fake TSS motif $lm_TSS\n"

	echo "[INFO] - Get interdistance for spacing to TSS..." 
	# setting smax ans smin value to a fixed one (filter around TSS here) : getting same enrichment for different zooming windows
	python3 $get_interdistances_TFBSTSS -pos $outdir/${name}/scores_TF/${name}.fasta.scores -neg $outdir/${name}/scores_TF/${name}.fasta.scores -pos2 $outdir/${name}/scores_TSS/${name}.fasta.scores -neg2 $outdir/${name}/scores_TSS/${name}.fasta.scores -ol $offset_left -or $offset_right -th ${thresholds_TF[@]} -smin 0 -smax $filter -lm $lm_TF -lm2 $lm_TSS -wi -o $outdir/${name}/out

	sed "1d" $outdir/${name}/out_spacing_pos21.tsv | awk -v OFS="\t" '{print $1,$2,$4,$3,$6,$5,$8,$7}' > $outdir/${name}/out_spacing_pos21_converted.tsv
	cat $outdir/${name}/out_spacing_pos12.tsv $outdir/${name}/out_spacing_pos21_converted.tsv |sed "1d" | sed "s/-/:/" | sed "s/:/\t/g" | sed "s/_/\t/" | awk -v OFS="\t" '{print $0,$3-$2+1}' | sed "1ichr\tStart\tEnd\tConf\tSpace\tScore1\tScore2\tmatricePosition1\tmatricePosition2\tcorrectedPosition1\tcorrectedPosition2\tSize" > $outdir/${name}/out_spacing_pos.tsv
	
	# ========================== CONFIGURATION COUNTS IN A SLIDING WINDOWS  =============================
	# Context : proportion test are performed on spacing_TSS result, but it is not sufficient to conclude 
	# anything about orientation enrichiment at specific distances. (Only taking global orientation  into 
	# account). From interdistnace result : pool every conformations met in a bp range of distance (10 for
	# example). (Relevant only for asymmetric matrices).

	if [[ "$matrix_type" == "ASYMMETRIC" ]]; then
		# sort file and defining an header to new file storing conformations counts
		awk 'NR == 1; NR > 1 {print $0 | "sort -k5 -n"}' $outdir/${name}/out_spacing_pos.tsv > $outdir/${name}/out_spacing_pos_sorted.tsv
		echo -e "start\tend\tcount_DRba\tcount_IR" > $outdir/${name}/counts_after_tss.tsv
		echo -e "start\tend\tcount_DRab\tcount_ER" > $outdir/${name}/counts_before_tss.tsv
		
		# sliding windows splitting into range 
		for ((i=0; i<=$filter; i++)); do
			start=$i
			end=$((i+$windows-1))
			count_DRab=0
			count_DRba=0
			count_ER=0
			count_IR=0

			awk -v OFS="\t" -v start="$start" -v end="$end" '
				NR > 1 && $5 >= start && $5 <= end { 
					if ($4 == "DR21") {
						count_DRba++
					}
					if ($4 == "IR") {
						count_IR++
					}
				}
				END {
					print start, end, count_DRba, count_IR
				}
			' $outdir/${name}/out_spacing_pos_sorted.tsv >> $outdir/${name}/counts_after_tss.tsv
			# before TSS = convert into negative distances
			awk -v OFS="\t" -v start_neg="$start" -v end_neg="$end" '
				NR > 1 && $5 >= start_neg && $5 <= end_neg { 
					if ($4 == "DR12") {
						count_DRab++
					}
					if ($4 == "ER") {
						count_ER++
					}
				}
				END {
					print -start_neg, -end_neg, count_DRab, count_ER
				}
			' $outdir/${name}/out_spacing_pos_sorted.tsv >> $outdir/${name}/counts_before_tss.tsv
		done
		
		# cleaning 
		rm $outdir/${name}/out_spacing_pos_sorted.tsv
	fi 

	echo $matrix_type ${thresholds_TF[@]} $outdir/${name} 

	# ----- Launching R script for statistics and plots -----
	echo "[INFO] - Launching Zacing_TFBSTSS.R to compute Z-score statistics and plots"
	$PATHRscript $Zacing_2TSS -f $outdir/${name}/out_spacing_pos.tsv -m $matrix_type -thtf ${thresholds_TF[@]} -thtss ${thresholds_TSS[@]} -n $outdir/${name} -len_thtf $length_thtf -len_thtss $length_thtss -max $maxSpace
}

#-------------------------------------------------------------------------------
add_coverage(){
	# TODO deprecated ? 
	norm="NA"
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-b)
				local bedgraphs=("${!2}")
				echo "-> bedgraphs to add coverage to the table: ${bedgraphs[@]}";shift 2;;
			-n)
				local names=("${!2}")
				echo "-> names of added exp to: ${names[@]}";shift 2;;
			-t)
				local table=$2
				echo "-> table ouput by initial_comp set to: ${2}";shift 2;;
			-od)
				local out_dir=$2
				echo "-> name of out directory set to: ${2}";shift 2;;
			-norm)
				local norm=$2
				echo "-> normalisation choosen: ${2}"; shift 2;;
			-h)
				usage add_coverage; return;;
			--help)
				usage add_coverage; return;;
			*)
				echo "Error in arguments"
				echo $1; usage add_coverage; return;;
		esac
	done

	local Errors=0
	if [ -z $bedgraphs ]; then echo "ERROR: -b argument needed or needs to be a LIST";Errors+=1;fi
	if [ -z $names ]; then echo "ERROR: -n argument needed or needs to be a LIST";Errors+=1;fi
	if [ -z $table ]; then echo "ERROR: -t argument needed";Errors+=1;fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi

	if [ ${#names[@]} -ne ${#bedgraphs[@]} ]; then
		echo "ERROR: -n & -b have to be lists of same length"; Errors+=1
	fi

	if [ $Errors -gt 0 ]; then usage add_coverage; return 1; fi
	# no need for mkdir as this function is called after pairwise_comp which already created the out_dir
	local log=$out_dir/log.txt
	awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3,$4}NR>1{print $1,$2,$3,$4}' $table | sort -k1,1 -k2,2n > $out_dir/tmp_peaks.bed
	local peak_file=$out_dir/tmp_peaks.bed
	local tmp_table=$out_dir/tmp_table.bed
	cp $table $tmp_table

	local list_files=()
	local i=0

	for bdg in ${bedgraphs[@]};
	do
		echo $bdg
		echo ${names[$i]}
		bedtools intersect -a $peak_file -b $bdg -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/${names[$i]}.inter
		python3 $compute_coverage -i $out_dir/${names[$i]}.inter -m
		if [[ $norm == "RPKM" ]]; then
			echo "RPKM used"
			awk  '{print (1000*$5)/($3-$2)}' $out_dir/${names[$i]}.inter.cov | sed "1i${names[$i]}" > $out_dir/${names[$i]}.inter.cov.normed
		elif [[ $norm == "CPM" ]]; then
			echo "CPM used"
			awk  '{print $5/($3-$2)}' $out_dir/${names[$i]}.inter.cov | sed "1i${names[$i]}" > $out_dir/${names[$i]}.inter.cov.normed
		else
			echo "no norm"
			awk  '{print $5}' $out_dir/${names[$i]}.inter.cov | sed "1i${names[$i]}" > $out_dir/${names[$i]}.inter.cov.normed
		fi
		list_files+=("$out_dir/${names[$i]}.inter.cov.normed")
		local i=$i+1
	done
	# local files=`join_by " " "${list_files[@]}"`
	#paste $tmp_table $files > $table
	paste $tmp_table ${list_files[@]} > $out_dir/table.tsv
	rm $peak_file $tmp_table
	rm $out_dir/*.inter.cov.normed
	rm $out_dir/*.inter.cov
	rm $out_dir/*.inter

}

#-------------------------------------------------------------------------------
add_score(){
	# TODO deprecated ?
	local pocc=false
	while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-m)
				local matrices=("${!2}")
				echo "matrice to computes scores column(s) set to: ${matrices[@]}";shift 2;;
			-n)
				local name=("${!2}")
				echo "name for score column(s) set to: ${name[@]}";shift 2;;
			-t)
				local table=$2
				echo "table set to: ${2}";shift 2;;
			-od)
				local out_dir=$2
				echo "output directory set to: ${2}";shift 2;;
			-po)
				local pocc=true
				echo "pocc mode activated";shift 2;;
			-h)
				usage add_coverage; return;;
			--help)
				usage add_coverage; return;;
			*)
				echo "Error in arguments"
				echo $1; usage add_coverage; return;;
		esac
	done

	local Errors=0
	if [ -z $matrices ]; then echo "ERROR: -m argument needed or needs to be a LIST";Errors+=1;fi
	if [ -z $name ]; then echo "ERROR: -n argument needed or needs to be a LIST";Errors+=1;fi
	if [ -z $table ]; then echo "ERROR: -t argument needed";Errors+=1;fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi

	if [ ${#name[@]} -ne ${#matrices[@]} ]; then
		echo "ERROR: -n & -m have to be lists of same length"; Errors+=1
	fi

	if [ $Errors -gt 0 ]; then usage add_coverage; return 1; fi
	mkdir -p $out_dir
	# awk -v OFS="\t" 'NR!=1{print $1, int($2+(($3-$2)/2)-25), int($2+(($3-$2)/2)+25), $4}' $table > $out_dir/tmp_peaks.bed
	local peak_file=$out_dir/tmp_peaks.bed
	local tmp_table=$out_dir/tmp_table.bed
	local peak_fasta=$out_dir/tmp_fasta.fa
	# sort -k1,1 -k2,2n $table > $tmp_table
	cat $table > $tmp_table

	if [[ $matrix == *".xml"* ]]; then
		awk -v OFS="\t" 'NR!=1{print $1, $2-2, $3, $4}' $table > $out_dir/tmp_peaks.bed
	else
		awk -v OFS="\t" 'NR!=1{print $1, $2, $3, $4}' $table > $out_dir/tmp_peaks.bed
	fi
	bedtools getfasta -fi /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas -bed $peak_file -fo $peak_fasta

	local list_files=()
	local i=0
	for matrix in ${matrices[@]};
	do
		echo $matrix
		if [[ $matrix == *".xml"* ]]; then # TFFM scores computation
			$Python_TFFM $tffmscores -o ${out_dir}/tffm_scores.tsv -pos $peak_fasta -t ${matrix}
			sed 's/-/:/' ${out_dir}/tffm_scores.tsv | awk -v OFS="\t" -v FS="[:\t]" '$1!="None"{print $1,$2,$2+$4-2,$2+$5-1,$10;next}{print 0,0,0,0,0}' | awk -v OFS="\t" '{print $3,$4,$5}' | sed "1istartBS_${name[$i]}\tstopBS_${name[$i]}\t${name[$i]}" > ${out_dir}/${name[$i]}_scores.tsv 
			list_files+=("${out_dir}/${name[$i]}_scores.tsv")
		fi
		if [[ $matrix == *".txt"* ]]; then # K-mer scores computation
			echo "WIP"
		fi
		if [[ $matrix == *".pfm"* ]]; then  # PWM scores
			python3 $scores_prog -m ${matrix} -f $peak_fasta -o ${out_dir}/
	# 		echo "${name[$i]}" > "${out_dir}/tab_${name[$i]}.tsv"
			sort -u -k1,1 -k8,8nr ${out_dir}/$(basename $peak_fasta).scores | sort -u -k1,1 | sed 's/-/:/' | awk -v OFS="\t" -v FS="[:\t]" '{print $1,$2,$2+$4-1,$2+$5-1,$10}' | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $3,$4,$5}' | sed "1istartBS_${name[$i]}\tstopBS_${name[$i]}\t${name[$i]}" > "${out_dir}/tab_${name[$i]}.tsv"
			list_files+=("${out_dir}/tab_${name[$i]}.tsv")
			
			# NOTE Use this to add type (exon, prom etc) to your table
	# 		sort -u -k1,1 -k8,8nr ${out_dir}/$(basename $peak_fasta).scores | sort -u -k1,1 | sed 's/-/:/' | awk -v OFS="\t" -v FS="[:\t]" '{print $1,$2+$4,$2+$5,$1":"$2"-"$3,$10}' > ${out_dir}/tmp_${name[$i]}.tsv
	# 		bedtools intersect -a ${out_dir}/tmp_${name[$i]}.tsv -b /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.bed -wao |sort -u -k4,4 -k10,10nr | sort -u -k4,4 | sed 's/-/:/' | awk -v FS="[:\t]" -v OFS="\t" '{print $4,$5,$6,$7,$11,$2,$3}' | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $5,$6,$7}' | sed "1itype\tstartBS\tstopBS" > ${out_dir}/tmp_${name[$i]}2.tsv
	# 		list_files+=("${out_dir}/tmp_${name[$i]}2.tsv")
			
			if $pocc ; then # PWM Pocc
				python3 $pocc_pfm -s ${out_dir}/$(basename $peak_fasta).scores -o ${out_dir}/Pocc.pocc
				sort -nr ${out_dir}/Pocc.pocc > ${out_dir}/${name[$i]}_Pocc_scores.tsv
				list_files+=("${out_dir}${name[$i]}_Pocc_scores.tsv")
			fi
		fi
		local i=${i}+1
	done
	local files=$(join_by " " "${list_files[@]}")
	paste $tmp_table ${list_files[@]} > $out_dir/table.bed

	# rm $peak_file $tmp_table

}

#-------------------------------------------------------------------------------
cooking_meth(){
    # FUNCTION: With coverage comaprison data from DAP-seq and ampDAP-seq asscoiated with DNA-methylation map, 
	# this function compute the impact of DNA-methylation on TF binding. It is scale on three levels: 
	# 	- At the peaks data set level : compute the correlation between coverage fold change (DAP/ampDAP) 
	# 	  and methylation level (proportion of methylated cytosines in the peak)
	# 	- At the peaks level : compute the correlation between coverage fold change (DAP/ampDAP) 
	# 	  and methylation level (proportion of methylated cytosines in the peak)
	# 	- At the TFBS level : compute the correlation between coverage fold change (DAP/ampDAP) and methylation
	# 	  level (proportion of methylated cytosines in the TFBS)
    # 
	# USAGE: cooking_methylation -g [FILE] -p <FILE> -m <FILE> -sym [yes/no] -od <PATH> -c [INT] -c2 [INT] -thf [PATH]
	# 
	# ARGUMENTS: 
	# 	-g: 	FASTA of the genome. Default: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas
	#   -p:     Tab-seaparated file for bound region with chr, start, end, covDAP (DAP condition coverage)
	#           , covAMP (ampDAP condition coverage). Note that you can also use pairwize_compariosn function output 
	# 			file as an input (`*_RiL_RiP.tsv` file) if you have previously run it with DAP and ampDAP coverage files.
	# 	-m:     PFM matrix file of your TF of interest. 
    #   -n:     Name for output architecture. (Used for searching threshold in Scores_Distribution folder and for 
    #           naming output files). 
	#   -sym:   Either `yes` or `no` to specify if the PFM is symmetric or not. Default: no (asymmetric).
	#   -o:    Output directory where results will be stored.
	#   -thf:   Threshold to use for TFBS prediction. If not set, the function will search for it in the 
	# 			Scores_Distribution folder. Default: "null".
	# 	-c:     Manuall TFBS detection cutoff if you want to use a specific one instead of the one from Scores_Distribution
    #           folder. Used during PFM search under peaks. Default: 0 (not set).
	#   -c2:    Manuall TFBS detection cutoff if you want to use a specific one instead of the one from Scores_Distribution
	#           folder. Used during R statistics. Default: (0) (not set).

	echo '
	____ ____ ____ _  _ _ _  _ ____     _  _ ____ ___ _  _
	|    |  | |  | |_/  | |\ | | __     |\/| |___  |  |__|
	|___ |__| |__| | \_ | | \| |__] ___ |  | |___  |  |  |
	'

	local Errors=0; local thresholds_dir="null"; local cutoff=0; local cutoff2=(0)
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-p)
				local peaksCov=$2
				echo "-> Tab-seaparated file for bound region with chr, start, end, covDAP, covAMP: ${2}";shift 2;;
			-m)
				local matrix_pfm=$2
				echo "-> TF PFM matrix: ${2}";shift 2;;
			-l)
				local motifLength=$2
				echo "-> Motif length set to: ${2}";shift 2;;
			-c)
				local cutoff=$2
				echo "-> PFM score cutoff set to (used during TFBS detection under peaks): ${2}";shift 2;;
			-c2)
				local cutoff2=$2
				echo "-> PFM score cutoff set to (used during TFBS statistics): ${2}";shift 2;;
			-thf)
				local thresholds_dir=$2
				echo "-> Thresholds directory set to: ${2}"; shift 2;;
			-n)
				local name=$2
				echo "-> Sample name set to: ${2}"; shift 2;; 
			-sym)
				local symmetry=$2
				echo "-> [yes] or [no] to specifies wether or not the motif is symmetric: ${2}";shift 2;;
			-o)
				local outdir=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-h)
				usage cooking_meth; return;;
			--help)
				usage cooking_meth; return;;
			*)
				echo "Error in arguments"
				echo $1; usage cooking_meth; return;;
		esac
	done
	
	# Check arguments
	if [ -z $peaksCov ]; then echo "ERROR: -p argument needed";Errors+=1;fi
	if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
	if [ -z $matrix_pfm ]; then echo "ERROR: -m argument needed";Errors+=1;fi
	if [ -z $outdir ]; then echo "ERROR: -o argument needed";Errors+=1;fi
	if [ -z $symmetry ]; then echo "ERROR: -sym argument needed";Errors+=1;fi
	if [ $Errors -gt 0 ]; then echo "error somewhere"; return 1; fi

	# log 
	log=$outdir/$name/log.txt
	touch $outdir/$name/log.txt
	echo $Argline > $log

	# get motif length 
	echo "[INFO] - Get motif length..."
	tmp545=$(cat $matrix_pfm | wc -l)
	motifLength=$(expr $tmp545 - 2)
	echo "[REPORT] - Yout motif length is: $motifLength"

	# ----- Scanning PFM against all bound regions given 
	echo "[INFO] - Sacnning PFM on given peaks..."
	mkdir -p -m 774 $outdir/pfm_search

	re='^[+-]?[0-9]+([.][0-9]+)?$'
	if ! [[ $cutoff =~ $re ]] ; then
		if [ -f $cutoff ]; then
			local cutoff=$(awk '{print $3}' $cutoff2 )
			local cutoff2=$cutoff
		fi
	fi

	# Check if input file is a Comparison output file (generated by pairwize comparison for instance)
	if [[ $peaksCov == *"_RiL_RiP.tsv"* ]]; then 

		echo "[INFO] - Formatting peaks input coverage file ..."

		# get right column index (corresponding to ampDAP_RiL and DAP_RiL)
		read -ra header <<< $(awk 'NR==1{print tolower($0)}' $peaksCov) #extract column names into a liste
		# indexinitialisation
		ampdap_ril_index=-1
		dap_ril_index=-1
		# parse column name and get right index
		for i in "${!header[@]}"; do
			if [[ "${header[$i]}" == "ampdap_ril" ]]; then
				ampdap_ril_index=$i
			fi
			if [[ "${header[$i]}" == "dap_ril" ]]; then
				dap_ril_index=$i
			fi
		done

		echo "ampdap_ril index: $ampdap_ril_index"
		echo "dap_ril index: $dap_ril_index"

		# To use RiL coverage and filter out low coverage peaks
		# awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$8,$7}' $peaksCov > $outdir/tmp_peaks.tsv
		# awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$10,$9}' $peaksCov | awk '$4>3 && $4>3' > $outdir/tmp_peaks.tsv
		# awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$9,$7}' $peaksCov | awk '$4>3 && $4>3' > $outdir/tmp_peaks.tsv #AJ 16/06/2025 change column order due to modifications in comparison (pairwize instead of initial), so i select first DAP_RiL and then ampDAP_RiL
	 	
		## select first DAP_RiL and then ampDAP_RiL into a tmp peaks file used for methylation 
		awk -v ampdap_ril_idx="$ampdap_ril_index" -v dap_ril_idx="$dap_ril_index" -v OFS="\t" 'NR!=1 {print $1, $2, $3, $(dap_ril_idx+1), $(ampdap_ril_idx+1)}' $peaksCov | awk '$(dap_ril_idx+1)>3 || $(ampdap_ril_idx+1)>3' > $outdir/tmp_peaks.tsv

		local peaksCov=$outdir/tmp_peaks.tsv
	else 
		echo "[ERROR] - Wrong input file format (check columns names)"
	fi 

	echo "[REPORT] - Peaks formatted, scanning ..."
	head $peaksCov

	# TODO
	# 19/06/2024, RBM: add the line below to remove redudant peaks. in the bedtools line juste after $peaksCov is thus change for $outdir/tmp_peaks_redundancy_removed.tsv
	cat $peaksCov | cut -f 1,2,3 | sort -k1,1 -k2,2n -k3,3n | uniq > $outdir/tmp_peaks_redundancy_removed.tsv
	# get peaks FASTA
	bedtools getfasta -fi $genome -bed $outdir/tmp_peaks_redundancy_removed.tsv -fo $outdir/pfm_search/all_peaks.fasta
	
	# Compute scores
	echo "[INFO] - Computing PFM scores..."
	python3 $scores_prog -m $matrix_pfm -f $outdir/pfm_search/all_peaks.fasta -o $outdir/pfm_search
	pfmResults=$outdir/pfm_search/all_peaks.fasta.scores

	# ----- Detect significant TFBS 
	# searching distribution scores if th set to 0
	if [ "$cutoff" == "0" ]; then
		echo "[INFO] - User set thresholds not detected, searching for thresholds in Scores_Distribution folder"
		thresholds=() 
		for i in 50 70 90; do 
			file="scores_99_${i}_info.txt"
			echo ${thresholds_dir}/${name}/${file}
			if [ -f ${thresholds_dir}/${file} ]; then thresholds+=("$(cut -f 2 "${thresholds_dir}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${file} ]; then thresholds+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${file}")"); fi; 
			if [ -f ${thresholds_dir}/${name}/${file} ]; then thresholds+=("$(cut -f 2 "${thresholds_dir}/${name}/${file}")"); fi;
			if [ -f ${thresholds_dir}/Scores_Distribution/${name}/${file} ]; then thresholds+=("$(cut -f 2 "${thresholds_dir}/Scores_Distribution/${name}/${file}")"); fi; 
		done
		echo "[REPORT] - Found thresholds detected:"
		echo $thresholds
	
		local cutoff=$thresholds
		local cutoff2=$thresholds
	fi

	# ----- Compute methylation probabilities and statistics associated
	echo "[INFO] - Get methylation probabilities"
	python3 $full_methylation -p $peaksCov -g $genome -m $methMap -s $pfmResults -o $outdir -c $cutoff -l $motifLength
	
	echo "[INFO] - Compute methylation statistics and plots"
	$PATHRscript $plot_meth_full $outdir $motifLength $symmetry $cutoff2
	$PATHRscript $figs_meth_violin $outdir $matrix_pfm $symmetry
}

#-------------------------------------------------------------------------------
cons_scores(){
	local explicit_mode=false
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	local negfile="none" #to avoid leaving an unset variable
		case $1 in
			-f)
				local filein=$2
				echo "input peak file: ${2}";shift 2;;
			-o)
				local outdir=$2
				echo "output directory set to: ${2}";shift 2;;
			-m)
				local matrix=$2
				echo "TF matrix: ${2}"; shift 2;;
			-e)
				local extend=$2
				echo "Extension required on both sides: ${2}"; shift 2;;
			-n)
				local negfile=$2
				echo "Bed file containing negative gene positions: ${2}"; shift 2;;
			-nset)
				local negset=$2
				echo "Number of negative set required: ${2}"; shift 2;;
			-seed)
				local seed=$2
				echo "Specified seed: ${2}"; shift 2;;
			--explicit)
				local explicit_mode=true
				echo "explicit mode required: extension will be integrated in output files"; shift 1;;
			*)
				echo "Error in arguments"
				echo $1; usage comparison; exit;;
		esac
	done

	if [ -z $seed ]; then local seed=34562; fi
	if [[ $negfile == "none" ]]; then echo "no negative regions file given, shuffling will exclude input peak file sequences"; else echo -e "negative regions file: ${negfile}\nshuffling positive regions within negative regions given"; fi
	if [ -z $negset ]; then echo "number of negative sets non specified, using default (1)"; local negset=1; fi


	## usage
	# cons_scores -f /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/orthologs/conservation/tests/ChIP_DAP_DEG_peakscoord.bed -o /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/orthologs/conservation/LFY/ChIP_DAP_DEG -m /home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm -e 1000 -n /home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/neg_controls_DAP00001/all_nc_ATXG_positions.bed --explicit

	## debug files for LFY
	# outdir=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/orthologs/conservation/tests
	# filein=$outdir/ChIP_DAP_DEG_peakscoord.bed
	# matrix=/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm
	# extend=1000
	# negfile=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/results/neg_controls_DAP00001/all_nc_ATXG_positions.bed


	## Necessary paths and files:
	local A_thaliana_FASTA=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas
	# local scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/bin/scores.py
	local get_bestscore_prog=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/scripts/select_bestscore.py
	local genome_file=/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt
	local plot_scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/scripts/plot_cons_scores_v2.r
	local bwtool_path=/home/312.3-StrucDev/312.3.1-Commun/bin/sl6

	mkdir -p $outdir ## if output directory doesn't exist, create it



	## Extract filename from initial peaks file
	local basename=$(basename $filein)
	local file_extension="${basename##*.}"
	local filename="${basename%.*}"
	echo "filename: ${filename}"
	echo "file extension: ${file_extension}"



	## check if chrN coordinates in peak file are given correctly
	if [[ $(awk 'NR==1{print $1}' $filein) != "chr"[0-9] ]]; then
		if [[ $(awk 'NR==1{print $1}' $filein) == [0-9] ]]; then
			awk -v OFS="\t" '{print "chr"$1,$2,$3}' $filein | sort -k1,3 | uniq > $outdir/${filename}.bed
			local filein=$outdir/${filename}.bed
			echo "chromosome coord: \"N\""
		elif [[ $(awk 'NR==1{print $1}' $filein) == "Chr"[0-9] ]]; then
			sed 's/Chr/chr/g' $filein | sort -k1,3 | uniq > $outdir/${filename}.bed
			local filein=$outdir/${filename}.bed
			echo "chromosome coord: \"ChrN\""
		else
			echo "unknown chromosome formatting mode, exiting..."
			exit 0
		fi
	else
		echo "chromosome coordinates correctly formatted"
		if [[ ! -f $outdir/${filename}.bed ]]; then
			sort -k1,3 $filein | uniq > $outdir/${filename}.bed
		fi
		local filein=$outdir/${filename}.bed
	fi



	## If phylop and phastcons files are still bdg -> convert to bw
	if [[ ! -f /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP.bw ]]; then
		## Retrieve only chr1-5 from phastcons and phylop files
		grep -v "chrC" /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons_chrN.bedGraph | grep -v "chrM" > /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons_chrN.bedGraph
		
		grep -v "chrC" /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP_chrN.bedGraph | grep -v "chrM" > /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP_chrN.bedGraph
		
		## Convert bedgraph to bigwig for both files
		echo 'convert bedgraph files to bigwig'
		local bgtbw_path=/home/312.6-Flo_Re/312.6.1-Commun/Romain/ucscGenomeBrowser
		$bgtbw_path/bedGraphToBigWig /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons_chrN.bedGraph /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons.bw
		
		$bgtbw_path/bedGraphToBigWig /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP_chrN.bedGraph /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.txt /home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP.bw

	fi

	local phastcons=/home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhastCons.bw
	local phylop=/home/312.6-Flo_Re/312.6.1-Commun/data/Ath_PhyloP.bw



	## Get fasta sequences of input peak file, calculate scores and get best score position
	if [[ ! -f $outdir/${filename}_TFBS_coord.bed ]]; then
		bedtools getfasta -fi $A_thaliana_FASTA -bed $filein -fo $outdir/${filename}.fa
		
		echo "PWM scores calculation"
		python3 $scores_prog -m $matrix -f $outdir/${filename}.fa -o $outdir
		
		
		## Only keep best site per sequence (greatest score)
		# 	python3 $get_bestscore_prog -f $outdir/${filename}.fa.scores -o $outdir
		sort -k1,1 -k8,8nr $outdir/${filename}.fa.scores \
			| awk -v OFS="\t" '{print $2,$4,$8,$1}' \
			| uniq -f3 \
			| awk -v OFS="\t" '{print $4,$1,$3}' > $outdir/${filename}_bestscore.bed
		
		
		## add best score position info to original peak bed file
		sed 's/:/\t/g' $outdir/${filename}_bestscore.bed | sed 's/-/\t/' > $outdir/${filename}_bestscore_position.bed
		# old method after first sed: awk -v OFS="\t" '{gsub("-","\t",$2)}1' 
		
		
		## Cut coordinates around best TFBS 
		echo "retrieving TFBS coordinates..."
		if [[ $matrix == "/home/312.6-Flo_Re/312.6.1-Commun/data/LFY.pfm" ]]; then # calculating matrix length
			local matrixlength=19
		else
			local matrixlength=$(tail -n+3 $matrix | wc -l)
		fi

		## Retrieve TFBS coordinates
		awk -v OFS="\t" -v fullmatrixlength=$matrixlength '{print $1,$2+$4-1,$2+$4+fullmatrixlength-1}' $outdir/${filename}_bestscore_position.bed | sort -k1,1 -k2,3n -u > $outdir/${filename}_TFBS_coord.bed
		
		
	fi



	################ Extend and retrieve avg score per position ################ 
	####### METHOD FOR SHUFFLE AS IN PAPER #######
	echo "now using real shuffling method"

	if [[ ! -f $outdir/${filename}_shuffled_PhyloP.txt ]]; then
		echo "bwtool is going to be used now for positive controls"
	# 	echo "bwtool for TFBS coordinates"
		## Use these coordinates as input for bwtool
		$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phastcons $outdir/${filename}_shuffled_PhastCons.txt -firstbase

		$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phylop $outdir/${filename}_shuffled_PhyloP.txt -firstbase
	fi


	if [[ ! -f $outdir/${filename}_shuffled_PhyloP_medstdev.txt ]]; then
		echo "bwtool is going to be used to retrieve median stdev and number of values used, for all positive controls"
	# 	echo "bwtool for TFBS coordinates"
		## Use these coordinates as input for bwtool
		$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phastcons $outdir/${filename}_shuffled_PhastCons_medstdev.txt -expanded -firstbase

		$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord.bed $phylop $outdir/${filename}_shuffled_PhyloP_medstdev.txt -expanded -firstbase
	fi


	## if there are over 50k peaks in bed file, use only one negative set
	if [[ $(wc -l <$outdir/${filename}_TFBS_coord.bed) -gt 50000 ]]; then
		echo $(wc -l <$outdir/${filename}_TFBS_coord.bed)
		negset=1
		echo $negset
	fi

	# exit 0

	for ((i=1;i<=$negset;i++)); do
		if [[ ! -f $outdir/${filename}_shuffled_nset_all_PhyloP.txt ]]; then 
			## shuffle of these coordinates
			echo "shuffling..."
			((seed=seed*i))
			
			if [[ $negfile == "none" ]]; then
				bedtools shuffle -i $outdir/${filename}_TFBS_coord.bed -g $genome_file -chrom -seed ${seed} | sort -k1,3 | uniq > $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed
			else 
				bedtools shuffle -i $outdir/${filename}_TFBS_coord.bed -g $genome_file -incl $negfile -chrom -seed ${seed} | sort -k1,3 | uniq > $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed
			fi
			## bwtool of the coordinates
			echo "looking for phastcons and phylop coordinates for set ${i}"
			$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed $phastcons $outdir/${filename}_shuffled_nset${i}_PhastCons.txt -firstbase

			$bwtool_path/bwtool aggregate "$extend:$extend" $outdir/${filename}_TFBS_coord_shuffled_nset${i}.bed $phylop $outdir/${filename}_shuffled_nset${i}_PhyloP.txt -firstbase
		fi
	done


	if [[ ! -f $outdir/${filename}_shuffled_nset_all_PhyloP.txt ]]; then
		if [[ $negset != "1" ]]; then
			multijoin $outdir/${filename}_shuffled_nset_all_PhastCons.txt $outdir/${filename}_shuffled_nset[0-9]_PhastCons.txt
			ls $outdir | grep -P "^${filename}_shuffled_nset[0-9]+_PhastCons.txt" | awk -v dir=$outdir '{print dir"/"$1}' | xargs -d"\n" rm
			# rm $outdir/${filename}_shuffled_nset[0-9]_PhastCons.txt 
			
			multijoin $outdir/${filename}_shuffled_nset_all_PhyloP.txt $outdir/${filename}_shuffled_nset[0-9]_PhyloP.txt
			ls $outdir | grep -P "^${filename}_shuffled_nset[0-9]+_PhyloP.txt" | awk -v dir=$outdir '{print dir"/"$1}' | xargs -d"\n" rm
			# rm $outdir/${filename}_shuffled_nset[0-9]_PhyloP.txt 
			
			## remove negsets
			ls $outdir | grep -P "^${filename}_TFBS_coord_shuffled_nset[0-9]+.bed$" | awk -v dir=$outdir '{print dir"/"$1}' | xargs -d"\n" rm
		else
			mv $outdir/${filename}_shuffled_nset1_PhastCons.txt $outdir/${filename}_shuffled_nset_all_PhastCons.txt
			mv $outdir/${filename}_shuffled_nset1_PhyloP.txt $outdir/${filename}_shuffled_nset_all_PhyloP.txt
		fi
	fi

	echo "now plotting..."
	local R_cons_plots=/home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/LFYUFO_figs/bin/Rscript
	## plot
	if [[ ! -f $outdir/${filename}_shuffled_phastcons_phylop.pdf ]] || [[ ! -f $outdir/${filename}_shuffled_phastcons_phylop_${extend}.pdf ]]; then
		if [[ $explicit_mode == "true" ]]; then
			
			$R_cons_plots $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend --explicit
	# 		conda run -n LFYUFO_figs Rscript $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend --explicit
		else
			
			$R_cons_plots $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend
	# 		conda run -n LFYUFO_figs Rscript $plot_scores_prog -n ${filename}_shuffled -i $outdir -o $outdir -e $extend
		fi
	fi

}
