#################################################################################
#################################################################################
# ____ ____ _  _ ____ _ ____ _  _ ____ ____ ___ _ ____ _  _    ____ _ _    ____ #
# |    |  | |\ | |___ | | __ |  | |__/ |__|  |  | |  | |\ |    |___ | |    |___ #
# |___ |__| | \| |    | |__] |__| |  \ |  |  |  | |__| | \|    |    | |___ |___ #
#################################################################################
#################################################################################

# This file gather all installations and data needed for the pipeline:
# (You will need to complete it before using the TF_genomic_analysis/compil_functions.sh scripts)

 --> Scripts and tools developed in the team are avalaible in the /bin directory. 
        SCRIPT_DIR=/home/312.6-Flo_Re/312.6.1-Commun/TF_genomic_analysis
        PATH_TO_COMPIL="$SCRIPT_DIR"
        BIN_DIR="$SCRIPT_DIR/bin"

# --> Software used *via* CONDA environments can be installed in a new environment
        # with TF_genomic_analysis/conda.yml and TF_genomic_analysis/conda_python2.yml files (see tutorial).
        CONDA_ENV="/home/312.6-Flo_Re/312.6.1-Commun/Conda_forge/Flore/miniforge3/envs_dirs/TFgenomics"
        CONDA_MACS3_ENV="/home/312.6-Flo_Re/312.6.1-Commun/Conda_forge/Flore/miniforge3/envs_dirs/macs3"
        # Note that some tools (TFFM, related to MEME and DNAShape), use python 2.7, so you will
        # need to create a specific CONDA environment and using pyhton from this environment.
        Python_TFFM=/home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/py27_tffm/bin/python

        

#################################################################################
# PATH TO DATA
DATA_DIR="$SCRIPT_DIR/data"

#################################################################################
# GENERAL PATH SETTINGS
PATHRscript=path_to_your_R/bin/Rscript # change for your own path to Rscript (we used v4.4.0)
source $PATH_TO_COMPIL/compil_usages.sh #The TF_genomic_analysis/compil_usages.sh gather
                                        #all pipieline functions helper

#################################################################################
# PATH USED IN THE PIPELINE (BY FUNCTIONS FROM TF_genomic_analysis/compil_functions.sh)

# ---------- download_SRA
SRAtoolkit=/home/312.3-StrucDev/312.3.1-Commun/SRAtoolkit/sratoolkit.3.0.0-ubuntu64/bin/ # (v3.0.0)
export TMPDIR2=/nobackup/ # use a disk with some space (~20Go for safety) for SRAtoolkit

# ---------- mapping_Fastq_bowtie2
FASTQCmk=/home/312.3-StrucDev/312.3.1-Commun/FastQC/FastQC/fastqc # v0.11.7
NGmerge=/home/312.3-StrucDev/312.3.1-Commun/bin/sl7/NGmerge # v0.2_dev
PATH_TO_SAMBLASTER=/home/312.6-Flo_Re/312.6.1-Commun/Programs/samblaster # v0.1.26
PATH_TO_BOWTIE_INDEX="/home/312.3-StrucDev/312.3.1-Commun/bowtie/bowtie2-2.3.4.1-linux-x86_64/indexes/at"
PATH_TO_BOWTIE2=$CONDA_ENV/bin/ # v2.3.4.1
PATH_TO_SAMTOOLS=$CONDA_ENV/bin/ # v1.8 (using htslib v1.8)
PATH_TO_BEDTOOLS2=$CONDA_ENV/bin/ # v2.27.1

# ---------- peakcalling_MACS2
PATH_TO_MACS=$CONDA_ENV/bin/
mspc_mk=/home/312.6-Flo_Re/312.6.1-Commun/lib/MSPC/mspc/mspc # v5.4.0
bdg2bwig=$BIN_DIR/bedGraphToBigWig

# ---------- peakcalling_MACS3
PATH_TO_MACS3=$CONDA_MACS3_ENV/bin

# ---------- replicate_comp
Pairwise_comps=$BIN_DIR/Pairwise_comps.R

# ---------- pairwize_comparison
do_plot_quant=$BIN_DIR/Hist_cov_gen.r
merge_peaks=$BIN_DIR/merge_all_peaks.py
compute_coverage=$BIN_DIR/compute_coverage.py
comparisonPlotMaker=$BIN_DIR/ComparisonPLotMaker.R
StatEdgeR=$BIN_DIR/statsedgeR.R

# ---------- comparison
merge_peaks_Nsets=$BIN_DIR/merge_all_peaks_Nsets.py 

# ---------- analyze_decile
HeamapClassic=$BIN_DIR/Heatmap.R
HeatmapSymmetric=$BIN_DIR/HeatmapAR.R

# ---------- compute_rpkmrip_rpkmril
norma_cov_plots=$BIN_DIR/rpkmrip_rpkmril_perRep_per_Condition_covPlots.R
rpkmrip_rpkmril=$BIN_DIR/rpkmrip_rpkmril.R

# ---------- compute_motif
meme_prog=/home/prog/meme/meme_4.12.0/bin/meme-chip # v4.12.0
meme2meme=/home/prog/meme/meme_4.12.0/bin/meme2meme # v4.12.0
meme2pfm=$BIN_DIR/meme2pfm.sh
prepMEMEforPalTFFM=$BIN_DIR/prepMEMEforPalTFFM.py
pfmTOtffm=$BIN_DIR/get_tffm.py
transpose=$BIN_DIR/transpose_matrix.r

# ---------- compute_DNAShape
# All scripts came from the Git repository https://github.com/amathelier/DNAshapedTFBS
# Please clone this one in the desired directory and change the path below
PATH_TO_DNAshapedTFBS=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DNAshapedTFBS-master
ComputeDNAshaped=$PATH_TO_DNAshapedTFBS/DNAshapedTFBS.py
HeatmapDNAshape=$PATH_TO_DNAshapedTFBS/feature_importance_heatmap.py
araTha=/home/312.6-Flo_Re/312.6.1-Commun/data/DNAshape/A_thaliana # found in version in DATA_DIR
helt=$araTha/HeIT.bw;
mgw=$araTha/MGW.bw;
prot=$araTha/ProT.bw;
roll=$araTha/Roll.bw;
helt2=$araTha/HeIT2o.bw;
mgw2=$araTha/MGW2o.bw;
prot2=$araTha/ProT2o.bw;
roll2=$araTha/Roll2o.bw;

# ---------- compute_kmer
# All scripts came from the Git repository https://github.com/gifford-lab/GEM3
# Please clone this one in the desired directory and change the path below
GEM=/home/312.6-Flo_Re/312.6.1-Commun/Romain/gem/gem.jar

# ---------- compute_SeqConv
# All scripts came from the Git repository https://github.com/shenwei19/SeqConv
# Note that minor changes were added to the script to adapt it to our pipeline.
# We decided to keep the scripts in our /bin directory, not identical to the original ones.
SeqConv_train=$BIN_DIR/SeqConv_train.py
SeqConv_predict=$BIN_DIR/SeqConv_predict.py
SeqConv_motif=$BIN_DIR/SeqConv_motif.py
plot_logo=$BIN_DIR/plot_logo.R

# ---------- compute_NS
negative_set_script=$BIN_DIR/neg_set_gen.py

# ---------- compute_ROCS
# NB : ces scripts historiques vivent dans un repertoire separe du BIN_DIR actuel
ROCS_SCRIPTS_DIR=/home/312.6-Flo_Re/312.6.1-Commun/scripts/TFgenomicsAnalysis/bin
parse_KSM_scores=$ROCS_SCRIPTS_DIR/parse_KSM_scores.py
pocc_pfm=$ROCS_SCRIPTS_DIR/compute_POcc.py
tffmscores=$ROCS_SCRIPTS_DIR/get_best_score_tffm.py
scores_prog=$ROCS_SCRIPTS_DIR/scores.py 
plot_ROCS_prog=$BIN_DIR/plots_ROCS_multiple.py
plot_earlyROCS_prog=$BIN_DIR/plots_earlyROCS_multiple.py
dimer_builder=$ROCS_SCRIPTS_DIR/build_dimer_matrix.sh 

# ---------- compute_distribution
interactivedistrib=$BIN_DIR/InteractiveDistrib.r

# ---------- compute_space, compute_spacing2TFs and compute_spacing_TFBSTSS
spacing_mk=$BIN_DIR/get_interdistances.py
tffm_all_scores=$BIN_DIR/get_all_score_tffm.py
Zscore_spacing=$BIN_DIR/Zacing.R
get_interdistances_2TF=$BIN_DIR/get_interdistances_2TF.py
get_interdistances_TFBSTSS=$BIN_DIR/get_interdistances_TFBSTSS.py
Zacing_2TFs=$BIN_DIR/Zacing_2TF.R
builder_2matrix=$BIN_DIR/build_dimer_2matrix.sh
Zacing_2TSS=$BIN_DIR/Zacing_2TSS.R

# ---------- heatmap_reads #TODO to rm ? not used
heatmap_mk=$BIN_DIR/show_reads.py
bdg_to_bw=$BIN_DIR/bedGraphToBigWig
impactspacing_mk=$BIN_DIR/compute_spacingScorev2.r

# ---------- methylation
# Original publication PUMMEDID: 36782389
methMap=/home/312.6-Flo_Re/312.6.1-Commun/data/Methylation_maps/Zhu_lab_PNAS_2016/GSM1876327_Col.mC.bed # found in version in DATA_DIR
full_methylation=$BIN_DIR/full_methylation.py
plot_meth_full=$BIN_DIR/plot_meth_full.R
figs_meth_violin=$BIN_DIR/figs_meth_violin.R