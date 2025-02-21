#!/bin/bash

# Stop execution if we get an error along the way
#set -o errexit 

############################################
# Fit each accessibility model on sequences from all cell types
############################################

# Load required modules
ml system gcc

# Add executibles to path
export PATH=$HOME/work/lsgkm-svr/bin:$PATH

# Directory containing fit models
model_dir=/home/adufour/work/gskm/fit_models_1000bp

# Directory for prediction results
pred_dir=/home/adufour/work/gskm/model_predictions
if [ ! -d ${pred_dir} ]; then mkdir ${pred_dir}; fi

# Fasta files to predict
fasta_dir=/home/adufour/work/gskm/fastas_1000bp_randOnly
true_seq_files=(${fasta_dir}/*_true_seqs.fasta)

# Cell type headers for models
model_files=(${model_dir}/*.fold0.model.txt)
file_basenames=("${model_files[@]##*/}")
model_headers=("${file_basenames[@]/.fold0.model.txt/}")

# Define training/testing folds
all_chrs=($(seq 1 1 18))
all_chrs+=("X")
#all_chrs=("${all_chrs[@]/#/chr}")

# 10-fold cross validation
fold0=(1)
fold1=(2 18)
fold2=(3 13)
fold3=(6 X)
fold4=(5 16)
fold5=(4 15)
fold6=(7 14)
fold7=(11 17)
fold8=(9 12)
fold9=(8 10)
fold10=(Y)
fold11=(AEMK02000133.1 AEMK02000137.1 AEMK02000146.1)

#n_folds=($(seq 0 1 9))

n_folds=(0) # Only use one fold for testing cross-types
all_folds=("${n_folds[@]/#/fold}")

# Directory to hold testing data
test_dir=${fasta_dir}/testing_fastas
if [ ! -d ${test_dir} ]; then mkdir ${test_dir}; fi

########################################################################################
# Predict testing data for each fold with each cell type
########################################################################################

job_header=gkm_test

for model in "${model_headers[@]}"
do
    for fold in "${all_folds[@]}"
    do
        # Define fold chromosomes
        fchr=${fold}[@]
        test_chr=("${!fchr}")

        echo "Testing chromosomes for ${fold}..."

        # Model:
        model_file=${model_dir}/${model}.${fold}.model.txt

        # Output prefix
        outpre=${pred_dir}/${model}.${fold}

        for celltype in "${model_headers[@]}"
        do
            # Get true and null seq fasta files for testing
            true_seq_files=( $(printf "${fasta_dir}/${celltype}_%s_true_seqs.fasta " "${test_chr[@]}") )
            null_seq_files=( $(printf "${fasta_dir}/${celltype}_%s_null_seqs.fasta " "${test_chr[@]}") )

            true_result=${outpre}.pred.${celltype}.true.txt
            null_result=${outpre}.pred.${celltype}.null.txt

            # First true seqs:
            if [ ! -f ${true_result} ]
            then
                # Create files for combined testing data
                true_seqs=${test_dir}/${celltype}.${fold}.true_seqs.fasta
                if [ ! -f ${true_seqs} ]; then cat "${true_seq_files[@]}" > ${true_seqs}; fi

                # Predict accessibility
                true_log=${outpre}.pred.${celltype}.true.log
                echo "Predicting true ${celltype} accessibility using ${model}..."
                sbatch -t 12:00:00 --mem=10G --cpus-per-task=4 \
                --job-name=${job_header} --output=${true_log} --error=${true_log} \
                --wrap "gkmpredict -T 4 ${true_seqs} ${model_file} ${true_result}"
            else
                echo "File ${true_result} already exists! Skipping..."
            fi
            # Now null seqs:
            if [ ! -f ${null_result} ]
            then
                null_seqs=${test_dir}/${celltype}.${fold}.null_seqs.fasta
                if [ ! -f ${null_seqs} ]; then cat "${null_seq_files[@]}" > ${null_seqs}; fi

                # Predict accessibility
                null_log=${outpre}.pred.${celltype}.null.log
                echo "Predicting null ${celltype} accessibility using ${model}..."
                sbatch -t 12:00:00 --mem=10G --cpus-per-task=4 \
                --job-name=${job_header} --output=${null_log} --error=${null_log} \
                --wrap "gkmpredict -T 4 ${null_seqs} ${model_file} ${null_result}"
            else
                echo "File ${null_result} already exists! Skipping..."
            fi
        done
    done
done

########################################################################################
