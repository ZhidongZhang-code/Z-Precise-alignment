#!/bin/bash
#SBATCH -J cdhit
#SBATCH -N 1
#SBATCH -n 8

module load CentOS/7.9/Anaconda3/24.5.0
#source activate /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/miniconda/envs/bise-project/
#cd /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/project/Z-Precise-alignment/demo

#python3 /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/project/Z-Precise-alignment/src/main_alignment.py \
#--fasta_input 27_query_seed_urda.fas \
#--use_diamond  \
#--normailzationout urda.form \
#--config /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/project/Z-Precise-alignment/demo/config.yaml \
#--log_file 2025-6-2.log

python /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/project/Z-Precise-alignment/src/main_alignment.py  \
--fasta_input 27_query_seed_urda.fas  \
--normailzationout blast_normailzation.csv \
--extractcds blast_normailzation_excate_cds.fasta \
--extractprotein blast_normailzation_excate_protein.fasta \
--use_diamond \
--num_threads 1 --cpu 8 --log_file 2024_3_13.log \
--config ./config.yaml \
--target_sequence_file active_site_30aaQSYP30aa.fas \
--pfamkey FAD_binding_2 \
--clustalo_extract_out clu-extract-out.fasta \
--clustalo_out clu-out.fasta
#--keep_temp

#python /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/project/Z-Precise-alignment/src/main_excat.py \
#--input clu-extract-out.fasta \
#--position 33:Y \
#--extractprotein blast_normailzation_excate_protein.fasta \
#--extractcds blast_normailzation_excate_cds.fasta \
#--site_cds_output site-cds.fasta \
#--site_protein_output site-protein.fasta \
#--formatoutput urda-format \
#--config /cpfs01/projects-HDD/cfff-47998b01bebd_HDD/zzd_21250700018/project/Z-Precise-alignment/cfg/config.yaml \
#--log 2024_3_12.log 
#--keep_temp 
