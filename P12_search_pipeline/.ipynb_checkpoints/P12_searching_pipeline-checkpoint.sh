eval "$(conda shell.bash hook)"

# activate conda environment to do the search
conda deactivate
conda activate /scicore/home/schwede/soares0000/code/miniconda3/envs/astrochop

# OLD 1
# # # search the input P12 using 2 iterations over AFDB50 and ESMfold_Mgnify30 but only doing the second round for the hits with a probability higher than 80% > this already should bring OBfolds
# # python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in P12_MXstr_core.pdb -n_iterations 2 -db afdb50 mgnify_esm30 -iter_prob 0.8

# # # now search with the top dali matched domain with 1 iteration over AFDB50 and ESMfold_Mgnify30
# # python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in 8gme_matched_domain.pdb -n_iterations 1 -db afdb50 mgnify_esm30 

# OLD2
# # search the input P12 using 2 iterations over PDB100, BFVD, AFDB50 and ESMfold_Mgnify30 but only doing the second round for the hits with a probability higher than 80% > this already should bring OBfolds
# python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in P12_MXstr_core.pdb -n_iterations 2 -db pdb100 afdb50 mgnify_esm30 bfvd -iter_prob 0.8 

# # now search with the top dali matched domain with 1 iteration over PDB100, BFVD, AFDB50 and ESMfold_Mgnify30
# python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in 8gme_matched_domain.pdb -n_iterations 1 -db pdb100 afdb50 mgnify_esm30 

# # NEW WITH BFVD
# # search the input P12 using 2 iterations over PDB100, AFDB50 and ESMfold_Mgnify30 but only doing the second round for the hits with a probability higher than 80% > this already should bring OBfolds
# python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in P12_MXstr_core.pdb -n_iterations 2 -db pdb100 afdb50 mgnify_esm30 -iter_prob 0.8

# # now enrish with P12 matches in BFVD using the same parameters
# python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in P12_MXstr_core.pdb -n_iterations 2 -db bfvd -iter_prob 0.8

# # now search with the top dali matched domain with 1 iteration over PDB100, BFVD, AFDB50 and ESMfold_Mgnify30
# python /scicore/home/schwede/soares0000/code/foldseek_screen.py -in 8gme_matched_domain.pdb -n_iterations 1 -db pdb100 afdb50 mgnify_esm30 

# # now activate the protscape mamba environment
# conda deactivate

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"

mamba activate protscape

# # filter the pdbs in the AF2_models_selected_regions folder (the matched domains) to a max seqID of 50%
foldseek easy-cluster AF2_models_selected_regions AF2_models_selected_regions_foldseek_seqID50 /tmp --min-seq-id 0.50

# filter the seed and reference pdbs to a max seqID of 100% 
foldseek easy-cluster P12_MXstr_core.pdb 8gme_matched_domain.pdb OBfold_4MZ9.pdb reference_structures/8s4t_OB_singlechain.pdb reference_structures/5odl_OB_singlechain.pdb target_sequences /tmp --min-seq-id 1

# filter the Quinone P12 AF2 model pdbs to a max seqID of 100% (modelled with AF2-at-scicore using default parameters; only the ranked 0 model was then chosen)
foldseek easy-cluster AF2_models_best_ranking Quinone_P12_AF2_models /tmp --min-seq-id 1

# merge resulting reps with the sequences of the seeds 
cat target_sequences_rep_seq.fasta AF2_models_selected_regions_foldseek_seqID50_rep_seq.fasta Quinone_P12_AF2_models_rep_seq.fasta > AF2_models_selected_regions_foldseek_seqID50_rep_seq_with_seeds.fasta

# filter them to a max seqID of 100%
mmseqs easy-cluster AF2_models_selected_regions_foldseek_seqID50_rep_seq_with_seeds.fasta AF2_models_selected_regions_foldseek_seqID50_rep_seq_with_seeds /tmp --min-seq-id 1

# run protscape in mmseqs mode for the resulting reps fasta file and define clusters automatically
python3 /scicore/home/schwede/soares0000/code/protscape.py -infile AF2_models_selected_regions_foldseek_seqID50_rep_seq_with_seeds_rep_seq.fasta -mode mmseqs -find_clusters True  -n_iterations 2 -evalue_clustering 1e-4 -evalue_sigma -0.5

# copy the seeds, the references and the Quinone P12 structures to the selected regions folder and run protscape in foldseek mode for the resulting reps and define clusters automatically
cp *pdb AF2_models_selected_regions
cd AF2_models_best_ranking
cp *pdb ../AF2_models_selected_regions
cd ../reference_structures
cp *OB_singlechain.pdb ../AF2_models_selected_regions
cd ../

python3 /scicore/home/schwede/soares0000/code/protscape.py -infile AF2_models_selected_regions -mode foldseek -find_clusters True  -n_iterations 2 -evalue_clustering 1e-4 -evalue_sigma -0.5 -infasta AF2_models_selected_regions_foldseek_seqID50_rep_seq_with_seeds_rep_seq.fasta

# # submit clustering with clans to the cluster
# # mkdir -p slurm_output
# # sbatch sbatch_clans.sh AF2_models_selected_regions_foldseek_seqID90_rep_seq_with_seeds_mmseqs.clans

# # run protscape in mmseqs mode only with those sequences that have an EBAmin score better than 3.5 (SCOPe same fold level)
# python3 /scicore/home/schwede/soares0000/code/protscape.py -infile AF2_models_selected_regions_foldseek_seqID90_rep_seq_with_seeds_eba_better_than_3.5.fasta -mode mmseqs -find_clusters True  -n_iterations 2 -evalue_clustering 1e-4 -evalue_sigma -0.5