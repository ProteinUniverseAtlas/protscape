eval "$(conda shell.bash hook)"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"

# activate conda environment to do the search
mamba deactivate
mamba activate astrochop

# search the input P12 using 2 iterations over PDB100, AFDB50 and ESMfold_Mgnify30 but only doing the second round for the hits with a probability higher than 80% > this already should bring OBfolds
python ../scripts/foldseek_screen.py -in P12_MXstr_core.pdb -n_iterations 2 -db pdb100 afdb50 mgnify_esm30 -iter_prob 0.8

# now enrish with P12 matches in BFVD using the same parameters
python ../scripts/foldseek_screen.py -in P12_MXstr_core.pdb -n_iterations 2 -db bfvd -iter_prob 0.8

# now search with the top dali matched domain with 1 iteration over PDB100, BFVD, AFDB50 and ESMfold_Mgnify30
python ../scripts/foldseek_screen.py -in 8gme_matched_domain.pdb -n_iterations 1 -db pdb100 afdb50 mgnify_esm30 

# now activate the protscape mamba environment
mamba deactivate
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

# copy the seeds, the references and the Quinone P12 structures to the selected regions folder and run protscape in foldseek mode for the resulting reps and define clusters automatically
cp *pdb AF2_models_selected_regions
cd AF2_models_best_ranking
cp *pdb ../AF2_models_selected_regions
cd ../reference_structures
cp *OB_singlechain.pdb ../AF2_models_selected_regions
cd ../

python3 ../scripts/protscape.py -infile AF2_models_selected_regions -mode foldseek -find_clusters True  -n_iterations 2 -evalue_clustering 1e-4 -evalue_sigma -0.5 -infasta AF2_models_selected_regions_foldseek_seqID50_rep_seq_with_seeds_rep_seq.fasta
