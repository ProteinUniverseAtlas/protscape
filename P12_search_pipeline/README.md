# P12 search pipeline

This folder includes the computational workflow used for the study of the structural space of the P12 protein, as well as the metadata resulting from the analysis.

- `AF2_models_best_ranking`: includes the best model (model_0) generated with AlphaFold2 for the 55 redundant P12 proteins from (Quinones-Olvera et al)[https://www.nature.com/articles/s41467-024-47416-z]
- `AF2_models_selected_regions_summary.csv`: a table summarising the domains collected as part of the search pipeline, with the follwoing columns:
     - `file`: the internal file identifier of the domain
     - `start_residue`: the residue setting the start of the domain in the full-length structure
     - `end_residue`: the residue setting the end of the domain in the full-length structure
     - `acession_code`: the accession code of the full-length model in the corresponding source database
     - `source_db`: the source database

- `AF2_models_selected_regions_foldseek_seqID50_all_seqs.fasta`: the fasta file with the sequences of all domains in the `AF2_models_selected_regions_summary.csv` file
- `AF2_models_selected_regions_foldseek_seqID50_rep_seq.fasta`: the fasta file of the representatives of the dataset, filtered to a sequence identity of 50%
- `AF2_models_selected_regions_foldseek_seqID50_cluster.tsv`: the map between proteins and their representatives
- `P12_searching_pipeline.sh`: a shell script of the steps taken to generate the dataset and perform all-against-all comparisons

## Notes

1. To run the shell script, set up `protscape`, `astrochop` and `foldseek_screen`, as described in this repository.
2. The structure files of all the structure models used are not provided due to space, but can be generated from the `AF2_models_selected_regions_summary.csv` file. Just save them to a folder named `AF2_models_selected_regions` and run the `protscape` commands in the shell script