# ProtScape

ProtScape is a set of scripts used to create and visualise protein similarity networks.

Currently, 3 scripts are made available in the `scripts` folder:

- `astrochop.py`: it takes as input a PDB structure and chops it into domains using [Chainsaw](https://github.com/JudeWells/chainsaw). It is implemented so that it can also take a consensus approach if other tools are available.
- `foldseek_screen.py`: it does iterative structural searches over any Foldseek structure database. It uses the FoldSeek API for the searches and `astrochop.py` to find the protein domains matched.
- `protscape.py`: it takes as input a fasta file or a set of protein structures and carries out all-against-all mmseqs or foldseek alignments to generate a protein similarity network in [CLANS](https://academic.oup.com/bioinformatics/article/20/18/3702/202504) format.

## Installation

### Astrochop and Foldseek_screen

1. As `foldseek_screen` relies on `astrochop` and `astrochop` relies on `chainsaw`, first set up `chainsaw` as described in [https://github.com/JudeWells/chainsaw].
2. Change the path to `chainsaw` in the `astrochop.py` script
3. Setup the `astrochop` environment:

```
mamba -n astrochop -c bioconda -c salilab -c conda-forge -c anaconda --file astrochop_requirements.txt
```

### Protscape

1. Setup the `protscape` environment:
```
mamba -n protscape -c bioconda -c salilab -c conda-forge -c anaconda --file protscape_requirements.txt
```

## Usage

### Astrochop
```
mamba activate astrochop

python3 scripts/astrochop.py inpdb_structure.pdb
```
This activates the `astrochop` environment and then runs all provided choppers (so far, parsers are compatible with [Chainsaw](https://github.com/JudeWells/chainsaw) and [Merizo](https://github.com/psipred/Merizo)) and outputs a consensus list of domain regions.

### Foldseek_screen
```
mamba activate astrochop

python3 scripts/foldseek_screen.py -in inpdb_structure.pdb -n_iterations 2 -db pdb100 afdb50 mgnify_esm30 bfvd -iter_prob 0.8
```
This activates the `astrochop` environment and then searches for structure homologs in the pdb100, afdb50, mgnify_esm30 and bfvd FoldSeek databases using the FoldSeek API in 3Di mode. All matched proteins are then chopped into domains with `astrochop` and only the domains matched at a probability >80% are used as inputs for a second round of searches. Full-length protein matches are saved in a `AF2_models` folder and the coordinates of the selected domains in the `AF2_models_selected_regions` folder.

### ProtScape

Activate the environment
```
mamba activate protscape

# To generate a protein sequence similarity network from a fasta file and find clusters automatically
python3 scripts/protscape.py -infile infasta.fa -find_clusters True

# To generate a protein sequence similarity network from a set of structures and find clusters automatically
python3 scripts/protscape.py -infile directory_with_structures -mode mmseqs -find_clusters True

# To generate a protein structure similarity network
python3 scripts/protscape.py -infile directory_with_structures -mode foldseek -find_clusters True

```

## Notes

The `P12_search_pipeline` folder contains a shell script corresponding to the computational workflow used to collect and analyse the structure space of the P12 protein, as well as the metadata to visualise the network in Cosmograph and the structure models generated for the P12 homologs.