
To obtain the validation results from the paper for parallel CCs:
1. Take the target sequence from Validation.fasta
2. Replace cc_data.npz with the dataset filtered from homologs from this folder.

For antiparallel targets a limited version of the program is provided here (fold_cc_anti.py), 
which you will need to edit by commenting/uncommenting target sequences. 
A few extra alanine residues will need to be removed from the resulting model.
Its dataset (cc_data_anti.npz) was constructed from non-redundant PDB structures at 50% sequence identity 
(as there are plenty antiparallel fragments in the CC+ database), so it was not further filtered to exclude homologs.
