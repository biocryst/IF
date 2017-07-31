# IF
Coiled coils folding and intermediate filament models. Please see the 'Validation' folder if you want to reproduce the benchmarking results in our paper. For an example of how to generate the fragment libraries, look into the scripts in the 'Fragments' folder.

```
Usage: CCFold.py [-h] [-n 2|3] [-s] [-t one|termini|all] [-sp 100] [-bp 350] [-spp 50]
                 sequence [output_dir]


Positional arguments:
  sequence            Fasta file with the target sequence (hmodimer), 
                      or a pair of aligned sequences (heterodimer).
                     
  output_dir          Directory to output the model files.
  

Optional arguments:

  -h, --help          Show this help message and exit
 
  -n 2|3              Number of helices in a coiled coil. 2 (default) or 3 are
                      supported at the moment.

  -s                  Straighten the model. Prefer symmetric output at some
                      cost to smoothness. Number of fragments tried for each
                      segment is defined by the -spp parameter.
                      
  -t one|termini|all  Produce single model (default), vary termini or vary  
                      fragments in all windows. Number of fragments to vary is                      
                      defined by the -spp parameter.
                      
  -sp 100             Number of fragments to consider for initial 
                      segmentation.
                      
  -bp 350             Number of fragments to pick for mutual alignment.
  
  -spp 50             Number of fragments to use in shortest path outputs.
  ```
  
