# IF
Coiled coils folding and intermediate filament models
```
Usage: CCFold.py [-h] [-s] [-t one|termini|all] [-sp 100] [-bp 350] [-spp 50]
                 sequence [output_dir]


Positional arguments:
  sequence            Fasta file with the target sequence (hmodimer), 
                      or a pair of aligned sequences (heterodimer).
                     
  output_dir          Directory to output the model files.
  

Optional arguments:

  -h, --help          show this help message and exit
  
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
  
