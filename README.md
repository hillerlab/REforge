# REforge

REforge (Regulatory Element forward genomics) is a method to associate transciption factor binding site divergence in regulatory elements with phenotypic differences between species [1].

# Requirements:
* Python3 with Bio and rpy2 (requires R)
* Stubb_2.1 [2]
* tree_doctor [3]. A linux 64 bit binary is included. The source code is available at https://github.com/CshlSiepelLab/phast/ in src/util/

# Installation:

`git clone https://github.com/hillerlab/REforge`

Download Stubb (http://www.sinhalab.net/software) and newmat11 (http://www.robertnz.net/download.html)
```
cp /path/to/stubb_2.1.tar.gz /path/to/newmat11.tar.gz .
tar -xvf stubb_2.1.tar.gz
tar -xvf newmat11.tar.gz -C stubb_2.1/lib/newmat/

# REforge uses a slightly modified version of Stubb
cd stubb_2.1/
patch -p1 < ../REforge/stubb.patch
cd lib/newmat/
gmake -f nm_gnu.mak
cd ../../
make
export PATH=$PATH:`pwd`/bin
cd ../REforge/
export PATH=$PATH:`pwd`
```

# Example of screening 55 simulated CREs
```
# this directory contains a minimal example of simulated CREs
cd example

# Create a joblist file 'alljobs_simulation' containing the branch_scoring jobs
REforge.py data/tree_simulation.nwk data/motifs.wtmx data/species_lost_simulation.txt elements_simul.ls \
  --windowsize 200 --add_suffix _simulation -bg=background/

# Run job list as batch or in parallel 
bash alljobs_simulation > scores_simulation

# Run the association test
REforge_statistics.py data/tree_simulation.nwk data/motifs.wtmx data/species_lost_simulation.txt \
   elements_simul.ls --add_suffix _simulation
# This generates a file 'significant_elements_simulation' that alphabetically lists the elements, their P-value and the number of branches
```

# General workflow
## Input data
- species tree in newick format
- motif file in wtmx format (see example/data/motifs.wtmx as example)
- Phenotype-loss species list: one species per line (see example/data/species_lost_simulation.txt)
- CRE fastafiles: each file contains the sequence for every (ancestral or extant) species
- CRE list: path of fastafile of one CRE per line

## Step 1: Branch score computation
Generate the REforge branch_scoring commands for all CREs.
```
REforge.py <tree> <motiffile> <lost_species_list> <element_list>
```
This creates for every CRE a REforge_branch_scoring.py job. Each line in alljobs<suffix> consists a single job. Each job is completely independent of any other job, thus each job can be run in parallel to others.

Execute that alljobs file. Either sequentially by doing
```
chmod +x alljobs_simulation
./alljobs_simulation > scores_simulation
```
or run it in parallel by using a compute cluster.
Every job returns a line in the following format:
motif_file	CRE_file	(branch_start>branch_end:branch_score	)*	<
which should be concatenated into a file called "scores<suffix>".

## Step 2: Association test
```
REforge_statistics.py <tree> <motiffile> <lost_species_list> <element_list>
```
REforge_statistics.py classifies branches into trait-loss and trait-preserving and assesses the significance the association of this classification with the branch scores. 

## Common Parameters
#### REforge.py
```
--add_suffix <suffix>
Appends suffix to every generated file
--windowsize/-w <n>
Scoring window used in sequence scoring
--background/-bg <folder>
Background used for sequence scoring. Either a file or a folder with backgrounds for different GC contents
--scrCrrIter <n>
Number score correction iterations. 0 turns the score correction off
```
#### REforge_statistics.py
```
--add_suffix <suffix>
Appends suffix to every generated file
--filterspecies <comma separated list>
Exclude species from analyses
--elements <file>
Analyse only the elements specified <file>
```

## Special parameters
```
--verbose/-v
--debug/-d
```

#### REforge.py
```
--no_ancestral_filter
Turn of ancestral score filtering
--no_fixed_TP
Do not fix transition probabilites while computing branch scores
--filter_branch_threshold <x>
Filter branches if start and end node are below <x>
--filter_branches <file>
Like --filter_branch_threshold but with motif specific branch thresholds from <file>
--filter_GC_change <x>
Filter branches with a GC content change above <x>
--filter_length_change <x>
Filter branches with a relative length change above <x>
```
 
# References
[1] Langer BE, Roscito JG, Hiller M. REforge associates transcription factor binding site divergence in regulatory elements with phenotypic differences between species. submitted

[2] Sinha S, van Nimwegen E, Siggia ED. A Probabilistic Method to Detect Regulatory Modules. Bioinformatics, 19(S1), 2003

[3] Hubisz MJ, Pollard KS, Siepel A. PHAST and RPHAST: Phylogenetic Analysis with Space/Time Models. Briefings in Bioinfomatics 12(1):41-51, 2011.
