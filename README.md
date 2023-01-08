# A computational approach for the identification of distant homologs of bacterial riboswitches based on inverse RNA folding
Our computational method for searching distant homologs of bacterial riboswitches can be described as a pipeline that is separated into the following steps: (I) Retrieve the representative consensus structure of the target riboswitch class, (II) Mutations and design of RNA with similar structure, (III) Seed generation and building of initial covariance model, (IV) Covariance model expansion, and (V) Database search and Phylogenomic analysis. The approach was influenced by synthetic biology, in which sequences are designed to perform specific tasks. The first two steps are responsible for designing synthetic sequences which are structurally similar to target riboswitch class, while the last three steps, building covariance model and database search and phylogenomic analysis, are important for the identification of the orthologs of the predicted structural candidates in all the relevant genomes and their evolutionary connection with target riboswitch class. 

![Figure 1-pipeline](https://user-images.githubusercontent.com/26137763/210639624-bea590b6-2b6a-4388-8a0c-69616f45fc6f.png)

Altough each and every screening steps depicted in this pipeline are important, but, the most crucial steps are <b>(1) Mutations and design of RNA with similar structure of target riboswitches using RNAfbinv 2.0,</b> and <b>(2) Building covariance model and database search </b>. To use this pipeline we need to install the following softwares (i) RNAfbinv 2.0 (which also require Vienna RNA package), (ii) Infernal 1.1, (iii) BLAST+ software. Instead of installing BLAST+ user can also perform web-BLAST available in https://blast.ncbi.nlm.nih.gov/Blast.cgi using NCBI-RefSeq database. 

# 1. Mutations and design of RNA with similar structure of target riboswitches using RNAfbinv 2.0

RNAfbinv is a fragment based RNA design tool. It uses a simulated annealing process to optimize a 2D RNA structure.<br/>
RNAfbinv 2.0 can be easily installed as it is available on pypi (python 3 compatible). To install it simply run ```pip install rnafbinv```.

[Vienna RNA package](https://www.tbi.univie.ac.at/RNA/ "Vienna RNA home") is required for RNAfbinv to work. This must be installed separately.<br/>
Current version was tested with Vienna 2.4 and above. RNAfbinv will identify Vienna package if it's bin directory is in PATH.<br/>
If you wish to link a specific installation of Vienna set the VIENNA_PATH environment variable to the correct bin directory.

You can set Vienna location in python
```python
import os
os.environ['VIENNA_PATH'] = "VIENNA_BIN_DIR_PATH"
```

or directly via the vienna script
```python
from rnafbinv import vienna
vienna.set_vienna_path("VIENNA_BIN_DIR_PATH")
```

## Usage

The design process can be ran using the following code:
```python
from rnafbinv import RNAfbinvCL
RNAfbinvCL.main(command_line_arguments)
```

To generate a tree for a specific sequence / structure:<br/>
Structure is a dot bracket notation structure and sequence is an IUPAC string with the same length
```python
from rnafbinv import shapiro_tree_aligner
shapiro_tree_aligner.get_tree(sructure, sequence)
```

To compare two trees and score them:
alignment_rules has a default value and is optional
```python
from rnafbinv import shapiro_tree_aligner
shapiro_tree_aligner.align_trees(source_tree, tree_target, alignment_rules)
```

## GUI / Command line

You can download the RNAfbinv wrapper from main [RNAfbinv2.0 git repository](https://github.com/matandro/RNAsfbinv/). Here, I have attached the clone version of the RNAfbinv 2.0 for better explanation of this method to this pipeline context.<br/>
The main folder includes python code to run the GUI / command line and a configuration file:
* RNAfbinv.py - A GUI wrapper for RNAfbinv2.0
* RNAfbinvCL.py - A command line wrapper for RNAfbinv2.0
* **Required** varna_generator.py - Used to generate images based on [VARNA](http://varna.lri.fr/ "VARNA rna homepage")
* **Required** config.ini - Configuration file with paths to required software (information below).
* **Required** img folder with NoImage.png - used in GUI as a placeholder

If you remove the VARNA jar or do not have java installed, images will not be generated but the design process will proceed normally.<br/><br/>

To specify [vienna package](https://www.tbi.univie.ac.at/RNA/ "The ViennaRNA Package homepage") binary folder please update the 'VIENNA' parameter in config.ini (or set VIENNA_PATH environment variable)<br/>
To specify Java binary folder please update the 'JAVA' parameter in config.ini (or set JAVA_PATH environment variable)<br/>
To specify [VARNA](http://varna.lri.fr/ "VARNA rna homepage")'s jar file please update the 'VARNA' parameter in config.ini (or set VARNA_PATH environment variable)<br/>
Note that if the java or vienna package binaries are in your environment variables you may leave it empty.

Example to a valid config.ini file which has java installed and within the system's path:
```
[PATH]
VIENNA=~/ViennaRNA/bin/
#JAVA=
VARNA=~/VARNA/VARNAv3-93.jar
```

### Command line arguments:

```
usage: RNAfbinvCL.py [-h] [-l LOG_OUTPUT] [--verbose | --debug]
                     [-p {MFE,centroid}] [-i ITERATIONS] [--seed SEED]
                     [-t LOOK_AHEAD] [--reduced_bi REDUCED_BI] [-e]
                     [--seq_motif] [-m MOTIF_LIST] [-s STARTING_SEQUENCE | -r]
                     [--length LENGTH] [-f INPUT_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -l LOG_OUTPUT, --log_output LOG_OUTPUT
                        Path to output log file. (default: None)
  --verbose             Increase output verbosity. (default: False)
  --debug               Debug level logging. (default: False)
  -p {MFE,centroid}, --structure_type {MFE,centroid}
                        uses RNAfold centroid or MFE folding. (default: MFE)
  -i ITERATIONS, --iterations ITERATIONS
                        Sets the number of simulated annealing iterations.
                        (default: 100)
  --seed SEED           Random seed used in the random number generator.
                        (default: None)
  -t LOOK_AHEAD, --look_ahead LOOK_AHEAD
                        Number of look head mutation attempts for each
                        iteration. (default: 4)
  --reduced_bi REDUCED_BI
                        Remove extra penalty for removal or addition of bulges
                        and interior loops under the given size. Alignment
                        penalties still occur. (default: 0)
  -e, --circular        Designs a circular RNA. (default: False)
  --seq_motif           Enables increased penalty for insertion or deletions
                        within marked regions (lower case characters in
                        sequence constraint). The feature was added to control
                        multi base sequence constraints (sequence motifs).
                        Only valid within a specific structural motif.
                        (default: False)
  -m MOTIF_LIST, --motif_list MOTIF_LIST
                        A comma separated list of motifs that are targeted for
                        preservation with size.Single motif format: <motif
                        No>[M|H|E|I|S|B]<motif No of bases>. Use
                        rnafbinv.ListMotifs.list_motifs(structure) to retrieve
                        a list of legal motifs for a given structure.
                        (default: [])
  -s STARTING_SEQUENCE, --starting_sequence STARTING_SEQUENCE
                        The initial sequence for the simulated annealing
                        process in IUPAC nucleotide codes. (default: None)
  -r, --random_start    Start simulated annealing with a random sequence.
                        (default: False)
  --length LENGTH       Maximum variation in result length compared to target
                        structure. (default: 0)
  -f INPUT_FILE         Path of ini file that includes mandatory information.
                        Some options can also be set via file. command line
                        options take precedence. (default: None)
```

### Input file format (the '-f' parameter):

```
# mandatory
TARGET_STRUCTURE=<target structure>
TARGET_SEQUENCE=<target sequence>
# optional
TARGET_ENERGY=<target energy>
TARGET_MR=<target mutational robustness>
SEED=<random seed>
STARTING_SEQUENCE=<starting sequence>
ITERATION=<number of simulated annealing iterations>
```

# 2. Building covariance model and database search
The source code of the latest version of Infernal is available on http://eddylab.org/infernal/.

User can also install using the following package manager:

If you are using Debian, you can install with:
```
sudo apt-get install infernal infernal-doc
```
If you have conda, you can install with:
```
conda install -c bioconda infernal=1.1.2
```
With homebrew, you can install with
```
brew tap brewsci/bio
brew install infernal
```

# Usage
<img width="434" alt="image" src="https://user-images.githubusercontent.com/26137763/211174724-ad0990ba-801d-4d0e-bcbe-1cd22b4653d7.png">
Commands for building covariance model (example with building covariance model for bacterial purine riboswitches)
```
cmbuild -F purine.cm SEED_sequence_purine.fasta
cmcalibrate --mpi purine.cm
```
Commands for database searching (example with searching sequence database for bacterial purine riboswitches)
```
cmsearch --nohmmonly -T 30.00 -Z 742849.287494 purine.cm SEQDB.fasta > database_search_result.txt
```




