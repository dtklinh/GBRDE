# GBRDE
Graph-based Rigid Domains Estimator
# Compatibility
    Ubuntu 16.04 LTS or later
    Python 3.7
# Python packages dependencies:
    biopython (1.78 or later) from https://biopython.org/
    csb (1.2.5 or later) from https://csb.codeplex.com/
    louvain (0.6.1) from https://pypi.python.org/pypi/louvain/
    matplotlib (only to draw praph) (3.0.3 or later) from matplotlib.org
    numpy (1.16.4 or later) from https://numpy.org/
    python-igraph (0.7.1) from https://igraph.org/python/
    scikit-learn (0.21.3 or later) from https://scikit-learn.org/stable/
    scipy (1.3.1 or later) from https://www.scipy.org/
    MDAnalysis (0.20.1 or later) from https://www.mdanalysis.org
    cairocffi (1.1.0 or later)
# other package dependencies:
    clustal omega
    (sudo apt-get install -y clustalo)
    python3.7-dev
    (sudo apt install python3.7-dev)
# Usage
    Clone project from https://github.com/dtklinh/GBRDE.git
    (git clone https://github.com/dtklinh/GBRDE.git)
    
    Go to GBRDE folder and run one of those command lines:

    1. Run by a list of PDBID, for example
    ./venv/bin/python3.7 GBRDE.py -n 'adk' -t 'list' -i '1ake_A, 4ake_A' -a 7.5 -r 3.5 -g 1.0

    2. run by a xyz file, for example
    ./venv/bin/python3.7 GBRDE.py -n 'lys' -t 'xyz' -i '.test/data/lysozyme.xyz' -a 7.5 -r 3.5 -g 1.0

    3. Run by a pdb file, for example
    ./venv/bin/python3.7 GBRDE.py -n 'adk' -t 'pdb' -i './test/data/adk.pdb' -a 7.5 -r 3.5 -g 1.0
    
    GBDE.py parameters:

    -n <ProteinName> # any name
    -t <TypeOfFile>  # one in three {'list', 'pdb', 'xyz'}
    -i <inputfile>   # if TypeOfFile is 'pdb' or 'xyz', then it is a path to a file
                        if TypeOfFile is 'list', it is a list, e.g. '1ake_A, 4ake_A'
    -a <AA_cutoff>  # default: 7.5
    -m <init_membership>, # default: None
    -r <rigid_thres>   # default: 3.5
    -g <merging_thres> # default: 1.0
    -o <output_directory>' # default: current working directory