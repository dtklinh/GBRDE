##/usr/bin/python3.7
#!./venv/bin/python3.7

import os, sys, getopt
#import numpy as np
sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(),'venv/bin'))
sys.path.append(os.path.join(os.getcwd(),'venv/lib/python3.7/site-packages'))


#os.system('. ./venv/bin/activate')
#sys.exit(1)

from mainPackage.Functions import RigidDomainFinder
#from GraphPackage.Graph_Util_Funcs import partition_gievenNumPar
import igraph as ig
from GraphPackage.Graph_Config import List_Colors

#print('Number of arguments:', len(sys.argv), 'arguments.')
#print('Arguments list:', str(sys.argv))

try:
    opts, args = getopt.getopt(sys.argv[1:],"hn:t:i:o:a:m:r:g:",
                               ["name=","type=","input=","output_folder=","AA_cutoff=","init_membership=","rigidity_threshold=",
                                "merging_threshold="])
except getopt.GetoptError:
      print('GBDE.py -n <ProteinName> -t <TypeOfFile> -i <inputfile> -a <AA_cutoff> -m <init_membership>, -r <rigid_thres> -g <merging_thres>')
      sys.exit(2)


name = ''
typeOfFile = ''
file = ''
OutFolder = os.getcwd()
AA_cutoff = 7.5
init_mem = None
rigidity_thres = 3.5
merging_thres = 1.0

for opt, arg in opts:
    print(opt,arg,'\n')
    if opt in ('-h', '--help'):
        print('GBDE.py -n <ProteinName> -t <TypeOfFile> -i <inputfile> -a <AA_cutoff> -m <init_membership>, -r <rigid_thres> -g <merging_thres> -o <output_directory>')
        sys.exit()
    elif opt in ('-n', '--name'):
        name = str(arg)
    elif opt in ('-t','--type'):
        typeOfFile = str(arg)
    elif opt in ('-i', '-input'):
        tmp = [x.strip() for x in str(arg).split(',')]
        if len(tmp)==1:
            file = tmp[0]
        else:
            file = tmp
    elif opt in ('-a', '--AA_cutoff'):
        AA_cutoff = float(arg)
    elif opt in ('-m', '--init_membership'):

        tmp = [x.strip() for x in str(arg).split(',')]
        #tmp = [x for x in pattern.split(str(arg)) if x]
        if len(tmp)>1:
            init_mem = [int(i) for i in tmp]
    elif opt in ('-r','--rigidity_threshold'):
        rigidity_thres = float(arg)
    elif opt in ('-g', '--merging_threshold'):
        merging_thres = float(arg)
    elif opt in ('-o', '--output_folder'):
        OutFolder = str(arg)
    else:
        print('GBDE.py -n <ProteinName> -t <TypeOfFile> -i <inputfile> -a <AA_cutoff> -m <init_membership>, -r <rigid_thres> -g <merging_thres> -o <output_directory>')
        sys.exit()

PredLabels = []
ProtG = None
#####--------------------------------------
rdf = RigidDomainFinder(AA_cutoff_neighborhood = AA_cutoff, init_membership = init_mem,
                        merging_threshold=merging_thres, rigidity_threshold = rigidity_thres)
print('Argument List: name={}\n, type of file={}\n, file={}\n, AA_cutoff={}\n, init_mem={}\n, rigid_thres={}\n, merg_thres={}'.format(
    name,typeOfFile, str(file), str(AA_cutoff),str(init_mem), str(rigidity_thres), str(merging_thres)))
if typeOfFile.lower()=='list':
    #None
    PredLabels = rdf.segment_by_PDBIDs(file) # file: list of PDBIDs
elif typeOfFile.lower()=='xyz':
    #None
    PredLabels = rdf.segment_by_xyzFormat(file)
elif typeOfFile.lower()=='pdb':
    #None
    PredLabels = rdf.segment_by_PDBFile(file,'X','A')
else:
    print('Type of file has to be one of those three: list, xyz, pdb')
    sys.exit(1)
ProtG = rdf.get_protein_graph()
##---------------------------------------
ProtG.vs['color']= [List_Colors[idx] for idx in PredLabels]
ProtG.vs['Cluster'] = [[v.index] for v in ProtG.vs]
out = ig.plot(ProtG, target= os.path.join(OutFolder,name+'.png'), inline=True)
#out.save(os.path.join(OutFolder,name+'.png'))
labels = set(PredLabels)
Lines = []
Lines.append('Predicted Labels')
Lines.append(','.join(str(x) for x in PredLabels))
for l in labels:
    #count = sum([ 1 for idx, val in enumerate(PredLabels) if val==l])
    rigid = ProtG.rmsd([idx for idx, val in enumerate(PredLabels) if val==l])
    Lines.append('RMSD for rigid domain id:{}'.format(str(l)))
    Lines.append(str(rigid))
from Utils.MyIO import WriteList2File
WriteList2File(os.path.join(OutFolder,name+'_Result.txt'), Lines)
#import numpy as np
#np.savetxt(os.path.join(OutFolder,name+'.csv'), PredLabels, delimiter=",", fmt='%s')
