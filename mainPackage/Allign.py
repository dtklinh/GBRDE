import numpy as np
#from csb.bio.io.wwpdb import RemoteStructureProvider as PDB
#from csb.bio.io import StructureParser as SP

#from csb.bio.utils import rmsd
import os, sys
from Bio.PDB import *
#from scipy.spatial.distance import squareform
from Bio.PDB import Chain


def align_o(sequences):
    ''' Aligned function for given sequences
    :param sequences: Input sequences
    :return: aligned masks
    '''
    from mainPackage.clustalo import align
    align1 = align(sequences)
    L, N = align1.length, align1.size
    mask = np.zeros((N, L), 'i')

    for seq in align1:
        print(seq)
        i = int(seq.id[4:])
        m = np.array([x == '-' for x in seq.sequence])
        mask[i, :] = 1 - m
    return mask
	
def retrieve_struct(pdb_id, chain_id, load_local = False):
    '''Structure parser use mmCif supported more filetype
    :param pdb_id: PDB_ID protein
    :param chain_id: Alphabetic chain
    :return: structure and has_structured sequences
    '''
    from Bio.PDB import PDBList, FastMMCIFParser
    from Bio.PDB.Polypeptide import three_to_one
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, file_format='mmCif', pdir= './')
    parser = FastMMCIFParser()
    structure = parser.get_structure(pdb_id, pdb_id.lower()+'.cif')
    os.remove(pdb_id.lower()+'.cif')
    chain = structure[0][chain_id]
    coords = []
    structured_sequence = ''
    for residue in chain:
        if 'CA' in residue and residue['CA'].is_disordered() == 0:
            coords.append(residue['CA'].get_coord())
            structured_sequence += three_to_one(residue.resname)
        else:
            print((residue.is_disordered(), residue.id))
    return np.array(coords), str(structured_sequence)

def get_structure(input_text, load_local = False):
    '''
    Get structured and seuqences list from Webform input text data.
    :param input_text: string seperate by ,
    :param load_local: Optional loader for user upload PDB file.
    :return: aligned structures and sequence list
    '''
    if not load_local:
        pdb_id_ = input_text.split(',')
        structs =[]
        sequences = []
        for pdb_id in pdb_id_:
            pdb = pdb_id.split('_')[0]
            chain = pdb_id.split('_')[1]
            s1, seq = retrieve_struct(pdb_id=pdb, chain_id=chain)
            structs.append(s1)
            sequences.append(seq)
        mask = align_o(sequences=sequences)
        structs1 = np.zeros((mask.shape[0], mask.shape[1], 3))
        for i in range(mask.shape[0]):
            X1 = np.array(structs[i])
            for j, k in enumerate(np.where(mask[i] == 1)[0]):
                structs1[i, k] = X1[j]

    return np.array(structs1), mask

def get_coord_seq(C:Chain):
    from Bio.PDB.Polypeptide import three_to_one
    coords = []
    structured_sequence = ''
    for residue in C:
        if 'CA' in residue and residue['CA'].is_disordered() == 0:
            coords.append(residue['CA'].get_coord())
            structured_sequence += three_to_one(residue.resname)
        else:
            print((residue.is_disordered(), residue.id))
    return np.array(coords), str(structured_sequence)

def get_structure_local(Path2PDBFile, name, chainstr):

    structs = []
    sequences = []
    parser = PDBParser()
    structure = parser.get_structure(name, Path2PDBFile)
    for idx_model, model in enumerate(structure):
        chain = model[chainstr]
        coord, seq = get_coord_seq(chain)
        structs.append(coord)
        sequences.append(seq)
    arr = [len(s) for s in sequences]
    if min(arr) != max(arr):
        mask = align_o(sequences=sequences)
        structs1 = np.zeros((mask.shape[0], mask.shape[1], 3))
        for i in range(mask.shape[0]):
            X1 = np.array(structs[i])
            for j, k in enumerate(np.where(mask[i] == 1)[0]):
                structs1[i, k] = X1[j]
        return np.array(structs1), mask
    else:
        return np.array(structs), None




#Test
if __name__ == '__main__':
    structure_complex = np.empty((2,0,3))
    for c in ['A', 'B', 'C', 'D', 'E']:

        structure, mask = get_structure('{}_{},{}_{}'.format('4HFI',c,'4NPQ',c))

        desired_structured = np.compress((mask.mean(axis=0) == 1), structure, axis=1)
        structure_complex = np.concatenate((structure_complex,desired_structured), axis=1)

    np.savetxt('XYZ_1.txt',structure_complex[0,:,:])
    np.savetxt('XYZ_2.txt',structure_complex[1,:,:])


