from Bio.PDB.PDBParser import PDBParser
from GraphPackage.Graph_Util_Funcs import ConstructRealClusterGraph
from GraphPackage.Graph_Config import G_ConstructType, CutOffContact, List_Colors
from GraphPackage.AssistantObjects import DynDomEntry
import numpy as np
from scipy import spatial
from Bio.PDB import *
import sys
# def calDistanceMatrix(PDBID, ChainID, PDBFile):
#     parser = PDBParser(PERMISSIVE=1)
#     structure_id = PDBID
#     filename = PDBFile
#     structure = parser.get_structure(structure_id, filename)
#     for Model in structure:
#         for

def PreProcess_Local(Path2PDBFile:str, name:str, chainID:str):
    from mainPackage.Allign import get_structure_local
    # File contains all protein conformational states
    #Path2PDBFile = '../MyDataSet/adk/adk.pdb'
    structs, mask = get_structure_local(Path2PDBFile,name,chainID)
    if mask is not None:
        return np.compress((mask.mean(axis=0) == 1), structs, axis=1)
    else:
        return structs

def PreProcess_PDBIDs(LstPDB_Chains:list):
    joined_str = ','.join(LstPDB_Chains)
    from mainPackage.Allign import get_structure
    strucs, mask = get_structure(joined_str)
    desired_structured = np.compress((mask.mean(axis=0) == 1), strucs, axis=1)
    return desired_structured

def calc_DisMxs(XYZ):
    M = np.zeros((XYZ.shape[0],XYZ.shape[1], XYZ.shape[1]))
    for idx, xyz in enumerate(XYZ):
        M[idx, :, :] = spatial.distance_matrix(xyz, xyz)
    return M


def run_Alg(XYZ, Serial = 'ADK', cutoff_neighborhood = 7.5, init_membership = None, merging_threshold=1.0,
            rigidity_threshold = 3.5):
    # DisMatrices: m x L x L --> m distance matrices of L x L
    DisMatrices = calc_DisMxs(XYZ)


    Mem = [0]*DisMatrices.shape[1]
    Entry = DynDomEntry(None, Mem, DisMatrices, XYZ)
    #Entry.write_xyz_format('/home/linh/PycharmProjects/Git_Respo/Protein-Rigid-Domains-Estimation/CCR5.txt')
    #print('Success!!!')
    #sys.exit(1)
    #from GraphPackage.Graph_Config import NumberOfClusterInClusteringGraph
    from GraphPackage.Graph_Util_Funcs import ConstructGraph
    ProtG = ConstructGraph(Entry.DistanceMatrices, G_ConstructType, cutoff_neighborhood)
    ProtG['DynDomEntry'] = Entry
    G = ConstructRealClusterGraph(Entry.DistanceMatrices, Entry.Membership,init_membership = init_membership,
                                  Construct_Type=G_ConstructType, cut_off_threshold=cutoff_neighborhood)
    SquareMatFeature = G.calc_squareMatFeature(Entry.DistanceMatrices)
    G.vs['color'] = [List_Colors[v['TrueLabel']] for v in G.vs]
    G['DynDomEntry'] = Entry
    G['SquareMatFeature'] = SquareMatFeature
    G.vs['OriginalIndex'] = [v.index for v in G.vs]
    G.es['OriginalIndex'] = [e.index for e in G.es]
    G_Org_Indexs = [i for v in G.vs for i in v['Cluster']]
    G['serial'] = Serial
    Membership = Mem
    delete_indexs = [i for i in range(len(Membership)) if i not in G_Org_Indexs]
    # do iteration
    Arr = G.do_work_iteration_2(rmsd_thres=rigidity_threshold)
    Arr = G.do_merge(thres=merging_threshold, Arr_G=Arr)

    PredLabels = [-1] * len(Membership)
    for idx, c in enumerate(Arr):
        for v in c.vs:
            for i in v['Cluster']:
                PredLabels[i] = idx
    PredLabels = [i for j, i in enumerate(PredLabels) if j not in delete_indexs]
    return PredLabels, ProtG

class RigidDomainFinder(object):
    def __init__(self, AA_cutoff_neighborhood = 7.5, init_membership = None, merging_threshold=1.0, rigidity_threshold = 3.5):
        self.AA_cutoff_neighborhood = AA_cutoff_neighborhood
        self.init_membership = init_membership
        self.merging_threshold = merging_threshold
        self.rigidity_threshold = rigidity_threshold
        self.ProteinGraph = None
    def segment_by_PDBIDs(self, Lst_PDBs:list):
        struct = PreProcess_PDBIDs(Lst_PDBs)
        PredLabels, ProtG = run_Alg(struct,'Protein_name',self.AA_cutoff_neighborhood,self.init_membership,
                             self.merging_threshold, self.rigidity_threshold)
        self.ProteinGraph = ProtG
        return PredLabels
    def segment_by_PDBFile(self,Path2PDBFile:str, PDBID:str,ChainID:str):
        struct = PreProcess_Local(Path2PDBFile,PDBID,ChainID)
        PredLabels, ProtG = run_Alg(struct,PDBID, self.AA_cutoff_neighborhood,self.init_membership,
                             self.merging_threshold, self.rigidity_threshold)
        self.ProteinGraph = ProtG
        return PredLabels
    def segment_by_xyzFormat(self, Path_to_xyz_format:str):
        from MDAnalysis.coordinates.XYZ import XYZReader
        rd = XYZReader(Path_to_xyz_format)
        Structs = np.zeros((rd.n_frames,rd.n_atoms,3))
        for idx, ts in enumerate(rd.trajectory):
            tmp = ts.positions
            Structs[idx,:,:] = tmp
        PredLabels, ProtG = run_Alg(Structs,'Name',self.AA_cutoff_neighborhood,self.init_membership,
                             self.merging_threshold, self.rigidity_threshold)
        self.ProteinGraph = ProtG
        return PredLabels
    def get_protein_graph(self):
        return self.ProteinGraph


if __name__=="__main__":
    RDF = RigidDomainFinder()
    #PredLabels = RDF.segment_by_xyzFormat('../test/data/lysozyme.xyz')
    PredLabels = RDF.segment_by_PDBIDs(['3lww_A', '3lww_C'])
    print(PredLabels)
    sys.exit(1)
    Path2DisMx1 = '../ADK/dist_1.txt'
    Path2DisMx2 = '../ADK/dist_2.txt'
    Path2XYZ1   = '../ADK/xyz_1.txt'
    Path2XYZ2   = '../ADK/xyz_2.txt'
    #M1 = np.loadtxt(Path2DisMx1)
    #M2 = np.loadtxt(Path2DisMx2)
    #M = np.zeros((2,M1.shape[0], M1.shape[1]))
    #M[0,:,:] = M1
    #M[1,:,:] = M2
    XYZ_1 = np.loadtxt(Path2XYZ1)
    XYZ_2 = np.loadtxt(Path2XYZ2)
    XYZ = np.zeros((2,XYZ_1.shape[0], XYZ_1.shape[1]))
    XYZ[0,:,:] = XYZ_1
    XYZ[1,:,:] = XYZ_2

    PrdLabls = run_Alg(XYZ)
    print(PrdLabls)


