from mainPackage.Functions import RigidDomainFinder
import igraph as ig
from GraphPackage.Graph_Config import List_Colors
if __name__ == '__main__':
    rf = RigidDomainFinder(rigidity_threshold = 3.5)
    PredictedLabels = rf.segment_by_PDBIDs(['4ake_A', '1ake_A'])
    #PredictedLabels = rf.segment_by_PDBIDs(['4MBS_A', '6AKY_A'])
    #PredictedLabels = rf.segment_by_PDBIDs(['6aky_A', '4mbs_A', '6akx_A', '5uiw_A'])
    print(PredictedLabels)
    ProtG = rf.get_protein_graph()
    ProtG.vs['color'] = [List_Colors[i] for i in PredictedLabels]
    ig.plot(ProtG)

