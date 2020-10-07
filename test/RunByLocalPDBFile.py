from mainPackage.Functions import RigidDomainFinder
from GraphPackage.Graph_Config import List_Colors
from igraph import plot, Plot
import cairocffi
if __name__ == '__main__':
    rf = RigidDomainFinder()
    PredictedLabels = rf.segment_by_PDBFile('./data/adk.pdb','adk','A')
    print(PredictedLabels)
    #ProtG = rf.get_protein_graph()
    #ProtG.vs['color'] = [List_Colors[i] for i in PredictedLabels]
    #plot(ProtG)

    