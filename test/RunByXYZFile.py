from mainPackage.Functions import RigidDomainFinder
from GraphPackage.Graph_Config import List_Colors
import igraph as ig
if __name__ == '__main__':
    rf = RigidDomainFinder(AA_cutoff_neighborhood = 7.5, init_membership = None, merging_threshold=1.0, rigidity_threshold = 3.5)
    # VMD likewise xyz format
    PredictedLabels = rf.segment_by_xyzFormat('./data/CCR5.xyz')
    print(PredictedLabels)
    ProtG = rf.get_protein_graph()
    ProtG.vs['Cluster'] = [[v.index] for v in ProtG.vs]
    ProtG.vs['color'] = [List_Colors[idx] for idx in PredictedLabels]
    out = ig.plot(ProtG)
    labels = set(PredictedLabels)
    for l in labels:
        count = sum([ 1 for idx, val in enumerate(PredictedLabels) if val==l])
        rigid = ProtG.rmsd([idx for idx, val in enumerate(PredictedLabels) if val==l])
        print('{} -> Rigid: {}'.format(str(count), str(rigid)))
    #out.save(os.path.join(OutFolder, name + '.png'))