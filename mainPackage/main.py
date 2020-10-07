from mainPackage.Functions import RigidDomainFinder

if __name__=='__main__':
    # Initilize object RigidDomainFinder
    RDF = RigidDomainFinder()

    # Run by pdbIDs
    PredLabels = RDF.segment_by_PDBIDs(['1ake_A', '4ake_A', '1ank_A'])
    print(PredLabels)

    # Run by local PDB file
    # Path2PDBFile = 'adk.pdb'
    # PredLabels = RDF.segment_by_PDBFile(Path2PDBFile,'ADK','A')
    # print(PredLabels)
