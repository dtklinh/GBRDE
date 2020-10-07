'''Construct Graph from multiple Distance Matrices'''
import GraphPackage.Graph_Config as cf
import numpy as np
import igraph as ig
import math, louvain, sys
from igraph import Graph
from GraphPackage.GraphLibrary import D_Graph

def ConstructGraph(DisMxs, Construct_Type=None, cut_off_threshold = None, Membership=None):
    # print shape(ListDisMx)
    # print ListDisMx
    #  vertex_attrs={'Cluster':Idx_Of_Cluster, "ClusterTupleLabel":ClusterTupleLabel, 'TrueLabel':Ver_Label},
    #                    edge_attrs={'Connectivity':Ed_Conn, 'TrueLabel':EdgeLabel}

    num_Vertex = DisMxs.shape[1]
    # construct edges
    alpha = 1.0
    ed = list()
    w = list()
    Vertex_TrueLabel = []
    Edge_TrueLabel = []
    if Construct_Type is None:
        Construct_Type = cf.G_ConstructType
    if cut_off_threshold is None:
        cut_off_threshold = cf.CutOffContact
    for i in range(0, num_Vertex - 1):
        for j in range(i + 1, num_Vertex):
            vec = [M[i][j] for M in DisMxs]
            mi = np.amin(vec)
            ma = np.amax(vec)
            va = np.var(vec)
            wei = math.exp(-alpha * va)
            if Construct_Type == 0:
                wei = 0
                is_edge = False
                for val in vec:
                    if val <= cut_off_threshold:
                        wei += 1
                        is_edge = True
                if is_edge:
                    ed.append((i, j))
                    w.append(wei)
            elif Construct_Type == 1:
                if mi <= cut_off_threshold:
                    ed.append((i, j))
                    w.append(wei)
            elif Construct_Type == 2:
                if ma <= cut_off_threshold:
                    ed.append((i, j))
                    w.append(wei)
            else:
                print ('Unknown graph construction type!!! ' + str(Construct_Type))
    G = D_Graph(n = num_Vertex, edges = ed, edge_attrs={"weight":w})
    if Membership is not None:
        for V in G.vs:
            V['TrueLabel'] = Membership[V.index]
            V['Cluster'] = [V.index]
        for E in G.es:
            E['TrueLabel'] = 1 if G.vs[E.source]['TrueLabel'] == G.vs[E.target]['TrueLabel'] else -1
            E['Connectivity'] = 1
    return G

#-------------------------------------------------------------------------------------
def ConstructProteinGraph(DisMxs, Membership, Construct_Type=None, cut_off_threshold = None):
    G = ConstructGraph(DisMxs, Construct_Type, cut_off_threshold)
    #partitions, _ = partition_gievenNumPar(G, 20)
    partitions = [[v.index] for v in G.vs]
    Memmership_Idx = set(Membership)
    Memmership_Idx = list(Memmership_Idx)
    Memmership_Idx.sort()
    # vertex attributes: Idx of cluster, label of cluster vertex
    Idx_Of_Cluster = []
    Ver_Label = []
    ClusterTupleLabel = []
    for p in partitions:
        membership_of_p = [Membership[i] for i in p]
        Idx_Of_Cluster.append(tuple(p))
        ClusterTupleLabel.append([membership_of_p.count(i) for i in Memmership_Idx])
        Ver_Label.append(np.bincount(membership_of_p).argmax())
    # construct edge of Graph
    eds = []
    Ed_Conn = []
    EdgeLabel = []
    for i, p1 in enumerate(partitions):
        for j, p2 in enumerate(partitions):
            if i < j:
                conn = G.checkConnectivity(p1, p2)
                if conn > 0:
                    eds.append((i, j))
                    Ed_Conn.append(conn)
                    EdgeLabel.append(1 if Ver_Label[i] == Ver_Label[j] else -1)
    return D_Graph(len(Ver_Label), eds, vertex_attrs={'Cluster': Idx_Of_Cluster, "ClusterTupleLabel": ClusterTupleLabel,
                                                      'TrueLabel': Ver_Label},
                   edge_attrs={'Connectivity': Ed_Conn, 'TrueLabel': EdgeLabel})

def ConstructPerfectGraph(DisMxs, Membership, Construct_Type=None, cut_off_threshold = None):
    G = ConstructGraph(DisMxs, Construct_Type, cut_off_threshold)
    NumVerPerCluster = len(G.vs)/20
    G.vs['OriginalIndex'] = [V.index for V in G.vs]
    Map_Cluster_Label = dict()
    Set_Labels = list(set(Membership))
    for idx_Cluster, lab in enumerate(Set_Labels):
        DomainIdx = tuple([idx for idx, val in enumerate(Membership) if val == lab])
        Sub_Vertex = G.vs[DomainIdx]
        Sub_Graph = G.subgraph(Sub_Vertex)
        nVertex_SubGraph = len(Sub_Graph.vs)/NumVerPerCluster
        if nVertex_SubGraph > 1:
            partitions, _ = partition_gievenNumPar(Sub_Graph, nVertex_SubGraph)
            for p in partitions:
                OrgSubIdx = tuple([Sub_Graph.vs[i]['OriginalIndex'] for i in p])
                Map_Cluster_Label.update({OrgSubIdx:idx_Cluster})
        else:
            Map_Cluster_Label.update({DomainIdx:idx_Cluster})

    '''Construct edges of graph'''
    eds = []
    Ed_Connectivity = []
    Idx_Of_Cluster = []
    Ver_Label = []
    for idx1, d1 in enumerate(Map_Cluster_Label.iteritems()):
        Clst1, lab1 = d1[0], d1[1]
        Idx_Of_Cluster.append(tuple(Clst1))
        Ver_Label.append(lab1)
        for idx2, d2 in enumerate(Map_Cluster_Label.iteritems()):
            Clst2, lab2 = d2[0], d2[1]
            if idx2 > idx1 and set(Clst1) != set(Clst2):
                if G.checkConnectivity(Clst1, Clst2) > 0:
                    eds.append((idx1, idx2))
                    Ed_Connectivity.append(G.checkConnectivity(Clst1, Clst2))
    G2 = D_Graph(len(Idx_Of_Cluster), eds, vertex_attrs ={'Cluster':Idx_Of_Cluster, 'TrueLabel':Ver_Label}, edge_attrs ={'Connectivity':Ed_Connectivity})
    return G2



def ConstructRealClusterGraph(DisMxs, Membership, Construct_Type=None, cut_off_threshold = None, init_membership = None):
    from GraphPackage.Graph_Config import NumberOfClusterInClusteringGraph
    G = ConstructGraph(DisMxs, Construct_Type, cut_off_threshold)

    val = 0.75
    wei = [1]*len(G.es)
    if init_membership is not None:
        wei = [1 if init_membership[e.source] == init_membership[e.target] else val for e in G.es]

    partitions, _ = partition_gievenNumPar(G, NumPar=NumberOfClusterInClusteringGraph,edge_weight_factors = wei)
    Memmership_Idx = set(Membership)
    Memmership_Idx = list(Memmership_Idx)
    Memmership_Idx.sort()
    # vertex attributes: Idx of cluster, label of cluster vertex
    Idx_Of_Cluster = []
    Ver_Label      = []
    ClusterTupleLabel = []
    for p in partitions:
        membership_of_p = [Membership[i] for i in p]
        Idx_Of_Cluster.append(tuple(p))
        ClusterTupleLabel.append([membership_of_p.count(i) for i in Memmership_Idx])
        Ver_Label.append(np.bincount(membership_of_p).argmax())
    # construct edge of Graph
    eds = []
    Ed_Conn = []
    EdgeLabel = []
    for i, p1 in enumerate(partitions):
        for j, p2 in enumerate(partitions):
            if i < j:
                conn = G.checkConnectivity(p1, p2)
                if conn > 0:
                    eds.append((i,j))
                    Ed_Conn.append(conn)
                    EdgeLabel.append(1 if Ver_Label[i]==Ver_Label[j] else -1)
    return D_Graph(len(Ver_Label), eds, vertex_attrs={'Cluster':Idx_Of_Cluster, "ClusterTupleLabel":ClusterTupleLabel, 'TrueLabel':Ver_Label},
                   edge_attrs={'Connectivity':Ed_Conn, 'TrueLabel':EdgeLabel})






def partition_gievenNumPar(G, NumPar= None, edge_weight_factors = None):
    if NumPar is None:
        if 15 <= len(G.vs)/10 <= 30:
            NumPar = len(G.vs)/10
        else:
            NumPar = 20

    low = 0.001
    high = 0.75
    count = 0
    thres = None
    w = G.es['weight']
    if edge_weight_factors is not None:
        w = [a * b for a, b in zip(w, edge_weight_factors)]
    partitions = None
    if NumPar <=0 or NumPar > len(G.vs):
        print ("Numpar {} is wrong number".format(str(NumPar)))
        sys.exit()
    while True:
        thres = (low + high)/2
        partitions = louvain.find_partition(G, partition_type=louvain.CPMVertexPartition,
                weights=w,resolution_parameter = thres)
        count += 1
        if np.abs(len(partitions) - NumPar) ==0 or count > 30:
            break
        elif len(partitions) > NumPar:
            high = thres
            thres = (low + high)/2
        else:
            low = thres
            thres = (low + high)/2
    return (partitions, thres)

def writeArrayG2File(Arr_G, Path2File):
    from Utils.MyIO import WriteList2File
    L = []
    for G in Arr_G:
        line = '\t'.join([str(v['OriginalIndex']) for v in G.vs])
        L.append(line)
    WriteList2File(Path2File, L)


