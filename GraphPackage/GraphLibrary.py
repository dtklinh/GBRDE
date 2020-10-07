import igraph as ig
from igraph import Graph
import numpy as np
import sys, os
from GraphPackage.AssistantObjects import Feature_Workhouse
from GraphPackage.AssistantObjects import DynDomEntry
from csb.bio.utils import rmsd
from mainPackage.PathAndDir import Dir2TmpFile
#from GraphPackage.Graph_Config import Path2ViterbiJar
from mainPackage.PathAndDir import Path2ViterbiJar
from Utils.MyIO import ReadViterbiOutFile

class D_Graph(Graph, object):
    def __init__(self, n, edges, directed=False, graph_attrs={}, vertex_attrs={}, edge_attrs={}):
        super(D_Graph, self).__init__(n=n, edges = edges, directed=directed,
                                      graph_attrs=graph_attrs, vertex_attrs = vertex_attrs, edge_attrs=edge_attrs)

    def constructLineGraph(self, ConsiderTwoEndNodesConn=False):
        N = len(self.es)
        eds = []
        for E1 in self.es:
            e1_v1 = E1.source
            e1_v2 = E1.target
            S1 = set((e1_v1, e1_v2))
            for E2 in self.es:
                if E1.index < E2.index:
                    e2_v1 = E2.source
                    e2_v2 = E2.target
                    S2 = set((e2_v1, e2_v2))
                    v_midd_idx = list(S1.intersection(S2))
                    if len(S1.intersection(S2)) == 1:
                        v_midd_idx = list(S1.intersection(S2))[0]
                        v_left, v_right = S1.union(S2) - S1.intersection(S2)
                        if ConsiderTwoEndNodesConn==True:
                            if self.get_eid(v_left, v_right, directed=False, error=False)<0:
                                eds.append((E1.index, E2.index))
                        else:
                            eds.append((E1.index, E2.index))
        return (N, eds)


    def printscreen(self):
        print (', '.join('({}, {})'.format(str(E.source), str(E.target)) for E in self.es))
    def writeSynthesisFeature(self, SeedNum, Path2SynthesisFeature):
        return None
        '''Lines = []
        np.random.seed(SeedNum)
        alpha_wth, beta_wth = 2.0, 0.5
        alpha_acr, beta_acr = 2.0, 1.0
        alpha_endnodes, beta_endnodes = 1.0, 1.0
        rv_wth = InverseGamma(alpha=alpha_wth, beta=beta_wth)
        rv_acr = InverseGamma(alpha=alpha_acr, beta=beta_acr)
        rv_endnodes = InverseGamma(alpha=alpha_endnodes, beta=beta_endnodes)
        for V1 in self.vs:
            for V2 in self.vs:
                if V1.index < V2.index:
                    EdID = self.get_eid(V1.index, V2.index, directed=False, error=False)
                    if EdID >=0:
                        if V1['ClusterID'] == V2['ClusterID']:
                            Lines.append([V1.index, V2.index, 1, rv_wth.random(1)[0]])
                        else:
                            Lines.append([V1.index, V2.index, -1, rv_endnodes.random(1)[0]])
                    else:
                        Set1 = set(G.neighbors(V1))
                        Set2 = set(G.neighbors(V2))
                        if len(Set1.intersection(Set2)) > 0:
                            if V1['ClusterID']==V2['ClusterID']:
                                Lines.append([V1.index, V2.index, 2, rv_wth.random(1)[0]])
                            else:
                                Lines.append([V1.index, V2.index, -2, rv_endnodes.random(1)[0]])
        from MyIO import WriteNumericalList2File
        WriteNumericalList2File(Path2SynthesisFeature, Lines)'''
    def checkConnectivity(self, Cluster1, Cluster2):
        eids = [(i,j) for i in Cluster1 for j in Cluster2]
        ArrIdxs = self.get_eids(eids, directed=False, error=False)
        return sum(1 for i in ArrIdxs if i >=0)
    def get_TwoEndNode(self, Ed1_idx, Ed2_idx):
        Set1 = set([self.es[Ed1_idx].source, self.es[Ed1_idx].target])
        Set2 = set([self.es[Ed2_idx].source, self.es[Ed2_idx].target])
        Set_tmp = Set1.union(Set2) - Set1.intersection(Set2)
        Lst = list(Set_tmp)
        Lst.sort()
        if len(Lst) == 2:
            return Lst[0], Lst[1]
        else:
            return None, None
    def isEdgeTouch(self, ed_id1, ed_id2):
        Set1 = set([self.es[ed_id1].source, self.es[ed_id1].target])
        Set2 = set([self.es[ed_id2].source, self.es[ed_id2].target])
        Set_tmp = Set1.union(Set2) - Set1.intersection(Set2)
        Set_tmp = list(Set_tmp)
        if len(Set_tmp) != 2:
            return 0 # two edges do not touch
        elif self.get_eid(Set_tmp[0], Set_tmp[1], directed=False, error=False) >= 0:
            return 1 # touch but two end nodes is edge
        else:
            return 2 # touch and not an edge
    def calc_squareMatFeature(self, DistanceMx):
        ''' vertex has attribute Cluster '''
        Feature_Generator = Feature_Workhouse(DistanceMx)
        N = len(self.vs)
        Mx = np.zeros((2,N,N)) # Mx[0]: mean of variance, Mx[1]: median of variance
        for E in self.es:
            fea = Feature_Generator.getMeanVar(self.vs[E.source]['Cluster'], self.vs[E.target]['Cluster'])
            fea_2 = Feature_Generator.getMedianVar(self.vs[E.source]['Cluster'], self.vs[E.target]['Cluster'])
            Mx[0, E.source, E.target] = fea
            Mx[0, E.target, E.source] = fea
            Mx[1, E.source, E.target] = fea_2
            Mx[1, E.target, E.source] = fea_2
            for E2 in self.es:
                if E.index < E2.index:
                    i1, i2 = self.get_TwoEndNode(E.index, E2.index)
                    if i1 is None or i2 is None: continue
                    if self.get_eid(i1, i2, directed=False, error=False) < 0:
                        fea_extra = Feature_Generator.getMeanVar(self.vs[i1]['Cluster'], self.vs[i2]['Cluster'])
                        fea_extra_2 = Feature_Generator.getMedianVar(self.vs[i1]['Cluster'], self.vs[i2]['Cluster'])
                        Mx[0, i1, i2] = fea_extra
                        Mx[0, i2, i1] = fea_extra
                        Mx[1, i1, i2] = fea_extra_2
                        Mx[1, i2, i1] = fea_extra_2
        return Mx
    #------------------------------------------------
    def do_work(self):
        from PathAndDir import Dir2TmpFile
        # create Line Graph and run viterbi to detect cross-boder edges
        #Entry = self['DynDomEntry']
        #SquareMatFeature = self.calc_squareMatFeature(Entry.DistanceMatrices)
        SquareMatFeature = self['SquareMatFeature']
        N_LineG, eds_LineG = self.constructLineGraph(ConsiderTwoEndNodesConn=True)
        LineG = D_LineGraph(N_LineG, eds_LineG, graph_attrs={'OriginalGraph': self, 'SquareMatFeature': SquareMatFeature})
        LineG.setFeature_Sim(SquareMatFeature, self['T'][0], self['T'][1])

        Path2ViterbiFeature = os.path.join(Dir2TmpFile, '{}_vi_fea.txt'.format(str(self['serial'])))
        LineG.WriteFeature(Path2ViterbiFeature)
        Path2ViterbiOutFile = os.path.join(Dir2TmpFile, '{}_out.txt'.format(str(self['serial'])))
        LineG.runViterbi(Path2ViterbiJar, Path2ViterbiFeature, Path2ViterbiOutFile)
        PredLabels_Edges = ReadViterbiOutFile(Path2ViterbiOutFile)
        tmpG = self.copy()
        Pair_idx_ToBeRemoved = [(self.es[idx].source, self.es[idx].target) for idx, val in enumerate(PredLabels_Edges) if val == 0]
        tmpG.delete_edges(Pair_idx_ToBeRemoved)
        Clsts = tmpG.clusters()
        Clsts = list(Clsts)
        Arr = []
        for C in Clsts:
            Arr.append(tmpG.subgraph(C))
        return Arr
    def do_work_2(self):
        from mainPackage.PathAndDir import Dir2TmpFile
        # create Line Graph and run viterbi to detect cross-boder edges
        # Entry = self['DynDomEntry']
        # SquareMatFeature = self.calc_squareMatFeature(Entry.DistanceMatrices)
        SquareMatFeature = self['SquareMatFeature']
        N_LineG, eds_LineG = self.constructLineGraph(ConsiderTwoEndNodesConn=True)
        LineG = D_LineGraph(N_LineG, eds_LineG,
                            graph_attrs={'OriginalGraph': self, 'SquareMatFeature': SquareMatFeature})
        LineG.setFeature_Sim_2(SquareMatFeature, self['P'][0], self['P'][1])

        Path2ViterbiFeature = os.path.join(Dir2TmpFile, '{}_vi_fea.txt'.format(str(self['serial'])))
        LineG.WriteFeature(Path2ViterbiFeature)
        Path2ViterbiOutFile = os.path.join(Dir2TmpFile, '{}_out.txt'.format(str(self['serial'])))
        LineG.runViterbi(Path2ViterbiJar, Path2ViterbiFeature, Path2ViterbiOutFile)
        PredLabels_Edges = ReadViterbiOutFile(Path2ViterbiOutFile)
        tmpG = self.copy()
        #tmpG = self.as_undirected()
        #tmpG = self.deepcopy()
        Pair_idx_ToBeRemoved = [(self.es[idx].source, self.es[idx].target) for idx, val in enumerate(PredLabels_Edges)
                                if val == 0]
        tmpG.delete_edges(Pair_idx_ToBeRemoved)
        Clsts = tmpG.clusters()
        Clsts = list(Clsts)
        Arr = []
        for C in Clsts:
            Arr.append(tmpG.subgraph(C))
        return Arr
    #-------------------------------------
    def do_work_iteration(self, rmsd_thres=5.0):
        t1, t2 = 4, 8
        Arr_G = []
        #Entry = self['DynDomEntry']
        self['T'] = [t1,t2]
        from collections import deque
        stack = deque()
        init_arr_G = []
        Cls = self.clusters()
        Cls = list(Cls)
        Entry = self['DynDomEntry']
        if len(Cls)==1:
            stack.append(self)
        else:
            for C in Cls:
                G_i = self.subgraph(C)
                G_i['SquareMatFeature'] = G_i.calc_squareMatFeature(Entry.DistanceMatrices)
                stack.append(G_i)
        while len(stack) != 0:
            tmpG = stack.pop()
            if tmpG.rmsd() < rmsd_thres or len(tmpG.es)==0:
                Arr_G.append(tmpG)
                #print ('My T: {}'.format(str(tmpG['T'])))
            else:
                tmp_Arr = tmpG.do_work()
                if len(tmp_Arr) == 1 and len(tmp_Arr[0].es) == len(tmpG.es):
                    if tmpG['T'][0]==0.1 and tmpG['T'][1] ==0.1:
                        Arr_G.extend(tmp_Arr)
                        #print ('ADD XXXX')
                        continue
                    else:
                        t1, t2 = tmpG['T']
                        tmp_Arr[0]['T'] = [max(0.1, t1-0.1), max(0.1, t2-0.1)]

                for g in tmp_Arr:
                    g['SquareMatFeature'] = g.calc_squareMatFeature(Entry.DistanceMatrices)
                    stack.append(g)
        return Arr_G
    #--------------------------------------------
    def do_work_iteration_2(self, rmsd_thres=3.5):
        p1, p2 = 0.0, 0.0
        step = 0.05
        Arr_G = []
        # Entry = self['DynDomEntry']
        self['P'] = [p1, p2]
        from collections import deque
        stack = deque()
        init_arr_G = []
        Cls = self.clusters()
        Cls = list(Cls)
        Entry = self['DynDomEntry']
        if len(Cls) == 1:
            stack.append(self)
        else:
            for C in Cls:
                G_i = self.subgraph(C)
                G_i['SquareMatFeature'] = G_i.calc_squareMatFeature(Entry.DistanceMatrices)
                stack.append(G_i)
        while len(stack) != 0:
            tmpG = stack.pop()
            if tmpG.rmsd() < rmsd_thres or len(tmpG.es) == 0:
                Arr_G.append(tmpG)
                #print ('My P: {}'.format(str(tmpG['P'])))
            else:
                tmp_Arr = tmpG.do_work_2()
                if len(tmp_Arr) == 1 and len(tmp_Arr[0].es) == len(tmpG.es):
                    if tmpG['P'][0] == 0.9 and tmpG['P'][1] == 0.9:
                        Arr_G.extend(tmp_Arr)
                        print ('ADD XXXX')
                        continue
                    else:
                        t1, t2 = tmpG['P']
                        tmp_Arr[0]['P'] = [t1 + step, t2 + step]    #[max(0.1, t1 - 0.1), max(0.1, t2 - 0.1)]

                for g in tmp_Arr:
                    g['SquareMatFeature'] = g.calc_squareMatFeature(Entry.DistanceMatrices)
                    stack.append(g)
        return Arr_G
    #---------------------------------
    def rmsd(self, Vertex_IDs = None):
        if Vertex_IDs is None:
            Vertex_IDs = range(len(self.vs))
        Entry = self['DynDomEntry']
        X = Entry.X
        Prot_IDs = [idx for i in Vertex_IDs for idx in self.vs[i]['Cluster']]
        X_Part = X[:, Prot_IDs, :]
        N = X_Part.shape[0]
        Arr = []
        for idx1, x1 in enumerate(X_Part):
            for idx2, x2 in enumerate(X_Part):
                if idx1 < idx2:
                    Arr.append(rmsd(x1, x2))

        return np.max(Arr)
    def find2SubGraph(self, thres, Clsts):
        ''' search 2 sub graph which are most likely to be in 1 rigid domain '''
    #---------------------------
        #Clsts = self.clusters()
        # -----------------------
        Entry = self['DynDomEntry']
        X = Entry.X
        # -----------------------
        if len(Clsts) <= 2:
            print ('Num of Par is too small to be merged!\n Return original graph')
            return (None, None, -1)
        Clsts_ProtG = []
        val_max = -1
        c1_max, c2_max = [], []
        for idx1, C1 in enumerate(Clsts):
            for idx2, C2 in enumerate(Clsts):
                if idx1 < idx2:
                    ProtG_C1 = [idx for i in C1 for idx in self.vs[i]['Cluster']]
                    ProtG_C2 = [idx for i in C2 for idx in self.vs[i]['Cluster']]
                    ProtG_C12 = [idx for i in C1 + C2 for idx in self.vs[i]['Cluster']]
                    X_1 = X[:, ProtG_C1, :]
                    X_2 = X[:, ProtG_C2, :]
                    X_12 = X[:, ProtG_C12, :]
                    #if rmsd(X_1[0], X_1[1]) > rmsd(X_12[0], X_12[1]) or rmsd(X_2[0], X_2[1]) > rmsd(X_12[0], X_12[1]):
                     #   continue
                    val = (rmsd(X_1[0], X_1[1]) + rmsd(X_2[0], X_2[1])) / rmsd(X_12[0], X_12[1])
                    if val > val_max:
                        val_max = val
                        c1_max = C1
                        c2_max = C2
        if val_max > thres:
            return (c1_max, c2_max, val_max)
        else:
            return (None, None, val_max)
        #--------------------------------------------
    def merge(self, thres, Pair_VerIdx_ToBeRemoved):
        # clone graph to create cluster
        '''
        tmpG = self.copy()
        tmpG.delete_edges(Pair_VerIdx_ToBeRemoved)
        Clsts = tmpG.clusters()
        #-----------------------
        Entry = self['DynDomEntry']
        X = Entry.X
        #-----------------------
        if len(Clsts) <= 2:
            print 'Num of Par is too small to be merged!\n Return original graph'
            return self.copy()
        Clsts_ProtG = []
        val_max = -1
        c1_max, c2_max = [], []
        for idx1, C1 in enumerate(Clsts):
            for idx2, C2 in enumerate(Clsts):
                if idx1 < idx2 and np.any([True if (c1,c2) in Pair_VerIdx_ToBeRemoved or (c2, c1) in Pair_VerIdx_ToBeRemoved else False for c1 in C1 for c2 in C2]):
                    ProtG_C1 = [idx for i in C1 for idx in tmpG.vs[i]['Cluster']]
                    ProtG_C2 = [idx for i in C2 for idx in tmpG.vs[i]['Cluster']]
                    ProtG_C12 = [idx for i in C1+C2 for idx in tmpG.vs[i]['Cluster']]
                    X_1 = X[:, ProtG_C1, :]
                    X_2 = X[:, ProtG_C2, :]
                    X_12 = X[:, ProtG_C12, :]
                    val = (rmsd(X_1[0], X_1[1]) + rmsd(X_2[0], X_2[1]))/rmsd(X_12)
                    if val > val_max:
                        val_max = val
                        c1_max = C1
                        c2_max = C2
                        '''

        tmpG = self.copy()
        tmpG.delete_edges(Pair_VerIdx_ToBeRemoved)
        Clsts = tmpG.clusters()
        Clsts = list(Clsts)
        #tmpG.delete_edges(Pair_VerIdx_ToBeRemoved)
        C1, C2, thres_max = tmpG.find2SubGraph(thres, Clsts)
        Pair_VerIdx_ToBeRemoved_Copy = Pair_VerIdx_ToBeRemoved[:]
        while thres_max > thres and C1 is not None:
            '''
            arr = [(c1, c2) if c1 < c2 else (c2,c1) for c1 in C1 for c2 in C2]
            for item in arr:
                if item in Pair_VerIdx_ToBeRemoved_Copy:
                    Pair_VerIdx_ToBeRemoved_Copy.remove(item) '''
            Clsts.remove(C1)
            Clsts.remove(C2)
            Clsts.append(C1+C2)
            #tmpG = self.copy()
            #tmpG.delete_edges(Pair_VerIdx_ToBeRemoved_Copy)
            C1, C2, thres_max = self.find2SubGraph(thres, Clsts)
        return Clsts
    #--------------------------------------------


    def do_merge(self,thres, Arr_G):
        tmpG = self.copy()
        Cls = []
        for g in Arr_G:
            #Cls.append([v['OriginalIndex'] for v in g.vs])
            Cls.append([v2.index for v in g.vs for v2 in self.vs if v['OriginalIndex']==v2['OriginalIndex']])
        C1, C2, thres_max = tmpG.find2SubGraph(thres, Cls)
        #Pair_VerIdx_ToBeRemoved_Copy = Pair_VerIdx_ToBeRemoved[:]
        while thres_max > thres and C1 is not None:
            '''
            arr = [(c1, c2) if c1 < c2 else (c2,c1) for c1 in C1 for c2 in C2]
            for item in arr:
                if item in Pair_VerIdx_ToBeRemoved_Copy:
                    Pair_VerIdx_ToBeRemoved_Copy.remove(item) '''
            Cls.remove(C1)
            Cls.remove(C2)
            Cls.append(C1 + C2)
            # tmpG = self.copy()
            # tmpG.delete_edges(Pair_VerIdx_ToBeRemoved_Copy)
            C1, C2, thres_max = self.find2SubGraph(thres, Cls)
        Arr_G_New = []
        for C in Cls:
            try:
                Arr_G_New.append(self.subgraph(C))
            except:
                print ('$$$$')
                print (','.join([str(v.index) for v in self.vs]))
                print (Cls)
                print (C)
        return Arr_G_New
    #---------------------------------------------
    def do_merge_2(self, thres, Arr_G):
        tmpG = self.copy()
        Cls = []
        for g in Arr_G:
            # Cls.append([v['OriginalIndex'] for v in g.vs])
            Cls.append([v2.index for v in g.vs for v2 in self.vs if v['OriginalIndex'] == v2['OriginalIndex']])
        if len(Arr_G)<=2:
            return Arr_G
        min_rmsd = 1000.0
        c1_new, c2_new = None, None
        for idx1, c1 in enumerate(Cls):
            for idx2, c2 in enumerate(Cls):
                if idx1 < idx2:
                    C12_rmsd = self.rmsd(c1+c2)
                    if min_rmsd > C12_rmsd:
                        min_rmsd = C12_rmsd
                        c1_new, c2_new = c1, c2
        Cls.remove(c1_new)
        Cls.remove(c2_new)
        Cls.append(c1_new + c2_new)
        Arr_G_New = []
        for C in Cls:
            try:
                Arr_G_New.append(self.subgraph(C))
            except:
                print ('$$$$')
                print (','.join([str(v.index) for v in self.vs]))
                print (Cls)
                print (C)
        return Arr_G_New
    #-------------------------------------------------
    def split(self, Vertex_Idxs= None, Edge_Idxs = None):
        Entry = self['DynDomEntry']
        if Vertex_Idxs is None:
            Vertex_Idxs = [v.index for v in self.vs]
        if Edge_Idxs is None:
            Edge_Idxs = [e.index for e in self.es]


        SquareMatFeature = self['SquareMatFeature'] #.calc_squareMatFeature(Entry.DistanceMatrices)
        N_LineG, eds_LineG = self.constructLineGraph()
        LineG = D_LineGraph(N_LineG, eds_LineG, graph_attrs={'OriginalGraph': self, 'SquareMatFeature': SquareMatFeature})
        LineG.setFeature(SquareMatFeature)
        Path2ViterbiFeature = os.path.join(Dir2TmpFile, '{}.viterbi_fea.txt'.format(str(0)))
        LineG.WriteFeature(Path2ViterbiFeature)
        Path2ViterbiOutFile = os.path.join(Dir2TmpFile, '{}.viterbi_out.txt'.format(str(0)))
        LineG.runViterbi(Path2ViterbiJar, Path2ViterbiFeature, Path2ViterbiOutFile)
        PredLabels_Edges = ReadViterbiOutFile(Path2ViterbiOutFile)
        tmpG = self.copy()
        Pair_idx_ToBeRemoved_IndexAdjust = [(tmpG.vs[tmpG.es[idx].source]['OriginalIndex'], tmpG.vs[tmpG.es[idx].target]['OriginalIndex'])
                                for idx, val in enumerate(PredLabels_Edges) if val == 0]
        Pair_idx_ToBeRemoved = [(tmpG.es[idx].source, tmpG.es[idx].target) for idx, val in enumerate(PredLabels_Edges) if val == 0]
        tmpG.delete_edges(Pair_idx_ToBeRemoved)
        Cls = tmpG.clusters()
        Cls = list(Cls)
        return Cls, Pair_idx_ToBeRemoved_IndexAdjust

    def split_complete(self, Pair_Idx_ToBeRemoved):
        Total_RMSD = self.rmsd()
        tmpG = self.copy()
        tmpG.delete_edges(Pair_Idx_ToBeRemoved)
        Cls = tmpG.clusters()
        Cls = list(Cls)
        for C in Cls:
            local_rmsd = self.rmsd(C)
            if local_rmsd > 3.5 and local_rmsd/Total_RMSD > 0.35:
                sub_G = tmpG.subgraph(C)
                _, PP = sub_G.split()
                if len(PP) != 0:
                    Pair_Idx_ToBeRemoved.extend(PP)
        return Pair_Idx_ToBeRemoved
    #------------------------------------
    def num_mistake(self):
        if "ClusterTupleLabel" in self.vs.attribute_names():
            return sum([sum(V["ClusterTupleLabel"]) - max(V["ClusterTupleLabel"]) for V in self.vs])
        else:
            print ('There is no ClusterTupleLabel')
            return None


class D_LineGraph(D_Graph, object):
    def __init__(self, n, edges, directed = False, graph_attrs={}, vertex_attrs={}, edge_attrs={}):
        super(D_LineGraph, self).__init__(n= n, edges= edges, directed=directed, graph_attrs=graph_attrs, vertex_attrs=vertex_attrs,
                                          edge_attrs=edge_attrs)
        #self.OriginalGraph = OriginalGraph
        self.initial()
    def initial(self):
        EdgeWeights = []
        Edge_TwoEndNodeConnected = []
        Edge_MidleNode = []
        Edge_TwoEndNodes = []
        if 'TRUELABEL' in ( s.upper() for s in self['OriginalGraph'].es.attribute_names()):
            Vertexlabels = [V["TrueLabel"] for V in self['OriginalGraph'].es]
        elif 'TRUELABEL' in (s.upper() for s in self['OriginalGraph'].vs.attribute_names()):
            Vertexlabels = [1 if self['OriginalGraph'].vs[e.source]['TrueLabel'] == self['OriginalGraph'].vs[e.target]['TrueLabel']
                            else -1 for e in self['OriginalGraph'].es]
        for E in self.es:
            e1_org, e2_org = self['OriginalGraph'].es[E.source], self['OriginalGraph'].es[E.target]
            S1 = set((e1_org.source, e1_org.target))
            S2 = set((e2_org.source, e2_org.target))
            L_intersect = list(S1.intersection(S2))
            if len(L_intersect)==1:
                Edge_MidleNode.append(L_intersect[0])
                L_TwoEndNodes = list(S1.union(S2) - S1.intersection(S2))
                L_TwoEndNodes.sort()
                if len(L_TwoEndNodes) != 2:
                    print ("Sth is wrong"); sys.exit(1)
                E
                Edge_TwoEndNodes.append((L_TwoEndNodes[0], L_TwoEndNodes[1]))
                Edge_TwoEndNodeConnected.append(True if self['OriginalGraph'].get_eid(L_TwoEndNodes[0], L_TwoEndNodes[1],
                                                                                   directed=False, error=False)>=0 else False)
                EdgeWeights.append(1.0/(self['OriginalGraph'].degree(self['OriginalGraph'].vs[L_intersect[0]])-1))
        self.vs["TrueLabel"] = Vertexlabels
        self.es["weight"] = EdgeWeights
        self.es["TwoEndNodeConnected"] = Edge_TwoEndNodeConnected
        self.es["MidNode_ID"] = Edge_MidleNode
        self.es["TwoEndNodes_ID"] = Edge_TwoEndNodes
        self.alpha = 10.0
    def refine_parameters(self):
        Edge_TwoEndNodeConnected = []
        Edge_MidleNode = []
        Edge_TwoEndNodes = []
        EdgeWeights = []
        for E in self.es:
            e1_org, e2_org = self['OriginalGraph'].es[E.source], self['OriginalGraph'].es[E.target]
            S1 = set((e1_org.source, e1_org.target))
            S2 = set((e2_org.source, e2_org.target))
            L_intersect = list(S1.intersection(S2))
            if len(L_intersect)==1:
                Edge_MidleNode.append(L_intersect[0])
                L_TwoEndNodes = list(S1.union(S2) - S1.intersection(S2))
                L_TwoEndNodes.sort()
                if len(L_TwoEndNodes) != 2:
                    print ("Sth is wrong"); sys.exit(1)
                E
                Edge_TwoEndNodes.append((L_TwoEndNodes[0], L_TwoEndNodes[1]))
                Edge_TwoEndNodeConnected.append(True if self['OriginalGraph'].get_eid(L_TwoEndNodes[0], L_TwoEndNodes[1],
                                                                                   directed=False, error=False)>=0 else False)
                EdgeWeights.append(1.0/(self['OriginalGraph'].degree(self['OriginalGraph'].vs[L_intersect[0]])-1))
        self.es["weight"] = EdgeWeights
        self.es["TwoEndNodeConnected"] = Edge_TwoEndNodeConnected
        self.es["MidNode_ID"] = Edge_MidleNode
        self.es["TwoEndNodes_ID"] = Edge_TwoEndNodes

    def setFeature_Sim(self, SquareMatFeature, t1 = 1, t2 = 2):
        import heapq
        G = self['OriginalGraph']
        Arr_Ver = [SquareMatFeature[0,e.source, e.target] for e in G.es]
        Arr_Edge = [SquareMatFeature[0,e["TwoEndNodes_ID"][0], e["TwoEndNodes_ID"][1]] for e in self.es if e["TwoEndNodeConnected"]==False]
        # Arr_Ver = [SquareMatFeature[0,G.vs[e.source]['OriginalIndex'], G.vs[e.target]['OriginalIndex']] for e in G.es]
        # Arr_Edge = []
        # for e in self.es:
        #     if e["TwoEndNodeConnected"]==False:
        #         idx1, idx2 = G.vs[e["TwoEndNodes_ID"][0]]['OriginalIndex'], G.vs[e["TwoEndNodes_ID"][1]]['OriginalIndex']
        #         Arr_Edge.append(SquareMatFeature[0, idx1, idx2])

        Arr_Ver = np.array(Arr_Ver)
        Arr_Edge = np.array(Arr_Edge)
        from Outlier_Detection.OutlierLib import is_outlier
        from Outlier_Detection.Config import Edge_outlier_Thres, Vertex_outlier_Thres

        ol_ver_idx = is_outlier(Arr_Ver, Vertex_outlier_Thres)
        outlier_Ver = Arr_Ver[ol_ver_idx]
        normal_Ver = Arr_Ver[np.logical_not(ol_ver_idx)]
        additional_Ver = heapq.nlargest(int(max(1,len(outlier_Ver))/t1), normal_Ver)
        outlier_Ver = np.append(outlier_Ver, additional_Ver)

        ol_edge_idx = is_outlier(Arr_Edge, Edge_outlier_Thres)
        outlier_Edge = Arr_Edge[ol_edge_idx]
        normal_Edge = Arr_Edge[np.logical_not(ol_edge_idx)]
        additional_Edge = heapq.nlargest(int(max(1,len(outlier_Edge))/t2), normal_Edge)
        outlier_Edge = np.append(outlier_Edge, additional_Edge)

        for V in self.vs:
            v1_idx, v2_idx = G.es[V.index].source, G.es[V.index].target
            if SquareMatFeature[0, v1_idx, v2_idx] in outlier_Ver:
                V['score'] = [1, -1]
            else:
                V['score'] = [-1, 1]
        for E in self.es:
            w = E['weight']
            w = 1.0
            #w = 1
            if E["TwoEndNodeConnected"] == True:
                E['score'] = [0]*4
            else:
                v1_end_idx, v2_end_idx = E["TwoEndNodes_ID"]
                v_mid_idx = E["MidNode_ID"]
                V = SquareMatFeature[0,v1_end_idx, v2_end_idx]
                V1 = SquareMatFeature[0,v1_end_idx, v_mid_idx]
                V2 = SquareMatFeature[0,v2_end_idx, v_mid_idx]
                if V in outlier_Edge: # V is Big
                    if V1 in outlier_Ver and V2 in outlier_Ver:
                        E['score'] = [w, -w, -w, -w]
                    elif V1 in outlier_Ver and V2 not in outlier_Ver:
                        E['score'] = [-w, w, -w, -w]
                    elif V1 not in outlier_Ver and V2 in outlier_Ver:
                        E['score'] = [-w, -w, w, -w]
                    else:
                        if V1 < V2:
                            E['score'] = [-w, -w, w, -w]
                        else:
                            E['score'] = [-w, w, -w, -w]
                else: #V is small
                    if V1 not in outlier_Ver and V2 not in outlier_Ver:
                        E['score'] = [-w, -w, -w, w]
                    elif V1 in outlier_Ver and V2 in outlier_Ver:
                        E['score'] = [w, -w, -w, -w]
                    else:
                        E['score'] = [0] * 4

        #-------------------------------------------
    def setFeature_Sim_2(self, SquareMatFeature, percent_1=0.05, percent_2=0.05):
        import heapq, math
        G = self['OriginalGraph']
        Arr_Ver = [SquareMatFeature[0, e.source, e.target] for e in G.es]
        Arr_Edge = [SquareMatFeature[0, e["TwoEndNodes_ID"][0], e["TwoEndNodes_ID"][1]] for e in self.es if
                    e["TwoEndNodeConnected"] == False]
        # Arr_Ver = [SquareMatFeature[0,G.vs[e.source]['OriginalIndex'], G.vs[e.target]['OriginalIndex']] for e in G.es]
        # Arr_Edge = []
        # for e in self.es:
        #     if e["TwoEndNodeConnected"]==False:
        #         idx1, idx2 = G.vs[e["TwoEndNodes_ID"][0]]['OriginalIndex'], G.vs[e["TwoEndNodes_ID"][1]]['OriginalIndex']
        #         Arr_Edge.append(SquareMatFeature[0, idx1, idx2])

        Arr_Ver = np.array(Arr_Ver)
        Arr_Edge = np.array(Arr_Edge)
        from Outlier_Detection.OutlierLib import is_outlier
        from Outlier_Detection.Config import Edge_outlier_Thres, Vertex_outlier_Thres

        ol_ver_idx = is_outlier(Arr_Ver, Vertex_outlier_Thres)
        outlier_Ver = Arr_Ver[ol_ver_idx]
        normal_Ver = Arr_Ver[np.logical_not(ol_ver_idx)]
        additional_Ver = heapq.nlargest(int(math.floor(len(normal_Ver)*percent_1)), normal_Ver)
        outlier_Ver = np.append(outlier_Ver, additional_Ver)

        ol_edge_idx = is_outlier(Arr_Edge, Edge_outlier_Thres)
        outlier_Edge = Arr_Edge[ol_edge_idx]
        normal_Edge = Arr_Edge[np.logical_not(ol_edge_idx)]
        additional_Edge = heapq.nlargest(int(math.floor(len(normal_Edge)*percent_2)), normal_Edge)
        outlier_Edge = np.append(outlier_Edge, additional_Edge)

        for V in self.vs:
            v1_idx, v2_idx = G.es[V.index].source, G.es[V.index].target
            if SquareMatFeature[0, v1_idx, v2_idx] in outlier_Ver:
                V['score'] = [1, -1]
            else:
                V['score'] = [-1, 1]
        for E in self.es:
            w = E['weight']
            w = 1.0
            # w = 1
            if E["TwoEndNodeConnected"] == True:
                E['score'] = [0] * 4
            else:
                v1_end_idx, v2_end_idx = E["TwoEndNodes_ID"]
                v_mid_idx = E["MidNode_ID"]
                V = SquareMatFeature[0, v1_end_idx, v2_end_idx]
                V1 = SquareMatFeature[0, v1_end_idx, v_mid_idx]
                V2 = SquareMatFeature[0, v2_end_idx, v_mid_idx]
                if V in outlier_Edge:  # V is Big
                    if V1 in outlier_Ver and V2 in outlier_Ver:
                        E['score'] = [w, -w, -w, -w]
                    elif V1 in outlier_Ver and V2 not in outlier_Ver:
                        E['score'] = [-w, w, -w, -w]
                    elif V1 not in outlier_Ver and V2 in outlier_Ver:
                        E['score'] = [-w, -w, w, -w]
                    else:
                        if V1 < V2:
                            E['score'] = [-w, -w, w, -w]
                        else:
                            E['score'] = [-w, w, -w, -w]
                else:  # V is small
                    if V1 not in outlier_Ver and V2 not in outlier_Ver:
                        E['score'] = [-w, -w, -w, w]
                    elif V1 in outlier_Ver and V2 in outlier_Ver:
                        E['score'] = [w, -w, -w, -w]
                    else:
                        E['score'] = [0] * 4


##--------------------------------------------
    def setFeature(self, squareMatFeature):
        Vec_LG_Edges_TwoEndNotConn = [squareMatFeature[0, E['TwoEndNodes_ID'][0], E['TwoEndNodes_ID'][1]] for E in self.es if
             E["TwoEndNodeConnected"]==False]
        Vec_LG_Edges_TwoEndNotConn = np.matrix(Vec_LG_Edges_TwoEndNotConn).transpose()
        Vec_Diff = [squareMatFeature[0, E["MidNode_ID"], E['TwoEndNodes_ID'][0]] - squareMatFeature[0, E["MidNode_ID"], E['TwoEndNodes_ID'][1]]
                     for E in self.es if E["TwoEndNodeConnected"]==False]
        Vec_Diff = np.matrix(Vec_Diff).transpose()
        from sklearn.preprocessing import QuantileTransformer, MinMaxScaler
        QuanTileScaler = QuantileTransformer()
        QuanTileScaler.fit(Vec_LG_Edges_TwoEndNotConn)
        MMScaler = MinMaxScaler(feature_range=(-2,2))
        MMScaler.fit(Vec_Diff)
        from FuzzyLogic import estimateFuzzyFeature_Node, estimateFuzzyFeature_Edge2
        for V in self.vs:
            idx1, idx2 = self['OriginalGraph'].es[V.index].source, self['OriginalGraph'].es[V.index].target
            v = squareMatFeature[0, idx1, idx2]
            val1, val2 = estimateFuzzyFeature_Node(v, -1), estimateFuzzyFeature_Node(v, 1)
            V['score'] = [val1, val2]
            V['score'] = [0, 0]
        for E in self.es:
            w = E['weight']
            idx1, idx2 = E["TwoEndNodes_ID"]
            idx_mid = E["MidNode_ID"]
            v = squareMatFeature[0, idx1, idx2]
            v_diff = squareMatFeature[0, idx_mid, idx1] - squareMatFeature[0, idx_mid, idx2]
            val00, val01, val10, val11 = None, None, None, None
            if self['OriginalGraph'].get_eid(idx1, idx2, directed=False, error=False) >= 0:
                val00, val01, val10, val11 = 0, 0, 0, 0
                print ("Skip edge {}".format(str(E.index)))
            else:
                v = QuanTileScaler.transform(np.matrix(v))[0, 0]
                v_diff = MMScaler.transform(np.matrix(v_diff))[0, 0]
                val00, val01 = estimateFuzzyFeature_Edge2(v, v_diff, -1, -1), estimateFuzzyFeature_Edge2(v, v_diff, -1, 1)
                val10, val11 = estimateFuzzyFeature_Edge2(v, v_diff, 1, -1), estimateFuzzyFeature_Edge2(v, v_diff, 1, 1)
            E['score'] = [w*val00, w*val01, w*val10, w*val11]
    def WriteFeature(self, Path2OutFile):
        '''matFeature [2 x len(LG_Vertex) x len(LG_Vertex)] is scaled'''
        #print matFeature
        def getNodeFeature():
            fea = np.zeros((0,3))
            for V in self.vs:
                tmp = [V.index]
                tmp.extend(np.exp(V['score']))
                fea = np.append(fea, [tmp], 0)
            return fea
        def getEdgeFeature():
            fea = np.zeros((0,6))
            for E in self.es:
                tmp = [E.source, E.target]
                tmp.extend(np.exp(E['score']))
                fea = np.append(fea, [tmp], 0)
            return fea
        fea_node = getNodeFeature()
        fea_edge = getEdgeFeature()
        Lines = []
        Lines.append("numlabels:")
        Lines.append("2")
        Lines.append("nodes:")
        for i in range(fea_node.shape[0]):
            Lines.append('\t'.join(str(val) for val in fea_node[i,:]))
        Lines.append("edges:")
        for i in range(fea_edge.shape[0]):
            Lines.append('\t'.join(str(val) for val in fea_edge[i,:]))
        from Utils.MyIO import WriteList2File
        WriteList2File(Path2OutFile, Lines)
    def runViterbi(self, Path2ViterbiJar, Path2FeatureFile, Path2OutFile):
        import subprocess
        subprocess.call(['java', '-jar', Path2ViterbiJar, Path2FeatureFile, Path2OutFile, '10000'])
        print ("Finish running Viterbi!")
    def calculateLogScore(self, Labels):
        score = 0
        for V in self.vs:
            lab = Labels[V.index]

            score += V['score'][lab]
        for E in self.es:
            idx1, idx2 = E.source, E.target
            lab1, lab2 = Labels[idx1], Labels[idx2]
            idx = -1
            if lab1 != 1 and lab2 !=1:
                idx = 0
            elif lab1 != 1 and lab2 ==1:
                idx = 1
            elif lab1 == 1 and lab2 != 1:
                idx = 2
            else:
                idx = 3
            score += E['score'][idx]
        return score
    def Gibbs_Sampling(self, n_blocking, init_labels = None, num_run=10000):
        def do_work(label):
            return self.calculateLogScore(label)
        import random, itertools, concurrent.futures
        from scipy.special import logsumexp
        labels = init_labels if init_labels is not None else [1]*len(self.vs) #0 as -1
        Arr_Scores = []
        for i in range(num_run):
            selected_idxs = random.sample(range(len(labels)), n_blocking)
            lst = list(itertools.product([0, 1], repeat=len(selected_idxs)))
            Arr_Labels = []
            Arr_Score = []
            for Labs in lst:
                tmp = labels[:]
                for idx, val in zip(selected_idxs, Labs):
                    tmp[idx] = val
                Arr_Labels.append(tmp)
            for labs in Arr_Labels:
                Arr_Score.append(self.calculateLogScore(labs))
            # with concurrent.futures.ProcessPoolExecutor(8) as executor:
            #     futures_to_work = [executor.submit(do_work, Assign_label) for Assign_label in Arr_Labels]
            #     concurrent.futures.wait(futures_to_work)
            # for future in concurrent.futures.as_completed(futures_to_work):
            #     try:
            #         val = future.result()
            #         Arr_Score.append(val)
            #         print 'Finish serial: {}'.format(str(val))
            #     except:
            #         print "Error!!!"
            #         print future.result()
            logsum = logsumexp(Arr_Score)
            Arr_prob = [np.exp(s - logsum) for s in Arr_Score]
            Selected_Label = Arr_Labels[np.random.choice(range(len(Arr_prob)), 1, p=Arr_prob)[0]]
            Arr_Scores.append(self.calculateLogScore(Selected_Label))
            if i % 500 == 0:
                print (Selected_Label)
                print (Arr_Scores[-1])
                print ('---------------------')
        return Arr_Scores





#if __name__=='__main__':
    # idx_num = 2
    # Path2GraphStructure = "../TextFolder/Graph_{}.txt".format(str(idx_num))
    # Path2GraphCluster = "../TextFolder/Graph_{}_Cluster.txt".format(str(idx_num))
    # Path2SynthesisFeature = "../TextFolder/Graph_{}_Feature.txt".format(str(idx_num))
    # Path2ViterbiFeature = "../TextFolder/Graph_{}_Feature_Viterbi.txt".format(str(idx_num))
    # Path2ViterbiOut = "../TextFolder/Graph_{}_Feature_Out.txt".format(str(idx_num))
    # from GraphAssistFunc import loadD_Graph
    # G = loadD_Graph(Path2GraphCluster, Path2GraphStructure)
    # n, eds = G.constructLineGraph()
    # LineG = D_LineGraph(G, n, eds)
    #
    # matFeature = np.loadtxt(Path2SynthesisFeature)
    # ver_idx = np.abs(matFeature[:,2])==1
    # ed_idx = np.abs(matFeature[:,2])==2
    # mat_ver = matFeature[ver_idx, :]
    # mat_ed = matFeature[ed_idx, :]
    # from sklearn.preprocessing import QuantileTransformer
    # scaler = QuantileTransformer()
    # scaler.fit(np.matrix(mat_ver[:,3]).transpose())
    # val_ver = scaler.transform(np.matrix(mat_ver[:,3]).transpose())
    # matFeature[ver_idx,3] = val_ver.transpose()
    # scaler = QuantileTransformer()
    # scaler.fit(np.matrix(mat_ed[:,3]).transpose())
    # val_ed = scaler.transform(np.matrix(mat_ed[:,3]).transpose())
    # matFeature[ed_idx,3] = val_ed.transpose()
    #
    # squareMatFeature = np.ones((n,n))*-1
    # for vec in matFeature.tolist():
    #     squareMatFeature[int(vec[0]), int(vec[1])] = vec[3]
    #     squareMatFeature[int(vec[1]), int(vec[0])] = vec[3]
    #
    # LineG.setFeature(squareMatFeature)
    # LineG.WriteFeature(squareMatFeature, Path2ViterbiFeature)
    # LineG.runViterbi("../Script/ViterbiJar/ViterbiAlgorithm.jar", Path2ViterbiFeature, Path2ViterbiOut)
    #
    # from MyIO import ReadViterbiOutFile
    # ViterbiLabel = ReadViterbiOutFile(Path2ViterbiOut)
    # ViterbiLabel = [i if i >0 else -1 for i in ViterbiLabel]
    # print ("Viterbi Label\n", ViterbiLabel)
    # print (LineG.calculateLogScore(ViterbiLabel))
    # TrueLabels = LineG.vs['TrueLabel']
    # print ("True Label\n", TrueLabels)
    # print (LineG.calculateLogScore(TrueLabels))




