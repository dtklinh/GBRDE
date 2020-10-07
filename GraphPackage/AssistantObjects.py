import numpy as np

class Object_1(object):
    def __init__(self, cut_off_threshold, graph_construction_type, labels, resolutionPara):
        self.Cut_Off_Threshold = cut_off_threshold
        self.Graph_Construction_Type = graph_construction_type
        self.Labels = labels
        self.Resolution_Parameter = resolutionPara

#-------------------------------------------------------------------------
class Object_2(object):
    def __init__(self, serial,CommunityDetection_Type , Arr_Obj1, Membership):
        self.Serial_Number = serial
        self.Community_Detection_Type = CommunityDetection_Type
        #self.Resolution_Parameter = resolutionPara
        self.Array_Object_1 = Arr_Obj1
        self.Membership = Membership
#--------------------------------------------------------

class Feature_Workhouse(object):
    def __init__(self, Mx_3D):
        self.ThreeDim_Matrix = Mx_3D
    def getMeanVar(self, Cluster1, Cluster2):
        assert len(set(Cluster1).intersection(set(Cluster2))) == 0
        Arr = np.array(())
        for i in Cluster1:
            for j in Cluster2:
                tmp = [M[i,j] for M in self.ThreeDim_Matrix]
                Arr = np.append(Arr, np.var(tmp))
        return np.mean(Arr)
    def getMedianVar(self, Cluster1, Cluster2):
        Arr = np.array(())
        for i in Cluster1:
            for j in Cluster2:
                tmp = [M[i, j] for M in self.ThreeDim_Matrix]
                Arr = np.append(Arr, np.var(tmp))
        return np.median(Arr)

class DynDomEntry(object):
    def __init__(self, serial, membership, DistanceMxs, X):
        self.Serial = serial
        self.Membership = membership
        self.DistanceMatrices = DistanceMxs
        self.X = X
    def copy(self):
        import copy
        return DynDomEntry(self.Serial, copy.deepcopy(self.Membership), copy.deepcopy(self.DistanceMatrices), copy.deepcopy(self.X))
    def trim(self, Cluster_To_Be_Deleted):
        Membership = self.Membership
        Membership = np.delete(Membership, Cluster_To_Be_Deleted)
        DistanceMatrices = self.DistanceMatrices
        DistanceMatrices = np.delete(DistanceMatrices, Cluster_To_Be_Deleted, 1)
        DistanceMatrices = np.delete(DistanceMatrices, Cluster_To_Be_Deleted, 2)
        X = self.X
        X = np.delete(X, Cluster_To_Be_Deleted, 1)
        return DynDomEntry(self.Serial, Membership, DistanceMatrices, X)
    def keep(self, Cluster_To_Be_Kept):
        All_idx = range(len(self.Membership))
        Cluster_To_Be_Deleted = [i for i in All_idx if i not in Cluster_To_Be_Kept]
        return self.trim(Cluster_To_Be_Deleted)
    def write_xyz_format(self, Path2File):
        Lines = []
        for frame in self.X:
            Lines.append('{}'.format(str(self.X.shape[1])))
            Lines.append('Remark: XXX')
            for l in frame:
                Lines.append('C\t{}\t{}\t{}'.format(str(l[0]), str(l[1]), str(l[2])))
        with open(Path2File, mode='wt') as output:
            #output.write('\n'.join(Lines))
            for l in Lines:
                output.write("%s\n" % l)



