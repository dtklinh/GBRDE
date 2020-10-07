import numpy as np

def __error(TrueLabels, PredLabels):
    #D1 = np.equal.outer(TrueLabels, TrueLabels).astype('i')
    #D2 = np.equal.outer(PredLabels, PredLabels).astype('i')
    ''' TrueLabels, PredLabels are squared matrix'''
    num_err = sum(sum(np.logical_xor(TrueLabels,PredLabels).astype('i')))

    L = TrueLabels.shape[0]
    if L == 1:
        return 0
    return float(num_err)/(L*(L-1))

def __Obsolute_Num_Error(TrueLabels_Portion):
    L = len(TrueLabels_Portion)
    NumMajorityIdx = TrueLabels_Portion.count(max(set(TrueLabels_Portion), key=TrueLabels_Portion.count))
    return L - NumMajorityIdx


def eval_error(Membership, Partition_Labels):
    D1 = np.equal.outer(Membership, Membership).astype('i')
    D2 = np.equal.outer(Partition_Labels, Partition_Labels).astype('i')
    Set_Labels = set(Partition_Labels)
    N = len(Membership)
    Error = 0
    for label in Set_Labels:
        Idx_Of_Labels = [idx for idx, val in enumerate(Partition_Labels) if val == label]
        M1 = D1[np.ix_(Idx_Of_Labels, Idx_Of_Labels)]
        M2 = D2[np.ix_(Idx_Of_Labels, Idx_Of_Labels)]
        Error += __error(M1, M2)*len(Idx_Of_Labels)/N
    return Error

def eval_majority_error(Membership, Partition_Labels):
    Error = 0
    N = len(Membership)
    Set_Labels = set(Partition_Labels)
    for label in Set_Labels:
        Idx_Of_Labels = [idx for idx, val in enumerate(Partition_Labels) if val == label]
        Vec = [Membership[idx] for idx in Idx_Of_Labels]
        NumError = __Obsolute_Num_Error(Vec)
        Error += float(NumError)/N
    return Error

def random_Labels(L, num_partitions, seed):
    import random as rnd
    Res = [0]*L
    Lst = [int(i) for i in range(L)]
    #rnd.Random(seed).shuffle(Lst)
    count = 0
    k, m = divmod(L, num_partitions)
    for i in xrange(num_partitions):
        start = i * k + min(i, m)
        end =  (i + 1) * k + min(i + 1, m)
        for j in Lst[start:end]:
            Res[j] = count
        count += 1
    return Res

def error_MemMershipMatrix(TrueLabels, PredLabels):
    '''TrueLabels and PredLabels are arrays which have the same length'''
    assert len(TrueLabels) == len(PredLabels)
    N = len(TrueLabels)
    D1 = np.equal.outer(TrueLabels, TrueLabels).astype('i')
    D2 = np.equal.outer(PredLabels, PredLabels).astype('i')
    err = sum(sum(np.logical_xor(D1,D2).astype("i")))
    return float(err)/(N*(N-1))

def overlap(a, b, return_assignment=False):
    """
    Compare two segmentations using a linear assignment approach.
    PyPi page: https://pypi.python.org/pypi/munkres/
    """
    from munkres import Munkres, print_matrix

    labels_a = list(set(a))
    labels_a.sort()

    labels_b = list(set(b))
    labels_b.sort()

    overlap = np.zeros((len(labels_a),len(labels_b)),'i')

    for i, A in enumerate(labels_a):
        x = (np.array(a) == A).astype('i')
        for j, B in enumerate(labels_b):
            y = (np.array(b) == B).astype('i')
            overlap[i,j] = np.dot(x,y)

    if overlap.shape[0] > 1 and overlap.shape[1] > 1:
        m = Munkres()
        if overlap.shape[0] > overlap.shape[1]:
            i, j = zip(*m.compute(-overlap.T))
            assignments = zip(j,i)
        else:
            assignments = m.compute(-overlap)
    elif overlap.shape[0] == 1 and overlap.shape[1] > 1:
        assignments = [(0,overlap[0].argmax())]
    elif overlap.shape[0] > 1 and overlap.shape[1] == 1:
        assignments = [(overlap[:,0].argmax(), 0)]
    else:
        assignments = [(0,0)]

    x = np.zeros(len(a),'i')
    y = np.zeros(len(b),'i')

    for k, (i, j) in enumerate(assignments):
        x += (np.array(a) == labels_a[i]).astype('i') * k
        y += (np.array(b) == labels_b[j]).astype('i') * k

    i, j = np.transpose(assignments)

    if not return_assignment:
        return overlap[i,j].sum() / float(len(a))
    else:
        return overlap[i,j].sum() / float(len(a)), assignments





