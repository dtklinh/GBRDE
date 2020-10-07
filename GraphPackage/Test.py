import numpy as np
import cPickle as pickle
import os, sys
import igraph as ig
from FindNonRedundantScript import DynDomEntry
from Graph_Util_Funcs import ConstructGraph, partition_gievenNumPar
from Graph_Config import List_Colors
from AssistantObjects import Object_1, Object_2
from multiprocessing import Pool
import concurrent.futures

def do_work(Entry, cutoff_thres, graph_construction_type):
    #Entry = pickle.load(open(os.path.join(Dir2SelectedDynDomEntry, '{}.pkl'.format(str(int(serial)))), 'rb'))
    #Membership = Entry.Membership
    G = ConstructGraph(Entry.DistanceMatrices, Construct_Type=graph_construction_type, cut_off_threshold=cutoff_thres)
    partitions, resolution_para = partition_gievenNumPar(G, 20)
    Labels = [-1]*len(G.vs)
    count = 0
    for p in partitions:
        for i in p:
            Labels[i] = count
        count += 1
    print ('Finish thres: {}, G_type: {}'.format(str(cutoff_thres), str(graph_construction_type)))
    return Object_1(cutoff_thres, graph_construction_type,Labels, resolution_para)
    #return (serial, cutoff_thres, graph_construction_type, Labels, resolution_para)

if __name__=='__main__':

    Dir2SelectedDynDomEntry = '../MyDataSet/DynDom/SelectedDynDomEntry'
    Path2SelectedSerial = '../MyDataSet/DynDom/SerialList_1.txt'
    Dir2SelectedEntryFinal = '../MyDataSet/DynDom/SelectedEntryFinal'
    L = np.loadtxt(Path2SelectedSerial).astype('i')
    #L = list(L)
    print (L)


    G_Construct_Types = [0, 2]
    #CutOffThresholds = np.arange(7,15, 0.5)
    CutOffThresholds = [7.5, 10.5, 13.5]
    
    for serial in L:
        print ('In serial: {}'.format(str(serial)))
        Entry = pickle.load(open(os.path.join(Dir2SelectedDynDomEntry, '{}.pkl'.format(str(int(serial)))), 'rb'))
        Membership = Entry.Membership
        Arr = []
        with concurrent.futures.ProcessPoolExecutor(6) as executor:
            futures_to_work = [executor.submit(do_work,Entry, thres, G_type) for thres in CutOffThresholds for G_type in G_Construct_Types]
            concurrent.futures.wait(futures_to_work)
        for future in concurrent.futures.as_completed(futures_to_work):
            # entry = futures_to_work[future]
            try:
                val = future.result()
                if val is not None:
                    Arr.append(val)
            except:
                print ("Error")
        Obj2 = Object_2(int(serial), 'CPM', Arr, Membership)
        pickle.dump(Obj2, open(os.path.join(Dir2SelectedEntryFinal, '{}.pkl'.format(str(serial))), 'wb'), pickle.HIGHEST_PROTOCOL)
        print ('Finish serial: {}'.format(str(serial)))

    '''
    Entry = pickle.load(open(os.path.join(Dir2SelectedDynDomEntry, '{}.pkl'.format(str(int(L[33])))), 'rb'))
    o1 = do_work(Entry, 7.5, 0)
    sys.exit()
    Membership = Entry.Membership
    G = ConstructGraph(Entry.DistanceMatrices, Construct_Type=1, cut_off_threshold = 7.5)
    #G.vs['color'] = [List_Colors[Membership[v.index]] for v in G.vs]
    #Community detection
    partitions, thres = partition_gievenNumPar(G, 20)
    print "Num par: {}".format(str(len(partitions)))
    ig.plot(partitions)
    Par = [-1]*len(G.vs)
    count = 0
    for p in partitions:
        for i in p:
            Par[i] = count
        count += 1
    Obj1 = Object_1(10.5, 2, Par, thres)
    Obj2 = Object_2(15, 'CPM', [Obj1], Membership)
    pickle.dump(Obj2, open('../MyDataSet/DynDom/Obj2.txt', 'wb'), pickle.HIGHEST_PROTOCOL)
    tmp = pickle.load(open('../MyDataSet/DynDom/Obj2.txt', 'rb'))
    print tmp.Array_Object_1[0].Labels '''