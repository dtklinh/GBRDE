from mainPackage.PathAndDir import Dir2LineGraph
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
#import seaborn as sns
import os
from sklearn.preprocessing import MinMaxScaler
from Outlier_Detection.Config import Vertex_outlier_Thres, Edge_outlier_Thres

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh
def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def percentile_based_outlier(data, threshold=95):
    diff = (100 - threshold) / 2.0
    minval, maxval = np.percentile(data, [diff, 100 - diff])
    return (data < minval) | (data > maxval)

if __name__=='__main__':
    L = [2661]
    for serial in L:
        LineG = ig.Graph().Read_Picklez(os.path.join(Dir2LineGraph, '{}.pkl'.format(str(serial))))
        G = LineG['OriginalGraph']
        SquareMatFeature = LineG['SquareMatFeature']

        pred = [SquareMatFeature[0,e["TwoEndNodes_ID"][0], e["TwoEndNodes_ID"][1]] for e in LineG.es if e["TwoEndNodeConnected"]==False]
        #pred = [SquareMatFeature[0,e.source, e.target] for e in G.es]

        pred = np.array(pred)
        #pred = np.matrix(pred).transpose()
        #pred = MinMaxScaler().fit_transform(pred)


        #print pred
        is_out = is_outlier(pred, Vertex_outlier_Thres)

        #outliers = pred[is_out]
        #normal = pred[np.logical_not(is_out)]
        #thres = (np.max(normal) + np.min(outliers))/2
        #print "Thres: {}".format(str(thres))
        print (np.min(pred[is_out]))
        print (pred[is_out])
        print (len(pred[is_out]))
        #pred = pred.tolist()

        y = [1 if G.vs[e["TwoEndNodes_ID"][0]]['TrueLabel']==G.vs[e["TwoEndNodes_ID"][1]]['TrueLabel'] else -1 for e in LineG.es if e["TwoEndNodeConnected"]==False]
        #y = [1 if G.vs[e.source]['TrueLabel']==G.vs[e.target]['TrueLabel'] else -1 for e in G.es]
        print(len(y))

        mm = np.max(pred) #

        print (pred[is_out])
        print (np.zeros_like(pred[is_out]))

        #sns.distplot(pred, rug=True, hist=False)
        plt.plot(pred[is_out], np.zeros_like(pred[is_out]), 'ro', clip_on=False)

        #pred = pred.transpose().tolist()[0]
        #print len(pred)
        plt.hist([val for idx, val in enumerate(pred) if y[idx]==1],bins=np.arange(0,mm+0.1,mm/100), color='red', alpha=0.55)
        plt.hist([val for idx, val in enumerate(pred) if y[idx] == -1],bins=np.arange(0,mm+0.1,mm/100), color='blue', alpha=0.55)
        plt.show()

