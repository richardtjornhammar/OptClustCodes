import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import typing

from impetuous.special import counts,zvals
from impetuous.convert import read_xyz

def nearest_neighbor_graph_matrix ( distm:np.array, nn:int , bCheckSolution=False )->np.array:
    desc__ = """
     SLOW METHOD FOR CONSTRUCTING A NEAREST NEIGHBOR GRAPH
     REPRESENTATION OF A DISTANCE MATRIX GIVEN THAT EACH
     NODE SHOULD ONLY BE CONNECTED TO NN NEIHGBOURS
     A ROW CORRESPONDS TO A SPECIFIC NODES NEIGHBOURS
     NOTE : THIS CREATES A NON-SYMMETRIC DISTANCE MATRIX
    """
    from scipy.spatial.distance import squareform
    if len(np.shape(distm)) == 1 :
        distm = squareform ( distm )
    nn_distm = []
    global_cval = -np.inf
    for row in distm :
        cval = sorted(row)[nn-1]
        if cval > global_cval :
            global_cval = cval
        nrow = np.array([ rval if rval<=cval else np.inf for rval in row ])
        nn_distm.append( nrow )
    if bCheckSolution :
        print (desc__ , '\n' , np.sum(np.array(nn_distm)<np.inf,1)  )
    return ( np.array(nn_distm), global_cval )

if __name__ == '__main__' :

    XYZ = read_xyz( '../data/h20-32_h3o-1.xyz' )
    sample_info_sep = '.'
    print ( XYZ , len(XYZ) )
    df_ = pd.DataFrame( [crd[1] for crd in XYZ] ,columns=['D'+sample_info_sep+str(i) for i in range(3)] )
    df_ .loc[:,'Type'] = [l[0] for l in XYZ]
    print ( df_ )

    from scipy.spatial.distance import pdist,squareform
    #
    print ( pdist(df_.loc[:,[c for c in df_ if '.' in c]]) )
    nn_distm, cutoff = nearest_neighbor_graph_matrix(  pdist(df_.loc[:,[c for c in df_ if '.' in c]] ) ,
				 nn = 2 , bCheckSolution = True )
    print ( cutoff )
    nn_distm, cutoff = nearest_neighbor_graph_matrix(  pdist(df_.loc[:,[c for c in df_ if '.' in c]] ) ,
                                 nn = 3 , bCheckSolution = True )
    print ( cutoff )
    from impetuous.clustering import connectivity
    # NOW YOU CAN EMPLOY THE CONNECTIVITY ALGORITHM TO FIND
    # ALL THE CONNECTED CLUSTERS
    print ( connectivity( nn_distm , cutoff ) , cutoff )
    print ( squareform ( pdist(df_.loc[:,[c for c in df_ if '.' in c]]) ) )
