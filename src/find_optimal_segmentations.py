import pandas as pd
import numpy as np
import umap

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
#matplotlib.use('TkAgg',force=True)
import typing

from impetuous.optimisation import GRS
from impetuous.clustering import CCA_DBSCAN, connectivity, cluster_connections
from impetuous.special import counts,zvals,geomav
from impetuous.convert import read_xyz

def rename_order(A:list) ->list :
    S = sorted(list(set(A)))
    a0 = None
    R = dict()
    for a in A :
        if a0 is None :
            R[a] = 0
            a0   = 0
            continue
        if not a in R :
            R[a] = a0+1
            a0 = a0+1
    return ( [R[a] for a in A] )

def transform_min( xi:float, R:list , Sfunction=lambda x:np.min(x) ) -> np.array :
        return( transform ( xi=xi, R=R , Sfunction = Sfunction ) )

def transform_max( xi:float, R:list , Sfunction=lambda x:np.max(x) ) -> np.array :
        return( transform ( xi=xi, R=R , Sfunction = Sfunction ) )

def transform_mean( xi:float, R:list , Sfunction=lambda x:np.mean(x) ) -> np.array :
        return( transform ( xi=xi, R=R , Sfunction = Sfunction ) )

def transform ( xi:float, R:list , Sfunction=lambda x:np.mean(x) ) -> np.array :
    n,m = np.shape( R[0] )
    if len(R) == 1 :
        if n == m and np.sum(np.diag(R[0])) < 1E-6 and False :
            print ( 'CONNECTIVITY STEP' )	# USES DISTANC MATRIX
            interior , exterior = connectivity(R[0],xi)
        else :
            print ( 'SKLEARN DBSCAN STEP' )	# USES CORRDINATES
            labels = CCA_DBSCAN ( xi , R[0] )
            interior = counts( labels )
        retarr = np.array([ Sfunction( interior) , len(interior) ])
    else :
        print ( 'LCA STEP' )
        cls = cluster_connections( R[0] , xi , Z=R[1] )
        retarr = np.array([ Sfunction(counts(cls.tolist())) , len(set(cls)) ] )
    return ( retarr )


if __name__ == '__main__' :
    #
    print ( "CREATE UMAP COORDINATES" )
    #
    results_dir           = '../results/'
    data_dir              = '../data/'
    data_filenames = {    'pima'  : ['pima_analytes_df.tsv'       ,'pima_journal_df.tsv'  ]       ,
                          'water' : ['h20-32_h3o-1.xyz']                                          ,
                          'mnist' : ['mnist_data.tsv'             ,'mnist_target.tsv'     ]       }
    bRecalc         = True
    crd_desc        = { True:'umap' , False:'real' }
    #
    # RESTRICT TO AT MAX N = 50k ( 32GIG of RAM )
    Nrestrict	= 35000
    optims	= dict()
    #
    for what_data in [ 'water' , 'pima' , 'mnist' ] :
        sample_info_sep = {'water':'.','pima':'GSM','mnist':'pixel'}[what_data]
        bUseUmap        = {'water':False, 'pima':True  ,'mnist':True  }[what_data]
        DF = []
        for fn in data_filenames[what_data] :
            if '.csv' in fn or '.tsv' in fn :
                DF.append( pd.read_csv( data_dir + fn , sep={'csv':',','tsv':'\t'}[fn.split('.')[-1] ] )  )
            if '.xyz' in fn :
                XYZ = read_xyz( data_dir+fn )
                df_ = pd.DataFrame( [crd[1] for crd in XYZ] ,columns=['D'+sample_info_sep+str(i) for i in range(3)] )
                df_ .loc[:,'Type'] = [l[0] for l in XYZ]
                DF  .append( df_.copy() )
        #
        analyte_df      = DF[0].loc[:,[c for c in DF[0].columns.values if sample_info_sep in c ] ].iloc[:Nrestrict,:]
        #
        # THIS CLUSTERING FILE COTAIN THE KNOWN RESULTS
        clusters_df = pd.read_csv( results_dir + '.'+crd_desc[bUseUmap]+''+what_data+'_clusters.tsv'    , sep='\t' , index_col=0 )
        #
        # WE READ THE DISTANC MATRIX. IS NEEDED FOR SOME OF THE OTHER METHODS AS INPUT
        distdf	= pd.read_csv( results_dir + '.distm'+what_data+'_'+crd_desc[bUseUmap]+'.tsv' , sep='\t', index_col=0 )
        distm	= distdf.values
        sord	= sorted( distm[0] )
        extremes = [ sord[1]*0.1,sord[-1] ]
        #
        # WE READ THE COORDINATES
        if bUseUmap :
            axis_labels  = ['UMAP ' + str(i) for i in range(2) ]
            udf = pd.read_csv( results_dir + '.umap'+what_data+'_crd.tsv' , sep='\t' , index_col = 0 )
        else :
            udf = analyte_df
            axis_labels  = [ c for c in analyte_df.columns.values if sample_info_sep in c ]
        udf .loc[:,'Target'] = [ 'analyte' for i in range(len(udf)) ]
        udf .loc[:,'Colors'] = [ 'black'   for i in range(len(udf)) ]
        #
        # CCA_DBSCAN USES THE COORDINATES DIRECTLY SO WE NEED TO SUPPLY BOTH
        # COORDINATES AND DISTANCE MATRIX EXTREMUMS SO THAT THE OPTIMISATION
        # WILL WORK
        use_vals	= udf.loc[ : , axis_labels ].values
        N		= len ( use_vals )
        for doing_label,transform_func in zip( [ what_data+',min', what_data+',mean', what_data + ',max' ],[transform_min,transform_mean,transform_max]) :
            print ( 'DOING\t:\t', doing_label )
            optims [ doing_label ]  = GRS ( data = use_vals , extremes=extremes, transform=transform_func  , unimodal_function = lambda x:geomav(x)-N/(N-1) )
            print ( 'DONE \t:\t', doing_label )

    optimums_df = pd.DataFrame( np.array(list(optims.values())).reshape(-1,3) ,
			columns = [ k.split(',')[-1] for k in optims.keys() ][:3],
			index = [ k.split(',')[0] for k,i in zip(optims.keys(),range(len(optims.items()))) if i%3==0 ] )
    print ( optimums_df )
    optimums_df.to_csv( results_dir+'.optimums_df.tsv' , sep='\t' )
