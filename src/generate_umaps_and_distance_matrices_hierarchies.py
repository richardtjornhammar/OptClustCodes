import pandas as pd
import numpy as np
import umap

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
#matplotlib.use('TkAgg',force=True)
import typing

from impetuous.special import counts,zvals
from impetuous.convert import read_xyz

def get_analyte_df ( data_filenames , what_data ,sample_info = {'water':'.','pima':'GSM','mnist':'pixel'}, Nrestrict=35000 ) :
    from impetuous.convert import read_xyz
    sample_info_sep = sample_info[what_data]
    DF = []
    for fn in data_filenames[what_data] :
        if '.csv' in fn or '.tsv' in fn :
            DF.append( pd.read_csv( data_dir + fn , sep={'csv':',','tsv':'\t'}[fn.split('.')[-1] ] )  )
        if '.xyz' in fn :
            XYZ = read_xyz( data_dir+fn )
            df_ = pd.DataFrame( [crd[1] for crd in XYZ] ,columns=['D'+sample_info_sep+str(i) for i in range(3)] )
            df_ .loc[:,'Type'] = [l[0] for l in XYZ]
            DF  .append( df_.copy() )
    analyte_df  = DF[0].loc[:,[c for c in DF[0].columns.values if sample_info_sep in c ] ].iloc[:Nrestrict,:]
    return ( analyte_df )

if __name__ == '__main__' :
    #
    from impetuous.clustering import connectivity , connectedness , CCA_DBSCAN
    print ( "CREATE UMAP COORDINATES" )
    #
    results_dir           = '../results/'
    data_dir              = '../data/'
    data_filenames = {    'pima'  : ['pima_analytes_df.tsv'       ,'pima_journal_df.tsv'  ]       ,
                          'water' : ['h20-32_h3o-1.xyz']                                          ,
                          'mnist' : ['mnist_data.tsv'             ,'mnist_target.tsv'     ]       }
    sample_info     = {'water':'.','pima':'GSM','mnist':'pixel'}
    useUmap         = {'water':False, 'pima':True ,'mnist':True  }
    bRecalc         = True
    crd_desc        = { True:'umap' , False:'real' }
    lconn           = { 'water':1.0 , 'pima':1.0 , 'mnist' : 1.0  } # CURRENT DEFAULT SETTING
    nneig           = { 'water':15  , 'pima':15  , 'mnist' : 15   } # CURRENT DEFAULT SETTING
    mdist           = { 'water':0.1 , 'pima':0.1 , 'mnist' : 0.1  } # CURRENT DEFAULT SETTING
    #
    # RESTRICT TO AT MAX N = 50k ( 32GIG of RAM )
    Nrestrict = 35000
    #
    for what_data in [ 'water' , 'pima' , 'mnist' ] :
        sample_info_sep = sample_info[ what_data ]
        analyte_df = get_analyte_df ( data_filenames , what_data ,sample_info, Nrestrict )
        #
        # CREATE THE UMAP
        umappad      = umap.UMAP( random_state=42 , local_connectivity = lconn[what_data] ,
					min_dist=mdist[what_data] , n_neighbors=nneig[what_data] )
        if what_data == 'water' :
            umap_in_data = analyte_df.values
        else :
            umap_in_data = zvals(analyte_df)['z']
        #
        umap_crds    = umappad.fit_transform( umap_in_data )
        axis_labels  = ['UMAP ' + str(i) for i in range(2) ]
        udf          = pd.DataFrame( umap_crds ,
				columns = axis_labels ,
				index = analyte_df.index.values )
        udf.loc[:,'Colors'] = [ 'black' for i in range(len(analyte_df)) ]
        udf.loc[:,'Target'] = analyte_df.index.values.tolist()
        udf.to_csv( results_dir + '.umap'+what_data+'_crd.tsv' , sep='\t' )
        #
        # PRODUCE DISTANCE MATRIX
        distm       = pd.DataFrame([[None]])
        distdf      = None
        for bUseUmap in [False,True]:
            from scipy.spatial.distance import squareform,pdist
            if not bUseUmap :
                R = umap_in_data			# PRE UMAP DATA
            else :
                R = udf.loc[ : , axis_labels ].values	# POST UMAP DATA
            distm  = squareform( pdist( R ) )
            distdf = pd.DataFrame(distm,index=udf.index.values,columns=udf.index.values)
            distdf .to_csv( results_dir + '.distm'+what_data+'_'+crd_desc[bUseUmap]+'.tsv' , sep='\t')
            distm  = distdf.values
            #
            # CALCULATE HIERARCHY
            from impetuous.clustering import sclinkages
            linkages = sclinkages ( distm )
            links    = linkages['CL'] #( distm )
            links_df = pd.DataFrame( [ [k,v,len(k.split('.'))] for (k,v) in links.items() ] ,
                      columns = ['indices','distance','size'] )
            links_df .to_csv ( results_dir + '.'+crd_desc[bUseUmap]+''+what_data+'_links.tsv' , sep='\t' )
            clusters_df = pd.DataFrame ( [  [d,'.'.join([str(v_) for v_ in v]) ] for (d,v) in linkages['F'].items() ] )
            clusters_df.to_csv( results_dir + '.'+crd_desc[bUseUmap]+''+what_data+'_clusters.tsv' , sep='\t' )
            print ( clusters_df )

