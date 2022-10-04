import pandas as pd
import numpy as np
import umap

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import typing

from impetuous.special import counts,zvals,geomav
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
    data_filenames  = {    'pima'  : ['pima_analytes_df.tsv'       ,'pima_journal_df.tsv'  ]       ,
                           'water' : ['h20-32_h3o-1.xyz']                                          ,
                           'mnist' : ['mnist_data.tsv'             ,'mnist_target.tsv'     ]       }
    sample_info     = {'water':'.','pima':'GSM','mnist':'pixel'}
    useUmap         = {'water':False, 'pima':True  ,'mnist':True  }
    bRecalc         = True
    crd_desc        = { True:'umap' , False:'real' }
    bReadDistm      = True
    #
    # RESTRICT TO AT MAX N = 50k (32GIG of RAM)
    Nrestrict = 35000 # 50000
    #
    for what_data in [ 'water' , 'pima' , 'mnist' ] :
        bUseUmap = useUmap [ what_data ]
        sample_info_sep = sample_info[ what_data ]
        analyte_df = get_analyte_df ( data_filenames , what_data ,sample_info, Nrestrict )
        #
        N = len ( analyte_df.index.values )
        #
        if bUseUmap :
            axis_labels  = ['UMAP ' + str(i) for i in range(2) ]
            udf = pd.read_csv( results_dir + '.umap'+what_data+'_crd.tsv' , sep='\t' , index_col = 0 )
        else :
            udf			= analyte_df
            axis_labels		= [ c for c in analyte_df.columns.values if sample_info_sep in c ]
            udf.loc[:,'Target']	= ['analyte' for i in range(len(udf)) ]
            udf.loc[:,'Colors']	= ['black' for i in range(len(udf)) ]
        #
        # PRODUCE DISTANCE MATRIX
        distdf	= pd.read_csv( results_dir + '.distm'+what_data+'_'+crd_desc[bUseUmap]+'.tsv',sep='\t',index_col=0 )
        distm	= distdf.values
        print ( distm )
        #
        alf = 0.15
        if True :
            udf .loc[ :,'alpha' ] = alf
        else :
            udf .loc[ :,'alpha' ] = [ alf for i in range(len(udf)) ]
        #
        clusters_df			= pd.read_csv(	results_dir + '.' + crd_desc[bUseUmap] + '' + what_data + '_clusters.tsv' ,
						sep='\t' , index_col=0 )
        clusters_df.columns		= [ 'Distance','Cluster IDs' ]
        clusters_df.loc[:,'M']	= [ len(set(v.split('.')))  for v in clusters_df.loc[:,'Cluster IDs'].values ]
        clusters_df.loc[:,'S,mean']	= [ np.mean( counts( v.split('.') )) for v in clusters_df.loc[:,'Cluster IDs'].values ]
        clusters_df.loc[:,'S,min' ]	= [  np.min( counts( v.split('.') )) for v in clusters_df.loc[:,'Cluster IDs'].values ]
        clusters_df.loc[:,'S,max' ]	= [  np.max( counts( v.split('.') )) for v in clusters_df.loc[:,'Cluster IDs'].values ]
        #
        if True :
            #
            # FIND THE OPTIMAL INFORMATION QUOTA
            # SIZE IS THE AVERAGE _NEW_ CLUSTER SIZE
            # SIZE IS THE AVERAGE CLUSTER SIZE
            # MEAN X IS THE AVERAGE NUMBER OF NEIGHBORS
            #
            print ( clusters_df )
            #
            M = clusters_df.loc[:,'M'].values
            Se = clusters_df.loc[:,'S,mean'].values
            Ge = geomav( [Se , M] )
            #
            Si = clusters_df.loc[:,'S,min'].values
            Gi = geomav( [Si , M] )
            #
            Sa = clusters_df.loc[:,'S,max'].values
            Ga = geomav( [Sa , M] )
            #
            clusters_df.loc[:,'G,mean'] = Ge - N / (N-1)
            clusters_df.loc[:,'G,min']  = Gi - N / (N-1)
            clusters_df.loc[:,'G,max']  = Ga - N / (N-1)
            #
        #
        # HAIRY PLOTTER
        if True :
            pd.options.plotting.backend = 'holoviews'
            #
            import holoviews as hv
            from holoviews import opts
            import hvplot.pandas
            hv.extension('bokeh', logo=False )
            plottered = [ udf.hvplot.scatter( x=axis_labels[0] , y=axis_labels[1], c='Colors' ,
			hover_cols = ['Target'] , alpha='alpha' ,
			height=500, width=500 ) ,
			stat_df  .hvplot.line( x = 'Distance' , y = ['S,mean' ,'S,min', 'S,max','M'] , group_label = ['Type'] ,
							value_label = 'Heuristic Size', height = 500 , width = 500 ) ,
                        stat_df  .hvplot.line( x = 'Distance' , y = ['G,mean','G,min','G,max'] , group_label = ['Type'] ,
                                                        value_label = 'Unimodal Size', height = 500 , width = 500 ) ]
            #
            p1 , p2 = plottered[0] , plottered[1]
            p3 = plottered[2]
            #
            full_plot = p2 + p3
            full_plot .cols(2)
            #
            bokeh_p   = hv.render(full_plot)
            bokeh_p1  = hv.render(p2)
            bokeh_p2  = hv.render(p3)
            from bokeh.io import show, save, output_file
            from bokeh.layouts import column
            show ( bokeh_p )
            output_file('Fig'+what_data.upper()[0]+'_1.html')
            save ( bokeh_p1 )
            output_file('Fig'+what_data.upper()[0]+'_2.html')
            save ( bokeh_p2 )
