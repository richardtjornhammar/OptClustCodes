import pandas as pd
import numpy as np
import umap

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import typing

if __name__ == '__main__' :
    #
    print ( "HERE: OPTIMISATION TRIALS" )
    print ( "CLUSTERING ANALYSIS" )
    #
    results_dir           = '../results/'
    data_dir              = '../data/'
    sample_info           = {'water':'.','pima':'GSM','mnist':'pixel'}
    useUmap               = {'water':False, 'pima':True  ,'mnist':True  }
    bRecalc               = False
    bReadDistm            = True
    #
    for what_data in [ 'water','pima','mnist' ] :
        bUseUmap    = useUmap[what_data]
        bReadUmap   = bUseUmap
        crd_desc    = { True:'umap' , False:'real' }
        #
        # PRODUCE DISTANCE MATRIX
        #
        cdf                 = pd.read_csv ( results_dir + 'optclusts'+what_data+'_'+crd_desc[bReadUmap]+'.tsv' , sep='\t' , index_col=0 )
        epsilons            = { c.split(',')[0]:float(c.split(',')[1]) for c in cdf.columns if ',' in c }
        axis_labels         = [ c for c in cdf.columns if not ',' in c ][:2]
        new_labels          = [ c.split(',')[0] if ',' in c else c for c in cdf.columns ]
        cdf.columns         = new_labels
        cdf.loc[:, 'min' ]  = [str(int(i)) for i in cdf.loc[:,'min'].values  ]
        cdf.loc[:,'mean' ]  = [str(int(i)) for i in cdf.loc[:,'mean'].values ]
        cdf.loc[:, 'max' ]  = [str(int(i)) for i in cdf.loc[:,'max'].values  ]
        #
        print ( 'N_min  = ' , len(set(cdf.loc[:,'min' ].values)) )
        print ( 'N_mean = ' , len(set(cdf.loc[:,'mean'].values)) )
        print ( 'N_max  = ' , len(set(cdf.loc[:,'max' ].values)) )
        #
        from collections import Counter
        print ( Counter(cdf.loc[:,'min'].values) )
        print ( Counter(cdf.loc[:,'mean'].values) )
        print ( Counter(cdf.loc[:,'max'].values) )
        #
        if True :
            pd.options.plotting.backend = 'holoviews'
            #
            import holoviews as hv
            from holoviews import opts
            import hvplot.pandas
            hv.extension('bokeh', logo=False )
            plottered = [
                        cdf  .hvplot.scatter( x = axis_labels[0] , y =axis_labels[1] , color='min' ,
                                                height = 500 , width = 500, size = 1.13*30 ).opts(show_legend=False ) ,
                        cdf  .hvplot.scatter( x = axis_labels[0] , y =axis_labels[1] , color='mean' ,
                                                height = 500 , width = 500, size = 1.78*30 ).opts(show_legend=False ) ,
                        cdf  .hvplot.scatter( x = axis_labels[0] , y =axis_labels[1] , color='max' ,
                                                height = 500 , width = 500, size = 1.83*30 ).opts(show_legend=False ) , # legend_position='bottom') ,
            ]
            p1 , p2 = plottered[0] , plottered[1]
            p3 = plottered[2]
            #
            full_plot = p1  + p2 + p3
            full_plot .cols(3)
            #
            bokeh_p   = hv.render(full_plot)
            #
            from bokeh.io import show,save
            from bokeh.layouts import column
            show ( bokeh_p )
            save ( bokeh_p )
