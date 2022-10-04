import pandas as pd
import numpy  as np
#
# ../../Clustering/results/.umapwater_crd.tsv -> ../../Clustering/results/optclustswater_real.tsv
# ../../Clustering/results/.umappima_crd.tsv  -> ../../Clustering/results/optclustspima_umap.tsv
# ../../Clustering/results/.umapmnist_crd.tsv -> ../../Clustering/results/optclustsmnist_umap.tsv
#
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
    useUmap         = {'water':False, 'pima':True  ,'mnist':True  }
    bRecalc         = True
    crd_desc        = { True:'umap' , False:'real' }
    #
    # RESTRICT TO AT MAX N = 50k ( 32GIG of RAM )
    Nrestrict = 35000
    #
    for what_data in [ 'water' , 'pima' , 'mnist' ] :
        sample_info_sep = sample_info[what_data]
        bUseUmap        = useUmap[what_data]
	# PRODUCE THE RESULTS FRAME
        # READ UMAP
        crds_df = pd.read_csv( results_dir+'.umap'+what_data+'_crd.tsv' , sep='\t' , index_col=0 )
        #
        # ASSIGN OPTIMISATION COORDINATES
        # WE READ THE DISTANC MATRIX. IS NEEDED FOR SOME OF THE OTHER METHODS AS INPUT
        distdf  = pd.read_csv( results_dir + '.distm' + what_data + '_' + crd_desc[bUseUmap] + '.tsv' ,
				sep='\t' , index_col=0 )
        distm   = distdf.values
        print ( distm )
        #
        # READ OPTIMUM
        optims = pd.read_csv ( results_dir + '../results/.optimums_df.tsv', sep='\t', index_col=0 )
        print ( optims )
        #
        if False :
            #
            # CREATE CLUSTER LABELS
            # THESE TWO METHODS PRODUCE IDENTICAL RESULTS
            # THIS ONE IS FASTER BUT NEEDS COORDINATES INSTEAD
            # OF DISTANCE MATRIX AS INPUT
            #
            if bUseUmap :
                axis_labels  = ['UMAP ' + str(i) for i in range(2) ]
                udf = pd.read_csv( results_dir + '.umap'+what_data+'_crd.tsv' , sep='\t' , index_col = 0 )
            else :
                udf = get_analyte_df ( data_filenames , what_data ,sample_info, Nrestrict )
                axis_labels = [ c for c in udf.columns.values if sample_info_sep in c ]
                udf .loc[:,'Target'] = ['analyte' for i in range(len(udf)) ]
                udf .loc[:,'Colors'] = ['black' for i in range(len(udf)) ]

            optimal_clusters = [
                CCA_DBSCAN( optims.loc[what_data,'max' ] , udf.loc[ : , axis_labels ].values  ) ,
                CCA_DBSCAN( optims.loc[what_data,'mean'] , udf.loc[ : , axis_labels ].values  ) ,
                CCA_DBSCAN( optims.loc[what_data,'min' ] , udf.loc[ : , axis_labels ].values  )
            ]
        else :
            # THESE TWO METHODS PRODUCE IDENTICAL RESULTS
            optimal_clusters = [
                [ c[0] for c in connectivity ( distm , optims.loc[what_data, 'max'] )[1] ] ,
                [ c[0] for c in connectivity ( distm , optims.loc[what_data,'mean'] )[1] ] ,
                [ c[0] for c in connectivity ( distm , optims.loc[what_data, 'min'] )[1] ]
            ]
        print ( optimal_clusters )
        #
        cdf = pd.DataFrame(optimal_clusters, index = [	'max,' +str(optims.loc[what_data	,'max'	])	,
							'mean,'+str(optims.loc[what_data	,'mean'	])	,
							'min,' +str(optims.loc[what_data	,'min'	])	]).T
        cdf = pd.concat( [cdf.T,crds_df.T] ).T
        print ( cdf , what_data )
        #
        # THIS PART OF THE NAMING IS STRANGE THE UMAP DESIGNATION MEANS IF UMAP WAS USED FOR
        # ANNOTATIONS. THE COORDINATES IN THE FILE ARE THE ALL UMAP COORDINATES
        #
        cdf .to_csv( results_dir + 'optclusts'+what_data+'_'+crd_desc[bUseUmap]+'.tsv' , sep='\t' )
