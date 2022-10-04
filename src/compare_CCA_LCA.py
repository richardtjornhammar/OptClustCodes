import pandas as pd
import numpy as np
import umap

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
#matplotlib.use('TkAgg',force=True)
import typing
#
from impetuous.quantification   import invert_dict
from impetuous.clustering	import connectivity, connectedness, sclinkages
from impetuous.clustering	import cluster_connections as LCA_connectivity, CCA_DBSCAN
from impetuous.optimisation	import GRS
from impetuous.special		import geomav, unpack, counts, bagscore, fuzzy_locate, lint2lstr
from impetuous.convert		import read_xyz


def transform ( xi:float, R:list , bExact=True , Sfunction=lambda x:np.mean(x) ) -> np.array :
    #
    # NOTE THAT LENGTH OF INTERIOR EQUALS MAX OF EXTERIOR
    if bExact :
        interior,exterior = connectivity(R[0],xi)
        print ( connectivity ( R[0],xi ) )
        print ( connectedness( R[0],xi ) )
        exit(1)
        retarr = np.array([ Sfunction( interior) , len(interior) ])
    else :
        if len(R)>1 :
            cls = cluster_connections( R[0] , xi , Z=R[1] )
        else :
            cls = cluster_connections( R[0] , xi )
        retarr = np.array([ Sfunction(counts(cls.tolist())) , len(set(cls)) ] )
    return ( retarr )

def settest(listset1,listset2,logic=lambda a,b:len(a^b)==0 ):
    return ( [s1 for s1 in listset1 for s2 in listset2 if logic(s1,s2)] )

if __name__=='__main__':
    #
    print ( "HERE: OPTIMISATION TRIALS: LCA CCA COMPARISON : WATER" )
    #
    results_dir		= '../results/'
    data_dir		= '../data/'
    what_data		= 'water'
    #
    data_filenames = {  'pima'  : ['pima_analytes_df.tsv'       ,'pima_journal_df.tsv'  ]       ,
                        'water' : ['h20-32_h3o-1.xyz']                                          ,
                        'mnist' : ['mnist_data.tsv'             ,'mnist_target.tsv'     ]       }
    #
    # RESTRICT TO AT MAX N = 50k (32GIG of RAM)
    Nrestrict		= 35000
    sample_info_sep	= {'water':'.','pima':'GSM','mnist':'pixel'}[what_data]
    DF			= []
    for fn in data_filenames[ what_data ] :
        if '.csv' in fn or '.tsv' in fn :
            DF.append( pd.read_csv( data_dir + fn , sep={'csv':',','tsv':'\t'}[fn.split('.')[-1] ] )  )
        if '.xyz' in fn :
            XYZ = read_xyz( data_dir+fn )
            df_ = pd.DataFrame( [crd[1] for crd in XYZ] ,columns=['D'+sample_info_sep+str(i) for i in range(3)] )
            df_ .loc[:,'Type'] = [l[0] for l in XYZ]
            DF  .append( df_.copy() )
    #
    # RESTRICT
    analyte_df  = DF[0].loc[:,[c for c in DF[0].columns.values if sample_info_sep in c ] ].iloc[:Nrestrict,:]
    #
    bReadDistm  = True
    bUseUmap    = {'water':False,'pima':False,'mnist':False}[what_data]
    crd_desc    = { True:'umap' , False:'real' }
    #
    udf = analyte_df
    axis_labels  = [ c for c in analyte_df.columns.values if sample_info_sep in c ]
    udf .loc[:,'Target'] = ['analyte' for i in range(len(udf)) ]
    udf .loc[:,'Colors'] = ['black' for i in range(len(udf)) ]
    #
    # RETRIEVE DISTANCE MATRIX
    distdf = pd.read_csv( results_dir + '.distm'+what_data+'_'+crd_desc[bUseUmap]+'.tsv',sep='\t',index_col=0 )
    distm  = distdf.values
    #
    distances    = pd.read_csv( results_dir + '.'+crd_desc[bUseUmap]+''+what_data+'_clusters.tsv'    , sep='\t' , index_col=0 ).iloc[:,0].values
    me , mi , ma = np.mean(np.diff(distances)),np.min(distances),np.max(distances)
    ds = 2
    sample_distances =  np.array([ s*me+mi for s in range(int(np.floor((ma-mi)/me)))[::ds] ])
    #
    print ( sample_distances , crd_desc[bUseUmap] + '' + what_data )
    #
    # sanity checks
    print ( distm )
    d0	= distm.copy()
    d1	= distm.copy()
    d2	= distm.copy()
    DH	= []
    for epsilon in sample_distances :
        #
        # METHOD 1 CCA
        res0 = connectivity  ( d0,epsilon )
        A = [set([]) for i in range(len(res0[0]))]
        for r in res0[1] :
            A[r[0]] = A[r[0]]|set([r[1]])
        res0 = A
        #
        # METHOD 2 CCA
        res3  = connectedness( d1,epsilon )
        res10 = CCA_DBSCAN ( epsilon, udf.loc[:,axis_labels] )
        res1  = list(invert_dict({i:n for i,n in zip(range(len(res10)),res10)}).values())
        #
        # METHOD 3 LCA
        res2  = LCA_connectivity ( d2 , epsilon )
        #
        nmax = np.max( res2 )
        A = [set([]) for i in range(nmax)]
        for i,j in zip( range(len(res2)),res2) :
            A[j-1] = A[j-1]|set([i])
        res2  = A
        #
        print ('\nR0>\t', res0,'\nR1>\t',res1,'\nR2>\t',res2 )
        #
        M0  = len ( res0 )
        R0  = [ len ( r ) for r in res0 ]
        Sm0 = np.mean ( R0 )
        Si0 = np.min( R0 )
        Sa0 = np.max( R0 )
        #
        M1  = len ( res1 )
        R1  = [ len ( r ) for r in res1 ]
        Sm1 = np.mean ( R1 )
        Si1 = np.min( R1 )
        Sa1 = np.max( R1 )
        #
        M2  = len ( res2 )
        R2  = [ len ( r ) for r in res2 ]
        Sm2 = np.mean ( R2 )
        Si2 = np.min( R2 )
        Sa2 = np.max( R2 )
        #
        D = []
        D = [*D, *[ M0,Sm0,geomav([M0,Sm0]),Si0,geomav([M0,Si0]),Sa0,geomav([M0,Sa0]) ] ]
        D = [*D, *[ M1,Sm1,geomav([M1,Sm1]),Si1,geomav([M1,Si1]),Sa1,geomav([M1,Sa1]) ] ]
        D = [*D, *[ M2,Sm2,geomav([M2,Sm2]),Si2,geomav([M2,Si2]),Sa2,geomav([M2,Sa2]) ] ]
        print ( res0,len(res0),'\n',res1,len(res1),'\n', res2,len(res2) )
        print ( [ epsilon, M0 , M2 , len(settest(res0,res2)) ] )
        DH.append( D )
    heur_df = pd.DataFrame( DH,	columns = [	'M,CCA','S,mean,CCA','G,mean,CCA','S,min,CCA','G,min,CCA','S,max,CCA','G,max,CCA',
						'MX','SmX','GmX','SiX','GiX','SaX','GaX',
						'M,LCA','S,mean,LCA','G,mean,LCA','S,min,LCA','G,min,LCA','S,max,LCA','G,max,LCA'],
				index = sample_distances  )
    heur_df.loc[:,'Distance'] = sample_distances

    if True :
        pd.options.plotting.backend = 'holoviews'
        #
        import holoviews as hv
        from holoviews import opts
        import hvplot.pandas
        hv.extension('bokeh', logo=False )
        plottered = [
                        heur_df  .hvplot.line( x = 'Distance' , y = ['S,mean,CCA' ,'S,mean,LCA','M,CCA','M,LCA'] ,
							group_label = ['Type'] ,
                                                        value_label = 'Heuristic Size', height = 500 , width = 500 ) ,
                        heur_df  .hvplot.line( x = 'Distance' , y = ['S,min,CCA','S,min,LCA','S,max,CCA','S,max,LCA'] ,
							group_label = ['Type'] ,
                                                        value_label = 'Hueristic Size', height = 500 , width = 500 ) ,
                        heur_df  .hvplot.line( x = 'Distance' , y = ['G,mean,CCA','G,mean,LCA','G,min,CCA','G,min,LCA','G,max,CCA','G,max,LCA'] ,
							group_label = ['Type'] ,
                                                        value_label = 'Hueristic Size', height = 500 , width = 500 ) ]

        p1,p2 = plottered[0] , plottered[1]
        p3 = plottered[2]

        full_plot = p1 + p2 + p3
        full_plot .cols(3)
        #
        bokeh_p   = hv.render(full_plot)

        from bokeh.io import show,save
        from bokeh.layouts import column
        show ( bokeh_p )
        save ( bokeh_p )
    print ( heur_df.loc[:,[c for c in heur_df.columns if 'G' in c ]] )

