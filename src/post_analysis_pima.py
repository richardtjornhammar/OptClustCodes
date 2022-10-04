import numpy as np
import pandas as pd
import impetuous.quantification as impq
from impetuous.special import unpack,zvals


if __name__=='__main__' :

    print ( 'PIMA' )

    results_dir           = '../results/'
    data_dir              = '../data/'

    afunc = lambda x:x

    afn , jfn = 'pima_analytes_df.tsv' , 'pima_journal_df.tsv'
    analyte_df = pd.read_csv( data_dir + afn , sep='\t' , index_col=0 )
    journal_df = pd.read_csv( data_dir + jfn , sep='\t' , index_col=0 )

    data = zvals(analyte_df)['z']
    targets = journal_df.loc[['C(Type0)','C(Type1)'],data.columns.values].copy()
    targets = targets.rename( index={'C(Type0)':'Class','C(Type1)':'Sex'} )
    targets.loc['Sex']   = [ int(v=='male') for v in targets.loc[ 'Sex' ].values ]
    targets.loc['Class'] = [ int(v=='lean') for v in targets.loc['Class'].values ]
    #
    # RETRIEVE CLUSTERING SOLUTION
    clusters_df = pd.read_csv( results_dir + 'optclustspima_umap.tsv' , sep='\t', index_col=0)
    optims 	= [ v.split(',')[1] for v in clusters_df.columns.values if ',' in v ]
    clusters_df .columns = [ c.split(',')[0] if ',' in c else c for c in clusters_df.columns.values ]
    #
    print ( clusters_df )
    print ( journal_df )
    #
    if True :
        formula = 'Group ~ C(Class) + C(Sex)'
        #
        print ( 'CONDUCTING AN ANOVA TEST FOR EACH ANALYTE' )
        print ( 'USING FORMULA: ', formula )
        #
        all_analytes_df = impq.quantify_analytes( analyte_df = data , journal_df = targets , formula = formula )
        all_analytes_df .to_csv( results_dir + 'pima_ANOVA.tsv' , sep='\t' )
    else :
        all_analytes_df = pd.read_csv(  results_dir + 'pima_ANOVA.tsv' , sep='\t' , index_col=0 )

    which_solution = 'max'

    print ( all_analytes_df )

    clusters_df.loc[:,'index'] = clusters_df.index.values
    C = clusters_df.groupby(which_solution).apply(lambda x:','.join(str(y) for y in x.index.values))
    C .name = 'K'
    print ( C )

    results = []
    for K in C :
        #
        # HERE WE TEST FOR SIGNIFCANT ENRICHMENT OF THE CLUSTERS GIVEN
        # THE PREVIOUSLY DETERMINED SIGNIFCANT ANALYTES
        #
        subset = all_analytes_df.iloc[ [ int(i) for i in K.split(',')] , [-2,-1] ]
        SigAnalytes = all_analytes_df.iloc[ all_analytes_df.loc[:,'C(Class),q'].values<0.05 , : ]
        gs = impq.group_significance( subset ,
                        AllAnalytes = set(all_analytes_df.index.values) ,
			SigAnalytes = set(SigAnalytes.index.values) ,
                        alternative = 'greater' )
        results.append( [ gs[0],len(K.split(',')),','.join(set(subset.index.values)),len(set(subset.index.values)) ] )

    cluster_significance_df =  pd.DataFrame( results , columns=['cluster p value','cluster size','constituents', 'S' ] ).sort_values('cluster p value')
    cluster_significance_df.loc[:,'cluster q value'] = [q[0] for q in impq.qvalues(cluster_significance_df.loc[:,'cluster p value'].values) ]
    cluster_significance_df.to_csv(  results_dir + 'pima_'+which_solution+'_clusters_significance.tsv' , sep='\t'  )
    print ( cluster_significance_df , np.sum( cluster_significance_df.loc[:,'cluster size'].values ) )
    print ( np.sum( cluster_significance_df.loc[:,'cluster q value']<0.05) )
    #
    of=open(results_dir + 'pima_cluster1.txt','w')
    print ( len( cluster_significance_df.loc[:,'constituents'].iloc[0].split(',') ) )
    for name in  unpack([ v.split('///') if '///' in v else v for v in cluster_significance_df.loc[:,'constituents'].iloc[0].split(',') ]):
        print ( name , file=of )
