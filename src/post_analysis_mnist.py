import numpy as np
import pandas as pd
import impetuous.quantification as impq
from collections import Counter

if __name__ == '__main__' :
    #
    print ( 'MNIST' )
    #
    results_dir         = '../results/'
    data_dir            = '../data/'
    afunc		= lambda x:x
    #
    afn , jfn	= 'mnist_data.tsv' , 'mnist_target.tsv'
    N_restrict	= 35000
    #
    journal_df  = pd.read_csv( data_dir + jfn , sep='\t' , index_col=0 ).iloc[:N_restrict,:]
    targets     = journal_df.copy().iloc[:N_restrict]
    clusters_df = pd.read_csv( results_dir + 'optclustsmnist_umap.tsv' , sep='\t', index_col=0)
    optims      = [ v.split(',')[1] for v in clusters_df.columns.values if ',' in v ]
    clusters_df .columns = [ c.split(',')[0] if ',' in c else c for c in clusters_df.columns.values ]
    #
    print ( clusters_df )
    print ( journal_df  )
    #
    #which_solution = 'min'
    #which_solution = 'max'
    which_solution = 'mean'
    #
    clusters_df.loc[:,'index'] = clusters_df.index.values
    C = clusters_df.groupby(which_solution).apply(lambda x:','.join(str(y) for y in x.index.values))
    C .name = 'K'
    print ( C )
    #
    results = []
    for K in C :
        parts   = [ int(i) for i in K.split(',') ]
        content = [ targets.iloc[:,0].values[p] for p in parts ]
        ccnr    = Counter(content)
        L       = len(content)
        print ( content , ccnr , len(ccnr) , len(content) )
        R = []
        f0,f = 0,0
        for item in ccnr.items():
            f1 = item[1]/L
            if f1>f0:
                f0=f1
            R .append( str(item[0])+','+str(int(np.floor(item[1]/L*100))) + '%' )
        res = [ len(ccnr), L , f0*(1-np.exp(-(len(ccnr)-1))) , '|'.join(R) ]
        results.append(res)
        print ( res )
    #
    rdf = pd.DataFrame(results ,columns = [ 'cluster cases' ,'cluster size', 'specificity' , 'composition' ])
    #
    print ( rdf.sort_values('cluster size').iloc[-40:,:] )
    print ( rdf.sort_values('specificity').iloc[-40:,:] )
    N = np.sum( rdf.loc[:,'cluster size'].values )
    print ( N )
    print ( np.sum(rdf.iloc[ [v<0.5 for v in rdf.loc[:,'specificity'].values],: ].loc[ : ,'cluster size']) )
    print ( np.sum( rdf.iloc[ [v<0.5 for v in rdf.loc[:,'specificity'].values],: ] . loc[ : ,'cluster size'] )/N *100 ,'%')
