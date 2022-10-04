import pandas as pd
import numpy  as np

if __name__ == '__main__' :
    #
    fnames = [ '../../Clustering/results/pima_cluster1.txt' , '../results/pima_cluster1.txt' ]
    sets = { fname:set() for fname in fnames }
    for name in fnames :
        with open ( name,'r' ) as input :
            for line in input :
                nline = line.replace('\n','')
                sets[ name ] = sets[name] | set([nline])
    lsv = list( sets.values() )
    print ( lsv[0]&lsv[1] )
    print ( lsv[0]^lsv[1] )
