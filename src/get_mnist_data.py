import pandas as pd

if __name__=='__main__':
    from sklearn.datasets import fetch_openml
    mnist = fetch_openml("mnist_784", version=1)
    pd.DataFrame( mnist.target ).to_csv('../data/mnist_target.tsv'  ,sep='\t' )
    pd.DataFrame( mnist.data   ).to_csv('../data/mnist_data.tsv'    ,sep='\t' )

