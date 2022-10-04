# OptClustCodes
Codes needed to reproduce the https://arxiv.org/abs/2208.04720 paper

# To reproduce qualitatively equivalent results
Take these steps
```
cd src
```

# Enter an enviroment
To enter the work environment, issue:
```
nix-shell env39tf.nix
```
To enter the work environment and execute all code, issue:
```
env39tf_and_run_all.nix
```
To enter the work environment, show umaps and conduct the post analysis, issue :
```
env39tf_show_article_clusters_do_post_analysis.nix
```

# If that doesnt work
Just install impetuous-gfa using:
```
pip install impetuous-gfa
```
and manually run the codes in sequence

# Retrieve mnist data
```
python3 get_mnist_data.py
```

# Retrive water system data
Note that the supplied data is post CPMD simulation while the fetched data is pre simulation
```
python3 get_water_data.py
```

# Check the symmetry broken solution
For water this corresponds to building SNN with nearest neighbours=2 and 3
```
pythno3 distm_graph.py
```

# Retrive the obesity data
```
python3 get_pima_data.py
```

# Generate all the intermediate results data
```
python3 generate_umaps_and_distance_matrices_hierarchies.py
```

# Compare CCA and LCA
```
python3 compare_CCA_LCA.py
```

# Conduct the full hierarchy (LCA) solutions
```
python3 LCA_clustering.py
```

# Calculate the information function for water 
```
python3 find_optimal_segmentations.py
```

# Produce annotated UMAPS with optimal clustering solution
```
python3 produce_optclust_annotated_umap.py
```

# Show three UMAP solutions : Water,Pima,Mnist
```
python3 show_umap_solution.py
```

# Conduct the post analysis of the Pima data
```
python3 post_analysis_pima.py
```

# Conduct the compositional analysis of the MNIST data
```
python3 post_analysis_mnist.py
```
