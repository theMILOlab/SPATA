# Import ------------------------------------------------------------------

import sys
print(sys.version)

import importlib
spam_loader = importlib.find_loader('scvelo')
found = spam_loader is not None
print("scvelo found: ",found)



import scvelo as scv
import os 
import sys
import argparse
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 200
import csv
import pandas as pd
#%config InlineBackend.figure_format = 'retina'



# Get arguments -----------------------------------------------------------

parser = argparse.ArgumentParser()

#Arguments
parser.add_argument('-I', '--h5ad', required=True)
parser.add_argument('-U', '--UMAP', required=True)
parser.add_argument('-F', '--Folder', required=True)


args = parser.parse_args()
input_h5ad = args.h5ad
Dir= args.Folder
UMAP = args.UMAP

#Start Analysis
print("Start Analysis..")

os.chdir(Dir)
wd=os.getcwd()
print("Dir: ", wd)

print(input_h5ad)
adata = scv.read(input_h5ad)


#Run analysis

print("----- Model 3 --------")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


scv.tl.recover_dynamics(adata, max_iter=200, return_model=True)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

#implement UMAP of coordinates
umap_s=pd.read_csv(UMAP, sep=",")
umap_s=umap_s[['UMAP_1','UMAP_2']].to_numpy(dtype='float32')
adata.obsm['X_umap']=umap_s
#
#Plot stream

scv.pl.velocity_embedding_stream(adata, show=False, basis="umap", color="seurat_clusters",legend_loc='right margin',save='Stream_cluster.svg',palette='Dark2')


scv.tl.velocity(adata, mode='dynamical')
scv.tl.latent_time(adata)

#Outputs
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
scv.DataFrame(adata.uns['rank_velocity_genes']['names']).to_csv("rank_velocity_genes.csv")

adata.obs.to_csv("observations_metadata.csv")

scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names').to_csv("rank_dynamical_genes.csv")


print("Done..")





