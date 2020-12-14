# Import ------------------------------------------------------------------

import sys
print(sys.version)

import importlib
spam_loader = importlib.find_loader('scvelo')
found = spam_loader is not None
print("scvelo found: ",found)



import scvelo as scv
import scanpy as sc
import cellrank as cr
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
parser.add_argument('-S', '--nr_states', required=True)


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

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


scv.tl.recover_dynamics(adata, max_iter=200)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.tl.latent_time(adata)


#implement UMAP of coordinates
umap_s=pd.read_csv(UMAP, sep=",")
umap_s=umap_s[['UMAP_1','UMAP_2']].to_numpy(dtype='float32')
adata.obsm['X_umap']=umap_s
#
#Plot stream

scv.pl.velocity_embedding_stream(adata, show=False, basis="umap", color="seurat_clusters",legend_loc='right margin',save='Stream_cluster.svg',palette='Dark2')


#Outputs
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
scv.DataFrame(adata.uns['rank_velocity_genes']['names']).to_csv("rank_velocity_genes.csv")

scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names').to_csv("rank_dynamical_genes.csv")


print("Done.. with sc velo, start with cell rank")

if args.nr_states=="No":
    cr.tl.terminal_states(adata, cluster_key='seurat_clusters', weight_connectivities=0.2)
    print("No defined States")
else:
    cr.tl.terminal_states(adata, cluster_key='seurat_clusters', weight_connectivities=0.2, nr_states=int(args.nr_states))
    print("defined States")

cr.pl.terminal_states(adata, save='Terminal_states.svg')


cr.tl.initial_states(adata, discrete=True, cluster_key='seurat_clusters')
cr.pl.initial_states(adata, discrete=True, save='Initial_states.svg')

cr.tl.lineages(adata)
scv.tl.recover_latent_time(adata, root_key='initial_states_probs', end_key='terminal_states_probs')
scv.tl.paga(adata, groups='seurat_clusters', root_key='initial_states_probs', end_key='terminal_states_probs',
use_time_prior='velocity_pseudotime')

cr.pl.cluster_fates(adata, mode="paga_pie", cluster_key="seurat_clusters", basis='umap', legend_kwargs={'loc': 'top right out'}, legend_loc='top left out',
node_size_scale=5, edge_width_scale=1, max_edge_width=4, title='directed PAGA', save='Directed PAGA.svg')

cr.tl.lineage_drivers(adata).to_csv("lineage_driver.csv")
adata.obs.to_csv("observations_metadata.csv")
pd.DataFrame(adata.obsm["macrostates_bwd"], index=adata.obs.index.values).to_csv("Initial_states.csv")
pd.DataFrame(adata.obsm["macrostates_fwd"], index=adata.obs.index.values).to_csv("Terminal_states.csv")
pd.DataFrame(adata.obsm["to_terminal_states"], index=adata.obs.index.values).to_csv("2_Terminal_states.csv")

print("Done... with sc velo, start with cell rank")
