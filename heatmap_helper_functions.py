import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import seaborn as sns
import scanpy as sc
import anndata
from collections import defaultdict
import random
import math
all_colors = list(sns.xkcd_rgb.values())
from matplotlib import colors
joses_cmap=sns.color_palette(["#A020F0", "#9C1FEB", "#991EE6", "#961EE1", "#931DDC", "#8F1CD7", "#8C1CD2", "#891BCE", "#861AC9", "#821AC4", "#7F19BF", "#7C18BA", "#7918B5", "#7517B0", "#7216AC",
  "#6F16A7", "#6C15A2", "#69159D", "#651498", "#621393", "#5F138F", "#5C128A", "#581185", "#551180", "#52107B", "#4F0F76", "#4B0F71", "#480E6D", "#450D68", "#420D63",
  "#3F0C5E", "#3B0B59", "#380B54", "#350A4F", "#320A4B", "#2E0946", "#2B0841", "#28083C", "#250737", "#210632", "#1E062E", "#1B0529", "#180424", "#15041F", "#11031A",
  "#0E0215", "#0B0210", "#08010C", "#040007", "#010002", "#020200", "#070700", "#0C0C00", "#121200", "#171700", "#1C1C00", "#212100", "#262600", "#2B2B00", "#303000",
  "#363600", "#3B3B00", "#404000", "#454500", "#4A4A00", "#4F4F00", "#555500", "#5A5A00", "#5F5F00", "#646400", "#696900", "#6E6E00", "#737300", "#797900", "#7E7E00",
  "#838300", "#888800", "#8D8D00", "#929200", "#979700", "#9D9D00", "#A2A200", "#A7A700", "#ACAC00", "#B1B100", "#B6B600", "#BCBC00", "#C1C100", "#C6C600", "#CBCB00",
  "#D0D000", "#D5D500", "#DADA00", "#E0E000", "#E5E500", "#EAEA00", "#EFEF00", "#F4F400", "#F9F900", "#FFFF00"])
shalek_colors = colors.ListedColormap(joses_cmap.as_hex())
def get_colors(n_cols):
    out = []
    for i in range(n_cols):
        c= random.choice(all_colors)
        out += [c]
        all_colors.remove(c)
    return out
def reorder_from_labels(labels, index):
    # order based on labels:
    clust_to_sample = defaultdict(list)
    for i,s in enumerate(index):
        clust_to_sample[labels[i]] += [s]
    
    new_order = []
    for clust,samps in clust_to_sample.items():
        new_order += samps
    return new_order

def reorder_from_multiple_labels(labels_df,index,labels_order):
    clust_to_sample = defaultdict(list)
    cur_label = labels_order[0]
    labels_order = labels_order[1:]
    
    for i,s in enumerate(index):
        
        clust_to_sample[labels_df.loc[s,cur_label]] += [s]
    
    new_order = []
    # impose an order on the samples
    clusts = sorted(clust_to_sample.keys())
    for clust in clusts:
        samps = clust_to_sample[clust]
        if len(labels_order) == 0: # base case, just reordering on one label
            new_order += samps
        else:
            new_order += reorder_from_multiple_labels(labels_df, samps,labels_order)
    return new_order


def remove_ribo_and_mito(genes):
    mito_genes = [name for name in genes if name.startswith('MT')]
    ribo_genes = [name for name in genes if name.startswith("RP")]
    for g in mito_genes+ribo_genes:
        genes.remove(g)
def heatmap_from_gene_list(adata, gene_list, cell_ordering_keys=None,  gene_cluster=False, cells="all", legend=True, cell_cluster = False, vmax=4,vmin=-2,figpath="", gene_colors=None,z_score=0):
    '''
    creates a heatmap from a gene list ordered by groups if cell ordering keys are not none 
    '''
    if type(cells) ==str:
        cells = adata.obs_names
    if cell_cluster == False : 
        if type(cell_ordering_keys) != list:
            print("cell ordering keys must be specified as a list if cell cluster is false")
        row_order=reorder_from_multiple_labels(adata.obs.loc[cells], cells, cell_ordering_keys)
    else:
        # if we are not reordering, just use the existing order as an abritrary order to make labels match and clustering happens on the clustermap call
        row_order = adata.obs_names
    if type(cell_ordering_keys) == list:
        row_colors_dict = {}
        if legend:
            leg_colors = {}
        for o in cell_ordering_keys:
            if o+"_colors" in adata.uns.keys():
                lut = dict(zip(adata.obs[o].cat.categories, adata.uns[o+"_colors"]))
            else:
                lut=dict(zip(adata.obs.loc[cells,o].unique(), get_colors(len(adata.obs.loc[cells,o].unique()))))
            if legend and len(adata.obs.loc[cells,o].unique())<=7:
                leg_colors.update({str(o)+"_"+str(i):j for i,j in lut.items()})
            
            row_colors_dict[o]=adata.obs.loc[row_order,o].map(lut)


        row_colors = pd.DataFrame(row_colors_dict)
    
    
        #labels=sklearn.cluster.KMeans(n_clusters=7).fit_predict(curr_mat.T)
        #print(labels.shape)

        #new_order=reorder_from_labels(labels, curr_mat.columns)
    raw_adata = anndata.AnnData(adata.raw.X)
    raw_adata.obs_names = adata.obs_names
    raw_adata.var_names = adata.raw.var_names 
    #sc.pp.normalize_per_cell(raw_adata)
    import scipy.sparse
    if scipy.sparse.issparse(raw_adata.X):
        x_mat = raw_adata.X.todense()
    else:
        x_mat = raw_adata.X
    gene_df = pd.DataFrame(x_mat, columns=adata.raw.var_names, index=adata.obs_names)
    g=sns.clustermap(gene_df.loc[row_order,gene_list].T, method='ward',xticklabels=False, yticklabels=True,z_score=z_score,cmap=shalek_colors, col_colors=row_colors, col_cluster=cell_cluster,row_cluster=gene_cluster, figsize=(20,20),vmax=vmax,vmin=vmin,row_colors=gene_colors)
    if legend:
        for label,color in leg_colors.items():
            g.ax_col_dendrogram.bar(0, 0, color=color, label=label, linewidth=0)

        g.ax_col_dendrogram.legend(loc="center", ncol=4)
    if figpath != "":
        plt.savefig(figpath)
    return g

def heatmap_from_differential_genes(adata, diff_group,cell_ordering_keys, n_genes_per_group=30, remove_RM = True,genes=None, min_gene_count=50,vmax=4,vmin=-2,diffex_method="wilcoxon",figpath="",exclude_null=False, null_val="0", groups=None, legend=True):
    
    # differential expression only works on categorical groups, make this group categorical
    adata.obs[diff_group]=adata.obs[diff_group].astype('category')
    #print(adata.obs[diff_group].isnull().sum() )
    if exclude_null and type(groups) != type(None):
        print("Don't use exclude_null and group selection in the same run!")
        return
    if type(groups) != type(None):
        cells = adata.obs[adata.obs[diff_group].isin(groups)].index
        groups = list(groups)
    elif exclude_null: #and adata.obs[diff_group].isnull().sum() > 0:
        # remove the 'null' category from the group
        if  adata.obs[adata.obs[diff_group]==null_val].shape[0] == 0:
            print("No values have the null value")

        cells = adata.obs[adata.obs[diff_group]!=null_val].index
        groups = list(adata.obs.loc[cells,diff_group].unique())
    else:
        cells = adata.obs.index
        groups = 'all'
    
    print(groups)        
    # now look at the differential genes in these groups
    sc.tl.rank_genes_groups(adata, groupby=diff_group,method=diffex_method,only_positive=True, groups=groups )
    
    # get all the genes of interest
    if genes== None:
        genes = set()
        sc.pl.rank_genes_groups(adata, n_genes=n_genes_per_group, show=False)
        for i in [set(i) for i in list(adata.uns["rank_genes_groups"]["names"])[:n_genes_per_group]]:
            genes = genes.union(i)
        genes=list(genes)
    
    
    if remove_RM:
        remove_ribo_and_mito(genes)

         
    gene_df = pd.DataFrame(adata.raw.X, columns=adata.raw.var_names, index=adata.obs_names)

    gene_df = gene_df.loc[cells]
    
    counts = gene_df[gene_df > 0].count()

    genes = counts[genes].loc[counts[genes]>min_gene_count].index
    
    heatmap_from_gene_list(adata, genes, cell_ordering_keys=cell_ordering_keys,  gene_cluster=True, cells=cells, legend=legend, cell_cluster = False, vmax=vmax,vmin=vmin,figpath=figpath)
        

def ranking_to_csvs(adata,path):
    for i in adata.uns["rank_genes_groups"]["names"].dtype.names:
    
        genes = adata.uns["rank_genes_groups"]["names"][i]
        scores = adata.uns["rank_genes_groups"]["scores"][i]
        pd.DataFrame({"genes":genes,"scores":scores}).to_csv(path+"_"+i+".csv", sep=",")



def pc_heatmap(adata, pc_use, n_cells=500, n_genes=30):
    '''
    This is basically a copy of the seurat pc heatmap function
    it makes a heatmap of cells with high and low weights in that pc as columns and genes with high and low weights as rows

    pc_use can be a number or a list representing which pcs to use - it is 1 indexed
    '''
    

    # check that PCA has been run on adata
    try:
        gene_loadings = adata.varm["PCs"] # rows are genes columns are PCs - gene_loadings[:,0] is PC 1
    except ValueError:
        sc.pl.pca(adata)
        gene_loadings = adata.varm["PCs"]

    cell_loadings = adata.obsm["X_pca"] # rows are cells columns are PCs - cell_loadings[:,0] is PC 1

    # get which pcs to use
    if type(pc_use) != list:
        pc_use= list(pc_use)
    if len(pc_use) >4:
        n_rows = math.ceil(len(pc_use)/4)
        n_col = 4
    else:
        n_rows = 1
        n_col = len(pc_use)
    plt.subplots(n_rows,n_col,figsize=(10*n_rows, 7*n_col))
    
    for plt_loc,pc in enumerate(pc_use):
        # for the PC of interest, get the cells to include and the genes to include

        # top and bottom cells:
        sorted_cells = np.argsort(cell_loadings[:,pc-1])
        cell_indices = list(sorted_cells[:int(n_cells/2)])+ list(sorted_cells[-1*int(n_cells/2):])
        cell_names = adata.obs_names[cell_indices]

        # top and bottom genes
        sorted_genes = np.argsort(gene_loadings[:,pc-1])
        gene_indices = list(sorted_genes[:n_genes])+ list(sorted_genes[-1*n_genes:])
        gene_names = adata.var_names[gene_indices]


        final_df = pd.DataFrame(adata[cell_names][:,gene_names].X.T, index = gene_names)
        plt.subplot(n_rows,n_col,plt_loc+1)
        sns.heatmap(final_df,xticklabels=False,yticklabels=True,cmap=shalek_colors)
        plt.title("PC"+str(pc))
    plt.tight_layout()


def differential_genes_to_file(adata,filename,n_genes=50):
    f = open(filename,"w")
    nclusts = len(adata.uns["rank_genes_groups"]["names"][0])
    for n,i in enumerate([[i[j] for i in adata.uns["rank_genes_groups"]["names"][:n_genes]] for j in range(0,nclusts)]):
        f.write("cluster "+ str(n) + "\n")
        for j in i:
            f.write(j+"\n")
