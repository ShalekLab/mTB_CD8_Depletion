from scipy.stats import zscore
import scanpy as sc
import anndata
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import scipy.sparse
import seaborn as sns
import heatmap_helper_functions as hh

def grouped_dotplot(x_condition,y_condition,genes,adata,ordered_y_condition=[],within_gene_sep = .6,
                   between_gene_sep = .8,ordered_x_condition=[],pct_expressing_celltype=True):
    tmp_obs = adata.obs

    for G in genes:

        tmp_obs[G] = adata.raw[:,G].X
        tmp_obs[G+" on"] = adata.raw[:,G].X > 0

    means = tmp_obs.groupby([x_condition,y_condition])[genes].mean().stack().reset_index()
    if pct_expressing_celltype:
        pcts = (tmp_obs.groupby([x_condition,y_condition])[[g+" on" for g in genes]].sum()/tmp_obs.groupby([x_condition,y_condition])[[g+" on" for g in genes]].count())#.stack()
        pcts.columns = genes
        
    else: # do percent expressing across entire x condition category
        pcts = (tmp_obs.groupby([x_condition,y_condition])[[g+" on" for g in genes]].sum()/tmp_obs.groupby([x_condition])[[g+" on" for g in genes]].count())#.stack()
        pcts.columns = genes
    means["pcts"]=pcts.stack().reset_index()[0]
    means.columns = [x_condition,y_condition, "gene","mean","pcts"]
    #zscore the means
    means["zmeans"] =means.groupby("gene").transform(lambda x: zscore(x,ddof=1))["mean"]
    means["x_label_name"]= [means.loc[i,"gene"]+means.loc[i,x_condition] for i in means.index]
    x_coords = []#list(range(len(genes)*len(means[x_condition].unique())))
    
    if len(ordered_x_condition)==0:
        ordered_x_condition = means[x_condition].unique()
    x_labelnames = []
    x_coord_value = 0
    linepositions = []
    gene_label_locs = []
    for g in genes:
        x_labelnames += [g+l for l in ordered_x_condition]
        x_coords += [x_coord_value + between_gene_sep,]+ [x_coord_value+between_gene_sep + (l+1)*within_gene_sep for l in range(len(ordered_x_condition)-1)]
        added_space = between_gene_sep+(within_gene_sep*(len(ordered_x_condition)-1))
        gene_label_locs+=[x_coord_value + between_gene_sep+(within_gene_sep*((len(ordered_x_condition)-1)/2.0))]
        x_coord_value+= added_space
        linepositions += [x_coord_value + (between_gene_sep/2.0)]

    x_coord_map = dict(zip(x_labelnames,x_coords))
    means["xcoord"]= means["x_label_name"].map(x_coord_map)
    if len(ordered_y_condition) == 0:
        ordered_y_condition = means[y_condition].unique()
    y_coords = range(len(ordered_y_condition))
    y_coord_map =dict(zip(ordered_y_condition, y_coords))
    means["ycoord"] = means[y_condition].map(y_coord_map)
    figheight=len(y_coords)*.38
    figwidth=len(x_coords)*.4
    plt.figure(figsize=(figwidth,figheight))
    kwargs  =   {'edgecolor':"k", # for edge color
             'linewidth':1, # line width of spot
             'linestyle':'-', # line style of spot
            }
    ax=sns.scatterplot(data=means, x= "xcoord",y="ycoord",hue="zmeans",size="pcts",palette="Blues",sizes=(0, 250),**kwargs)
    ax.set_xticks(x_coords)
    ax.set_xticklabels(list(ordered_x_condition)*len(genes))
    ax.set_yticks(y_coords)
    ax.set_yticklabels(ordered_y_condition)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation = 90)
    ax.set_ylim((ax.get_ylim()[0]-.3,ax.get_ylim()[1]+.3))
    #ax.set_xlim((ax.get_xlim()[0]+(between_gene_sep-within_gene_sep),ax.get_xlim()[1]-(between_gene_sep-within_gene_sep)))
    for i,g in enumerate(genes):
        plt.text(gene_label_locs[i],ax.get_ylim()[1]+.3,g,horizontalalignment='center',multialignment="center")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),title="group expression means, \nzscored by gene")
    for xc in linepositions[:-1]:
        plt.axvline(x=xc, color='grey')






LOCs = {"LOC101925857":"TRAC",
             "LOC102144039":"Mamu-DRA",
                  "LOC102122418":"CYBB*",
                        "LOC102136468":"MafaHLA-DPB1",
                              "LOC102140945":"FCGR3",
                                    "LOC102147203":"GBP1",
       "LOC102142617":"HLA1-B14A",
                "LOC102129434":"RPS11-like",
                "LOC107129205":"TRGV108B",
                "LOC102115168":"TRDC",
                "LOC102143603":"PWWP-pseudo",
                "LOC107130355":"TRAV-HPB-MLT",
                "LOC102115251":"KIR3DL2",
                "LOC107131025":"IGLV1-BL2",
                "LOC102144764":"KIR3DL3",
                "LOC102116611":"KIR2DL3",
                "LOC102139636":"HLA-DB1",
                "LOC102132533":"KLRC2",
                "LOC102140945":"FCGR3A",
                "LOC102132169":"KLRC2-like",
                "LOC102122509":"KLRC4",
                "LOC102115805":"IGHV-MPC11",
                "LOC102133674":"KLRC2-like.2",
                "LOC102116023":"MAP3K8",
                "LOC102134776":"CYB5A",
                "LOC102145001":"SCART1",
                "LOC102124348":"uncharacterized gene",
                "LOC102128672":"TRG2C",
                "LOC102133485":"IFITM3-pseudo",
                "LOC102142127":"uncharacterized gene.2",
                "LOC102145938":"IFITM3-like",
                "LOC102134129":"IFI27-pseudo",
                "LOC102123623":"HIRA",
                "LOC102139778":"uncharacterized gene.3",
                "LOC102136964":"LAIR1",
                "LOC102135811":"BOLA2-pseudo",
                "LOC102143485":"uncharacterized gene.4",
                "LOC102116308":"IGHM-like",
                "LOC102132314":"H2AX",
                "LOC102142866":"FANCD2",
                "LOC102141578":"H3C1-like",
                "LOC102119226":"H1-5",
                "LOC107128728":"uncharacterized gene.5",
                "LOC107126576":"uncharacterized gene.6",
                "LOC102135858":"HGMB2-pseudo",
                "LOC102121131":"H1-3",
                "LOC102123264":"H1-4",
                "LOC102128802":"TUBA1C",
                "LOC102138750":"H1-2",
                "LOC102137365":"H2AC12",
                "LOC102136402":"APITD1",
                "LOC102127357":"APOBEC3D",
                "LOC102146491":"H2AZP1",
                "LOC102134721":"ANP32E-pseudo",
                "LOC107129707":"uncharacterized gene.7",
                "LOC102118091":"H2AC14",
                "LOC102124147":"SPIB",
                "LOC102131515":"LRR1",
                "LOC102144925":"uncharacterized gene.8",
                "LOC102134772":"DEK-pseudo",
                "LOC102142158":"AKR7A2",
                "LOC102119092":"H4",
                "LOC102133926":"HGMB1-like",
                "LOC102134966":"COA4-h",
                "LOC102139881":"MAP3K20",
                "LOC102129658":"uncharacterized gene.9",
                "LOC102128706":"IGHA2-like",
                "LOC102142558":"IGHG1-like",
                "LOC102122028":"RPL40-pseudo",
                "LOC102120885":"RPLP0-pseudo",
                "LOC102121775":"Gogo-B*0103A-like",
                "LOC102131664":"HLA-A-24-like",
                "LOC102137590":"RPLP0-pseudo.2",
                "LOC102136805":"CCL4",
                "LOC102122037":"GOLIM4-like",
                "LOC102142222":"HLA-B-15-like",
                "LOC101867070":"uncharacterized gene.10",
                "LOC102136862":"HLA-DPA1-like",
                "LOC102139771":"RPL14-pseudo",
                "LOC102123006":"ARL6IP1-pseudo",
                "LOC102141835":"HLA-B7-like",
                "LOC102116527":"Gogo-B*0103A-like.2",
                "LOC102140002":"Gogo-B*0103A-like.3",
                "LOC102132418":"EEF1A1-pseudo",
                "LOC102143352":"RPL40-S18-pseudo",
                "LOC102137002":"uncharactarized gene.11",
                "LOC101865988":"uncharactarized gene.12",
                "LOC107128653":"uncharactarized gene.13",
                "LOC102122467":"uncharactarized gene.14",
                "LOC107129085":"uncharactarized gene.15",
                "LOC101925240":"HSPA1A",
                "LOC102119366":"CCL4-like",
                "LOC102125474":"PKM2-pseudo",
                "LOC102116897":"MafaHLA-B7",
                "LOC102128579":"ARIH2-like",
                "LOC102116151":"HLA-B-37-like",
                "LOC102141029":"RPL21",
                "LOC107126798":"uncharactarized gene.16",
                "LOC102143405":"Mafa-B-PATR-B-like",
                "LOC102142031":"ZNF33A-like",
                "LOC102131261":"BIN1",
                "LOC107130579":"uncharactarized gene.17",
                "LOC102133096":"SET-like",
                "LOC102139567":"LRRC37A2",
                "LOC102142617":"HLA1-B14A",
                "LOC102129434":"40S S11",
                "LOC107129205":"TRGV108B",
                "LOC102136846":"HBA1/2/3",
                           "LOC102136192":"HBA1/2/3-like",
                           "LOC102145413":"Mafa-DOB",
                           "LOC102134487":"CYB561A3",
                           "LOC102133296":"IGKV1-39-like",
                           "LOC102140918":"SH3BP5",
                           "LOC102119202":"PPP1R16B",
                           "LOC102126451":"TNFRSF13B",
                           "LOC102133073":"AKAP2",
                           "LOC107130289":"IGHG3-like",
                           "LOC102144414":"Mafa-DRB1-13-like",
                           "LOC107128164":"uncharactarized gene.18",
                           "LOC102126820":"GRAP",
                           "LOC102141176":"Mafa-DQA1-like",
                           "LOC102144867":"RPS11-like",
                           "LOC107128718":"uncharactarized gene.19",
                           "LOC102122749":"NAPSA-like",
                           "LOC102144782":"Mafa-DQB1-like",
                           "LOC102129412":"LIMD1",
                           "LOC102137481":"Mafa-DOA",
                           "LOC107126411":"IGLC6-like",
                           "LOC102144295":"HLA-A-11-like",
                           "LOC102127115":"RPS23-pseudo",
                           "LOC102133073":"AKAP2",
                           "LOC102136305":"uncharactarized gene.20",
                           "LOC102146935":"IFITM3-like",
                           "LOC102130734":"A2M",
                           "LOC102142030":"TM4SF1-like",
                           "LOC102119459":"PPFIBP1",
                           "LOC102140918":"SH2BP5",
                           "LOC102146822":"PTPRJ", "LOC102124330":"CXCL5","LOC102120556":"STAT5A",
                           "LOC102130733":"IFNL3","LOC102138968":"IFNL2","LOC102147203":"GBP1","LOC102146847":"GBP3",
                            "LOC102134512":"CCL23","LOC102134112":"CCL15","LOC102125111":"PF4","LOC102116354":"ZNF675",
                           }
def replace_rankgenesgroups_dotplot_locs(ax_out):
    lbls = []
    for label in ax_out["mainplot_ax"].get_yticklabels():
        if label.get_text() in LOCs:
            label.set_text(LOCs[label.get_text()]+"*")
        lbls.append(label)
    ax_out["mainplot_ax"].set_yticklabels(lbls,fontdict={'fontsize':10})
    lbls = []
    for label in ax_out["mainplot_ax"].get_xticklabels():
        if label.get_text() in LOCs:
            label.set_text(LOCs[label.get_text()]+"*")
        lbls.append(label)
    ax_out["mainplot_ax"].set_xticklabels(lbls,fontdict={'fontsize':10})



def save_filtered_rankgenesgroups(adata_no_doublets,directory):
    names = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]["names"])
    logchnges = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]["logfoldchanges"])
    pvals_adj = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pvals_adj'])
    pvals = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pvals'])
    pts = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pts'])
    #pts_rest = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pts_rest'])
    method = adata_no_doublets.uns["rank_genes_groups"]['params']['method']
    for i in names.columns:
        names_include = names[i].notna()
        #print(names_include)
        n_df = pd.DataFrame(index=names[i][names_include], columns=["logfoldchanges","pvals_adj","pvals","pt_expressing"])#,"pt_other_expressing"])
        n_df["logfoldchanges"] = logchnges[i].values[names_include]
        n_df["pvals_adj"] = pvals_adj[i].values[names_include]
        n_df["pvals"] = pvals[i].values[names_include]
        n_df["pt_expressing"]=pts.loc[names[i][names_include],i].values
        #n_df["pt_other_expressing"]=pts_rest.loc[names[i][names_include],i].values
        if "/" in i:
            i=i.replace("/","_")
        
        n_df.to_csv(directory+"filtered_"+method+"_all_cells_"+i+".csv")



def boxplot_sample_proportions(adata, x_value, color_value,hue="treatment",figsize=(10,5), plottype="box",order=None,hue_order=None,edgecolor=False,swap=False):
    tmp = adata.obs.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0)

    m=tmp.divide(tmp.sum(axis=1), axis=0)
    props = []

    i=0
    if hue+"_colors" in adata.uns and not swap:
        color_dict = dict(zip(adata.obs[hue].cat.categories,adata.uns[hue+"_colors"]))
    elif color_value+"_colors" in adata.uns and swap:
        color_dict = dict(zip(adata.obs[color_value].cat.categories,adata.uns[color_value+"_colors"]))
    else:
        color_dict=None
    for sample in m.index:

        for celltype in m.columns:
            vals = [sample,m.loc[sample,celltype],celltype,adata.obs.loc[adata.obs[x_value]==sample,hue].unique()[0]]
            props.append(vals)
            i+=1
    props_df = pd.DataFrame(props,columns=[x_value,x_value+"_proportion",color_value,hue])
    #sns.boxplot(x="celltype", y="sample_proportion", hue="treatment", data=tips)
    props_df[hue]=props_df[hue].astype("category")
    plt.figure(figsize=figsize)
    if swap:
        old_hue = hue
        hue=color_value
        color_value=old_hue
        old_hue_order=hue_order
        hue_order=order
        order=old_hue_order
    if plottype=="box":

        p=sns.boxplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df, palette=color_dict,hue_order=hue_order,linewidth=3)
        
        if edgecolor==True:
            for i,box in enumerate(p.artists):
                box.set_edgecolor(box.get_facecolor())
                r,g,b,a = box.get_facecolor()
                box.set_facecolor((r,g,b,.3))
            swarm_palette=color_dict
        else:
            swarm_palette={i:"white" for i in color_dict}
        if hue_order is None:
            hue_order = adata.obs[hue].cat.categories
        sns.swarmplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df,  dodge=True,palette=swarm_palette,edgecolor="white",linewidth=.7,size=3.2, hue_order=hue_order)
        plt.legend(p.artists,hue_order)
    if plottype=="bar":
        p=sns.barplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df,palette=color_dict,hue_order=hue_order)
    p.set_xticklabels(p.get_xticklabels(),rotation=90)

def new_adata_from_raw(adata, cluster_values, cluster_key="leiden"):
    adata_0 = anndata.AnnData(adata[adata.obs[cluster_key].isin(cluster_values)].raw.X)
    adata_0.obs_names = adata[adata.obs[cluster_key].isin(cluster_values)].obs_names
    adata_0.var_names = adata[adata.obs[cluster_key].isin(cluster_values)].raw.var_names
    adata_0.obs = adata[adata.obs[cluster_key].isin(cluster_values)].obs
    adata_0.raw = adata_0
    return adata_0


def split_umap_by_category(adata, groupby, colorby=None,gene=None, nrows=None, ncols=None, figsize_multiplier=4, markersize=.3, color_dict = None,other_color="lightgray",gene_cmap="Reds", dim_reduction="umap",cat_order=None):
    '''
    plot one umap by each category in groupby and color the plot by categories in color_by

    if nrows or ncols is none, set them to 1 row and num unique categories cols

    color_dict will be set to the adata color parameter if it is none
    '''


    if nrows is None or ncols is None:
        nrows = 1
        ncols = len(adata.obs[groupby].unique())
    if colorby is None:
        colorby = groupby
    

    if color_dict is None and gene is None:
        color_dict = dict(zip(adata.obs[colorby].cat.categories,adata.uns[colorby+"_colors"]))
    if gene is not None:
        if scipy.sparse.issparse(adata.raw[:,gene].X):
            genedf = pd.DataFrame(adata.raw[:,gene].X.todense(),index=adata.obs_names)
        else:
            genedf = pd.DataFrame(adata.raw[:,gene].X,index=adata.obs_names)
        gene_min = genedf.min()
        gene_max = genedf.max()
    if gene is None: 
        fig,ax = plt.subplots(nrows,ncols,   figsize=(figsize_multiplier*ncols,figsize_multiplier*nrows))
    else:
        fig,ax = plt.subplots(nrows,ncols,   figsize=(figsize_multiplier*(ncols+.5),figsize_multiplier*nrows))
    if cat_order is not None:
        vals=cat_order
    else:
        vals = sorted(list(adata.obs[groupby].unique()))
    val_to_coords = {}
    val_to_exp={}
    for ind,val in enumerate(vals):
        val_to_coords[val] = adata.obsm["X_"+dim_reduction][adata.obs[groupby]==val]
        if gene is not None:
            val_to_exp[val] = genedf.loc[adata.obs[groupby]==val].values

    for ind,val in enumerate(vals):
        other_vals = list(vals)
        other_vals.remove(val)
        row = int(ind/ncols)
        col = ind%ncols
        if nrows >1:
            this_ax = ax[row,col]
        else:
            this_ax = ax[col]
        
        for v in other_vals:
            this_ax.scatter(val_to_coords[v][:,0],val_to_coords[v][:,1],c=other_color,s=markersize,marker="o",alpha=1)
        
        adata_sub = adata[adata.obs[groupby]==val]
        vals_2 = list(adata_sub.obs[colorby].unique())
        val_to_coords2={}
        for ind2,val_2 in enumerate(vals_2):
            val_to_coords2[val_2] = adata_sub.obsm["X_"+dim_reduction][adata_sub.obs[colorby]==val_2]
            if gene is None:
                this_ax.scatter(val_to_coords2[val_2][:,0],val_to_coords2[val_2][:,1],s=markersize,marker="o",c =color_dict[val_2],alpha=1 )
            else:
                m=this_ax.scatter(val_to_coords2[val_2][:,0],val_to_coords2[val_2][:,1],s=markersize,marker=".",c = val_to_exp[val_2].T.flatten(),cmap=gene_cmap,vmin=gene_min,vmax=gene_max)
                if ind == len(adata.obs[groupby].unique())-1:
                    plt.colorbar(m,ax=ax)
        this_ax.set_title(val)
        this_ax.set_xticks([])
        this_ax.set_yticks([])
    if gene is not None:
        fig.suptitle(gene)
    return fig


def qcplots(gran_adata, groupby="leiden", gs4=None,fig=None, donor_colname = "M.Number",sample_colname="sample",include_stackedbars=True):
    import matplotlib.gridspec as gridspec
    from matplotlib import ticker
    if gs4 is None:
        if include_stackedbars:
            fig=plt.figure(figsize=(7,15))
            gs4 = gridspec.GridSpec(6,1)
        else:
            fig=plt.figure(figsize=(7,11))
            gs4 = gridspec.GridSpec(4,1)
    ax_tc = fig.add_subplot(gs4[0, 0])
    #else:
        #gs4 = ax.get_subplotspec()
        #ax_tc=ax
    sc.pl.violin(gran_adata, "total_counts",groupby=groupby,rotation=90,ax=ax_tc,show=False,stripplot=False)
    ax_tc.set_xlabel("")
    ax_tc.set_xticklabels([])
    ax_tc.set_xticks([])
    ax_tc.set_ylabel("n_UMI")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1))
    ax_tc.yaxis.set_major_formatter(formatter)
    ax_mito = fig.add_subplot(gs4[1, 0])
    sc.pl.violin(gran_adata, "percent_mito",groupby=groupby,rotation=90,ax=ax_mito,show=False, stripplot=False)
    ax_mito.set_xlabel("")
    ax_mito.set_xticklabels([])
    ax_mito.set_xticks([])
    ax_mito.set_ylabel("%mito")
    ax_genes = fig.add_subplot(gs4[2, 0])
    sc.pl.violin(gran_adata, "n_genes_by_counts",groupby=groupby,rotation=90,ax=ax_genes,show=False, stripplot=False)
    ax_genes.set_xlabel("")
    ax_genes.set_xticklabels([])
    ax_genes.set_xticks([])
    ax_genes.set_ylabel("n_genes")
    formatter_g = ticker.ScalarFormatter(useMathText=True)
    formatter_g.set_scientific(True) 
    formatter_g.set_powerlimits((-1,1))
    ax_genes.yaxis.set_major_formatter(formatter_g)
    ax_doublet = fig.add_subplot(gs4[3, 0])
    sc.pl.violin(gran_adata, "doublet_scores",groupby=groupby,rotation=90,ax=ax_doublet,show=False, stripplot=False)
    ax_doublet.set_ylabel("doublet\nscores")
    if include_stackedbars:
        ax_doublet.set_xlabel("")
        ax_doublet.set_xticklabels([])
        ax_doublet.set_xticks([])
        ax_sample = fig.add_subplot(gs4[4, 0])
        hh.normalized_stacked_bar_plot(gran_adata, groupby,sample_colname,ax=ax_sample,legend=False)
        ax_sample.set_xlabel("")
        ax_sample.set_xticklabels([])
        ax_sample.set_xticks([])
        ax_sample.set_ylabel("pct cells")
        ax_monkey = fig.add_subplot(gs4[5, 0])
        hh.normalized_stacked_bar_plot(gran_adata, groupby,donor_colname,ax=ax_monkey,legend=False)
        ax_monkey.set_ylabel("pct cells")
