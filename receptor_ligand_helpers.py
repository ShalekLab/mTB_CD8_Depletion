import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def reduce_receptor_ligand_list(receptor_ligand_list, genes_in_matrix):
    '''
    reduce a list of receptors and ligands (pandas dataframe with column for "Receptor" and column for "Ligand")
    to only include genes in the expression matrix of interest
    '''
    new_list = receptor_ligand_list[receptor_ligand_list["Receptor"].apply(lambda i: i in genes_in_matrix)]
    new_list=new_list[new_list["Ligand"].apply(lambda i:i in genes_in_matrix)]
    new_list.index = new_list["Ligand"]+"-"+new_list["Receptor"]
    return new_list


def calc_pct_expressing(new_adata, gran, celltype_col):
    '''
    Calculates the percent of cells expressing each gene in each celltype in 'celltype_col' for the sample 'gran' 
    
    returns a dataframe where columns are genes and rows are celltypes
    '''
    this_gran_adata = new_adata[new_adata.obs["sample"]==gran]
    this_gran_exp = pd.DataFrame(this_gran_adata.X, columns=this_gran_adata.var_names)
    gr_zero = this_gran_exp >0
    gr_zero["celltype"] = this_gran_adata.obs[celltype_col].values
    pct_exp = gr_zero.groupby("celltype").sum().divide(gr_zero.groupby("celltype").count())
    return pct_exp


def mean_values_in_group(samples_in_group, matrix_key, per_sample_results):
    all_celltype_pairs = set()
    for sample in samples_in_group:
         all_celltype_pairs= all_celltype_pairs.union(set(per_sample_results[sample][matrix_key].columns))
    data_sum = pd.DataFrame(data = np.zeros((per_sample_results[samples_in_group[0]][matrix_key].shape[0],len(all_celltype_pairs))),
                           index = per_sample_results[samples_in_group[0]][matrix_key].index,
                           columns = list(all_celltype_pairs))
    #print(data_sum)
    for sample in samples_in_group:
        #print(per_sample_results[sample][matrix_key])
        shared_celltype_pairs = list(set(per_sample_results[sample][matrix_key].columns).intersection(set(all_celltype_pairs)))
        #print(data_sum.loc[per_sample_results[sample][matrix_key].index,per_sample_results[sample][matrix_key].columns])
        data_sum.loc[per_sample_results[sample][matrix_key].index,shared_celltype_pairs]= data_sum.loc[per_sample_results[sample][matrix_key].index,shared_celltype_pairs] + per_sample_results[sample][matrix_key].fillna(0)[shared_celltype_pairs]
    
    return data_sum/len(samples_in_group)


def mean_values_in_group_filtered(samples_in_group, matrix_key,pval_key, per_sample_results):
    all_celltype_pairs = set()
    for sample in samples_in_group:
         all_celltype_pairs= all_celltype_pairs.union(set(per_sample_results[sample][matrix_key].columns))
    data_sum = pd.DataFrame(data = np.zeros((per_sample_results[samples_in_group[0]][matrix_key].shape[0],len(all_celltype_pairs))),
                           index = per_sample_results[samples_in_group[0]][matrix_key].index,
                           columns = list(all_celltype_pairs))
    #print(data_sum)
    for sample in samples_in_group:
        #print(per_sample_results[sample][matrix_key])
        masked_score = per_sample_results[sample][matrix_key].copy()
        masked_score[per_sample_results[sample][pval_key] > 0.05] = 0.0
        
        shared_celltype_pairs = list(set(masked_score.columns).intersection(set(all_celltype_pairs)))
        #print(data_sum.loc[per_sample_results[sample][matrix_key].index,per_sample_results[sample][matrix_key].columns])
        
        data_sum.loc[masked_score.index,shared_celltype_pairs]= data_sum.loc[masked_score.index,shared_celltype_pairs] + masked_score.fillna(0)[shared_celltype_pairs]
    
    return data_sum/len(samples_in_group)



import seaborn as sns
from scipy.stats import zscore

def grouped_dotplot_senderboxed(x_condition,y_condition,ligands_reduced,adata,ordered_y_condition=[],within_gene_sep = .6,
                   between_gene_sep = .8,ordered_x_condition=[],pct_expressing_celltype=True,pos_lfc=None,neg_lfc=None):
    if type(pos_lfc) == None or type(neg_lfc) == None:
        print("you need to specify the x condition variables for positive or negative logfoldchanges")
    tmp_obs = adata.obs
    genes = list(ligands_reduced["ligand"].unique())
    ligands_reduced["condition"] = pos_lfc
    ligands_reduced.loc[ligands_reduced["lfc"]<0, "condition"] = neg_lfc
    ligand_grouped = ligands_reduced.groupby("ligand")
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
    means["boxed"]=1
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
        conditions_senders = ligand_grouped.get_group(g).groupby("condition")
        for cond in list(conditions_senders["condition"].unique()):
            print(cond)
            sender_celltypes = conditions_senders.get_group(cond[0])["sender"].unique()
            means.loc[(means["gene"]==g)&(means[y_condition].isin(sender_celltypes))&(means[x_condition]==cond[0]),"boxed"] = 0
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
    kwargs_box = {'edgecolor':"k","linewidth":1,"linestyle":"-"}
    ax = sns.scatterplot(data=means, x="xcoord",y="ycoord",size="boxed",sizes=(0, 250), marker="s",color="w",**kwargs_box)
    sns.scatterplot(data=means, x= "xcoord",y="ycoord",hue="zmeans",size="pcts",palette="Blues",sizes=(0, 250),ax=ax,**kwargs)
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

        
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

# code borrowed from https://github.com/icbi-lab/scirpy/blob/master/scirpy/pl/_vdj_usage.py
def _gapped_ribbons(
    data,
    ax,
    xstart=1.2,
    gapfreq = 1.0,
    gapwidth = 0.4,
    ribcol = None,
    fun = lambda x: x[3] + (x[4] / (1 + np.exp(-((x[5] / x[2]) * (x[0] - x[1]))))),
    figsize = (3.44, 2.58),
    dpi=300):
    """Draws ribbons using `fill_between`
    Called by VDJ usage plot to connect bars.
    Parameters
    ----------
    data
        Breakpoints defining the ribbon as a 2D matrix. Each row is an x position,
        columns are the lower and upper extent of the ribbon at that position.
    ax
        Custom axis, almost always called with an axis supplied.
    xstart
        The midpoint of the first bar.
    gapfreq
        Frequency of bars. Normally a bar would be drawn at every integer x position,
        hence default is 1.
    gapwidth
        At every bar position, there will be a gap. The width of this gap is identical
         to bar widths, but could also be set to 0 if we need continous ribbons.
    ribcol
        Face color of the ribbon.
    fun
        A function defining the curve of each ribbon segment from breakpoint to
        breakpoint. By default, it is a sigmoid with 6 parameters:
         * range between x position of bars,
         * curve start on the x axis,
         * slope of curve,
         * curve start y position,
         * curve end y position,
         * compression factor of the sigmoid curve
    figsize
        Size of the resulting figure in inches.
    figresolution
        Resolution of the figure in dpi.
    Returns
    -------
    Axes object.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        if isinstance(ax, list):
            ax = ax[0]

    spread = 10
    xw = gapfreq - gapwidth
    slope = xw * 0.8
    x, y1, y2 = [], [], []
    for i in range(1, len(data)):
        xmin = xstart + (i - 1) * gapfreq
        tx = np.linspace(xmin, xmin + xw, 100)
        xshift = xmin + xw / 2
        p1, p2 = data[i - 1]
        p3, p4 = data[i]
        ty1 = fun((tx, xshift, slope, p1, p3 - p1, spread))
        ty2 = fun((tx, xshift, slope, p2, p4 - p2, spread))
        x += tx.tolist()
        y1 += ty1.tolist()
        y2 += ty2.tolist()
        x += np.linspace(xmin + xw, xstart + i * gapfreq, 10).tolist()
        y1 += np.zeros(10).tolist()
        y2 += np.zeros(10).tolist()
    if ribcol is None:
        ax.fill_between(x, y1, y2, alpha=0.6)
    else:
        ax.fill_between(x, y1, y2, color=ribcol, alpha=0.6)

    return ax

def calc_ligand_celltype_barheights(celltype_order, cellgrouped, ligand_order,second_celltype_order,second_celltype_category="sender", spacer=0.007):
    y_startstops_cells = []
    y_startstops_ligands = []
    y_val_cells = 0
    y_val_ligand = 0
    y_val_ribbon = 0
    ribbon_startstops = []
    for sender in celltype_order:
        thisgroup = cellgrouped.get_group(sender)

        height = thisgroup["magnitude"].sum()
        midpoint = y_val_cells + (height/2.0)
        y_startstops_cells.append((y_val_cells, height, sender, midpoint))
        y_val_cells += height
        for ligand in ligand_order:
            if ligand in thisgroup["ligand"].values:
                this_lig = thisgroup.loc[thisgroup["ligand"]==ligand]
                
                height = this_lig["magnitude"].sum()
                midpoint = y_val_ligand + (height/2.0)
                y_startstops_ligands.append((y_val_ligand, height, ligand,midpoint ))
                y_val_ligand += height
                for cell2 in second_celltype_order:
                    if cell2 in this_lig[second_celltype_category].values:
                        height = this_lig.loc[this_lig[second_celltype_category]==cell2,"magnitude"].sum()
                        ribbon_startstops.append([sender,cell2,ligand,y_val_ribbon, y_val_ribbon+height])
                        y_val_ribbon += height
        y_val_cells += spacer
        y_val_ligand += spacer
        y_val_ribbon += spacer
    return y_startstops_cells,y_startstops_ligands,ribbon_startstops,y_val_ligand

def plot_rec_ligand_ribbon(relevant_intr,cell_colors, sender_celltype_order = [],ligand_order = [], receiver_celltype_order=[],spacer=0.007,
                          widths = 1.5,ligand_width = 1.5,between_cell_ligand_spacing=0.15,between_sr_spacing=5,vmin = -6,vmax = 6,cmap=None):
    '''
    relevant_intr is a dataframe with columns sender, ligand, and receiver with the logfoldchange of the interaction in column "0"
    
    for sender, lignad, and receiver orders, if a value in the list does not exist in the corresponding relevant_intr column, it will be skipped,
    and if a value is not in the list but is in the corresponding relevant_intr column, it will also be skipped
    
    spacer dictates the vertical space between celltypes
    
    widths determines the horizontal space of each column
    
    between_cell_ligand_spacing dictates the amount of space between each cell and ligand column
    
    between_sr_spacing is the space between the sender ligands and the receiver ligands, which will be filled with the ribbons
    '''
    if len(ligand_order) == 0:
        ligand_order = relevant_intr["ligand"].unique()
    elif len(set(relevant_intr["ligand"].unique())-set(ligand_order)) > 0:
        relevant_intr = relevant_intr[relevant_intr["ligand"].isin(ligand_order)]
    relevant_intr["magnitude"] = np.abs(relevant_intr["0"])
    sender_ligand_sizes = relevant_intr.groupby(["sender","ligand","receiver"])["magnitude"].sum()
    sender_ligand_sizes_normed = (sender_ligand_sizes/sender_ligand_sizes.sum()).reset_index()
    #print(sender_ligand_sizes_normed)
    receiver_ligand_sizes = relevant_intr.groupby(["receiver","ligand","sender"])["magnitude"].sum()
    receiver_ligand_sizes_normed = (receiver_ligand_sizes/receiver_ligand_sizes.sum()).reset_index()

    if len(sender_celltype_order) == 0:
        celltype_order = sender_ligand_sizes_normed["sender"].unique()
    else:
        celltype_order = sender_celltype_order
    if len(receiver_celltype_order) == 0:
        receiver_celltype_order = receiver_ligand_sizes_normed["receiver"].unique()
    #ligand_order = relevant_intr["ligand"].unique()
    cellgrouped = sender_ligand_sizes_normed.groupby("sender")
    print(sender_ligand_sizes_normed)
    

    y_startstops_cells_senders,y_startstops_ligands_sender, sender_ribbon_startstops,sender_y_val_ligand = calc_ligand_celltype_barheights(celltype_order, cellgrouped, ligand_order, receiver_celltype_order,second_celltype_category="receiver",spacer=spacer)            
    ribbon_df = pd.DataFrame(sender_ribbon_startstops, columns = ["sender","receiver","ligand","sender_ystart","sender_yend"])
    ribbon_df.index = ribbon_df["sender"]+ribbon_df["receiver"]+ribbon_df["ligand"]
    
    

    
    cellgrouped = receiver_ligand_sizes_normed.groupby("receiver")
    
    
    sender_cell_start_x = 1 
    sender_ligand_start_x = sender_cell_start_x + widths+between_cell_ligand_spacing

    
    receiver_ligand_start_x = sender_ligand_start_x + ligand_width+between_sr_spacing
    receiver_cell_start_x = receiver_ligand_start_x + ligand_width + between_cell_ligand_spacing
    y_startstops_cells_receiver,y_startstops_ligands_receiver,receiver_ribbon_startstops,receiver_y_val_ligand = calc_ligand_celltype_barheights(receiver_celltype_order, cellgrouped, ligand_order,celltype_order,second_celltype_category="sender",spacer= spacer)            
    ribbon_df_r = pd.DataFrame(receiver_ribbon_startstops, columns = ["receiver","sender","ligand","receiver_ystart","receiver_yend"])
    ribbon_df_r.index = ribbon_df_r["sender"]+ribbon_df_r["receiver"]+ribbon_df_r["ligand"]

    ribbons_combined = pd.concat([ribbon_df, ribbon_df_r[["receiver_ystart","receiver_yend"]]],axis=1)

    width=10
    height=15
    dpi=300
    ligand_colors = dict(zip(ligand_order, list(mcolors.XKCD_COLORS.values())[20:20+len(ligand_order)]))
    fig, ax = plt.subplots(figsize=(width , height ))
    for start_y,height,label,midpoint in y_startstops_cells_senders:

        ax.add_patch(
                Rectangle(xy=(sender_cell_start_x, start_y), width=widths,
                          height=height, facecolor=cell_colors[label])
            )
        plt.text(sender_cell_start_x,midpoint, label,va='center',ha="right")
    for start_y,height,label,midpoint in y_startstops_ligands_sender:

        ax.add_patch(
                Rectangle(xy=(sender_ligand_start_x, start_y), width=ligand_width,
                          height=height, facecolor=ligand_colors[label])
            )
        plt.text(sender_ligand_start_x+(widths/2.0),midpoint, label,va='center',ha="center")
    for start_y,height,label,midpoint in y_startstops_cells_receiver:

        ax.add_patch(
                Rectangle(xy=(receiver_cell_start_x, start_y), width=widths,
                          height=height, facecolor=cell_colors[label],)
            )
        plt.text(receiver_cell_start_x+widths,midpoint, label,va='center',ha="left")
    for start_y,height,label,midpoint in y_startstops_ligands_receiver:

        ax.add_patch(
                Rectangle(xy=(receiver_ligand_start_x, start_y), width=ligand_width,
                          height=height, facecolor=ligand_colors[label])
            )
        plt.text(receiver_ligand_start_x+(widths/2.0),midpoint, label,va='center',ha="center")

    

    # when 
    if cmap is None:
        cmap = plt.cm.RdGy_r
        cmap.set_under("k")
        cmap.set_over("r")
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    relevant_intr.index = relevant_intr["sender"]+relevant_intr["receiver"]+relevant_intr["ligand"]
    print(ribbons_combined)
    for i in ribbons_combined.index:
        cval = relevant_intr.loc[i,"0"]
        #print(cval)
        if np.abs(cval) == vmax:
            # out of bounds!
            cval = 1000*cval
        elif cval > vmax-0.1:
            cval = vmax-0.1
        elif cval < vmin+0.1:
            cval = vmin+0.1
        # ^ all of that is so we can use a specific color for the very max and the very min 
        _gapped_ribbons([[ribbons_combined.loc[i,"sender_ystart"],ribbons_combined.loc[i,"sender_yend"]],
             [ribbons_combined.loc[i,"receiver_ystart"],ribbons_combined.loc[i,"receiver_yend"]]],ax,
            xstart=sender_ligand_start_x + (ligand_width),
            gapfreq = between_sr_spacing,
            gapwidth = 0,
            ribcol = cmap(norm(cval)),
            figsize = (width , height ))    
    ax.set_xlim(0, receiver_cell_start_x +widths+1)
    ax.set_ylim(0,np.max([sender_y_val_ligand,receiver_y_val_ligand]))
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap= cmap),alpha=0.6)
    
    
import networkx as nx
def intr_plot_sender_receiver(group_results,median_celltype_props_this_group,cell_colors,width_multiplier=40,node_order=[]):
    sums = group_results.sum()
    sums.name = "score"
    summed_interactions = pd.DataFrame(sums)
    if len(node_order) == 0:
        node_order=median_celltype_props_this_group.index
    
    summed_interactions["sender"] = [i.split("_to_")[0]+"~s" for i in summed_interactions.index]
    summed_interactions["receiver"] = [i.split("_to_")[1]+"~r" for i in summed_interactions.index]
    #ummed_interactions["counts"] = group_results_dict["CD8a"]["filtered_mean_weighted_scores"].count()
    G = nx.MultiDiGraph()
    
    if len(node_order) > 0:
        G.add_nodes_from([i+"~s" for i in node_order],bipartite=0)
        G.add_nodes_from([i+"~r" for i in node_order],bipartite=1)
    for i in summed_interactions.index:
        G.add_weighted_edges_from(ebunch_to_add=[(summed_interactions.loc[i,"sender"],summed_interactions.loc[i,"receiver"],summed_interactions.loc[i,"score"])])
    #l, r = nx.bipartite.sets(G)
    
    sep = (1-.1)/len(node_order)
    
    pos = {}
    pos.update({i+"~s":[-1,.05+m*sep] for m,i in enumerate(node_order)})
    pos.update({i+"~r":[1,.05+m*sep] for m,i in enumerate(node_order)})
    plt.figure(figsize=(15,15))
    #pos = nx.bipartite_layout(G,[i+"~s" for i in node_order])
    print(pos)
    nx.draw_networkx_nodes(G,pos,node_color = [cell_colors[i.split("~")[0]] for i in G.nodes],
                    node_size = [median_celltype_props_this_group[i.split("~")[0]]*7000 for i in G.nodes])
    
    labels = nx.draw_networkx_labels(G,pos )
    group1_edges = []
    group1_widths= []
    group1_colors = []
    group2_edges = []
    group2_widths= []
    group2_colors=[]



    for e in G.edges:
        if e not in group1_edges and e not in group2_edges: #and (e[1],e[0],0) in G.edges or e[2]>0:


            group1_edges.append(e)
            group1_widths.append(e[0].split("~")[0]+"_to_"+e[1].split("~")[0])
            group1_colors.append(cell_colors[e[0].split("~")[0]])
            if (e[1],e[0],0) in G.edges:
                group2_edges.append((e[1],e[0],0))
                group2_widths.append(e[1].split("~")[0]+"_to_"+e[0].split("~")[0])
                group2_colors.append(cell_colors[e[1].split("~")[0]])
            #if e[2] > 0:
            #   group2_edges.append(e)
    nx.draw_networkx_edges(G, pos, group1_edges, width=np.abs(summed_interactions.loc[group1_widths,"score"]*width_multiplier),edge_color=group1_colors,arrowstyle="->")
    nx.draw_networkx_edges(G, pos, group2_edges, width=np.abs(summed_interactions.loc[group2_widths,"score"]*width_multiplier),edge_color=group2_colors,
                          connectionstyle='arc3,rad=0.2',arrowstyle="->")
    plt.margins(x=0.2,y=0.2)
    return node_order


def celltype_diff_intr_plot(fold_changes,width_multiplier=40,node_order=[]):
    summed_interactions = fold_changes.copy()
    summed_interactions["sender"] = [i.split("_to_")[0] for i in summed_interactions.index]
    summed_interactions["receiver"] = [i.split("_to_")[1] for i in summed_interactions.index]
    #ummed_interactions["counts"] = group_results_dict["CD8a"]["filtered_mean_weighted_scores"].count()
    G = nx.MultiDiGraph()
    if len(node_order) > 0:
        G.add_nodes_from(node_order)
    for i in summed_interactions.index:
        G.add_weighted_edges_from(ebunch_to_add=[(summed_interactions.loc[i,"sender"],summed_interactions.loc[i,"receiver"],summed_interactions.loc[i,"0"])])
    
    node_order = G.nodes

    plt.figure(figsize=(15,15))
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G,pos,node_color = [cell_colors[i] for i in G.nodes])
                    #node_size = [median_celltype_props_this_group[i]*7000 for i in G.nodes])
    
    labels = nx.draw_networkx_labels(G,pos )
    group1_edges = []
    group1_widths= []
    group1_colors = []
    group2_edges = []
    group2_widths= []
    group2_colors=[]



    for e in G.edges:
        if e not in group1_edges and e not in group2_edges: #and (e[1],e[0],0) in G.edges or e[2]>0:


            group1_edges.append(e)
            group1_widths.append(e[0]+"_to_"+e[1])
            #group1_colors.append(cell_colors[e[0]])
            if (e[1],e[0],0) in G.edges:
                group2_edges.append((e[1],e[0],0))
                group2_widths.append(e[1]+"_to_"+e[0])
                #group2_colors.append(cell_colors[e[1]])
            #if e[2] > 0:
            #   group2_edges.append(e)
            
    #to match these edge color scales use edge_vmin,edge_vmax
    nx.draw_networkx_edges(G, pos, group1_edges, width=np.abs(summed_interactions.loc[group1_widths,"0"])*width_multiplier,edge_color=summed_interactions.loc[group1_widths,"0"],edge_cmap=plt.cm.bwr,arrowstyle="->",edge_vmin=-8,edge_vmax=8)
    nx.draw_networkx_edges(G, pos, group2_edges, width=np.abs(summed_interactions.loc[group2_widths,"0"])*width_multiplier,edge_color=summed_interactions.loc[group1_widths,"0"],edge_cmap=plt.cm.bwr,
                          connectionstyle='arc3,rad=0.2',arrowstyle="->",edge_vmin=-8,edge_vmax=4)
    plt.margins(x=0.2,y=0.2)
    return node_order




def diff_intr_plot_sender_receiver(fold_changes,width_multiplier=40,node_order=[]):
    summed_interactions = fold_changes.copy()
    
    if len(node_order) == 0:
        node_order=[l.split("_to_")[0] for l in summed_interactions.index]
    
    summed_interactions["sender"] = [i.split("_to_")[0]+"~s" for i in summed_interactions.index]
    summed_interactions["receiver"] = [i.split("_to_")[1]+"~r" for i in summed_interactions.index]
    summed_interactions["absscore"] = np.abs(summed_interactions["0"])
    node_size_mult = 10
    sender_sizes = summed_interactions.groupby("sender")["absscore"].sum()*node_size_mult
    receiver_sizes = summed_interactions.groupby("receiver")["absscore"].sum()*node_size_mult

    
    #ummed_interactions["counts"] = group_results_dict["CD8a"]["filtered_mean_weighted_scores"].count()
    G = nx.MultiDiGraph()
    
    if len(node_order) > 0:
        sender_node_order = [i+"~s" for i in node_order if i+"~s" in sender_sizes.index]
        G.add_nodes_from(sender_node_order,bipartite=0)
        receiver_node_order = [i+"~r" for i in node_order if i+"~r" in receiver_sizes.index]
        G.add_nodes_from(receiver_node_order,bipartite=1)
    for i in summed_interactions.index:
        G.add_weighted_edges_from(ebunch_to_add=[(summed_interactions.loc[i,"sender"],summed_interactions.loc[i,"receiver"],summed_interactions.loc[i,"0"])])
    #l, r = nx.bipartite.sets(G)
    
    #calc_sender_receiver_node_sizes
        sep = (1-.1)/len(node_order)
    
    pos = {}
    sender_counter = 0
    receiver_counter = 0
    for m,i in enumerate(node_order):
        if i+"~s" in sender_sizes.index:
            sender_counter +=  np.sqrt(sender_sizes.loc[i+"~s"])
            pos[i+"~s"] = [-1, sender_counter]
            sender_counter += sep + np.sqrt(sender_sizes.loc[i+"~s"])
        
            
        
        
        
        if i+"~r" in receiver_sizes.index:
            receiver_counter += np.sqrt(receiver_sizes.loc[i+"~r"])
            pos[i+"~r"] = [1, receiver_counter]
            receiver_counter += sep + np.sqrt(receiver_sizes.loc[i+"~r"])
        
    #pos.update({i+"~s":[-1,.05+m*sep] for m,i in enumerate(node_order)})
    #pos.update({i+"~r":[1,.05+m*sep] for m,i in enumerate(node_order)})
    plt.figure(figsize=(15,15))
    #pos = nx.bipartite_layout(G,[i+"~s" for i in node_order])
    print(G.nodes)
    print(len(G.nodes))
    print(len(sender_node_order))
    print(len(receiver_node_order))
    nx.draw_networkx_nodes(G,pos,node_color = [cell_colors[i.split("~")[0]] for i in G.nodes], node_size =[sender_sizes.loc[i] for i in sender_node_order]+ [receiver_sizes.loc[i] for i in receiver_node_order] )
    
    labels = nx.draw_networkx_labels(G,pos )
    group1_edges = []
    group1_widths= []
    group1_colors = []
    group2_edges = []
    group2_widths= []
    group2_colors=[]



    for e in G.edges:
        if e not in group1_edges and e not in group2_edges: #and (e[1],e[0],0) in G.edges or e[2]>0:


            group1_edges.append(e)
            group1_widths.append(e[0].split("~")[0]+"_to_"+e[1].split("~")[0])
            #group1_colors.append(cell_colors[e[0].split("~")[0]])
            #if (e[1],e[0],0) in G.edges:
                #group2_edges.append((e[1],e[0],0))
                #group2_widths.append(e[1].split("~")[0]+"_to_"+e[0].split("~")[0])
                #group2_colors.append(cell_colors[e[1].split("~")[0]])
            #if e[2] > 0:
            #   group2_edges.append(e)
    nx.draw_networkx_edges(G, pos, group1_edges, width=np.abs(summed_interactions.loc[group1_widths,"0"]*width_multiplier),edge_color=summed_interactions.loc[group1_widths,"0"],edge_cmap=plt.cm.bwr,arrowstyle="->",edge_vmin=-8,edge_vmax=8)
    #nx.draw_networkx_edges(G, pos, group2_edges, width=summed_interactions.loc[group2_widths,"0"]*width_multiplier,edge_color=summed_interactions.loc[group1_widths,"0"],edge_cmap=plt.cm.bwr,arrowstyle="->")
    plt.margins(x=0.2,y=0.2)
    return node_order