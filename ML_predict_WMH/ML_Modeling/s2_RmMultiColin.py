

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Utility.Training_Utilities import *
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from collections import defaultdict
pd.options.mode.chained_assignment = None  # default='warn'

def find_cluster_best_f(my_X, f_idx_lst, imp_f, f_id, imp_method):
    f_df = pd.DataFrame({'Phenotypes':  my_X.columns[f_idx_lst]})
    merged_df = pd.merge(f_df, imp_f, how='inner', on=[f_id])
    merged_df.sort_values(by = imp_method, ascending=False)
    return merged_df[f_id][0]

def get_imp_analy(Imp_df, top_prop, imp_method = 'TotalGain_cv'):
    imp_score, iter = 0, 0
    while imp_score < top_prop:
        imp_score += Imp_df[imp_method][iter]
        iter+=1
    return iter+1

dpath = '/Volumes/JasonWork/Projects/OtherProjects/WLB/WMH/'
outimg = dpath + 'Results/RM_BL_BD/s2_RmMultiColin.png'
outfile = dpath + 'Results/RM_BL_BD/s2_RmMultiColin.csv'
cov_df = pd.read_csv(dpath + 'Data/WMH_Modifiable_bl_data.csv')
target_df = pd.read_csv(dpath + 'Data/WMH_target_rm_bl_bd.csv', usecols = ['eid', 'target_y', 'Region_code'])

mydf = pd.merge(target_df, cov_df, how = 'inner', on = ['eid'])

my_f_df = pd.read_csv(dpath + 'Results/RM_BL_BD/s1_FeaImportance.csv')
my_f_df.sort_values(by = 'TotalCover_cv', ascending=False, inplace = True)
top_nb = get_imp_analy(my_f_df, top_prop = 0.9, imp_method = 'TotalCover_cv')
my_f_df = my_f_df.iloc[:top_nb,:]
#my_f_df = my_f_df.iloc[:50,:]
my_f_lst = my_f_df.Phenotypes.tolist()

my_X = mydf[my_f_lst]
y = mydf.target_y
my_label = my_f_df.Phenotypes.tolist()

corr = np.array(my_X.corr(method='spearman'))
#corr = np.array(my_X.corr(method='pearson'))
corr = np.nan_to_num(corr)
corr = (corr + corr.T) / 2
np.fill_diagonal(corr, 1)
distance_matrix = 1 - np.abs(corr)
dist_linkage = hierarchy.ward(squareform(distance_matrix))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 10))
dendro = hierarchy.dendrogram(dist_linkage, labels=my_f_lst, ax=ax2)
ax2.set_xticklabels(dendro["ivl"], rotation=60, fontsize=10, horizontalalignment='right')
ax2.axhline(y = 0.5, color = 'r', linewidth=2, linestyle = '--')
dendro_idx = np.arange(0, len(dendro["ivl"]))
ax1.imshow(corr[dendro["leaves"], :][:, dendro["leaves"]])
ax1.set_xticks(dendro_idx)
ax1.set_yticks(dendro_idx)
ax1.set_xticklabels(dendro["ivl"], rotation=60, fontsize=10, horizontalalignment='right')
ax1.set_yticklabels(dendro["ivl"], fontsize=10)
fig.tight_layout()
plt.show()
plt.savefig(outimg)

cluster_ids = hierarchy.fcluster(dist_linkage, .5, criterion="distance")
cluster_id_to_feature_ids = defaultdict(list)
for idx, cluster_id in enumerate(cluster_ids):
    cluster_id_to_feature_ids[cluster_id].append(idx)

my_f_df['Cluster_ids'] = cluster_ids
selected_f = [find_cluster_best_f(my_X, v, my_f_df, f_id = 'Phenotypes', imp_method='TotalCover_cv') for v in cluster_id_to_feature_ids.values()]
select = ['*' if item in selected_f else '' for item in my_f_df.Phenotypes]
my_f_df['Selected'] = select

myout_df = my_f_df[['Phenotypes', 'ShapValues_cv', 'TotalGain_cv', 'TotalCover_cv', 'Ensemble', 'Cluster_ids', 'Selected']]
myout_df.to_csv(outfile, index = False)

print('Finished')
