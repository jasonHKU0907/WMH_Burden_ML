
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
pd.options.mode.chained_assignment = None  # default='warn'

def get_nb_f(mydf, loss_col):
    my_loss = mydf[loss_col].tolist()
    i = 0
    while((my_loss[i+1]<my_loss[i])|(my_loss[i+2]<my_loss[i])|(my_loss[i+3]<my_loss[i])):
        i+=1
    return i+1

dpath = '/Volumes/JasonWork/Projects/OtherProjects/WLB/WMH/'
outfile = dpath + 'Results/RM_BL_BD/s5_Pred.csv'
my_f_df = pd.read_csv(dpath + 'Results/RM_BL_BD/s4_SFS.csv')
nb_f = get_nb_f(my_f_df, loss_col = 'MSE_mean_cv')
nb_f = 9
my_f_lst = my_f_df.Phenotypes.tolist()[:nb_f]

cov_df = pd.read_csv(dpath + 'Data/WMH_Modifiable_bl_data.csv')
target_df = pd.read_csv(dpath + 'Data/WMH_target_rm_bl_bd.csv', usecols = ['eid', 'target_y', 'Region_code'])
mydf = pd.merge(target_df, cov_df, how = 'inner', on = ['eid'])

fold_id_lst = [i for i in range(10)]

y_test_full = np.zeros(shape = [1,1])
for fold_id in fold_id_lst:
    test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
    y_test_full = np.concatenate([y_test_full, np.expand_dims(mydf.iloc[test_idx].target_y, -1)])

y_pred_full_prev = y_test_full
y_pred_full = np.zeros(shape = [1,1])
eid_lst = np.zeros(shape = [1,1])

for fold_id in fold_id_lst:
    train_idx = mydf['Region_code'].index[mydf['Region_code'] != fold_id]
    test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
    X_train, X_test = mydf.iloc[train_idx][my_f_lst], mydf.iloc[test_idx][my_f_lst]
    y_train, y_test = mydf.iloc[train_idx].target_y, mydf.iloc[test_idx].target_y
    my_lgb = LGBMRegressor(objective='regression', metric='l2', verbosity=0, seed=2023)
    my_lgb.set_params(**my_params)
    my_lgb.fit(X_train, y_train)
    y_pred = my_lgb.predict(X_test)
    y_pred_full = np.concatenate([y_pred_full, np.expand_dims(y_pred, -1)])
    eid_lst = np.concatenate([eid_lst, np.expand_dims(mydf.iloc[test_idx].eid, -1)])

pred_df = pd.DataFrame({'eid':eid_lst[:,0],'y_test':y_test_full[:,0], 'y_pred':y_pred_full[:,0]})
pred_df.to_csv(outfile, index = False)



'''
pred_df[['y_test', 'y_pred']].corr()

corr_matrix = np.corrcoef(pred_df.y_pred, pred_df.y_test)
corr = corr_matrix[0,1]
corr**2

from sklearn.metrics import r2_score
r2_score(pred_df.y_pred, pred_df.y_test)

import matplotlib.pyplot as plt

plt.scatter(pred_df.y_pred, pred_df.y_test)


from sklearn.linear_model import LinearRegression
model = LinearRegression()
X, y = mydf[my_f_lst], mydf.target_y
column_means = X.mean()
X.fillna(column_means, inplace = True)
model.fit(X, y)
model.score(X, y)
target_df = pd.read_csv(dpath + 'Data/WMH_target.csv')
target_df['a'] = np.log(target_df['whm_raw'])
'''


