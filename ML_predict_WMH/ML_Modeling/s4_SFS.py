
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'

dpath = '/Volumes/JasonWork/Projects/OtherProjects/WLB/WMH/'
outfile = dpath + 'Results/RM_BL_BD/s4_SFS.csv'
cov_df = pd.read_csv(dpath + 'Data/WMH_Modifiable_bl_data.csv')
prs_df = pd.read_csv(dpath + 'Data/Standard_PRS_data.csv')
target_df = pd.read_csv(dpath + 'Data/WMH_target_rm_bl_bd.csv', usecols = ['eid', 'target_y', 'Region_code'])

mydf = pd.merge(target_df, cov_df, how = 'inner', on = ['eid'])

my_f_df = pd.read_csv(dpath + 'Results/RM_BL_BD/s3_FeaImportance.csv')
my_f_df.sort_values(by = 'TotalCover_cv', ascending=False, inplace = True)
my_f_lst = my_f_df.Phenotypes.tolist()

fold_id_lst = [i for i in range(10)]

y_test_full = np.zeros(shape = [1,1])
for fold_id in fold_id_lst:
    test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
    y_test_full = np.concatenate([y_test_full, np.expand_dims(mydf.iloc[test_idx].target_y, -1)])

y_pred_full_prev = y_test_full
tmp_f, MSE_all_lst, MSE_mean_lst, MSE_sd_lst = [], [], [], []
MAE_all_lst, MAE_mean_lst, MAE_sd_lst = [], [], []


for f in tqdm(my_f_lst):
    tmp_f.append(f)
    my_X = mydf[tmp_f]
    MSE_cv, MAE_cv = [], []
    y_pred_full = np.zeros(shape = [1,1])
    for fold_id in fold_id_lst:
        train_idx = mydf['Region_code'].index[mydf['Region_code'] != fold_id]
        test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
        X_train, X_test = mydf.iloc[train_idx][tmp_f], mydf.iloc[test_idx][tmp_f]
        y_train, y_test = mydf.iloc[train_idx].target_y, mydf.iloc[test_idx].target_y
        my_lgb = LGBMRegressor(objective='regression', metric='l2', verbosity=0, seed=2023)
        my_lgb.set_params(**my_params)
        my_lgb.fit(X_train, y_train)
        y_pred = my_lgb.predict(X_test)
        MSE_cv.append(np.round(mean_squared_error(y_test, y_pred), 3))
        MAE_cv.append(np.round(mean_absolute_error(y_test, y_pred), 3))
        y_pred_full = np.concatenate([y_pred_full, np.expand_dims(y_pred, -1)])
    my_mse, my_mae = mean_squared_error(y_test_full, y_pred_full), mean_absolute_error(y_test_full, y_pred_full)
    MSE_all_lst.append(my_mse)
    MSE_mean_lst.append(np.mean(MSE_cv))
    MSE_sd_lst.append(np.std(MSE_cv))
    MAE_all_lst.append(my_mae)
    MAE_mean_lst.append(np.mean(MAE_cv))
    MAE_sd_lst.append(np.std(MAE_cv))

out_df = pd.DataFrame({'Phenotypes': tmp_f,
                       'MSE_mean_cv': MSE_mean_lst, 'MSE_sd_cv':MSE_sd_lst, 'MSE_all':MSE_all_lst,
                       'MAE_mean_cv': MAE_mean_lst, 'MAE_sd_cv':MAE_sd_lst, 'MAE_all':MAE_all_lst})

out_df.to_csv(outfile, index = False)

print('finished')

