
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
outfile = dpath + 'Results/RM_BL_BD/s5_Pred_nonIMG.csv'
my_f_df = pd.read_csv(dpath + 'Results/RM_BL_BD/s4_SFS.csv')
nb_f = get_nb_f(my_f_df, loss_col = 'MSE_mean_cv')
nb_f = 9
my_f_lst = my_f_df.Phenotypes.tolist()[:nb_f]

cov_df = pd.read_csv(dpath + 'Data/Full_data.csv')
target_df = pd.read_csv(dpath + 'Data/WMH_target_rm_bl_bd.csv', usecols = ['eid', 'target_y'])
mydf_train = pd.merge(target_df, cov_df, how = 'inner', on = ['eid'])

X_train, y_train = mydf_train[my_f_lst], mydf_train.target_y
my_lgb = LGBMRegressor(objective='regression', metric='l2', verbosity=0, seed=2023)
my_lgb.set_params(**my_params)
my_lgb.fit(X_train, y_train)

mydf_test = cov_df.loc[cov_df.train_test_id == 'nonIMG']
X_test = mydf_test[my_f_lst]
y_pred = my_lgb.predict(X_test)
pred_df = pd.DataFrame({'eid':mydf_test.eid.tolist(), 'y_pred':y_pred})
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
