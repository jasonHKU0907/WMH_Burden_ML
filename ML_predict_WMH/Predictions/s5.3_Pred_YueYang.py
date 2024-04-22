
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from lightgbm import LGBMRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
pd.options.mode.chained_assignment = None  # default='warn'

dpath = '/Volumes/JasonWork/Projects/OtherProjects/WLB/WMH/'
outfile = dpath + 'Results/RM_BL_BD/s5_Pred_YueYang.csv'
my_f_lst = ['Age', 'Weight', 'Systolic_blood_pressure', 'Body_mass_index', 'Smoking_status', 'ALP', 'Creatinine', 'HbA1c']
ukb_cov_df = pd.read_csv(dpath + 'Data/WMH_Modifiable_bl_data.csv', usecols = ['eid'] + my_f_lst + ['Pulse_rate'])
ukb_target_df = pd.read_csv(dpath + 'Data/WMH_target_rm_bl_bd.csv', usecols = ['eid', 'target_y', 'Region_code'])
ukb_df = pd.merge(ukb_target_df, ukb_cov_df, how = 'inner', on = ['eid'])

yy_df = pd.read_csv(dpath + 'Data/YueyangData.csv', usecols = ['eid'] + my_f_lst + ['target_y'])
yy_df['HbA1c'] = yy_df['HbA1c']*10.93-23.5
yy_df['Smoking_status'] = yy_df['Smoking_status']*2
yy_df['Pulse_rate'] = ukb_df.Pulse_rate.median()

my_f_lst = my_f_lst + ['Pulse_rate']

my_lgb = LGBMRegressor(objective='regression', metric='l2', verbosity=0, seed=2023)
my_lgb.set_params(**my_params)
my_lgb.fit(ukb_df[my_f_lst], ukb_df.target_y)
y_pred = my_lgb.predict(yy_df[my_f_lst])

pred_df = pd.DataFrame({'eid':yy_df.eid,'y_test':yy_df.target_y, 'y_pred':y_pred})
pred_df[['y_test', 'y_pred']].corr()
pred_df.to_csv(outfile, index = False)

