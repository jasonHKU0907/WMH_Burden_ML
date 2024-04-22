
import numpy as np
import pandas as pd
from tqdm import tqdm
from Utility.Training_Utilities import *
from sklearn.metrics import mean_squared_error, mean_absolute_error
pd.options.mode.chained_assignment = None  # default='warn'

def get_avg_output(mydf, nb_iters):
    idx_lst = [ele for ele in range(len(mydf))]
    mae_lst, mse_lst, cor_lst = [], [], []
    for i in tqdm(range(nb_iters)):
        random.seed(i)
        bt_idx = [random.choice(idx_lst) for _ in range(len(idx_lst))]
        mydf_bt = mydf.copy()
        mydf_bt = mydf_bt.iloc[bt_idx, :]
        mydf_bt[['y_test', 'y_pred']].corr()
        mae_lst.append(mean_absolute_error(mydf_bt.y_test, mydf_bt.y_pred))
        mse_lst.append(mean_squared_error(mydf_bt.y_test, mydf_bt.y_pred))
        cor_lst.append(np.corrcoef(mydf_bt.y_test, mydf_bt.y_pred)[0, 1])
    results_df = pd.DataFrame([mae_lst, mse_lst, cor_lst]).T
    results_df.columns = ['MAE', 'MSE', 'COR']
    mae_out = [np.round(results_df.MAE.median(),5),
               np.round(results_df.MAE.quantile(0.025),5),
               np.round(results_df.MAE.quantile(0.975),5),
               '{:.3f}'.format(results_df.MAE.median()) + ' [' +
               '{:.3f}'.format(results_df.MAE.quantile(0.025)) + ' - ' +
               '{:.3f}'.format(results_df.MAE.quantile(0.975)) + ']']
    mse_out = [np.round(results_df.MSE.median(), 5),
               np.round(results_df.MSE.quantile(0.025), 5),
               np.round(results_df.MSE.quantile(0.975), 5),
               '{:.3f}'.format(results_df.MSE.median()) + ' [' +
               '{:.3f}'.format(results_df.MSE.quantile(0.025)) + ' - ' +
               '{:.3f}'.format(results_df.MSE.quantile(0.975)) + ']']
    cor_out = [np.round(results_df.COR.median(), 5),
               np.round(results_df.COR.quantile(0.025), 5),
               np.round(results_df.COR.quantile(0.975), 5),
               '{:.3f}'.format(results_df.COR.median()) + ' [' +
               '{:.3f}'.format(results_df.COR.quantile(0.025)) + ' - ' +
               '{:.3f}'.format(results_df.COR.quantile(0.975)) + ']']
    myout_df = pd.DataFrame([mae_out, mse_out, cor_out]).T
    myout_df.columns = ['MAE', 'MSE', 'COR']
    myout_df.index = ['Median', 'LBD', 'UBD', 'OUT']
    return (results_df, myout_df)

dpath = '/Volumes/JasonWork/Projects/OtherProjects/WLB/WMH/'
outfile1 = dpath + 'Results/RM_BL_BD/Evaluations/RegEval/s6_Eval_ukb_24485.csv'
outfile2 = dpath + 'Results/RM_BL_BD/Evaluations/RegEval/s6_Eval_ukb_24485_summary.csv'
mydf = pd.read_csv(dpath + 'Results/RM_BL_BD/s5_Pred_x24485_2_0.csv')

results_df, myout_df = get_avg_output(mydf, nb_iters=1000)

results_df.to_csv(outfile1, index = False)
myout_df.to_csv(outfile2)



outfile11 = dpath + 'Results/RM_BL_BD/Evaluations/RegEval/s6_Eval_raw_ukb_24485.csv'
outfile21 = dpath + 'Results/RM_BL_BD/Evaluations/RegEval/s6_Eval_raw_ukb_24485_summary.csv'
mydf['y_test'] = np.exp(mydf['y_test'])/1000
mydf['y_pred'] = np.exp(mydf['y_pred'])/1000
results_df, myout_df = get_avg_output(mydf, nb_iters=1000)

results_df.to_csv(outfile11, index = False)
myout_df.to_csv(outfile21)


