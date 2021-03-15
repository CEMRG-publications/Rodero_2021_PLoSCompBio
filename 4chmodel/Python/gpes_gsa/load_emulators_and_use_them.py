from gpe import GPEmul
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import random
import torch
from utils.metrics import R2Score, MSE, MAPE


np.set_printoptions(formatter={'all':lambda x: '{:.4f}'.format(x)})


def main():
	seed = 8
	np.random.seed(seed)
	random.seed(seed)
	torch.manual_seed(seed)

	#----------------------------------------------------------------
	ylab = ['EDV_LV','Myo_vol_LV','ESV_LV','SV_LV','V1_LV','EF1_LV','dPdtmax_LV','dPdtmin_LV','PeakP_LV','tpeak_LV','QRS_LV','AT1090_LV','AT_LV','EDV_RV','Myo_vol_RV','ESV_RV','SV_RV','EF_RV','V1_RV','EF1_RV','dPdtmax_RV','dPdtmin_RV','PeakP_RV','tpeak_RV','IRT_RV']
	# X_train, X_val, y_train, y_val = load train and val datasets

	# n_restarts = 3
	# lr = 0.1
	# max_epochs = 1000
	# patience = 20
	# device = torch.device('cpu')

	# #--------
	# # TRAINING

	# emul = GPEmul(X_train, y_train, device, learn_noise=False, scale_data=True)
	# emul.train(X_test, y_test, n_restarts, lr, max_epochs, patience, savepath='./checkpoint.pth', watch_metric='R2Score')
	# emul.print_stats()

	# #--------
	# # SAVING TRAINED GPE

	# model = emul.best_model
	# filename = 'gpe.pth'
	# torch.save(model, filename)

	# #--------
	# # LOADING (you need the training dataset and validation dataset)

	# emul = GPEmul(X_train, y_train, device, learn_noise=False, scale_data=True)
	# model = 'gpe.pth'
	# emul.model.load_state_dict(torch.load(model, map_location=device))
	# emul.model.to(device)

	#----------------
	# SMART LOADING

	tag = 'tov7_final'
	idx_feature = 0

	path = '/data/PCA_GP/' + tag 
	mse_score_cv = np.loadtxt(path + '/' + str(idx_feature) + '/' + 'mse_split_test_scores.txt', dtype=float)
	best_split = np.argmin(mse_score_cv)

	emul = GPEmulUtils.load_emulator(path, best_split, watch_metric='MSE')

	#------------------------------------------------------------------
	# LOAD NEW POINTS

	X_new = np.loadtxt(path + '/' + 'x_new.txt',dtype = float)
	y_original = np.loadtxt(path + '/' + str(idx_feature) + '/' + 'y_new.txt',dtype = float)
	y_new = np.copy(y_original[:, idx_feature])

	#----------------------------------------------------------------
	# PLOT PREDICTIONS VS REAL

	# height = 9.36111
	# width = 5.91667
	# fig, axis = plt.subplots(1, 1, figsize=(1.5*width/2, 1.5*height/3))

	mean, std = emul.predict(X_new)
	mse = MSE(emul.tensorize(y_new), emul.tensorize(mean))
	print('MSE is {:.4f}'.format(mse))
	l = np.argsort(mean)

	CI = 3
	axis.scatter(np.arange(len(l)), y_test[l], facecolors='none', edgecolors='C0')
	axis.scatter(np.arange(len(l)), mean[l], facecolors='C0', s=16)
	axis.errorbar(np.arange(len(l)), mean[l], yerr=CI*std[l], c='C0', ls='none', lw=0.5)
	axis.set_ylabel(ylab[idx_feature], fontsize=12)
	axis.set_title('MSE = {:.4f}'.format(mse))
	axis.set_xticks([])
	axis.set_xticklabels([])
	plt.figlegend(['observed', 'predicted', 'uncertainty ({} STD)'.format(CI)], loc='upper center')
	# plt.show()
	plt.savefig('/home/crg17/Pictures/testingto2save')

if __name__ == '__main__':
	main()