import diversipy
from gpe import GPEmul
import numpy as np
import random
from sklearn.model_selection import KFold, StratifiedKFold, train_test_split
import sys
import timeit
import torch
from utils.concurrent import execute_task_in_parallel


FOLD = 5
LEARNING_RATE = 0.1
MAX_EPOCHS = 1000
N_RESTARTS = 10
PATIENCE = 20
SEED = 8


def cv(X_train, y_train, X_val, y_val, device, n_restarts, lr, max_epochs, patience, split, path):
	print('Split {}...'.format(split))

	np.savetxt(path + '{}_X_train.txt'.format(split), X_train, fmt='%.6f')
	np.savetxt(path + '{}_y_train.txt'.format(split), y_train, fmt='%.6f')
	np.savetxt(path + '{}_X_val.txt'.format(split), X_val, fmt='%.6f')
	np.savetxt(path + '{}_y_val.txt'.format(split), y_val, fmt='%.6f')

	emul = GPEmul(X_train, y_train, device, learn_noise=False, scale_data=True)
	emul.train(X_val, y_val, n_restarts, lr, max_epochs, patience, savepath=path + '{}_checkpoint.pth'.format(split))
	emul.print_stats()

	return emul.metric_score, emul.best_model


def main():
	seed = SEED
	np.random.seed(seed)
	random.seed(seed)
	torch.manual_seed(seed)

	start_time = timeit.default_timer()

	fold = FOLD
	lr = LEARNING_RATE
	max_epochs = MAX_EPOCHS
	n_restarts = N_RESTARTS
	patience = PATIENCE

	tag = sys.argv[1]
	path_in = 'data/' + tag + '/'
	
	X = np.loadtxt(path_in + 'X.txt', dtype=float)
	Y = np.loadtxt(path_in + 'Y.txt', dtype=float)

	##===============================================================
	## FOR TOV
	##===============================================================
	# idx_sep = 17

	# X_ct = np.copy(X[:idx_sep, :])
	# X_art = np.copy(X[idx_sep:, :])

	# Y_ct = np.copy(Y[:idx_sep, :])
	# Y_art = np.copy(Y[idx_sep:, :])
	
	# z = np.concatenate((np.zeros(X_ct.shape[0], dtype=int), np.ones(X_art.shape[0], dtype=int)))

	# X_ct, Y_ct = shuffle(X_ct, Y_ct, random_state=seed)
	# X_art, Y_art = shuffle(X_art, Y_art, random_state=seed)

	# X = np.vstack((X_ct, X_art))
	# Y = np.vstack((Y_ct, Y_art))

	# idx_feature = int(sys.argv[2])
	# y = np.copy(Y[:, idx_feature])

	# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

	# path_out = 'tmp/' + tag + '/' + str(idx_feature) + '/'

	# skf = StratifiedKFold(n_splits=fold)
	##===============================================================

	idx_feature = int(sys.argv[2])
	y = np.copy(Y[:, idx_feature])

	device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

	path_out = 'tmp/' + tag + '/' + str(idx_feature) + '/'

	kf = KFold(n_splits=fold, shuffle=True, random_state=seed)

	inputs = {
		i: (X[idx_train], y[idx_train], X[idx_val], y[idx_val], device, n_restarts, lr, max_epochs, patience, i, path_out)
		for i, (idx_train, idx_val) in enumerate(kf.split(X))
	}
	results = execute_task_in_parallel(cv, inputs)

	r2_score_cv = [results[i][0] for i in range(fold)]
	model_state_cv = [results[i][1] for i in range(fold)]
	
	for i in range(fold):
		filename = '{}_gpe.pth'.format(i)
		torch.save(model_state_cv[i], path_out + filename)

	np.savetxt(path_out + 'r2_split_test_scores.txt', np.array(r2_score_cv), fmt='%.6f')

	print('\n')
	print('R2 split test scores: {}'.format(np.around(r2_score_cv, decimals=4)))
	print('R2 mean cross-validation test score: [{:.4f}]'.format(np.array(r2_score_cv).mean()))
	
	print('Elapsed Time: {:.4f} sec'.format(timeit.default_timer() - start_time))

if __name__ == '__main__':
	main()