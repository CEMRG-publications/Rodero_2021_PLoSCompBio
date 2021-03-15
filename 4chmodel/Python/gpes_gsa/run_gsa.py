from gpe import GPEmul, GPEmulUtils
from itertools import combinations
import matplotlib.gridspec as grsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from SALib.sample import saltelli
from SALib.analyze import sobol
from scipy.special import binom
import seaborn as sns
import sys
import timeit
import torch
from utils.design import get_minmax
from utils.design import read_labels


SEED = 8


def f(emul, X, n_draws):
	y_samples = emul.sample(X, n_draws=n_draws)
	return y_samples


def main():
	seed = SEED
	random.seed(seed)
	np.random.seed(seed)
	torch.manual_seed(seed)

	#========================
	# GPE loading
	#========================
	tag = sys.argv[1]
	idx_feature = sys.argv[2]

	path = 'tmp/' + tag + '/' + idx_feature + '/'
	r2_score_cv = np.loadtxt(path + 'r2_split_test_scores.txt', dtype=float)
	best_split = np.argmax(r2_score_cv)

	emul = GPEmulUtils.load_emulator(path, best_split)

	#========================
	# GSA - SALIB
	#========================
	start_time = timeit.default_timer()
	label = read_labels('data/'+ tag +'/ylabels.txt')	

	N = 1000
	n_draws = 1000

	X_train = np.loadtxt(path + '{}_X_train.txt'.format(best_split), dtype=float)
	I = get_minmax(X_train)

	index_i = read_labels('data/'+ tag +'/xlabels.txt')
	index_ij = ['({}, {})'.format(c[0], c[1]) for c in combinations(index_i, 2)]

	D = len(index_i)
	problem = {
		'num_vars': D,
		'names': index_i,
		'bounds': I
	}

	X_sobol = saltelli.sample(problem, N, calc_second_order=True) # N x (2D + 2) | if calc_second_order=False --> N x (D + 2)
	Y = f(emul, X_sobol, n_draws)

	ST = np.zeros((0, D), dtype=float)
	S1 = np.zeros((0, D), dtype=float)
	S2 = np.zeros((0, int(binom(D, 2))), dtype=float)
	for i in range(n_draws):
		S = sobol.analyze(problem, Y[i], calc_second_order=True, parallel=True, n_processors=24, seed=seed)
		total_order, first_order, (_, second_order) = sobol.Si_to_pandas_dict(S)
		ST = np.vstack((ST, total_order['ST'].reshape(1, -1)))
		S1 = np.vstack((S1, first_order['S1'].reshape(1, -1)))
		S2 = np.vstack((S2, np.array(second_order['S2']).reshape(1, -1)))

	print('GSA - Elapsed time: {:.4f} sec'.format(timeit.default_timer() - start_time))

	np.savetxt(path + 'STi_salib_{}.txt'.format(N), ST, fmt='%.6f')
	np.savetxt(path + 'Si_salib_{}.txt'.format(N), S1, fmt='%.6f')
	np.savetxt(path + 'Sij_salib_{}.txt'.format(N), S2, fmt='%.6f')

	df_ST = pd.DataFrame(data=ST, columns=index_i)
	df_S1 = pd.DataFrame(data=S1, columns=index_i)
	df_S2 = pd.DataFrame(data=S2, columns=index_ij)

	plt.style.use('seaborn')
	gs = grsp.GridSpec(2, 2)
	fig = plt.figure(figsize=(2*8.27, 4*11.69/3))
	ax0 = fig.add_subplot(gs[0, 0])
	ax1 = fig.add_subplot(gs[0, 1])
	ax2 = fig.add_subplot(gs[1, :])
	sns.boxplot(ax=ax0, data=df_S1)
	sns.boxplot(ax=ax1, data=df_ST)
	sns.boxplot(ax=ax2, data=df_S2)
	ax0.set_ylim(0, 1)
	ax0.set_title('First-order effect', fontweight='bold', fontsize=12)
	ax0.set_xticklabels(ax0.get_xticklabels(), rotation=45, horizontalalignment='right')
	ax1.set_ylim(0, 1)
	ax1.set_title('Total effect', fontweight='bold', fontsize=12)
	ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, horizontalalignment='right')
	ax2.set_ylim(0, 1)
	ax2.set_title('Second-order effect', fontweight='bold', fontsize=12)
	ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, horizontalalignment='right')
	plt.savefig(path + idx_feature + '_' + label[int(idx_feature)] + '_si_distr_salib_{}.pdf'.format(N))

if __name__ == '__main__':
	main()