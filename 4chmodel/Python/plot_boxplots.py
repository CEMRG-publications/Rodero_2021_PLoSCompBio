#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
import pandas._libs

normalised_residuals = pandas.read_csv("/data/normalised_residuals.txt",sep="\t",header=None)
# print(normalised_residuals)

# plt.boxplot(normalised_residuals)
# plt.show()