import scipy
import numpy as np
import sklearn.datasets
import spkmeans

def print_matrix(matrix):
    for dp in matrix:
        for j in range(len(dp) - 1):
            print('{:.3f},'.format(dp[j]), end='')
        print('{:.3f}'.format(dp[-1]))

def print_array(array):
    for i in range(len(array) - 1):
        print('{:.3f},'.format(array[i]), end='')
    print('{:.3f}'.format(array[-1]))

def createWightedAdjacencyMatrix(observations, n):
     W = np.full((n, n), 0, dtype=np.float64)
     for i in range(n):
         W[i] = np.linalg.norm(observations - observations[i], axis=1) / 2
     W = np.exp(-W)
     np.fill_diagonal(W, val=0)
     return W
