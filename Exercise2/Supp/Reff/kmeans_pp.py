import argparse

import numpy as np
import pandas as pd
import mykmeanssp as km

parser = argparse.ArgumentParser()
parser.add_argument('K', type=int)
parser.add_argument('N', type=int)
parser.add_argument('d', type=int)
parser.add_argument('MAX_ITER', type=int)
parser.add_argument('filename', type=str)
args = parser.parse_args()

if args.K >= args.N or args.MAX_ITER < 1 or args.d < 1 or args.K < 1 or args.filename == "":
    print("The arguments you have given are invalid. Please make sure that: " + '\n')
    print("1. k <= n." + '\n')
    print("2. All arguments except n and filename are at least 1, and natural numbers only." + '\n')
    print("3. Make sure you entered a valid filename.")
    exit(0)

observations = pd.read_csv(args.filename, header=None)
observationsArr = observations.to_numpy()
indexes_to_print = np.zeros(args.K, np.intc)


def k_means_pp(K, N, d, MAX_ITER, observationsArr):
    np.random.seed(0)
    centroids = np.zeros((K, d))
    ind = np.random.choice(N)
    indexes_to_print[0] = ind
    centroids[0] = observationsArr[ind]
    distances = np.zeros((K, N))
    distances[0] = np.power((observationsArr - centroids[0]), 2).sum(axis=1)
    for i in range(1, K):
        probs = np.min(distances[:i,], axis=0)
        probsSum = probs / probs.sum()
        ind = np.random.choice(N, p = probsSum)
        centroids[i] = observationsArr[ind]
        indexes_to_print[i] = ind
        distances[i] = np.power((observationsArr - centroids[i]), 2).sum(axis=1)

    return centroids


centroids = k_means_pp(args.K, args.N, args.d, args.MAX_ITER, observationsArr)
observationsArr = observationsArr.tolist()
centroids = centroids.tolist()
indexes_to_print = indexes_to_print.tolist()
to_print = km.kmeanspp(observationsArr, centroids, indexes_to_print, args.MAX_ITER, args.K, args.N, args.d)




