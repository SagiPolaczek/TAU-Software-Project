import sys
import pandas as pd
import numpy as np
import mykmeanssp as km


def combine_inputs(input_1, input_2):
    # Read data
    data_1 = pd.read_csv(input_1, header=None)
    data_2 = pd.read_csv(input_2, header=None)

    # Merge
    data = pd.merge(data_1, data_2, on=0)
    data.set_index(0, inplace=True)

    # Return the data sorted.אבא
    return data.sort_index()


# Initialize centroids as described in the kmeanspp algorithm
def init_centroids(data, K):
    N, d = data.shape
    # Convert the df into an numpy array (speed)
    dataArr = data.to_numpy()

    centroids = np.zeros((K, d))
    centroids_indices = [0 for i in range(K)]
    
    # Generate "Random" seed
    np.random.seed(0)

    # Select m_1 randomly from x_1 -> x_N
    random_idx = np.random.choice(N)

    # Set the 0 centroid
    centroids_indices[0] = random_idx
    centroids[0] = dataArr[random_idx]

    # Init Distance & Probabilities arrays
    D = np.zeros(N)
    P = np.zeros(N)

    Z = 1
    while Z < K:
        # For each vector x_i, compute the minimum distance from a centroid
        for i in range(N):
            x_i = dataArr[i]
            min_dis = float('inf')
            for j in range(Z):
                m_j  = centroids[j]
                # Computing distance between x_i and m_j
                curr_dis = np.power((x_i - m_j), 2).sum()
                if (curr_dis < min_dis):
                    min_dis = curr_dis

            D[i] = min_dis

        # P = Normalized D. Probabilities.
        P = D / D.sum()
        
        new_index = np.random.choice(N, p=P)
        centroids[Z] = dataArr[new_index]
        centroids_indices[Z] = new_index


        Z += 1

    return centroids.tolist() , centroids_indices


# Reading user's Command Line arguments 
inputs = sys.argv

assert len(inputs) == 4 or len(inputs) == 5, "The program can have only 4 or 5 arguments!"

try:
    K = int(inputs[1])
except:
    assert False, "Arguments must be a positive numbers!"

assert K > 0, "K must be positive!"

max_iter = 300

input_1 = inputs[2]
input_2 = inputs[3]

if len(inputs) == 5:       
    try:
        max_iter = int(inputs[2])
    except:
        assert False, "Arguments must be a positive numbers!"

    assert max_iter >= 0, "Max Iterations must be non-negative!"

    input_1 = inputs[3]
    input_2 = inputs[4]


# Combine both input files by INNER JOIN using the first column in each file as a key.
data = combine_inputs(input_1, input_2)
N, d = data.shape
data_points = data.values.tolist()

assert K < N, "K must be smaller than N!"

initial_centroids, centroids_indices = init_centroids(data, K)


final_centroids = km.fit(initial_centroids, data_points, centroids_indices, N, d, K, max_iter)

# Print indices centroids
for i in range(len(centroids_indices)):
    if (i == len(centroids_indices) -1):
        print(centroids_indices[i])
    else:
        print(centroids_indices[i], end=",")

    

# Print centroids
for i in range(len(final_centroids)):
    for j in range(len(final_centroids[i])):
        if(j == len(final_centroids[i])-1):
            # print(format(final_centroids[i][j], '.4f'))
            print(np.round(final_centroids[i][j], 4))
        else:
            print(np.round(final_centroids[i][j], 4), end=",")
