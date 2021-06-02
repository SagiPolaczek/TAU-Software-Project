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


    return data



def init_centroids(data, K):
    N, d = data.shape

    centroids = [[0 for i in range(d)] for i in range(K)]
    centroids_indices = [0 for i in range(K)]
    
    # Select m_1 randomly from x_1 -> x_N
    np.random.seed(0)
    random_idx = np.random.choice(N)
    centroids_indices[0] = random_idx

    m_1 = data.loc[random_idx].tolist() #####
    centroids[0] = m_1

    D = [0 for i in range(N)] ##### maybe numpy?
    P = [0 for i in range(N)] ##### maybe numpy?
    Z = 1
    while Z < K:
        # For each x_i
        for i in range(N):
            x_i = data.loc[i].tolist()
            min_dis = float('inf')
            for j in range(1, Z+1):
                m_j  = centroids[j-1]
                curr_dis = compute_distance(x_i, m_j)
                min_dis = min(curr_dis, min_dis)

            D[i] = min_dis


        sum_D = sum(D)

        for i in range(N):
            P[i] = (D[i] / sum_D)
        
        new_index = np.random.choice(N, p=P)
        centroids[Z] = data.loc[new_index].tolist()
        centroids_indices[Z] = new_index


        Z += 1

        

    return centroids, centroids_indices


# Compute distance (Norm)^2 between a two points in R^d
def compute_distance(u, v):
    sum = 0
    for i in range(len(u)):
        sum += (u[i]-v[i])**2
    return sum




# Reading user's Command Line arguments <3 
inputs = sys.argv

assert len(inputs) == 4 or len(inputs) == 5, "The program can have only 4 or 5 arguments!" # SCEPTIC

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


# STAGE 2 - Combine both input files by INNER JOIN using the first column in each file as a key.
data = combine_inputs(input_1, input_2)
N, d = data.shape
data_points = data.values.tolist()

assert K < N, "K must be smaller than N!"

initial_centroids, centroids_indices = init_centroids(data, K)


final_centroids = km.fit(initial_centroids, data_points, centroids_indices, N, d, K, max_iter)

# Print indices centroids !!!
for i in range(len(centroids_indices)):
    if (i == len(centroids_indices) -1):
        print(centroids_indices[i])
    else:
        print(centroids_indices[i], end=",")

    

# Print centroids
for i in range(len(final_centroids)):
    for j in range(len(final_centroids[i])):
        if(j == len(final_centroids[i])-1):
            print(format(final_centroids[i][j], '.4f'))
        else:
            print(format(final_centroids[i][j], '.4f'), end=",")
