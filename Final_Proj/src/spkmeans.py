import sys
import pandas as pd
import numpy as np
import myspkmeans as spk
from enum import Enum

# Goals (enum) definition
class Goals(Enum):
    SPK = "spk"
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"

    def has_value(item):
        return item in [v.value for v in Goals.__members__.values()]


# Read data from input file
def read_data(file_name):

    input_file = open(file_name, 'r') 
    data_points = []

    while True:
        line = input_file.readline()
        if not line:
            break
        data_point = [float(x) for x in line.split(",")]
        data_points.append(data_point)

    input_file.close()
    
    return data_points


# Initialize centroids as described in the kmeanspp algorithm (HW2)
def init_centroids(data_points, K):
    N, d = data.shape
    # Convert the df into an numpy array (speed)
    
    # dataArr = data.to_numpy() -> from HW2. since now we receive python list we to convert differently -> DELETE LATER
    dataArr = np.array(data_points)

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


def main():
    # Reading user's Command Line arguments
    inputs = sys.argv

    assert len(inputs) == 4 or len(inputs) == 5, "The program can have only 4 or 5 arguments!"

    assert inputs[1].isnumeric(), "K must be an integer"

    K = int(inputs[1])
    assert K >= 0, "K must be non-negative!"

    max_iter = 300

    goal = inputs[2]
    file_name = inputs[3]

    if len(inputs) == 5:
        assert inputs[2].isnumeric(), "max_iter must be an integer"
        assert max_iter >= 0, "max_iter must be non-negative!"

        max_iter = int(inputs[2])
        goal = inputs[3]
        file_name = inputs[4]

    assert not Goals.has_value(goal), "Invalid goal input!"

    data_points = read_data(file_name)

    # TODO: can we assume unempty input or validate?
    N = len(data_points)
    d = len(data_points[0])

    assert K < N, "K must be smaller than N!"

    initial_centroids, centroids_indices = init_centroids(data_points, K)

    result = km.fit(initial_centroids, data_points, centroids_indices, N, d, K, max_iter, goal) # added goal here

    if goal == "spk":
        # Print indices centroids
        for i in range(len(centroids_indices)):
            if (i == len(centroids_indices) -1):
                print(centroids_indices[i])
            else:
                print(centroids_indices[i], end=",")

    # Print centroids \ required matrix
    for i in range(len(result)):
        for j in range(len(result[i])):
            if(j == len(result[i])-1):
                print(np.round(result[i][j], 4))
            else:
                print(np.round(result[i][j], 4), end=",")
    
