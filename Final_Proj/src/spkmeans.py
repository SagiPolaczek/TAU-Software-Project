import sys
import numpy as np
import myspkmeans as spk

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
def init_centroids(data_points, N, d, K):
    
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

data_points = read_data(file_name)

N = len(data_points)
assert N > 0, "N must be greater than zero!"
d = len(data_points[0])
assert d > 0, "d must be greater than zero!"

assert K < N, "K must be smaller than N!"

if goal == "spk":
    data_points = spk.fit_init_spk()
    initial_centroids, centroids_indices = init_centroids(data_points, N, d, d) # K = d
    spk.fit_finish_spk(initial_centroids, data_points, centroids_indices, N, d, d, max_iter) # K = d
elif goal == "wam":
    spk.fit_general(data_points, N, d, K, max_iter, ord('w'))
elif goal == "ddg":
    spk.fit_general(data_points, N, d, K, max_iter, ord('d'))
elif goal == "lnorm":
    spk.fit_general(data_points, N, d, K, max_iter, ord('l'))
elif goal == "jacobi":
    spk.fit_general(data_points, N, d, K, max_iter, ord('j'))
else:
    assert False, "Invalid goal input!"