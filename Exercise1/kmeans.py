import sys

# Main function
def k_means(K, max_iter=200):

    # Read data from input to a list of datapoints
    data_points = read_data()

    assert K < len(data_points), "K must be smaller than N!"

    # Initialize clusters
    clusters = [[] for i in range(K)] # will create K empty clusters

    # Initialize centroids
    centroids = []
    for i in range(K):
        centroids.append(data_points[i])

    # Main Loop
    seen_changes = True
    iter_count = 0

    while (seen_changes and iter_count < max_iter): # If one of them is false then we should stop
        seen_changes = False
        # Initialize temporary clusters
        temp_clusters = [[] for i in range(K)]
        # Loop through all data points
        for i in range(len(data_points)):
                # Retrieve data point
                data_point = data_points[i]
                cluster_index = find_closest_centroid(centroids, data_point)
                if (i != cluster_index or iter_count == 0): # If we changed data point's cluster
                    seen_changes = True
                temp_clusters[cluster_index].append(data_point)
        # Update centroids with the new clusters
        clusters=temp_clusters
        update_centroids(centroids, clusters)
        iter_count += 1

    # Final return
    return centroids

# Read data from input file
def read_data():

    data_points = []
    data_point = []

    while True:
        try: 
            data_point = input().split(sep=",")
            for i in range(len(data_point)):
                data_point[i] = float(data_point[i])
            data_points.append(data_point)
        except EOFError:
            break

    return data_points


# Find the closest centroid to a data point
def find_closest_centroid(centroids, data_point):
    min_distance = float('inf')
    min_index = 0
    for i in range(len(centroids)):
        curr_distance = compute_distance(centroids[i], data_point)
        if(curr_distance < min_distance):
            min_distance = curr_distance
            min_index = i
    return min_index

# Compute distance (Norm)^2 between a two points in R^d
def compute_distance(u, v):
    sum = 0
    for i in range(len(u)):
        sum += (u[i]-v[i])**2
    return sum

# Updated centroid with using new clusters
def update_centroids(centroids, clusters):
    for i in range(len(clusters)):
        # Init a ZERO vector
        sum = [0 for i in range(len(centroids[0]))]
        # Sum all the vectors in a cluster
        if(len(clusters[i]) != 0): # if the cluster is empty we won't do anything
            for j in range(len(clusters[i])):
                sum = vectors_sum(sum, clusters[i][j])
            # Divide the result vector by the size of the cluster
            for j in range(len(sum)):
                sum[j] = sum[j]/len(clusters[i])
            centroids[i] = sum

# Sum two vector (With changing the original)
def vectors_sum(v, u):
    res = []
    for i in range(len(v)):
        res.append(v[i] + u[i])
    return res

# Read user inputs and call k_means function
inputs = sys.argv

assert len(inputs) == 2 or len(inputs) == 3, "The program can have only 2 or 3 arguments!"

try:
    K = int(inputs[1])
except:
    assert False, "Arguments must be a positive numbers!"

assert K > 0, "K must be positive!"



result = [[]]
if len(inputs) == 2:          # call k_mean with default max_iter
    result = k_means(K)
elif len(inputs) == 3:        # call k_mean with input max_iter
    try:
        max_iter = int(inputs[2])
    except:
        assert False, "Arguments must be a positive numbers!"
    
    assert max_iter >= 0, "Max Iterations must be non-negative!"
    result = k_means(K, max_iter)

# print result
for i in range(len(result)):
    for j in range(len(result[i])):
        if(j == len(result[i])-1):
            print(format(result[i][j], '.4f'))
        else:
            print(format(result[i][j], '.4f'), end=",")