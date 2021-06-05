from sklearn.cluster import KMeans
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt

# Load data from sklearn API
data = load_iris().data

# Compute inertias in a pythonic way
inertias = [KMeans(n_clusters=k, init='k-means++',
                random_state=0).fit(data).inertia_ for k in range(1,11)]
# Define indices
k_indices = [i for i in range(1,11)]

# Plot the graph
plt.plot(k_indices, inertias)
plt.xlabel("K")
plt.ylabel("Inertia")
plt.title('Elbow Method for selection of optimal "K" clusters')
plt.annotate("Elbow",xy=(2,inertias[1]),xytext=(4, 200),arrowprops=dict(arrowstyle="->"))

# Save the graph into a png file
plt.savefig('elbow.png')

print("Done!")