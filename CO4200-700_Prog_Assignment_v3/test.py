from sklearn.cluster import AgglomerativeClustering
import numpy as np
# Protein sequences or identifiers
protein_sequences = ["Seq1", "Seq2", "Seq3", "Seq4", "Seq5"]

D = np.array([[0, 1, 2, 4, 4],
              [1, 0, 2, 4, 4],
              [2, 2, 0, 2, 2],
              [4, 4, 2, 0, 1],
              [4, 4, 2, 1, 0]])

model = AgglomerativeClustering(linkage='average', metric='precomputed', distance_threshold=None)

# Cluster the objects based on the distances
cluster = model.fit(D)
n_objects = len(cluster.labels_)  # This is just the number of sequences

# Cluster.children_ gives the hierarchical clustering
# Leaves are the original objects and labelled from 0 to 4 here
# The internal nodes in the cluster are labelled from 5 onwards
next_node = n_objects
print (n_objects)
for i, merge in enumerate(cluster.children_):

    if merge[0] < n_objects:
        seq1 = merge[0]
    else:
        f"Cluster_{merge[0]}"

    if merge[1] < n_objects:
        seq2 = merge[1]
    else:
        f"Cluster_{merge[1]}"

    print(f"Align {protein_sequences[seq1]} with {protein_sequences[seq2]} to give Cluster_{next_node}")
    next_node += 1

print("done")