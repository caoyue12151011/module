# Original code from: http://peekaboo-vision.blogspot.com/2012/02/
# simplistic-minimum-spanning-tree-in.html
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform


def MST(P,copy_X=True):
    """X are edge weights of fully connected graph"""
    X = squareform(pdist(P))
    if copy_X:
        X = X.copy()

    if X.shape[0] != X.shape[1]:
        raise ValueError("X needs to be square matrix of edge weights")
    n_vertices = X.shape[0]
    spanning_edges = []
    
    # initialize with node 0: 
    visited_vertices = [0]
    num_visited = 1
    # exclude self connections:
    diag_indices = np.arange(n_vertices)
    X[diag_indices, diag_indices] = np.inf
    
    while num_visited != n_vertices:
        new_edge = np.argmin(X[visited_vertices], axis=None)
        # 2d encoding of new_edge from flat, get correct indices
        new_edge = divmod(new_edge, n_vertices)
        new_edge = [visited_vertices[new_edge[0]], new_edge[1]] 
        # add edge to tree
        spanning_edges.append(new_edge)
        visited_vertices.append(new_edge[1])
        # remove all edges inside current tree
        X[visited_vertices, new_edge[1]] = np.inf
        X[new_edge[1], visited_vertices] = np.inf 
        num_visited += 1
    return np.vstack(spanning_edges)


def test_mst():
    P = np.random.uniform(size=(50, 2))
    edge_list = MST(P)
    plt.scatter(P[:, 0], P[:, 1])
    
    for edge in edge_list:
        i, j = edge
        plt.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], c='r')
        plt.text(P[i, 0], P[i, 1],str(edge[0]))
        plt.text(P[j, 0], P[j, 1],str(edge[1]))
    plt.show()
    print(edge_list)
































    
