import symnmf
import sys
import numpy as np
from sklearn.metrics import silhouette_score
import kmeans

def general_error(msg="An Error Has Occured"):
    """
    Prints an error message and exits the program.

    Args:
        msg (str, optional): Error message. Defaults to "An Error Has Occured".
    """
    print(msg)
    exit(1)

def calc_coeff(points_lst, lables):
    """
    Calculates the silhouette score for the clustering.

    Args:
        points_lst (list of list of float): Dataset.
        lables (list of int): Cluster labels.

    Returns:
        float: Silhouette score.
    """
    points_lst_np = np.array(points_lst)
    labels_np = np.array(lables)
    score = silhouette_score(points_lst_np, labels_np)
    return score

def h_matrix_to_lables(h_mat):
    """
    Converts an H matrix to cluster labels by assigning each point to the cluster with the highest value.

    Args:
        h_mat (list of list of float): H matrix.

    Returns:
        list of int: Cluster label for each point.
    """
    lables = [-1 for _ in range(len(h_mat))]
    for i in range(len(h_mat)):
        max_val = -1
        max_idx = -1
        for j in range(len(h_mat[i])):
            if h_mat[i][j] > max_val:
                max_idx = j
                max_val = h_mat[i][j]
        lables[i] = max_idx
    return lables

def main():
    """
    Main function to compare k-means and symmetric NMF clustering using silhouette score.
    """
    if len(sys.argv) != 3:
        general_error()
    k = int(sys.argv[1])  # we can assume k is a legal value
    path = sys.argv[2]

    # kmeans
    points_lst, kmeans_lables = kmeans.run_kmeans_alg(path, k)
    kmeans_score = calc_coeff(points_lst, kmeans_lables)

    # symnmf
    dataset = symnmf.txt_input_to_list(path)
    h_mat = symnmf.get_goal_matrix("symnmf", k, dataset)
    symnmf_lables = h_matrix_to_lables(h_mat)
    symnmf_score = calc_coeff(points_lst, symnmf_lables)

    # printing
    print(f"nmf: {symnmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")

if __name__ == "__main__":
    main()
