import symnmf
import sys
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score

np.random.seed(1234)
MAX_ITER = 300
EPSILON = 0.0001

'''point object'''
class Point:
    def __init__(self, cord, index = 1):
        self.cord = cord  # list of coordinates (of size dimension)
        self.dim = len(cord)  # dimension size of point
        self.index = index  # index in the original dataset 

    def euc_distance(self, p):    
        """
        calculates the eucledian distance between this point, and another one
        """
    
        dist = 0
        for i in range(self.dim):
            dist += (self.cord[i] - p.cord[i])**2
        
        return dist ** 0.5
    
'''cluster of points'''
class Cluster:
    def __init__(self):
        self.pointsLst = []
        self.centroid = Point([])

    ''' add new point to cluster'''
    def add_point(self, p):
        self.pointsLst.append(p)

    '''calculates new centroid
    returns true if there was a change greater than epsilon in the centroid, false if no change'''
    def calc_centroid(self, eps):
        if(self.pointsLst == []):
            return False

        centPoint = Point([0]*self.centroid.dim) #new centroid
        
        #calculate sum of all point for each coordinate
        for p in self.pointsLst:
            for i in range(p.dim):
                centPoint.cord[i] += p.cord[i]
        
        changed = False
        
        #calculate average for each coordinate
        #check the difference from current centroid
        for i in range(len(centPoint.cord)):
            centPoint.cord[i] = centPoint.cord[i]/len(self.pointsLst)

        if(centPoint.euc_distance(self.centroid) >= eps):
            changed = True
        
        self.centroid = centPoint
        return changed
        
    '''sets all points'''
    def restart_list(self):
        self.pointsLst = []

def general_error(msg = "An Error Has Occured"):
    print(msg)
    exit(1)

def add_point_to_closest_cluster(clstr_lst, p):
    """
    p - a point of type "Point"
    clstr_lst - a list of K clusters each one of type "Cluster"
    The function adds the point to the closest cluster
    """
    min_cluster = None
    min_dist = float('inf')
    for clstr in clstr_lst:
        dist = p.euc_distance(clstr.centroid)
        if dist < min_dist:
            min_cluster = clstr
            min_dist = dist
    
    min_cluster.add_point(p)

'''create k clusters, each one receives one of the first k points'''
def init_clusters(points_lst, k, epsilon):
    cluster_lst = []
    for i in range(0,k):
        clust = Cluster()
        clust.add_point(points_lst[i])
        clust.centroid = Point([0]*points_lst[i].dim)
        clust.calc_centroid(epsilon)
        cluster_lst.append(clust)

    return cluster_lst

def is_float(x):
    try:
        float(x)
        return True
    except:
        return False

def print_cluster_lst(lst):
    for c in lst:
        print(lst_to_str(c.centroid.cord))

def lst_to_str(lst):
    str = ""
    for x in lst:
        str += '%.4f'%x + ','
    return str [:-1]

"""
txt_input_to_list - converts the path given to a list of lists
@path: the path of the file
"""
def txt_input_to_list(path):
    try:
        with open(path, 'r') as f:
            lines = [line.strip() for line in f]
            lines = [[float(x) for x in line.split(',')] for line in lines]
    except:
        general_error()
    
    return lines

def list_to_points_list(lst):
    points_lst = [Point(lst[i], i) for i in range(len(lst))]
    return points_lst

def lineToStr(lst):
    """
    given a list of floats as strings, return a list of floats with same values 
    """ 
    return [float(x) for x in lst]

"""
matrix_to_str - converts a list of lists (float) to formatted str ready for printing
@mat: the matrix to convert
"""
def matrix_to_str(mat):
    rows = []
    for row in mat:
        row_str = ','.join([f'{x:.4f}' for x in row])
        rows.append(row_str)

    return '\n'.join(rows)

"""
creates a list of size N where each cell represents the point in the dataset in this index's cluster index
@cluster_lst: list of clusters
@N: amount of points in the dataset
"""
def cluster_lst_to_lables(cluster_lst, N):
    points_clstr_idx = [-1 for i in range(N)]  # will hold each points final cluster index
    for i in range(len(cluster_lst)):
        for pnt in cluster_lst[i].pointsLst:
            points_clstr_idx[pnt.index] = i
    
    return points_clstr_idx

def kmeans_lables(points_lst, k):
    cluster_lst = init_clusters(points_lst, k, EPSILON)

    # iterate over points and reclsuters them MAX_ITER times
    for i in range(MAX_ITER):
        changed = False

        # restart all lists of points in each cluster
        for c in cluster_lst:
            c.restart_list()

        # adds each point the new assigned cluster
        for p in points_lst:
            add_point_to_closest_cluster(cluster_lst, p)
        
        # recalculates centroid for each cluster
        for c in cluster_lst:
            if ((c.calc_centroid(EPSILON))):
                changed = True

        # if no cenrtroid changed in more than epsilon, stop iterating
        if (changed == False):
            break
    
    lables = cluster_lst_to_lables(cluster_lst, len(points_lst))
    return lables

def calc_coeff(points_lst, lables):
    points_lst_np = np.array(points_lst)
    labels_np = np.array(lables)

    score = silhouette_score(points_lst_np, labels_np)
    return score

def create_random_matrix(rows, cols, upper_bound, lower_bound = 0):
    return [[np.random.uniform(low=lower_bound, high=upper_bound) for _ in range(cols)] for _ in range(rows)]

"""
get_average_val_of_matrix - Calculates the average value of all of the cells in matrix
@mat: the matrix
"""
def get_average_val_of_matrix(mat):
    sum = 0
    for row in mat:
        for val in row:
            sum += val
    cnt = len(mat) * len(mat[0])
    return sum / cnt

def h_matrix_to_lables(h_mat):
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

def nmf_lables(points_lst, k):
    n = len(points_lst)
    w_mat = symnmf.norm(points_lst)
    w_mat_avg_val = get_average_val_of_matrix(w_mat)
    upper_val = 2 * ((w_mat_avg_val / k)**0.5)
    h_mat = create_random_matrix(n, k, upper_val)

    # call the symnmf() function in symnmf module
    res_mat = symnmf.symnmf(h_mat, w_mat, k)

    lables = h_matrix_to_lables(res_mat)
    return lables

def main():
    if len(sys.argv) != 3:
        general_error()
    
    k = int(sys.argv[1])
    path = sys.argv[2]

    points_lst_raw = txt_input_to_list(path)
    points_lst_objects = list_to_points_list(points_lst_raw) 

    # kmeans
    kmeans_score = calc_coeff(points_lst_raw, kmeans_lables(points_lst_objects, k))

    # symnmf
    symnmf_score = calc_coeff(points_lst_raw, nmf_lables(points_lst_raw, k))

    #printing
    print(f"nmf: {symnmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")


if __name__ == "__main__":
    main()
