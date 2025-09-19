ERR_MSG = "An Error Has Occured"
MAX_ITER = 300
EPSILON = 0.0001

def general_error(msg = ERR_MSG):
    print(msg)
    exit(1)

class Point:
    """
    Represents a point in a multi-dimensional space.

    Attributes:
        cord (list of float): Coordinates of the point.
        dim (int): Dimension of the point.
        index (int): Index of the point in the original dataset.
    """
    def __init__(self, cord, index = -1):
        """
        Initializes a Point object.

        Args:
            cord (list of float): Coordinates of the point.
            index (int, optional): Index in the original dataset. Defaults to -1.
        """
        self.cord = cord
        self.dim = len(cord)
        self.index = index

    def euc_distance(self, p):
        """
        Calculates the Euclidean distance between this point and another point.

        Args:
            p (Point): Another point.

        Returns:
            float: Euclidean distance.
        """
        dist = 0
        for i in range(self.dim):
            dist += (self.cord[i] - p.cord[i]) ** 2
        return dist ** 0.5
    
class Cluster:
    """
    Represents a cluster of points.

    Attributes:
        pointsLst (list of Point): Points in the cluster.
        centroid (Point): Centroid of the cluster.
    """
    def __init__(self):
        """
        Initializes a Cluster object.
        """
        self.pointsLst = []
        self.centroid = Point([])

    def add_point(self, p):
        """
        Adds a point to the cluster.

        Args:
            p (Point): Point to add.
        """
        self.pointsLst.append(p)

    def calc_centroid(self, eps):
        """
        Calculates the new centroid of the cluster.

        Args:
            eps (float): Epsilon threshold for centroid change.

        Returns:
            bool: True if centroid changed more than epsilon, False otherwise.
        """
        if not self.pointsLst:
            return False

        centPoint = Point([0] * self.centroid.dim)
        for p in self.pointsLst:
            for i in range(p.dim):
                centPoint.cord[i] += p.cord[i]

        for i in range(len(centPoint.cord)):
            centPoint.cord[i] = centPoint.cord[i] / len(self.pointsLst)

        changed = centPoint.euc_distance(self.centroid) >= eps
        self.centroid = centPoint
        return changed

    def restart_list(self):
        """
        Removes all points from the cluster.
        """
        self.pointsLst = []

def add_point_to_closest_cluster(clstr_lst, p):
    """
    adds a point to the closest cluster based on Eucledian distance.
    Args:
        clstr_lst (list of Cluster) : List of clusters.
        p (Point): Point to add
    """
    min_cluster = None
    min_dist = float('inf')
    for clstr in clstr_lst:
        dist = p.euc_distance(clstr.centroid)
        if dist < min_dist:
            min_cluster = clstr
            min_dist = dist
    min_cluster.add_point(p)

def make_points_lst(path):
    try:
        with open(path, 'r') as f:
            points_lst = [Point([float(x) for x in line.strip().split(',')], i) for i, line in enumerate(f)]
    except:
        general_error()
    return points_lst

def init_clusters(points_lst, k, epsilon):
    """
    Initializes k clusters, each with one of the first k points.

    Args:
        points_lst (list of Point): List of points.
        k (int): Number of clusters.
        epsilon (float): Epsilon threshold for centroid change.

    Returns:
        list of Cluster: Initialized clusters.
    """
    cluster_lst = []
    for i in range(k):
        clust = Cluster()
        clust.add_point(points_lst[i])
        clust.centroid = Point([0] * points_lst[i].dim)
        clust.calc_centroid(epsilon)
        cluster_lst.append(clust)
    return cluster_lst

def get_lables(cluster_lst, N):
    """
    Converts a list of clusters to a list of cluster labels for each point.

    Args:
        cluster_lst (list of Cluster): List of clusters.
        N (int): Number of points in the dataset.

    Returns:
        list of int: Cluster label for each point.
    """
    points_clstr_idx = [-1 for _ in range(N)]
    for i in range(len(cluster_lst)):
        for pnt in cluster_lst[i].pointsLst:
            points_clstr_idx[pnt.index] = i
    return points_clstr_idx

def get_cluster_list(points_lst, k):
    cluster_lst = init_clusters(points_lst, k, EPSILON)
    
    for _ in range(MAX_ITER):
        changed = False
        
        # 1. empty clusters from points
        for c in cluster_lst:
            c.restart_list()
        
        # 2. add each point to the closest cluster
        for p in points_lst:
            add_point_to_closest_cluster(cluster_lst, p)
        
        # 3. recalculate each cluster's centroid
        for c in cluster_lst:
            if c.calc_centroid(EPSILON):
                changed = True
        
        if not changed:  # the centroid's didnt change by epsilon, and so no need to continue the loop
            break

    return cluster_lst
        
def points_lst_to_regular_lst(points_lst):
    return [p.cord for p in points_lst]

def run_kmeans_alg(path, k):    
    points_lst = make_points_lst(path)
    points_lst_regular = points_lst_to_regular_lst(points_lst)

    # validate k
    if k >= len(points_lst):
        general_error()
    
    # get final cluster list
    cluster_lst = get_cluster_list(points_lst, k)

    # get lables
    lables = get_lables(cluster_lst, len(points_lst))

    # return tuple
    return points_lst_regular, lables