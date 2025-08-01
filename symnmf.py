import numpy as np
import sys
import symnmf

ERR_MSG = "An Error Has Occured"

def general_error(msg = ERR_MSG):
    print(msg)
    exit(1)

"""
is_float - Checks if a string x represents a floating number
@x: the string
"""
def is_float(x):
    try:
        float(x)
        return True
    except:
        return False

"""
set_data - checks cmd args valididty
@data: an array of cmd args
"""
def set_data(data):
    # 1. making sure the input is in correct length
    if len(data) != 4:
        general_error()

    # 2. read arguments from CMD
    k = int(sys.argv[1])  # int, number of clusters
    goal = sys.argv[2]  # 'symnf' / 'sym' / 'ddg' / 'norm'
    path = sys.argv[3]  # path to data set


    # 3. check validity
    if not is_float(k):
        general_error()
    
    k = float(k)
    if (k <= 1 or k % 1 != 0):
        general_error()

    return int(k), goal, path

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


def main():
    k, goal, path = set_data(sys.argv)
    dataset = txt_input_to_list(path)
    n = len(dataset)
    
    if k >= n:
        general_error()

    if goal == "symnf":
        pass #TODO
    elif goal == "sym":
        print(symnmf.sym(dataset))
    elif goal == "ddg":
        print(symnmf.ddg(dataset))
    elif goal == "norm":
        print(symnmf.norm(dataset))

if __name__ == "__main__":
    main()
