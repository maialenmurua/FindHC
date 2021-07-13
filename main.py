import numpy as np
import pandas as pd
import sys
import shlex


#################################################
#################################################
######### MAIN FUNCTION TO EXECUTE ##############
#################################################
#################################################



# Opens the file with the adjacency matrix 
# of the graph

def ReadFiles(name):
    raw = []
    with open('{0}.txt'.format(name),'r') as f:
        for line in f:
            raw.append(line.split())
    data = pd.DataFrame(raw)
    matrix=np.array(data.values,dtype=float)
    return matrix
    
    
# Converts the adjacency matrix in S matrix

    
def ConvertMatrix(matrix):
    N=matrix.shape[0]
    N_edges = int(np.max(np.sum(matrix,axis=1)))
    A=np.zeros((N,N_edges)).astype(int)
    for i in range(N):
        count = 0
        for j in range(N):                    
            if(matrix[i,j]!= 0):         
                A[i,count] = j+1              
                count = count + 1 
    return A 

# Opens the file with the permutation 
    
def ReadPerm(name):
    raw=[]
    with open('{0}.txt'.format(name),'r') as f:
        for line in f:
             raw.append(line.split())
    array=np.array(pd.DataFrame(raw).values,dtype=int)
    return list(np.concatenate(array)) 
    
    
# Calls to the BF

def Call_BF(file_name,perm,branching_method,ind_degree2,solver,collapsing):
    adj_mat=ReadFiles(file_name)
    if perm=='No':
        matrix=ConvertMatrix(adj_mat)
    else:
        permutation=ReadPerm(perm)
        aux=adj_mat[:,permutation]
        aux_matrix=aux[permutation,:]
        matrix=ConvertMatrix(aux_matrix)
    if collapsing=='No':
        from BF_classes_branching import BF_Algorithms
        algorithms=BF_Algorithms(matrix,0.99,0,branching_method,ind_degree2,solver)
    elif collapsing=='Yes':
        from BF_classes_collapse import BF_Algorithms
        algorithms=BF_Algorithms(matrix,0.99,0,branching_method,ind_degree2,solver)
    print(algorithms.Branch_and_Fix())
    
if __name__ == '__main__':
    
    graph = sys.argv[1]
    perm = sys.argv[2]
    branching= sys.argv[3]
    simp = sys.argv[4]
    solver = sys.argv[5]
    collapse = sys.argv[6]
    Call_BF(graph,perm,branching,simp,solver,collapse)


