import copy
import numpy as np
import aux_bf
from graph_class import HGraph

###########################################################
###########################################################
# AUXILIAR FUNCTIONS TO IMPLEMENT GLOBAL BRANCHING METHOD #
###########################################################
###########################################################





# Given a matrix of the current graph and a fixed arc (i,j)
# it eliminates some arcs: all arcs that emanate from i
# all arcs that reach j and (j,i)


def Update_Matrix(matrix,fixed_arc): 
    new_matrix=np.array(matrix) # copy the original matrix           
    new_matrix[fixed_arc[0]-1]=0 # eliminate all arcs emanating from i
    candidates=np.where(new_matrix==fixed_arc[1])  # arcs that reach j        
    for i in range(len(candidates[0])):
        if len(np.nonzero(new_matrix[candidates[0][i]]))>0:
            new_matrix[candidates[0][i],candidates[1][i]]=0 # eliminate all arcs that reach j
    new_matrix[fixed_arc[0]-1,0]=fixed_arc[1] # restablish conection (i,j)
    if len(np.nonzero(new_matrix[fixed_arc[1]-1]))>0:
        if len(np.where(new_matrix[fixed_arc[1]-1]==fixed_arc[0])[0])>0:
            new_matrix[fixed_arc[1]-1,np.where(new_matrix[fixed_arc[1]-1]==fixed_arc[0])[0]]=0 # eliminate (j,i)
    for i in range(new_matrix.shape[0]):
        new_matrix[i]= aux_bf.pushZerosToEnd(new_matrix[i]) # each row of the matrix, the zeros to the end
    idx = np.argwhere(np.all(new_matrix[..., :] == 0, axis=0)) # eliminate columns with all zeros
    new_matrix = np.delete(new_matrix, idx, axis=1)
    return new_matrix


# Given the updated matrix obtain the fixed arcs

    
def Update_arcs(new_matrix): 
    graph=HGraph(new_matrix)
    graph.FixedArcs()
    new_fixed_arcs=graph.new_fixed_arcs
    return new_fixed_arcs

# Given a matrix of the current graph and a set of fixed arcs U
# it eliminates some arcs as in Update_Matrix    

    
def Update_Matrix_From_List(matrix,list_arcs):
    new_matrix=np.array(matrix)           
    for fixed_arc in list_arcs:
        new_matrix[fixed_arc[0]-1]=0     
        candidates=np.where(new_matrix==fixed_arc[1])          
        for i in range(len(candidates[0])):
            if len(np.nonzero(new_matrix[candidates[0][i]]))>0:                      
                new_matrix[candidates[0][i],candidates[1][i]]=0
        new_matrix[fixed_arc[0]-1,0]=fixed_arc[1]
        if len(np.nonzero(new_matrix[fixed_arc[1]-1]))>0:
            if len(np.where(new_matrix[fixed_arc[1]-1]==fixed_arc[0])[0])>0:
                new_matrix[fixed_arc[1]-1,np.where(new_matrix[fixed_arc[1]-1]==fixed_arc[0])[0]]=0
        for i in range(new_matrix.shape[0]):
            new_matrix[i]= aux_bf.pushZerosToEnd(new_matrix[i])  
    idx = np.argwhere(np.all(new_matrix[..., :] == 0, axis=0))
    new_matrix = np.delete(new_matrix, idx, axis=1)
    return new_matrix  

# Given the current matrix a fixed arc (branching candidate) and current fixed arcs U
# it counts how many arcs is possible to fix with the branching candidate
    
      

def Arcs_Fixing(A,fixed_arc,fixed_arcs):
    count=0 
    new_matrix=Update_Matrix(A,fixed_arc)
    new_fixed_arcs=Update_arcs(new_matrix)
    last_fixed_arcs =  list(set(tuple(row) for row in new_fixed_arcs)-set(tuple(row) for row in fixed_arcs)-set(fixed_arc))     
    count=count+len(last_fixed_arcs)
    initial=copy.copy(new_fixed_arcs)
    while len(last_fixed_arcs)>0:
        new_matrix=Update_Matrix_From_List(new_matrix,last_fixed_arcs)
        current=Update_arcs(new_matrix)
        if initial==current:
            last_fixed_arcs=[]
        else:
            last_fixed_arcs =  list(set(tuple(row) for row in current)-set(tuple(row) for row in initial))
            count=count+len(last_fixed_arcs)
        initial=copy.copy(current)
    return count

