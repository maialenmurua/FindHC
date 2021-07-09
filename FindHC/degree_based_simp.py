import numpy as np
import aux_bf


##############################################
##############################################
## SIMPLIFICATION OF VERTICES OF DEGREE TWO ##
##############################################
##############################################




# It finds the vertices of 2-degree in a matrix
    
def Find_2degree_Neighbors(mat):    
    n,N = mat.shape     # n:Number of vertices,  N: max number of neighbors
    in_degrees = np.zeros((n)).astype(int)    # For each vertex, the in-degree
    out_degrees = np.zeros((n)).astype(int)   # For each vertex, the out-degree
    
    for i in  range(n):
        for j in range(N):
            if mat[i,j]>0:
                out_degrees[i] =  out_degrees[i] + 1
                in_degrees[mat[i,j]-1] =  in_degrees[mat[i,j]-1] + 1
        
    
    neighbors_degree_2=np.zeros(n)             # The vertices of degree 2 of vertex i are stored in  
                                      # the list neighbors_degree_2
    for i in  range(n):
        if (in_degrees[i]==1) and (out_degrees[i]==1):
            if (mat[mat[i,0]-1,0] != i+1):                # For being a 2-degree node
                neighbors_degree_2[i]=1              # if the incoming arc is (i,j)
                                                          # the outgoing arc can not be (j,i)
                                                          # i.e., doubly connected does not count
            
        elif (in_degrees[i]==1) and (out_degrees[i]==2):  # For being a 2-degree node
                                                          # one of the 2 out arcs should go back to i
            out1 = mat[i,0]-1    
            out2 = mat[i,1]-1 
            out1_goes_to_i = np.sum(mat[out1,:]==i+1)>0
            out2_goes_to_i = np.sum(mat[out2,:]==i+1)>0
            if out1_goes_to_i and out2_goes_to_i:           # We apply XOR operator (^)
                neighbors_degree_2[i]=1
                
            
        elif (in_degrees[i]==2) and (out_degrees[i]==2):  # For being a 2-degree node
                                                          # the two out arcs should go to nodes
                                                          # that have arcs that go i
            out1 = mat[i,0]-1    
            out2 = mat[i,1]-1 
            out1_goes_to_i = np.sum(mat[out1,:]==i+1)>0
            out2_goes_to_i = np.sum(mat[out2,:]==i+1)>0
            if out1_goes_to_i and  out2_goes_to_i:
                neighbors_degree_2[i]=1
                
        elif (in_degrees[i]==2) and (out_degrees[i]==1):  # For being a 2-degree node
                                                          # one of the two in arcs should come from the
                                                          # node that receives the arc from i
            out1 = mat[i,0]-1                
            out1_goes_to_i = np.sum(mat[out1,:]==i+1)>0            
            if out1_goes_to_i:
                neighbors_degree_2[i]=1
    
                
                
    return neighbors_degree_2  
    
# It simplifies the matrix taking into account the vertices of degree 2 
    
def Degree2_Vertices(A):        
    N=A.shape[0]
    while True:  
        Vertices_Updated = 1 # while matrix A is updated due to arc elimination vertices_updated=1

        feasible_solution = 1 # we assume that is feasible to find an HC in A, unless it is shown in the analysis.
        
        vertices_degree_2=Find_2degree_Neighbors(A) # All the vertices with deg 2 in the current A matrix.     
        while (Vertices_Updated):
            Vertices_Updated = 0
            for i in range(N):  # for each vertex i
                Vertices_Updated_i = 0
                all_neighbors_i = list(A[i,np.nonzero(A[i])[0]]-1) # all the vertices that emanate from i
                j = 0
                l_neighbors_degree_two = 0
                neighbors_degree_2 = []                 
                while (j<len(all_neighbors_i)) & (l_neighbors_degree_two<=2):  
                    value = all_neighbors_i[j]                    
                    if vertices_degree_2[i]==1:
                        neighbors_degree_2 = neighbors_degree_2 + [value]
                        l_neighbors_degree_two = l_neighbors_degree_two + 1    
                    j = j +1 	
                         
                if (l_neighbors_degree_two>2): # in this case is not possible to find an HC
                    feasible_solution = 0
                    return [], feasible_solution
                    
                    
                # all the vertices that emanate from i and are of degree 2        
                reached_by_i=np.where(A==i+1)[0] # all the vertices that arrive to i                                     
                j = 0
                l_reached_by_degree_two = 0
                reached_by_degree_2 = [] 
                while (j<len(reached_by_i)) & (l_reached_by_degree_two<=2):  
                    value = reached_by_i[j]
                    if vertices_degree_2[i]==1:
                        reached_by_degree_2 = reached_by_degree_2 + [value]
                        l_reached_by_degree_two = l_reached_by_degree_two + 1    
                    j = j +1 	
                if (l_reached_by_degree_two>2): # in this case is not possible to find an HC
                    feasible_solution = 0
                    return [], feasible_solution                    
                
                elif (l_neighbors_degree_two==2):  # if there are two vertices that emanate from i of degree 2   
                    pos1=np.where(A[i]==neighbors_degree_2[0]+1)[0] # identify the positions where the vertices
                    pos2=np.where(A[i]==neighbors_degree_2[1]+1)[0] # of deg 2 appear
                    A[i]=0  # Eliminate all the arcs that emanate from i 
                    A[i,pos1]=neighbors_degree_2[0]+1 # Restablish the connections of the vertices of degree 2 
                    A[i,pos2]=neighbors_degree_2[1]+1
                    A[i]=aux_bf.pushZerosToEnd(A[i]) # in matrix A, all zeros are at the end
                        
                    if len(all_neighbors_i)>2: # if there are more vertices that emanate from i apart from those of deg 2 
                        Vertices_Updated_i = 1 # at least one arc was eliminated and vertices_updated=1

                if (l_reached_by_degree_two==2): # if there are two vertices that arrive to i of degree 2     
                    pos1=np.where(A[reached_by_degree_2[0]]==i+1)[0] # identify the positions where the vertices of degree 2 
                    pos2=np.where(A[reached_by_degree_2[1]]==i+1)[0] # are reached by i. 
                    for k in range(len(reached_by_i)):
                        A[reached_by_i[k], np.where(A[reached_by_i[k]]==i+1)[0]]=0 #eliminate all the vertices that arrive  
                                                                                             # to i 

                    A[reached_by_degree_2[0],pos1]=i+1 # restablish the connections of vertices of degree 2. 
                    A[reached_by_degree_2[1],pos2]=i+1
                    A[reached_by_degree_2[0]]=aux_bf.pushZerosToEnd(A[reached_by_degree_2[0]])  # in matrix A, all zeros are 
                    A[reached_by_degree_2[1]]=aux_bf.pushZerosToEnd(A[reached_by_degree_2[1]])  #are at the end

                    for k in range(len(reached_by_i)):  # in matrix A, all zeros are 
                        A[reached_by_i[k]]=aux_bf.pushZerosToEnd(A[reached_by_i[k]]) #are at the end
                   
                            
                    if len(reached_by_i)>2: # if there are more vertices that arrive to i apart from those of deg 2
                        Vertices_Updated_i = 1 # at least one arc was eliminated and vertices_updated=1                           
                    
                if Vertices_Updated_i==1:    
                    Vertices_Updated = Vertices_Updated + Vertices_Updated_i                   
	               
                vertices_degree_2=Find_2degree_Neighbors(A) # update the list of vertices with deg 2 in the updated matrix
                        
        idx = np.argwhere(np.all(A[..., :] == 0, axis=0)) # identify if there is a column of A with all zeros
        A = np.delete(A, idx, axis=1) # if yes, eliminate that column
                                
        return A, feasible_solution     # return the updated matrix A and if it is feasible to find an HC.
