import numpy as np
import networkx as nx




####################################################
####################################################
## IMPLEMENTS FUNCTIONS FOR BF COLLAPSE ALGORITHM ##
####################################################
####################################################


# It finds the fixed paths 
   
def Find_Path_Candidates(mat,n):    
                  
    if mat.shape[1]>1:    
        
        fixed_indices=np.where(mat[:,1]==0)[0]
        unique,counts=np.unique(mat[fixed_indices,0],return_counts=True)
        if any(c>1 for c in list(counts))==True:
            unfeasible=True
            graph_is_circuit=False
            at_least_one_path=False
            n_intermediate_nodes=0
            path_candidates=np.zeros((n)).astype(int)
            return unfeasible,graph_is_circuit,at_least_one_path,n_intermediate_nodes,path_candidates

        else:
            candidates = []   # Fixed arcs are those for which there is only one connection
            [uniq_vals,counts] = np.unique(mat[:,0],  return_counts=True) 
            repeated = np.where(counts>1)[0]
        
            if len(repeated)>0:
                repeated = uniq_vals[repeated]           
            for i in range(n):
                if mat[i,0]>0 and mat[i,1]==0:
                    if not(mat[i,0] in repeated):
                        candidates.append(i)
       
            path_candidates = np.zeros((n)).astype(int) # path_candidates initialized as vector of zeros
            unfeasible=False
            graph_is_circuit = False                    # i.e., mat[i,1]==0
    else:              

        orig_mat=Convert_Matrix(mat)
        G=nx.Graph(orig_mat)
        if  ((len(np.unique(mat[:,0])))== n) & (len(np.where(mat==0)[0])==0) &(nx.is_connected(G)==True):
            graph_is_circuit = True
            unfeasible=False 
        else:
            unfeasible=True
            graph_is_circuit=False                    # If the matrix has only  one column, and assuming the graph
        at_least_one_path = False                     # is connected, we have a HC and there is not any path, 
        path_candidates = 2*np.ones((n)).astype(int)  # All vertices are "intermediate"  
        n_intermediate_nodes = n
        return unfeasible,graph_is_circuit,at_least_one_path,n_intermediate_nodes, path_candidates
            

    # To count the appearences of each vertex in a subpath, we visit all vertices
    # connected to a single vertex and increase the frequency of both
    for vertex in candidates:
        
        path_candidates[vertex] = path_candidates[vertex] + 1
        path_candidates[mat[vertex,0]-1] = path_candidates[mat[vertex,0]-1] + 1
       

    
    # The number of intermediate nodes is equal the number of fixed vertices with degree 2
    n_intermediate_nodes = np.sum(path_candidates==2)
    
    if n_intermediate_nodes==n:        # The graph is a  set of circuits (more than one)
        at_least_one_path = False
    else:                                                  # There is at least one path if at least 1 vertex
        at_least_one_path =  np.sum(path_candidates==1)>0  # have degree 1
    
    return unfeasible,graph_is_circuit,at_least_one_path,n_intermediate_nodes, path_candidates
    


def FindSubpaths(mat,path_candidates,n):
    
    # The candidate list is initialized with the indices of all vertices
    # that appear at least once in a fixed subpath
    candidates = [i for i in range(n) if path_candidates[i]>0]

    # The list of paths is initialized empty
    list_paths = []
    
    to_visit = True   # A variables that indicates whether the search for fixed subpaths  should proceed
    while to_visit:
        one_found = False          # Within this loop, initially, no  subpath has been found 
        for vertex in candidates:  # Each vertex in the list of candidates is inspected            
            if path_candidates[vertex]==1 and mat[vertex,0]>0 and mat[vertex,1]==0:   # A vertex of degree 1 is found 
                current_vertex = vertex                           # If mat[vertex,1]==0 it is necessarilly
                one_found = True                                  # a starting vertex                
                path = []                                         
                while len(candidates)>0 and  mat[current_vertex,0]>0 and mat[current_vertex,1]==0 and \
                not (path_candidates[mat[current_vertex,0]-1]<2  and mat[mat[current_vertex,0]-1,0]>0  \
                and mat[mat[current_vertex,0]-1,1]==0):
   
                                                                      # Starting from it we recover
                    path.append(current_vertex)                       # all vertices in the subpath
                    next_vertex = mat[current_vertex,0]-1             # Connected vertex                   
                    candidates.remove(current_vertex)                 # Once added to the subpath, it is
                    current_vertex=next_vertex                        # removed from candidates

                    
                path.append(current_vertex)                           # The last vertex of this subpath is added
                candidates.remove(current_vertex)
                list_paths.append(path)                               # Subpath is added to list of paths
        to_visit = len(candidates)>0 and one_found                    # We continue extracting subpaths from
                                                                      # candidates while it is not empty
                                                                      # and in the last pass a subpath was found
    return list_paths, candidates    
    
def create_reduced_matrix(mat,n,vertices_to_keep):        
    new_index = -1*np.ones(n).astype(int)       # new_index is initialized having all values -1
    for i in range(len(vertices_to_keep)):      # The vertices to keep are given indices from 0 to n-k
            new_index[vertices_to_keep[i]] = i  # in the same order as they appear in  vertices_to_keep         
    
    reduced_matrix = mat[vertices_to_keep,:]    # The reduced_mat is created using only rows in vertices_to_keep
    for i in range(reduced_matrix.shape[0]):    # Using the new_index vertices in reduced_mat are re-enumerated
        for j in range(reduced_matrix.shape[1]):
            if reduced_matrix[i,j]>0:
                reduced_matrix[i,j] = new_index[reduced_matrix[i,j]-1]+1               

    return new_index, reduced_matrix
    
def deflate_graph(mat,thresh_red_vertices):
    n = mat.shape[0]
    
    # The structures are initialized as empty
    list_fixed_arcs = []
    list_coded_paths = []
    vertices_to_keep = []
    reduced_mat = []
    cycle_vertices = []
        
    # We determine whether the graph is a circuit, whether there is at least one path, and which
    # vertices belong to the paths
    unfeasible,graph_is_circuit, at_least_one_path, n_intermediate_nodes, path_candidates = Find_Path_Candidates(mat,n)

    # For simplification there should be at least one path and the number of nodes that could be 
    # removed should be at least thresh_red_vertices
    if  ((at_least_one_path==True) and (n_intermediate_nodes>=thresh_red_vertices) and (unfeasible==False)):   
        
        # We compute all fixed subpaths and also determine if there are fixed cycles
        list_subpaths, cycle_vertices = FindSubpaths(mat,path_candidates,n)

        # If there is a fixed cycle we halt the simplification
        if len(cycle_vertices)>0: 
            return unfeasible,graph_is_circuit, cycle_vertices, vertices_to_keep, list_fixed_arcs, list_coded_paths, reduced_mat
      
        # The vertices that are kept in the new deflated graph are those that are not intermediate
        # Either they are not in any fixed path (degree 0), 
        # or they are the start and end vertices of a path (degree 1)
        vertices_to_keep = np.where(path_candidates<2)[0]
        
        # We create the new, simplified matrix, by removing intermediate nodes
        # and re-enumerating the remaining nodes (those in vertices_to_keep)
        # Together with the simplified matrix, we get which are the indices of the old
        # vertices in new matrix (new_index). 
        new_index, reduced_mat = create_reduced_matrix(mat,n,vertices_to_keep) 
        
        # We will update the structures needed to reconstruct the original graph
        # In list_fixed_arcs we store each single arc that replaces a subpath. The vertices
        # of the arcs are stored in the indexing corresponding to the new (deflated) graph
        # In list_coded_paths we store each original subpath
        # The paths are stored using the indexing corresponding to the original graph
        # list_fixed_arcs[i] has the arc that represents the subpath in list_coded_paths[i]
        for l in list_subpaths:        # For each path 
            l_size = len(l)         
            if l_size>2:                                     # If the subpath has at least three elements
                first_vertex = new_index[l[0]]               # New indices of first and last vertices in path
                second_vertex = new_index[l[-1]]
                arc = [first_vertex,second_vertex]            # The new fixed arc is formed by the two vertices
                reduced_mat[first_vertex,0] = second_vertex+1 # We add the new arc to the deflated matrix
                list_fixed_arcs.append(arc)                   # We add the fixed arc to the list
                list_coded_paths.append(l)                    # The path is added to the list                
                
    return unfeasible,graph_is_circuit, cycle_vertices, vertices_to_keep, list_fixed_arcs, list_coded_paths,reduced_mat
    


def inflate_circuit(circuit_for_reduced_mat,vertices_to_keep,list_fixed_arcs,list_coded_paths):
    circuit_length = len(circuit_for_reduced_mat)  # Length of the circuit in the reduced graph
    expanded_circuit = []            # We will reconstruct the original circuit here as a sequence
                                     # of subpaths and vertics
                                    
                               
    for i in range(circuit_length):    # We go over the circuit and for each fixed arc, we add the
        if i==circuit_length-1:        # corresponding subpath
            arc = [circuit_for_reduced_mat[i],circuit_for_reduced_mat[0]]  # Possible arc between last and first
        else:
            arc = [circuit_for_reduced_mat[i],circuit_for_reduced_mat[i+1]] # Arcs within the circuit 
        try:
            pos = list_fixed_arcs.index(arc)   # Index of the arc in our list. If the arc does not exist "Value error"
            expanded_circuit.append(list_coded_paths[pos])  # If the arc exists we add the corresponding subpath
        except ValueError:
            pos = -1
            expanded_circuit.append([vertices_to_keep[circuit_for_reduced_mat[i]]]) # If the arc does exist we 
            pass                                                                    # add the corresponding single
                                                                                    # vertex in the original graph
    
    # We remove possible repetitions of the same vertex that maybe the last of a path and first of the next
    circuit = []          
    for subpath in expanded_circuit:
        circuit = circuit+subpath
    final_circuit = [circuit[i] for i in range(len(circuit)) if circuit[i]!=circuit[i-1]]

    return final_circuit
    
    
# Converts the matrix A into the adjacency matrix
    
def Convert_Matrix(mat):
    n=mat.shape[0]
    adj_mat=np.matrix(np.zeros((n,n)))

    for i in range(n):
        for j in range(n):
            if any(c==j+1 for c in list(mat[i])):
                adj_mat[i,j]=1
                
    return adj_mat

    
def Get_HC(hg):
    orig_hg=[[hg[len(hg)-1]+1,hg[0]+1]]

    for i in range(len(hg)-1):
        orig_hg.append([hg[i]+1,hg[i+1]+1])
    orig_hg=sorted(orig_hg) 
    return orig_hg

# Converts an HC in the format of (2,3,4,1) 
# to [[1,2], [2,3], [3,4], [4,1]]

    
def Get_HC_from_array(hg):
    final_hc=[]
    for i in range(len(hg)):
        final_hc.append([i+1,hg[i][0]])
    return final_hc


def Get_HC_from_list(hg):
    final_hg=[hg[0][0]-1,hg[0][1]-1]
    for i in range(1,len(hg)-1):
        final_hg.append(hg[final_hg[i]][1]-1)
    return final_hg
