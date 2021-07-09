import numpy as np

###########################################################
###########################################################
############## AUXILIAR FUNCTIONS BF ######################
###########################################################
###########################################################


# Given a row of a matrix, the zeros are moved to the end of the row.

def pushZerosToEnd(arr): 
    n=len(arr)
    count = 0 
    for i in range(n): 
        if arr[i] != 0:               
            arr[count] = arr[i] 
            count+=1
    while count < n: 
        arr[count] = 0
        count += 1
    return arr
    
       
# Given a solution vector x, identifies which node
# has two non-zero entries
       
def Identify_Splitting_Node(x,tol,start_nodes):
    higher_tol = np.where(np.array(x)>tol)[0] 
    nodes_in_sol = [start_nodes[i] for i in higher_tol]
    unique_nodes, counts = np.unique(nodes_in_sol, return_counts=True)
    splitting_node_position = np.where(counts==2)[0]	
    if splitting_node_position !=None:
        splitting_node = unique_nodes[splitting_node_position][0]
    else:
        splitting_node=-1
    return splitting_node

# Given a solution vector x, identifies arcs that may create a subcycle    
       
def Inspect_Potencial_Subcycles(x,mu,beta,tol,d):
    # Check whether there are repeated entries that take value mu/(1-beta)
    val=mu/(1-beta)
    count=np.where((np.array(x)>tol)&(abs(np.array(x)-val)<tol))[0]    
    candidate_loop_vertices=[]          
    if len(count)>1: # Identify the arcs related with those variables
        for i in range(len(count)):
            candidate_loop_vertices.append(d[count[i]])
    return candidate_loop_vertices
