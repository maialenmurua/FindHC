import numpy as np    
import operator
import aux_branching6


####################################################
####################################################
## IMPLEMENTS THE SIX DIFFERENT BRANCHING METHODS ##
####################################################
####################################################



# Logical order of visiting arcs

def Branching_Method1(A,splitting_node): 
    arcs=[]    
    candidates=A[splitting_node][A[splitting_node]!=0] 
    for i in range(len(candidates)):
        arcs.append([splitting_node+1,int(candidates[i])])
    return arcs



# Method 2: 1st (i,j) 2nd (i,k) then the rest of the arcs in logical order
# It is assumed that xij<=xik 


def Branching_Method2(A,splitting_node,x,hg): 
    
    # The nodes that is visiting the splitting node related to x                                                 
    branching=[[item for item in hg if item[0]==splitting_node+1][0][1]-1,
              [item for item in hg if item[0]==splitting_node+1][1][1]-1]  
    x=np.array(x)
    vertex1=x[x>0][splitting_node] # The entries of x related to the arcs of branching
    vertex2=x[x>0][splitting_node+1]
    arc_1=[splitting_node+1,branching[0]+1] # The two arcs that emanate from the splitting node in x
    arc_2=[splitting_node+1,branching[1]+1]
    if vertex1<=vertex2: # The order depends on the value that have on x. 
        arcs=[arc_1,arc_2]
    else:
        arcs=[arc_2,arc_1]
    candidates=A[splitting_node][(A[splitting_node]!=branching[0]+1)& # the other arcs are added in logical order
    (A[splitting_node]!=branching[1]+1)&(A[splitting_node]!=0)]
    for i in range(len(candidates)):
        arcs.append([splitting_node+1,candidates[i]])
    return arcs # Returns a list with the branching candidates
    

# Method 3: 1st (i,k), 2nd (i,j) then the rest of the arcs in logical order
# It is assumed that xij<=xik 

def Branching_Method3(A,splitting_node,x,hg): 
     
    #The nodes that is visiting the splitting node related to x                 
    branching=[[item for item in hg if item[0]==splitting_node+1][0][1]-1,
              [item for item in hg if item[0]==splitting_node+1][1][1]-1]  
    x=np.array(x)
    vertex1=x[x>0][splitting_node] # The entries of x related to the arcs of branching
    vertex2=x[x>0][splitting_node+1]
    arc_1=[splitting_node+1,branching[0]+1] # The two arcs that emanate from the splitting node in x
    arc_2=[splitting_node+1,branching[1]+1]
    if vertex1<=vertex2: # Is the same as Method 2 but changing the order of (i,j) and (i,k)
        arcs=[arc_2,arc_1]
    else:
        arcs=[arc_1,arc_2]
    candidates=A[splitting_node][(A[splitting_node]!=branching[0]+1)&
    (A[splitting_node]!=branching[1]+1)&(A[splitting_node]!=0)] # the other arcs are added in logical order
    for i in range(len(candidates)):
        arcs.append([splitting_node+1,candidates[i]])
    return arcs # Returns a list with the branching candidates


# Method 4: 1st rest of the arcs in logical order and then (i,j), (i,k)
# (i,k) rest of the arcs in logical order and (i,j)


def Branching_Method4(A,splitting_node,x,hg): 
    
    #The nodes that is visiting the splitting node related to x
    branching=[[item for item in hg if item[0]==splitting_node+1][0][1]-1,  
              [item for item in  hg if item[0]==splitting_node+1][1][1]-1]     
    x=np.array(x)
    vertex1=x[x>0][splitting_node] # The entries of x related to the arcs of branching
    vertex2=x[x>0][splitting_node+1]
    arcs=[]
    candidates=A[splitting_node][(A[splitting_node]!=branching[0]+1)& # the other arcs are added in logical order
        (A[splitting_node]!=branching[1]+1)&(A[splitting_node]!=0)]
    for i in range(len(candidates)):
        arcs.append([splitting_node+1,candidates[i]])
    arc_1=[splitting_node+1,branching[0]+1]
    arc_2=[splitting_node+1,branching[1]+1] # The two arcs that emanate from the splitting node in x
    if vertex1<=vertex2:
        arcs.append(arc_1) # The order depends on the value that have on x. 
        arcs.append(arc_2)
    else:
        arcs.append(arc_2)
        arcs.append(arc_1)
    return arcs # Returns a list with the branching candidates
    

# Method 5: 1st (i,k) then the rest of the arcs in logical order and finally (i,j)
# (i,k) rest of the arcs in logical order and (i,j)

    
def Branching_Method5(A,splitting_node,x,hg): 

    #The nodes that is visiting the splitting node related to x
    branching=[[item for item in hg if item[0]==splitting_node+1][0][1]-1,
    [item for item in hg if item[0]==splitting_node+1][1][1]-1]
    x=np.array(x)
    vertex1=x[x>0][splitting_node] # The entries of x related to the arcs of branching
    vertex2=x[x>0][splitting_node+1]
    if vertex1<=vertex2: # The order depends on the value that have on x, (i,k)
        arcs=[[splitting_node+1,branching[0]+1]]
    else:
        arcs=[[splitting_node+1,branching[1]+1]]
    candidates=A[splitting_node][(A[splitting_node]!=branching[0]+1)&  # the other arcs are added in logical order
    (A[splitting_node]!=branching[1]+1)&(A[splitting_node]!=0)]
    for i in range(len(candidates)):
        arcs.append([splitting_node+1,candidates[i]])
    if vertex1<=vertex2: # The order depends on the value that have on x, (i,j)
        arcs.append([splitting_node+1,branching[1]+1])
    else:
        arcs.append([splitting_node+1,branching[0]+1])
    return arcs # Returns a list with the branching candidates


# Method 6: Global branching method, selects the arc that fixes most arcs. 
    
def Branching_Method6(A,splitting_node,fixed_arcs): 
    count=[]
    candidates=A[splitting_node][A[splitting_node]!=0]
    for i in range(len(candidates)):
        count.append(aux_branching6.Arcs_Fixing(A,[splitting_node+1,candidates[i]],fixed_arcs)) 
    # Count how many arcs can be fixed per each candidate
    d=dict(zip(candidates,count))
    sorted_d = sorted(d.items(), key=operator.itemgetter(1),reverse=True)
    arcs=[]
    for i in range(len(sorted_d)):
        arcs.append([splitting_node+1,sorted_d[i][0]])
    return arcs # Returns a list with the branching candidates

