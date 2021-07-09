import numpy as np


#######################################################
#######################################################
################# HGRAPH CLASS ########################
#######################################################
#######################################################

# This class contains information of the graph.

class HGraph:
    def __init__(self,A):
        self.A = A       # Original graph as a matrix
        self.N  = self.A.shape[0]  # Number of vertices
        self.N_edges = self.A.shape[1] # Number of maximum degree
        
    
    def RowSums(self): ### How many vertices are visited from each vertex i
        self.row_sums=[]
        for i in range(self.N): 
           for a in range(self.N_edges): 
               if self.A[i,a] == 0:
                   self.row_sums.append(a)
                   break     
               if a ==self.N_edges-1:
                   self.row_sums.append(a+1)      
    
    def StartNodes(self): # Vertices that emanate from vertex i
        self.RowSums()
        self.start_nodes = [i for i in range(self.N) for a in range(self.row_sums[i])] 
        self.start_nodes = np.array(self.start_nodes)

     
    def GoInto(self,fixed_arcs): # Vertices that are reached from fixed arcs
        self.index=[]
        for i in range (len(fixed_arcs)):
            self.index.append(fixed_arcs[i][1])    
    
    def FixedArcs(self): # Fixes arcs from the A matrix
        self.RowSums()
        self.new_fixed_arcs=[]   
        self.a=np.where(np.array(self.row_sums)==1)[0]
        for i in range(len(self.a)):
            self.new_fixed_arcs.append([self.a[i]+1,int(self.A[self.a[i],0])])
