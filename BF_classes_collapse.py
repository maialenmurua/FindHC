
import numpy as np
import copy
import time
import simplification as simp
from graph_class import HGraph
import branching_methods
import aux_bf
import degree_based_simp
import linear_programs



###############################################
###############################################
######## BRANCH-AND-FIX-COLLAPSE ##############
###############################################
###############################################



# Global variables
start = time.time()
calls=1     


 


# The BF_Algorithms class will implement the Branch and Fix algorithm and
# all the auxiliary algorithms required
# It needs one object of type BF_Constraints and one object of type BF_HPGraph


class BF_Algorithms:
    def __init__(self,A,beta,lev,branching_method,ind_degree2,solver):
        self.graph=HGraph(A)
        self.N = self.graph.N
        self.mu=1/(2*self.N) 
        self.A=self.graph.A
        self.graph.RowSums()
        self.row_sums=self.graph.row_sums
        self.graph.StartNodes()
        self.start_nodes=self.graph.start_nodes
        self.beta =beta
        self.Set_Bound()
        self.tol=1e-010  # Tolerance
        self.fixed_arcs=[]
        self.A=np.array(A)
        self.level=lev
        self.branching_method=int(branching_method)
        self.ind_degree2=ind_degree2
        self.solver=solver


    # Set the bound used to test Hamiltonicity
    def Set_Bound(self):
        self.bound = (1-(self.N-1)*self.mu)*(1-self.beta)+self.mu*(self.beta-self.beta**self.N)
        self.bound = self.bound/((1-self.beta)*(1-self.beta**self.N))

       
    def Arcs_In_Vector(self): 
        count=0
        self.d={}
        for i in range(self.N):
            for a in range(self.row_sums[i]):
                self.d.update({count:[i+1,self.A[i,a]]})
                count=count+1
    
    # Get an HC from a given solution vector x             
    def Get_Sol(self,x): 
        self.Arcs_In_Vector()
        self.hg=[]
        higher_tol = list(np.where(np.array(x)>0)[0])
        for i in range(len(higher_tol)):
            self.hg.append(self.d[higher_tol[i]])
            
   
    # Eliminate arcs in the current graph (matrix) given a fixed arc
    def Update_Matrix(self,fixed_arc):
        self.graph.RowSums()
        self.A[fixed_arc[0]-1]=0
        candidates=np.where(self.A==fixed_arc[1])
        for i in range(len(candidates[0])):
            if self.graph.row_sums[candidates[0][i]]>1:
                self.A[candidates[0][i],candidates[1][i]]=0               
        self.A[fixed_arc[0]-1,0]=fixed_arc[1]
        if len(np.where(self.A[fixed_arc[1]-1]==fixed_arc[0])[0])>0:
            if len(np.nonzero(self.A[fixed_arc[1]-1])[0])>1:
                self.A[fixed_arc[1]-1,np.where(self.A[fixed_arc[1]-1]==fixed_arc[0])[0]]=0
        for i in range(self.A.shape[0]):
            self.A[i]= aux_bf.pushZerosToEnd(self.A[i])
        idx = np.argwhere(np.all(self.A[..., :] == 0, axis=0))
        self.A = np.delete(self.A, idx, axis=1)
                
    
    # Fixed additional arcs from the updated matrix
    def Update_arcs(self):
        self.graph=HGraph(self.A)
        self.graph.FixedArcs()
        self.fixed_arcs.extend(x for x in self.graph.new_fixed_arcs if x not in self.fixed_arcs)
        
    
    #Update the constraints based on the updated graph  
    def Update_constraints(self):  
       self.A=self.graph.A
       self.graph.RowSums()
       self.row_sums=self.graph.row_sums
       
                              
    def Branch_and_Fix(self):
        
        global calls

        fixed_arcs=self.fixed_arcs # save the actual graph's fixed arcs
        
        if self.ind_degree2=='Yes':
            adj_mat,feasible_solution=degree_based_simp.Degree2_Vertices(self.A) # Update the adjacency matrix when possible 
            # Fixing arcs that are connected to two arcs of degree 2		
            if feasible_solution==0: # This means it is not possible to find a HC with the current graph
                return -1,[],self.level,calls
         
            else:
                self.__init__(adj_mat,self.beta,self.level,self.branching_method,self.ind_degree2,self.solver) # Redefined the graph with the new matrix

            self.fixed_arcs=fixed_arcs  #
        
        elif self.ind_degree2=='No':
            self.__init__(self.A,self.beta,self.level,self.branching_method,self.ind_degree2,self.solver)
        
        status, second_x, function_value =  linear_programs.Second_LP_Algorithm(self.graph,self.beta,self.fixed_arcs,self.solver) # Check feasibility with Second LP
   
        if status!=1 or function_value-self.bound>self.tol:
            return -1,[],self.level,calls      # Fathom the current branch

        if self.level!=0: # Inspect if there are subcycles when last arc was fixed
            self.Arcs_In_Vector()
            candidate_loop_vertices =  aux_bf.Inspect_Potencial_Subcycles(second_x,self.mu,self.beta,self.tol,self.d)          

            if len(candidate_loop_vertices)>0:
                if set(tuple(x) for x in candidate_loop_vertices).issubset(tuple(x) for x in self.fixed_arcs)==True:
                  
                    return -1,[] ,self.level,calls  # Fathom the current branch
               
        status, x, function_value = linear_programs.First_LP_Algorithm(self.graph,self.beta,self.fixed_arcs,self.solver) 
	# Obtain the solution with the First LP

       
        if status!=1: 
            return -1,[],self.level,calls        # Fathom the current branch
            
        splitting_node=  aux_bf.Identify_Splitting_Node(x,self.tol,self.start_nodes)   # Find splitting node

        self.Get_Sol(x)
    
        if splitting_node==-1: # An HC was found return code=1   

            return 1, self.hg,self.level,calls
            
        else:
             
             if self.branching_method==1: # Apply the corresponding branching method
                 mylist=branching_methods.Branching_Method1(self.A,splitting_node)
             elif self.branching_method==2: 
                 mylist=branching_methods.Branching_Method2(self.A,splitting_node,x,self.hg)
             elif self.branching_method==3: 
                 mylist=branching_methods.Branching_Method3(self.A,splitting_node,x,self.hg)
             elif self.branching_method==4: 
                 mylist=branching_methods.Branching_Method4(self.A,splitting_node,x,self.hg)
             elif self.branching_method==5:
                 mylist=branching_methods.Branching_Method5(self.A,splitting_node,x,self.hg)
             elif self.branching_method==6:
                 mylist=branching_methods.Branching_Method6(self.A,splitting_node,self.fixed_arcs)
                
             found=False
             i=0
             while (found==False) & (i<len(mylist)):
                
                fixed_arc=[splitting_node+1,int(mylist[i][1])] # Fix an arc
           

            # Calls the Branch and Fix algorithm recursively. To do that
            # we create a new Branch_and_Fix Object with new objects of graph
            # and constraints and apply the Branch and Fix algorithm on it 

                NewBF = BF_Algorithms(self.A, self.beta,self.level+1,self.branching_method,self.ind_degree2,self.solver)
               
                NewBF.Update_Matrix(fixed_arc) ### Update the graph, fix more arcs if possible
                NewBF.Update_arcs()
                initial=copy.copy(NewBF.fixed_arcs)
              
                if len(NewBF.fixed_arcs)==0: # Check whether new arcs were fixed or not
                    last_fixed_arcs=[]
                else:
                    last_fixed_arcs=list(set(tuple(row) for row in NewBF.fixed_arcs )-set(tuple(row) for row in self.fixed_arcs)-set(fixed_arc))  
             
    
                while len(last_fixed_arcs)>0: # If new arcs were fixed, eliminate and fixed arcs until more arcs cannot be fixed. 
                    for arc in last_fixed_arcs:
                        NewBF.Update_Matrix(arc)
                    NewBF.Update_arcs()
                    current=NewBF.fixed_arcs
        
                    if initial==current:
                        last_fixed_arcs=[]
                    else: 
                        last_fixed_arcs=list(set(tuple(row) for row in current)-set(tuple(row) for row in initial))  
                    initial=copy.copy(current)

                NewBF.Update_constraints() # Update new graphs's parameters 
                
                unfeasible,graph_is_circuit, cycle_vertices, vertices_to_keep,list_fixed_arcs, list_coded_paths,reduced_mat=simp.deflate_graph(NewBF.A,1)

                # Deflate the graph  
                
                if graph_is_circuit==True:
                    final_hg=simp.Get_HC_from_array(NewBF.A)

                    return 1, final_hg,self.level,calls
               
                elif (len(cycle_vertices)>0) | (unfeasible==True): # A subcycle was found is unfeasible
               
                    code=-1
                 
                elif len(reduced_mat)>0:      # The matrix was reduced
                                             # We will find the HCP in the reduced matrix and inflate it later 
                
                    reduced_BF = BF_Algorithms(reduced_mat,self.beta,self.level+1,self.branching_method,self.ind_degree2,self.solver)
                    code,reduced_hg,level,calls = reduced_BF.Branch_and_Fix() # Call the algorithm recursively
                    
                    if code==1:

                        reduced_hg2=simp.Get_HC_from_list(reduced_hg)
                        aux_hg=simp.inflate_circuit(reduced_hg2,vertices_to_keep,list_fixed_arcs,list_coded_paths)# Inflate the graph
                        
                        hg=simp.Get_HC(aux_hg) # Convert to the suitable format
                
 
                    else:
                        pass
                    
                else:
                    code,hg,level,calls = NewBF.Branch_and_Fix() # Call the algorithm recursively
                
                calls+=1
                found=(code==1)
                i=i+1

             if found==True:
                 return 1, hg,self.level,calls
             else:
                 return -1, [],self.level,calls





      

