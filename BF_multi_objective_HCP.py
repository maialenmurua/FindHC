
import numpy as np
import copy
import pareto_set as pareto_set
import linear_programs
import time
import aux_bf as aux_bf
import graph_class as graph
import aux_bf_mo as aux_bf_mo
import degree_based_simp as degree_based_simp
import branching_methods as branching_methods




###############################################
###############################################
######## MULTI-OBJECTIVE BRANCH-AND-FIX #######
###############################################
###############################################

start = time.time()
PERIOD_OF_TIME = 1
time_counter=0        

list_hg=[]
PS=pareto_set.NonDominatedSet()


# The BF_Algorithms class will implement the Multi-objective Branch-and-Fix algorithm and
# all the auxiliary algorithms required
# It needs one object of type BF_Constraints and one object of type BF_HPGraph


class BF_Algorithms:
    def __init__(self,A,beta,lev,branching_method,c1_file,c2_file,ind_degree2,solver,file_name,perm):
        self.graph=graph.HGraph(A)
        self.c1_file=c1_file
        self.c2_file=c2_file
        self.c1,self.c2=aux_bf_mo.CostMatrix(c1_file,c2_file)
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
        self.branching_method=branching_method
        self.ind_degree2=ind_degree2
        self.solver=solver
        self.file_name=file_name
        self.perm=perm

        
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
        higher_tol = list(np.where(np.array(x)>self.tol)[0])
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
        self.graph=graph.HGraph(self.A)
        self.graph.FixedArcs()
        self.fixed_arcs.extend(x for x in self.graph.new_fixed_arcs if x not in self.fixed_arcs)
        
    
    # Update the constraints based on the updated graph  
    def Update_constraints(self):
       self.A=self.graph.A
       self.graph.RowSums()
       self.row_sums=self.graph.row_sums  
       self.graph.StartNodes()
       self.start_nodes=self.graph.start_nodes               
              
    # Get partial solutions        
    def Partial_Solutions(self,time_counter):
        if self.perm=='No':
            with open('partial_output_{0}_{1}.txt'.format(self.file_name,time_counter), 'w') as f:
                print('Number of HCs:',len(list_hg),file=f)
                for i in range(len(PS.entries)):
                    print('Pareto Set:',PS.entries[i].objectives,file=f)
        else:
            with open('partial_output_{0}_{1}_{2}.txt'.format(self.file_name,self.perm,time_counter), 'w') as f:
                print('Number of HCs:',len(list_hg),file=f)
                for i in range(len(PS.entries)):
                    print('Pareto Set:',PS.entries[i].objectives,file=f)
            
            
    def Branch_and_Fix(self):
        global time_counter  
        global PERIOD_OF_TIME
        while True:
            fixed_arcs=self.fixed_arcs # save the actual graph's fixed arcs
            if self.ind_degree2=='Yes':
                adj_mat,feasible_solution=degree_based_simp.Degree2_Vertices(self.A) # Update the adjacency matrix when possible 
            # Fixing arcs that are connected to two arcs of degree 2		
                if feasible_solution==0: # This means it is not possible to find a HC with the current graph
                    return -1,[]
         
                else:
                    self.__init__(adj_mat,self.beta,self.level,self.branching_method,self.c1_file,self.c2_file,self.ind_degree2,self.solver,self.file_name,self.perm) 
                # Redefined the graph with the new matrix
                self.fixed_arcs=fixed_arcs  #
            elif self.ind_degree2=='No':
                self.__init__(self.A,self.beta,self.level,self.branching_method,self.c1_file,self.c2_file, self.ind_degree2,self.solver,self.file_name,self.perm)
                self.fixed_arcs=fixed_arcs  

            status, second_x, function_value =  linear_programs.Second_LP_Algorithm(self.graph,self.beta,self.fixed_arcs,self.solver)

            
            if status!=1 or function_value-self.bound>self.tol:
                return -1,[]      # This means that no HC was found

            if self.level!=0:
                self.Arcs_In_Vector()
                candidate_loop_vertices = aux_bf.Inspect_Potencial_Subcycles(second_x,self.mu,self.beta,self.tol,self.d)          
                if len(candidate_loop_vertices)>0:
                    if set(tuple(x) for x in candidate_loop_vertices).issubset(tuple(x) for x in self.fixed_arcs)==True:
                        return -1,[]
       
            if (self.fixed_arcs!=[])&(list_hg!=[]):
                val=aux_bf_mo.Compute_Fitness_Arcs(self.fixed_arcs,self.c1,self.c2)
                for i in range(len(PS.entries)):
                    if (val[0]>=PS.entries[i].objectives[0])&(val[1]>=PS.entries[i].objectives[1]):
                        return -1,[]
            
           
            status, x, function_value = linear_programs.First_LP_Algorithm(self.graph,self.beta,self.fixed_arcs,self.solver) 
            
            if status!=1:
                return -1,[]       # This means that no HC was found
            
            splitting_node= aux_bf.Identify_Splitting_Node(x,self.tol,self.start_nodes)   # Find splitting node
            self.Get_Sol(x)

            if splitting_node==-1: 
                 list_hg.append(self.hg)

                 sol=pareto_set.Solution(np.array([len(list_hg)]),np.array(aux_bf_mo.Compute_Fitness_Arcs(self.hg,self.c1,self.c2)))
                 PS.add(sol)                                  
                 return 1, []
                                                    
            else:
             
                if self.branching_method==1:
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
            
                i=0
                while (i<len(mylist)):
                
                    fixed_arc=[splitting_node+1,int(mylist[i][1])]

            # Calls the Branch and Fix algorithm recursively. To do that
            # we create a new Branch_and_Fix Object with new objects of graph
            # and constraints and apply the Branch and Fix algorithm on it 

                    NewBF = BF_Algorithms(self.A, self.beta,self.level+1,self.branching_method,self.c1_file,self.c2_file,
                                          self.ind_degree2,self.solver,self.file_name,self.perm)
               
                    NewBF.Update_Matrix(fixed_arc) ### Update the graph
                    NewBF.Update_arcs()
                    initial=copy.copy(NewBF.fixed_arcs)

                    if len(NewBF.fixed_arcs)==0:
                        last_fixed_arcs=[]
                    else:
                        last_fixed_arcs=list(set(tuple(row) for row in NewBF.fixed_arcs )-set(tuple(row) for row in self.fixed_arcs)-set(fixed_arc))  
             
                
                    while len(last_fixed_arcs)>0:
                        for arc in last_fixed_arcs:
                            NewBF.Update_Matrix(arc)
                        NewBF.Update_arcs()
                        current=NewBF.fixed_arcs

                        if initial==current:
                            last_fixed_arcs=[]
                        else: 
                            last_fixed_arcs=list(set(tuple(row) for row in current)-set(tuple(row) for row in initial))  
                        initial=copy.copy(current)
                    NewBF.Update_constraints() # Update new object's parameters  
                
                    code,hg = NewBF.Branch_and_Fix()
                    
                    i=i+1
                    if (time.time() > start + (PERIOD_OF_TIME/6))&(time_counter==0):
                        time_counter=time_counter+1
                        self.Partial_Solutions(time_counter)
                    elif (time.time() > start + 2*(PERIOD_OF_TIME/6))&(time_counter==1):
                        time_counter=time_counter+1
                        self.Partial_Solutions(time_counter)
                    elif (time.time() > start + 3*(PERIOD_OF_TIME/6))&(time_counter==2):
                        time_counter=time_counter+1
                        self.Partial_Solutions(time_counter)
                    elif (time.time() > start + 4*(PERIOD_OF_TIME/6))&(time_counter==3):
                        time_counter=time_counter+1
                        self.Partial_Solutions(time_counter)
                    elif (time.time() > start + 5*(PERIOD_OF_TIME/6))&(time_counter==4):
                        time_counter=time_counter+1
                        self.Partial_Solutions(time_counter)
                    elif (time.time() > start + PERIOD_OF_TIME)&(time_counter==5): 
                        time_counter=time_counter+1
                        self.Partial_Solutions(time_counter)
                    elif (time.time() > start + PERIOD_OF_TIME)&(time_counter==6):

                        break

                return 0, []
            



