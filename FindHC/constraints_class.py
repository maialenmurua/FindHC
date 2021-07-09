import numpy as np


#############################################
#############################################
############ CONSTRAINTS CLASS ##############
#############################################
#############################################


# It implements the structures need to build the constraints of the linear programs

class Constraints:
    def __init__(self,graph,beta):
        self.graph=graph  # Attributes of the HGraph class
        graph.RowSums()
        self.A=graph.A
        self.N =graph.N
        self.row_sums=graph.row_sums
        self.beta=beta # beta parameter
        self.mu=1/(2*self.N) # mu parameter
        self.number_constraints = np.sum(self.row_sums)  # Empty lists to initialize the variables
        self.zero_row=np.zeros(self.number_constraints)
        self.f=np.zeros(self.number_constraints)                      
        self.Aeq_constraint8=np.zeros((self.N,self.number_constraints))
        self.Aeq_constraint9=np.zeros(self.number_constraints)
        self.new_f=np.zeros(self.number_constraints)
        self.Aeq_constraint14=np.zeros(self.number_constraints)
        self.beq_constraint13=self.mu # Right side of constraint (2.13)
        self.beq_constraint14=self.mu*(self.beta-self.beta**self.N)/(1-self.beta) # Right side of constraint (2.14)
       
        
    def F(self,fixed_arcs): # Objective function 
        self.graph.GoInto(fixed_arcs)
        gointo=self.graph.index
        
        for i in range(len(fixed_arcs)):
            if gointo[i]!=1:                                        
                index_j=np.sum(self.row_sums[0:fixed_arcs[i][1]-1])
                for a in range(self.row_sums[fixed_arcs[i][1]-1]):
                    self.f[index_j+a]=self.f[index_j+a]+1
        
                index_i=int(np.sum(self.row_sums[0:fixed_arcs[i][0]-1])) 
                for b in range(self.row_sums[fixed_arcs[i][0]-1]):     
                    self.f[index_i+b]=self.f[index_i+b]-self.beta
            else:
                for a in range(self.row_sums[0]):
                    self.f[a]= self.f[a]-self.beta**self.N
                index=np.sum(self.row_sums[0:fixed_arcs[i][0]-1])
                self.f[index]=self.f[index]+self.beta


## Constraints numeration comes from 
## M. Haythorpe. Markov Chain Based Algorithms for the Hamiltonian Cycle Problem.
## PhD dissertation. University of South Australia


                                     
    def AeqConstraint8(self): # Left side of constraint (2.8) 
       for j in range(self.N):                                
           count = 1
           for i in range(self.N):
               for a in range(self.row_sums[i]):
                   self.Aeq_constraint8[j,count-1]=(int(i==j) - self.beta*int(self.A[i,a]==j+1))
                   count = count + 1
    
    def AeqConstraint9(self): # Left side of constraint (2.9)
        for a in range(self.row_sums[0]): 
            self.Aeq_constraint9[a] = 1

       
    def beqConstraint8(self): # Right side of constraint (2.8)
        self.beq_constraint8 = np.concatenate(([1-(self.N-1)*self.mu], 
                               np.linspace(self.mu,self.mu,self.N-1)))  

        
    def beqConstraint9(self): # Right side of constraint (2.9)
          self.beq_constraint9=((1-(self.N-1)*self.mu)*(1-self.beta) + 
          self.mu*(self.beta-self.beta**self.N))/((1-self.beta)*(1-self.beta**self.N))

    
    def FeasibilityF(self): # Objective function for the extra feasibility check
        for a in range(self.row_sums[0]): 
            self.new_f[a] = 1
  
    
    def AeqConstraint13(self,arc): # Left side of constraint (2.13) 
        self.Aeq_constraint13=np.zeros(int(self.number_constraints))
        index_ik1=np.sum(self.row_sums[0:arc[1]-1])   
        for a in range(self.row_sums[arc[1]-1]):
            self.Aeq_constraint13[index_ik1+a]=1       
        ik1_locate=np.where(self.A[arc[0]-1]==arc[1])[0][0]
        index_ik=np.sum(self.row_sums[0:arc[0]-1])
        self.Aeq_constraint13[int(index_ik+ik1_locate)]=-self.beta
             
    
    def AeqConstraint14(self,arc): # Left side of constraint (2.14)
        for a in range(self.row_sums[0]):
            self.Aeq_constraint14[a]= -self.beta**self.N
        index=np.sum(self.row_sums[0:arc[0]-1])
        self.Aeq_constraint14[index]=self.beta
    
      
    def AddConstraints(self,constraint1,constraint2): # Join two constraints
        return np.append(constraint1,constraint2,axis=0)
