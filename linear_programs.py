import numpy as np
from constraints_class import Constraints



##################################################
##################################################
######## IMPLEMENTS THE LINEAR PROGRAMS ##########
##################################################
##################################################




# This function adds constraints 2.8, 2.10, 2.14, and 2.13
# to the object of constraints. It is called before using
# the Second LP algorithm and First LP algorithm
    
    
    
## Constraints numeration comes from 
## M. Haythorpe. Markov Chain Based Algorithms for the Hamiltonian Cycle Problem.
## PhD dissertation. University of South Australia

def Add_Relevant_Constraint(graph,beta,fixed_arcs):
    
    constraints=Constraints(graph,beta)
    
    constraints.AeqConstraint8()
    constraints.beqConstraint8()
    
    constraints.AeqConstraint9()
    constraints.beqConstraint9()
        
    Aeq_8=constraints.Aeq_constraint8  # Add constraint 2.8 here
    beq_8=constraints.beq_constraint8
        
    Aeq_9=constraints.Aeq_constraint9  # Add constraint 2.9 here
    beq_9=constraints.beq_constraint9

    # For each arc in  fixed_arcs we add a constraint
        
    additional_Aeq=[]
    additional_beq=[]        
        
    for arc in fixed_arcs: 

        if arc[1]==1:     # If arc goes to node 1
            constraints.AeqConstraint14(arc)
            additional_Aeq.append(constraints.Aeq_constraint14)  # Add constraint 2.14
            additional_beq.append(constraints.beq_constraint14)
        else: 
            constraints.AeqConstraint13(arc)
            additional_Aeq.append(constraints.Aeq_constraint13)  # Add constraint 2.13
            additional_beq.append(constraints.beq_constraint13)
    return constraints,Aeq_8,beq_8,Aeq_9,beq_9, additional_Aeq, additional_beq
        
       
# Implements the Second LP Algorithm 
       
def Second_LP_Algorithm(graph,beta,fixed_arcs,solver): 
        
    # Adds constraints 2.8 and 2.14 or 2.13
    constraints,Aeq_8,beq_8,Aeq_9,beq_9, additional_Aeq, additional_beq=Add_Relevant_Constraint(graph,beta,fixed_arcs)        

    if len(additional_Aeq)>1:
        additional_Aeq=np.vstack(additional_Aeq)

    if len(fixed_arcs)==0:                                
        Aeq=Aeq_8     
        beq=beq_8     
    else:                                 
        Aeq=constraints.AddConstraints(Aeq_8,additional_Aeq)
        beq=constraints.AddConstraints(beq_8,additional_beq)
    constraints.FeasibilityF()    
    if solver=='cplex': # Call to the solver
        import cplex_solver
        status,x,function_value = cplex_solver.CPLEXPY(constraints.new_f,constraints.zero_row,Aeq,beq)
    elif solver=='gurobi':
        import gurobi
        status,x,function_value=gurobi.GUROBIPY(constraints.new_f,Aeq,beq)
    return status, x, function_value


# Implements the First LP Algorithm 
        
def First_LP_Algorithm(graph,beta,fixed_arcs,solver):
    
    # Adds constraints 2.8 and 2.9
    constraints,Aeq_8,beq_8,Aeq_9,beq_9, additional_Aeq, additional_beq=Add_Relevant_Constraint(graph,beta,fixed_arcs)  
    Aeq_new=constraints.AddConstraints(Aeq_8,[Aeq_9])
    beq_new=constraints.AddConstraints(beq_8,[beq_9])
    constraints.F(fixed_arcs)
    if solver=='cplex': # Call to the solver
        import cplex_solver
        status,x,function_value = cplex_solver.CPLEXPY(constraints.f,constraints.zero_row,Aeq_new,beq_new)
    elif solver=='gurobi':
        import gurobi
        status,x,function_value=gurobi.GUROBIPY(constraints.f,Aeq_new,beq_new) 
    return status,x, function_value
