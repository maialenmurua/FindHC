
import sys
sys.path.append('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/python/3.5/x86-64_linux')
import cplex
import numpy as np   
   

#############################################################
#############################################################
########### SOLVE THE LP WITH CPLEX SOLVER ##################
#############################################################
#############################################################



# Using the components of the linear program 
# creates all the data structure to call CPLEX
   
def CPLEXPY(f,zero_row,Aeq,beq):
    problem = cplex.Cplex() # Create an instance of a linear problem to solve 
    problem.objective.set_sense(problem.objective.sense.minimize) # minimum of the objective function      
    objective=f.tolist() # f is the objective function
    lower_bounds=zero_row.tolist() # Lower bounds are zeros 
    upper_bounds=np.repeat(np.inf,len(zero_row)).tolist() # There are no upper bounds
    problem.variables.add(obj = objective, lb = lower_bounds,ub=upper_bounds)
    Aeq=Aeq.tolist()
    index=np.arange(len(objective)).tolist()
    constraints = []
    for i in range(len(Aeq)):
        constraints.append([index,Aeq[i]]) ### build constraints with variable indices (left size)
    rhs=beq.tolist()  ### right side of the constraints
    constraint_senses = np.repeat('E',len(constraints)).tolist() ### equality constraints
    problem.linear_constraints.add(lin_expr = constraints,
                               senses = constraint_senses, # define the constraints
                               rhs = rhs)
    problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal) # we are solving the primal
    problem.parameters.simplex.tolerances.feasibility.set(1e-9) # the permitted minimum tolerance
    problem.parameters.preprocessing.aggregator.set(0)
    problem.parameters.preprocessing.reduce.set(0)
    problem.solve() # Solve the problem
    status = problem.solution.get_status() # get the status
    if status!=1: 
        return status, [], 0 # It is not feasible
    else:
        x=problem.solution.get_values() # Get the solution
        function_value=problem.solution.get_objective_value()
        return status,x,function_value
