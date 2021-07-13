import gurobipy as gp
from gurobipy import GRB



#############################################################
#############################################################
########### SOLVE THE LP WITH GUROBI SOLVER ##################
#############################################################
#############################################################



# Using the components of the linear program 
# creates all the data structure to call GUROBI
    
    
def GUROBIPY(f,Aeq,beq):
    m = gp.Model("lp")
    x = m.addMVar(shape=Aeq.shape[1], vtype=GRB.CONTINUOUS, name="x")
    m.setObjective(f @ x, GRB.MINIMIZE)
    m.addConstr(Aeq @ x == beq, name="c")
    m.optimize()
    if m.status!=2: # It is not feasible
        return 2,[],0
    else:
        return 1,list(x.X),m.objVal # Get the solution
