import numpy as np
import copy

class Solution:
   def __init__(self,x,objectives):
        self.x = x
        self.objectives = objectives
        

"""
 * This is an example implementation of a non-dominated set. It updates the set whenever new
 * solutions are added.
 *
"""
class NonDominatedSet:
    # entries of the non-dominated set
    def __init__(self):
        self.entries = []

    """
     * Add a solution to the non-dominated set
     * @param s The solution to be added.
     * @return true if the solution was indeed added. Otherwise false.
    """    
  
    def add(self,s):
        isAdded = True
        aux_list = copy.copy(self.entries)
        for other in self.entries:
            rel = self.getRelation(s,other)
            # if dominated by or equal in design space
            if (rel == -1 or (rel == 0 and  self.equalsInDesignSpace(s,other))):
                isAdded = False
                break
            elif (rel == 1):
                aux_list.remove(other)
        if (isAdded):
            aux_list.append(s)
        self.entries = aux_list
        return isAdded

        """
     * This is used for non-dominated sorting and returns the dominance relation
     * @param other solution to compare with
     * @return returns 1 if dominates, -1 if dominated and 0 if indifferent
    """
    def getRelation(self,sol,other): 
        val = 0        
        for i in range(len(sol.objectives)):
            #print(i,val,self.objectives[i])
            if (sol.objectives[i] < other.objectives[i]):                
                if (val == -1): return 0
                val = 1
            elif (sol.objectives[i] > other.objectives[i]):
                if (val == 1): return 0
                val = -1            
        return val

    
    """
     * @param other solution to compare with
     * @return True if tour and packing plan is equal
    """
    def equalsInDesignSpace(self,sol,other):
        return np.array_equal(sol.x,other.x) 


