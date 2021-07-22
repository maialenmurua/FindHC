import numpy as np    
import pandas as pd


        
def Compute_Fitness_Arcs(arcs,c1,c2):
    obj_func1=[]
    obj_func2=[]
    for i in range (len(arcs)): 
        obj_func1.append(c1[int(arcs[i][0])-1,int(arcs[i][1])-1])
        obj_func2.append(c2[int(arcs[i][0])-1,int(arcs[i][1])-1])
    return np.sum(obj_func1),np.sum(obj_func2)


def CostMatrix(c1_file,c2_file):
    data1=pd.read_csv('{0}.txt'.format(c1_file),sep=" ",header=None)
    matrix1=np.matrix(data1)
    data2=pd.read_csv('{0}.txt'.format(c2_file),sep=" ",header=None)
    matrix2=np.matrix(data2) 
    return matrix1,matrix2     
