


#FindHC

FindHC is a Python library to deal with the Hamiltonian cycle problem. 

# Installation

pip install -r conda_env/requirements.txt

CPLEX solver: https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290

chmod +x cplex_studio[edition].linux-x86-64.bin 
./cplex_studio[edition].linux-x86-64.bin

Installation path Ubuntu: /opt/ibm/ILOG/CPLEX_Studio[edition]

It also allows GUROBI solver: https://www.gurobi.com/downloads/gurobi-optimizer-eula/

# Usage

# Branch-and-Fix
python3 main.py file_name perm branching simplification solver collapse

file_name: the name of .txt file with the adjacency matrix
perm={'No','perm_file_name'}, perm_file_name: the name of .txt file with the permutation
branching={1,2,3,4,5,6}
simplification={'No','Yes'}
solver={'cplex','gurobi'}
collapse={'No','Yes'}

# Multi-objective Branch-and-Fix

python3 main_mo.py file_name perm branching c1 c2 simplification solver

file_name: the name of .txt file with the adjacency matrix
perm={'No','perm_file_name'}, perm_file_name: the name of .txt file with the permutation
branching={1,2,3,4,5,6}
c1: the name of .txt file with the weight matrix for objective 1
c2: the name of .txt file with the weight matrix for objective 2
simplification={'No','Yes'}
solver={'cplex','gurobi'}


# Authors

Maialen Murua 

maialen.murua@tecnalia.com

Roberto Santana

roberto.santana@ehu.eus
