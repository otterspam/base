from __future__ import division
import numpy as np
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from pyomo.environ import *
from pyomo.opt import SolverFactory

# Sudoku Solver

# Package requirements: Installation Guide
#  conda install -c https://conda.anaconda.org/conda-forge pyomo pyomo.extras
#  conda install glpk ipopt_bin -c cachemeorg
#  conda install pandas
#  conda install numpy

##### Mixed-Integer Program Formulation:
## Decision variables:
## x[ii,jj,kk] in {0,1}
## ii denotes row number, jj denotes column number, kk denotes value
## e.g. x[1,2,3] = 1 means that cell (1,2) has value 3


##### Inputs: can either be explicit or read in from excel
# explicit inputs in the form of [row,col,val]
# e.g. if cell (1,1) has value 5 then [1, 1, 5] should be in this list
 current_board = np.array([
     [1, 4, 5],
     [1, 6, 1],
     [2, 3, 1],
     [2, 7, 6],
     [3, 1, 8],
     [3, 2, 4],
     [3, 4, 9],
     [3, 7, 7],
     [3, 9, 5],
     [4, 6, 3],
     [4, 8, 4],
     [4, 9, 6],
     [5, 1, 6],
     [5, 6, 7],
     [5, 8, 9],
     [5, 9, 2],
     [6, 2, 5],
     [6, 6, 8],
     [7, 1, 9],
     [7, 7, 4],
     [8, 3, 6],
     [8, 9, 3],
     [9, 1, 4],
     [9, 2, 1],
     [9, 9, 8]])

# inputs from file ignores any integers not between 1 and 9
# puzzle_input = pd.read_excel('<Redacted file>',
#     sheet_name='PuzzleInput', header=None, dtype=pd.Int64Dtype())

current_board_list = []
for row in range(9):
    for col in range(9):
        if puzzle_input.iloc[row,col] in list(range(1,10)):
            current_board_list.append([row+1,col+1,puzzle_input.iloc[row,col]])

current_board = np.array(current_board_list)

##### Set up optimization problem
# Optimization solver options: uses glpk as MIP solver and a gap tolerance of 1%
# note that the problem is really a feasibility one and so the gap will rarely be used
opt = SolverFactory('glpk')
opt.options["mipgap"] = 0.01
model = ConcreteModel()

model.rowl = RangeSet(1,9)
model.coll = RangeSet(1,9)
model.numl = RangeSet(1,9)
model.x = Var(model.rowl, model.coll, model.numl, within = Binary)

# Read-in file requires that the initial cells be populated
starting_board = np.zeros((9,9), dtype=np.int)
for fixed_cell in range(current_board.shape[0]):
    ii = current_board[fixed_cell,0]
    jj = current_board[fixed_cell,1]
    kk = current_board[fixed_cell,2]
    model.x[ii,jj,kk].fix(1)
    starting_board[ii-1, jj-1]=kk

# constraint: one value per cell:
# for each ii,jj sum(x[ii,jj,kk] over kk) = 1
def CellOccupancy(model, ii, jj):
    return sum(model.x[ii,jj,kk] for kk in model.numl) == 1
model.c1 = Constraint(model.rowl, model.coll, rule=CellOccupancy)

# constraint: each number appears in each row:
# for each kk,jj sum(x[ii,jj,kk] over ii) = 1
def RowSumReq(model, jj, kk):
    return sum(model.x[ii,jj,kk] for ii in model.rowl) == 1
model.c2 = Constraint(model.coll, model.numl, rule=RowSumReq)

# constraint: each number appears in each column:
# for each ii,kk sum(x[ii,jj,kk] over jj) = 1
def ColSumReq(model, ii, kk):
    return sum(model.x[ii,jj,kk] for jj in model.coll) == 1
model.c3 = Constraint(model.rowl, model.numl, rule=ColSumReq)

# constraint: each number appears in each predefined box
# box 1 has ii=1:3, jj=1:3; box 2 has ii=1:3, jj=4:6 and so on
model.c4 = ConstraintList()
for kk in model.numl:
    for boxi in range(3):
        for boxj in range(3):
            model.c4.add(
                model.x[3*boxi+1,3*boxj+1,kk] +
                model.x[3*boxi+1,3*boxj+2,kk] +
                model.x[3*boxi+1,3*boxj+3,kk] +
                model.x[3*boxi+2,3*boxj+1,kk] +
                model.x[3*boxi+2,3*boxj+2,kk] +
                model.x[3*boxi+2,3*boxj+3,kk] +
                model.x[3*boxi+3,3*boxj+1,kk] +
                model.x[3*boxi+3,3*boxj+2,kk] +
                model.x[3*boxi+3,3*boxj+3,kk] == 1)

# objective: since this is really a feasibility problem, anything will do here
def ObjectiveRule(model):
    return sum( model.x[ii,jj,kk]*(9*ii+81*jj+kk) 
    for ii in model.rowl for jj in model.coll for kk in model.numl)
model.o = Objective(rule=ObjectiveRule, sense=maximize)

##### Solve the optimization problem
results = opt.solve(model, logfile="sudoku.log")#, tee=True)

##### Outputs:
solution_board = np.zeros((9,9), dtype=np.int)
model.solutions.load_from(results)
for v in model.component_objects(Var, active=True): 
    varobject = getattr(model, str(v))
    for index in varobject:
        if varobject[index].value == 1:
            solution_board[index[0]-1, index[1]-1] = index[2]

def print_board(arr):
    print("+-------+-------+-------+")
    for ii, row in enumerate(arr):
        print(
            "|", np.array2string(row[0:3]).lstrip('[').rstrip(']'), "|",
            np.array2string(row[3:6]).lstrip('[').rstrip(']'), "|",
            np.array2string(row[6:9]).lstrip('[').rstrip(']'), "|") 
        if ii % 3 == 2:
            print("+-------+-------+-------+")

print()
print()
print("<Redacted> Sudoku solver version 1.0")
print()
if (results.solver.termination_condition == TerminationCondition.infeasible):
    print("This puzzle cannot be solved")
    print("No refunds")
else:
    print("Your starting board:")
    print_board(starting_board)
    print()
    print("Your solution is:")
    print_board(solution_board)
print()
print("Thank you for using <Redacted> for your puzzle cheating needs")
