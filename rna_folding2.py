from pyomo.environ import ConcreteModel, Var, Objective, Constraint, SolverFactory, Set
from pyomo.environ import maximize, Binary, NonNegativeReals, RangeSet, Param, ConstraintList, PositiveIntegers
import  numpy as np


model = ConcreteModel()

#inserisco dei dati (fantocci)
data = [1,1,1,1,1,2,3,2,4,4,3,2,2,3,2,1,1,1,2,3]
D = {index+1 : value for index,value in enumerate(data)}
n = len(data)

#indices
model.I = RangeSet(1, n)
model.J = RangeSet(1, n)

#parameter
model.nucleot = Param(model.I, initialize = D)

#variable, aggiungo Q
model.P = Var(model.I, model.J, within = Binary, initialize = 0)
model.Q = Var(model.I, model.J, within = Binary, initialize = 0)

#costruisco la matrice dei costi C
def matrice_costi(data):
    n = len(data)
    C = np.ones((n,n))
    for i in range(n):
        for j in range(n):
          if not ((data[i] == 1 and data[j] == 4) or (data[i] == 4 and data[j] == 1)):
              if not ((data[i] == 2 and data[j] == 3) or (data[i] == 3 and data[j] == 2)):
                  C[i,j] = 0.7
    return C

C = matrice_costi(data)

#constraint
model.onebond = Constraint(model.I, rule = lambda model, i: sum(model.P[i,j] + model.P[j,i] for j in model.J) <= 1)

def stack_pairs1(model, i, j):
    return model.P[i,j] + model.P[i+1,j-1] - model.Q[i,j] <= 1

def stack_pairs2(model, i, j):
    return 2*model.Q[i,j] - model.P[i,j] - model.P[i+1,j-1] <= 0

model.stack_pairs = ConstraintList()
for i in range(1,n+1):
    for j in range(i+1,n+1):
        model.stack_pairs.add(stack_pairs1(model,i,j))
        model.stack_pairs.add(stack_pairs2(model,i,j))

def no_cross_rule(model, i, j, h, k):
    return model.P[i,j] + model.P[h,k] <= 1

model.no_cross_pair = ConstraintList()
for i in range(1,n+1):
    for j in range(i,n+1):
        for h in range(1,i):
            for k in range(i+1,j):
                model.no_cross_pair.add(no_cross_rule(model,i,j,h,k))
                model.no_cross_pair.add(no_cross_rule(model,j,i,h,k))
                model.no_cross_pair.add(no_cross_rule(model,i,j,k,h))
                model.no_cross_pair.add(no_cross_rule(model,j,i,k,h))

def no_consecutive_rule(model,i):
    return model.P[i,i+1] <= 0

model.no_consecutive_pair = ConstraintList()
for i in range(1,n):
    model.no_consecutive_pair.add(no_consecutive_rule(model,i))

#obj function
model.obj = Objective(expr = sum(C[i-1,j-1]*model.P[i,j] + model.Q[i,j] for i in model.I for j in model.J), sense = maximize)


sol = SolverFactory('glpk').solve(model)
sol_json = sol.json_repn()

if sol_json['Solver'][0]['Status'] != 'ok':
    print("Problem unsolved")
if sol_json['Solver'][0]['Termination condition'] != 'optimal':
    print("Problem unsolved")

for i in range(1,n+1):
    for j in range(1,n+1):
        if model.P[i,j] == 1:
            print(i,j)
