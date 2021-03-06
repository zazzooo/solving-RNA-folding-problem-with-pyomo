from pyomo.environ import ConcreteModel, Var, Objective, Constraint, SolverFactory, Set
from pyomo.environ import maximize, Binary, RangeSet, Param, ConstraintList, PositiveIntegers
from numpy import loadtxt

model = ConcreteModel()

# d = open('covid_sequence.txt', encoding="utf8")
# dat = d.read()
# d.close()
# dat = dat.replace('\n','')
# dat = dat.replace('A','1').replace('C','2').replace('G','3').replace('T','4')
# data = []
# for i in dat:
#     data.append(int(i))

data = [1,1,1,1,1,2,3,2,4,4,3,2,2,3,2,1,1,1,2,3] #,3,3,3,3,3,4,4,4,4,2,2,2,1,2,2,2,4,3,2,1,1,1,1,1,1,4
D = {index+1 : value for index,value in enumerate(data)}
n = len(data)
#indices
model.I = RangeSet(1, n)
model.J = RangeSet(1, n)

#parameter
model.nucleot = Param(model.I, initialize = D)

#variable
model.P = Var(model.I, model.J, within = Binary, initialize = 0)
model.Q = Var(model.I, model.J, within = Binary, initialize = 0)


#constraint
model.onebond = Constraint(model.I, rule = lambda model, i: sum(model.P[i,j] + model.P[j,i] for j in model.J) <= 1)

def paired_expr(model, i, j):
    if model.nucleot[i] == model.nucleot[j]:
        return model.P[i,j] <= 0
    elif model.nucleot[i] == 1 and model.nucleot[j] == 2:
        return  model.P[i,j] <= 0
    elif model.nucleot[i] == 2 and model.nucleot[j] == 1:
        return  model.P[i,j] <= 0
    elif model.nucleot[i] == 1 and model.nucleot[j] == 3:
        return model.P[i,j] <= 0
    elif model.nucleot[i] == 3 and model.nucleot[j] == 1:
        return model.P[i,j] <= 0
    elif model.nucleot[i] == 2 and model.nucleot[j] == 4:
        return  model.P[i,j] <= 0
    elif model.nucleot[i] == 4 and model.nucleot[j] == 2:
        return  model.P[i,j] <= 0
    elif model.nucleot[i] == 3 and model.nucleot[j] == 4:
        return model.P[i,j] <= 0
    elif model.nucleot[i] == 4 and model.nucleot[j] == 3:
        return model.P[i,j] <= 0
    return model.P[i,j] <= 1

model.bond_exprs = ConstraintList()
for i in range(1,n+1):
    for j in range(1,n+1):
        model.bond_exprs.add(paired_expr(model,i,j))

def no_cross_rule(model, i, j, h, k):
    return model.P[i,j] + model.P[h,k] <= 1

def stack_pairs1(model, i, j):
    return model.P[i,j] + model.P[i+1,j-1] - model.Q[i,j] <= 1

def stack_pairs2(model, i, j):
    return 2*model.Q[i,j] - model.P[i,j] - model.P[i+1,j-1] <= 0


model.no_cross_pair = ConstraintList()
for i in range(1,n+1):
    for j in range(i,n+1):
        for h in range(1,i):
            for k in range(i+1,j):
                model.no_cross_pair.add(no_cross_rule(model,i,j,h,k))
                model.no_cross_pair.add(no_cross_rule(model,j,i,h,k))
                model.no_cross_pair.add(no_cross_rule(model,i,j,k,h))
                model.no_cross_pair.add(no_cross_rule(model,j,i,k,h))

model.stack_pairs = ConstraintList()
for i in range(1,n+1):
    for j in range(i+1,n+1):
        model.no_cross_pair.add(stack_pairs1(model,i,j))
        model.no_cross_pair.add(stack_pairs2(model,i,j))

def no_consecutive_rule(model,i):
  return model.P[i,i+1] <= 0

model.no_consecutive_pair = ConstraintList()
for i in range(1,n):
  model.no_consecutive_pair.add(no_consecutive_rule(model,i))

#obj function
model.obj = Objective(expr = sum(model.P[i,j] + model.Q[i,j] for i in model.I for j in model.J), sense = maximize)

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


#try a plot
