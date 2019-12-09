############################################################################
# File Name: model12_mySimplex.py                                          #
# Author: Geonsik Yu, Purdue University, IE Dept                           #
# LP problem (Model 12: Alloy Blending) from:                              #
# https://sites.math.washington.edu/~burke/crs/407/models/m12.html         #
############################################################################
import MySimplex

## STEP 1. Set up what we need. -----------------------------------------------------------
## Declare variable names
variables = ["x1", "x2", "x3", "x4", "x5"]
## Delare a list of coefficients  of each variable in the objective function (same order)
obj_coeffs = [19.0, 17.0, 23.0, 21.0, 25.0]
## Delare a list of lowerbounds of each variable
lowerbounds = [0.0, 0.0, 0.0, 0.0, 0.0] 

## Declare contraint names
constraint_names = ["Tin(%)", "Zinc(%)", "Lead(%)"]
## Delare a list of RHS constants of each constraints
righthand = [40.0, 35.0, 25.0]
## Delare a list of inequality directions of each constraints
senses = ['E', 'E', 'E']
## Set coefficients of each variables in each constraints:
lin_expr = [[60.0, 25.0, 45.0, 20.0, 50.0],
	        [10.0, 15.0, 45.0, 50.0, 40.0],
	        [30.0, 60.0, 10.0, 30.0, 10.0]]

## STEP 2. Generate LP problem object -----------------------------------------------------
## Generate an LP problem framework
problem = MySimplex.SimplexProblem()
## Set objective as minimization
problem.setObjectiveDirection( Max=False )
## Set variables and objective function
problem.setVariables( Names=variables, ObjCoeffs=obj_coeffs, Lowerbounds=lowerbounds )
## Set constraints
for idx in range(len(lin_expr)):
	problem.addConstraint( Name = constraint_names[idx]
						, rowVec = lin_expr[idx]
						, ineq_dir = senses[idx]
						, RHS = righthand[idx] )

## STEP 3. Solve the problem --------------------------------------------------------------
problem.setup()
Tableau = problem.buildTableau()
Tableau = problem.solve(Tableau)

