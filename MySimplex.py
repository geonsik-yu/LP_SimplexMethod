from builtins import object
import numpy as np

"""
Programming Requirements:
(1) Design a routine to do the conversion to the standard form LP.
(2) Design a routine to detect if the LP at hand is feasible or not. (15 pt)
(3) Design a routine to detect if the coefficient matrix associated with 
   the standard form of the LP has full row rank; and if not, how to remove 
   the redundant constraints. (15 pt)
(4) Design a routine to check if a basic feasible solution (BFS) with the identity 
    matrix as the corresponding basis is readily available. 
    If yes, you can just go ahead to use the BFS to start your simplex method; 
    otherwise, use the methods we learned in class to start your simplex with 
    an identity matrix as the initial basis. (15 pt)
(5) Your code needs to implement one rule to prevent the simplex method from
	cycling. (15 pt)
(6) Termination: your code needs to be able to handle both cases at termination 
	(6-1) a finite optimal solution or (6-2) unbounded (15 pt)
"""

from numpy import dot, zeros
from numpy.linalg import matrix_rank, norm

### Global threshold value for checking zero value with floating point errors.
THRESHOLD = 1e-8


"""
Problem Standardization.

(1) Transform the objective direction into "Minimization" and
    negate the signs of objective coefficients.

(2) Transform the input constraint:
    2-1) [LHS <= RHS form]: pass.
    2-2) [LHS >= RHS form]: change to [-LHS <= -RHS].
    2-3) [LHS = RHS form]: split into [LHS <= RHS form] and [-LHS <= -RHS form].

(3) Check the resulting RHS's from (2). 
    IF all the RHS's are positive:
        Go add (+1 coeff) slack variables for each constraints.
        Check the problem type as "COMMON".
        Check the problem feasibility as "FEASIBLE"
    Else:
        Go add (+1 coeff) slack variables for constraints with positive RHS's (name S_i),
           add (-1 coeff) slack variables for constraints with negative RHS's (name S_i), and
           add (+1 coeff) artificial variables for constraints with negative RHS's (name A_i).
        Check the problem type as "TWOSTEP".

(3) Reduce problem by removing redundant row vectors from [A; b] matrix. 
    

(4) IF the problem type is "COMMON":
        Run tableau simplex steps until it meets one of the [Tableau Simplex's Termination Conditions]
    ELSE:
        (a) Put a secondary objective function with (+1 coeff) for all artificial variables.
        (b) Run tableau simplex steps until it meets one of the [Tableau Simplex's Termination Conditions]
        (c) IF the optimal objective value from (b) is zero:
                Then put the 
            ELSE:
                Check the problem feasibility as "FEASIBLE" and end procedure.

[Tableau Simplex's Termination Conditions]
(T1) If all the reduced costs are negative.
    - "Optimal corner point found"
    - Return the optimal solution and the optimal objective value.
(T2) If for all positive reduced cost columns, there are only non-positive elememts.
    - "Unboundedness found"
    - Return "-inf" for minimization problems and "+inf" for maximization problems.




3-1. If the standardized problem's initial tableau can be start from it, then:
	3-1-1) Run the simplex steps until the termination condition is satisfied.
	3-1-2-a)  
	3-1-2-b) Output 

3-2. Otherwise, generate "1st step subproblem" for [Two Step Method].
	3-2-a) Run the simplex steps until the termination condition is satisfied.

4-1. If the "1st step subproblem" 

4-2) Otherwise, 
	
	terminate the procedure.



"""


def findPivot(self):
    pass

def checkOptimality(self):



class SimplexProblem(object):
    """

    """

    def __init__(self):
        ############################################################################
        # Constructor: This function constructs an LP problem object               #
        ############################################################################

        ## Attribute that will be initialized by setObjectiveDirection() function.
        self._objDirection = None	 # Decides whether the given problem is maximization or minimization.

        ## Attributes that will be initialized by setVariables() function.
        self._objCoeffs = None		 # Coefficient values in the objective function for variables.
        self._varNames = None		 # Names for variables.
        self._varUpperbounds = None  # Upperbounds for non-slack / non-artificial variables. 
        self._varLowerbounds = None  # Lowerbounds for non-slack / non-artificial variables.
        self._varCount = None 		 # Number of non-slack, non-artificial objCoeffs.
        self._objSecondary = None    # Coefficient values in the secondary objective function
                                     # for non-slack / non-artificial variables.

        # Attributes that will be initialized by addConstraint function.
        self._conNames = None		 # List of names for all the constraints.
        self._AMatrix = None		 # Coefficient Matrix "A" 
        self._ineqDirs = None		 # 
        self._RHSs = None

        self._basisSet = []

   	def setObjectiveDirection(self, Max = True):
        ############################################################################
        # This function initializes the direction of optimization: If the problem  #
        # is on maximization, Max = True (Default), otherwise, Max = False.        #
        ############################################################################
        self._objDirection = Max

    def setVariables(self, Names, ObjCoeffs, Upperbounds=None, Lowerbounds=None):
        ############################################################################
        # This function initializes the decision variable set, their names, their  #
        # coefficients in the objective function, and their upper & lower bounds.  # 
        ############################################################################
        self._objCoeffs = None
        self._varNames = None
        self._varUpperbounds = None
        self._varLowerbounds = None        

   	def addConstraint(self, Name, rowVec, ineq_dir, RHS):
        ############################################################################
        # This function adds one constraint at a time. Each constraint requires to #
        # have coefficient vector with the same order with the variable vector, one#
        # inequality direction, and one right-hand side constant value.            #
        ############################################################################
        self._conNames.append(Name)
        self._AMatrix.append(rowVec)
        self._ineqDirs.append(ineq_dir)
        self._RHSs.append(RHS)

    def solve():
        ############################################################################

        ############################################################################



    def run_oneIteration():
        ############################################################################

        ############################################################################

    def check_optimality():
        ############################################################################

        ############################################################################

	def rref(self, ProblemMatrix):
        ############################################################################
        # This function returns (1) the reduced row echelon form (RREF) of the     #
        # input matrix and (2) the list of original row indices from the input     #
        # matrix of each row of RREF.    										   #
        # 																		   #
        # * ProblemMatrix	= 	[A; b]											   #
        # 																		   #
        # Reference Source: 													   #
        # https://stackoverflow.com/questions/7664246/python-built-in-function-to- #
        #  do-matrix-reduction/7665269#7665269               					   #
        ############################################################################
		Matrix = InputMatrix.copy()
		rows, cols = Matrix.shape
		r = 0
		pivots_pos = []
		row_exchanges = np.arange(rows)
		for c in range(cols):

			## Find the pivot row:
			pivot = np.argmax (np.abs (Matrix[r:rows,c])) + r
			m = np.abs(Matrix[pivot, c])
			if m <= THRESHOLD:
				## Skip column c, making sure the approximately zero terms are
				## actually zero.
				Matrix[r:rows, c] = np.zeros(rows-r)
			else:
				## keep track of bound variables
				pivots_pos.append((r,c))

				if pivot != r:
					## Swap current row and pivot row
					Matrix[[pivot, r], c:cols] = Matrix[[r, pivot], c:cols]
					row_exchanges[[pivot,r]] = row_exchanges[[r,pivot]]
					

				## Normalize pivot row
				Matrix[r, c:cols] = Matrix[r, c:cols] / Matrix[r, c];

				## Eliminate the current column
				v = Matrix[r, c:cols]
				## Above (before row r):
				if r > 0:
					ridx_above = np.arange(r)
					Matrix[ridx_above, c:cols] = Matrix[ridx_above, c:cols] - np.outer(v, Matrix[ridx_above, c]).T
				## Below (after row r):
				if r < rows-1:
					ridx_below = np.arange(r+1,rows)
					Matrix[ridx_below, c:cols] = Matrix[ridx_below, c:cols] - np.outer(v, Matrix[ridx_below, c]).T
				r += 1
			## Check if done
			if r == rows:
				break;
		return (Matrix, row_exchanges)

	def reduceProblemMatrix(self, ProblemMatrix):
        ############################################################################
        # This function returns the problem matrix [A;b] after removing row vector #
        # redundancy.    														   #
        ############################################################################
		rref_results = rref(ProblemMatrix)
		rowSums = rref_results[0].sum(axis=1, dtype='float')
		remainderIndices = []
		for idx in range(len(rowSums)):
			if rowSums[idx] < THRESHOLD:
				break;
			remainderIndices.append(rref_results[1][idx])
		return InputMatrix[remainderIndices,:]

    def loss(self, X, y=None):
        """
        Evaluate loss and gradient for the four-layer convolutional network.

        Input / output: Same API as TwoLayerNet in fc_net.py.
        """
        W1, b1 = self.params['W1'], self.params['b1']
        W2, b2 = self.params['W2'], self.params['b2']
        W3, b3 = self.params['W3'], self.params['b3']
        W4, b4 = self.params['W4'], self.params['b4']

        # pass conv_param to the forward pass for the convolutional layer
        # Padding and stride chosen to preserve the input spatial size
        filter_size = W1.shape[2]
        conv_param = {'stride': 1, 'pad': (filter_size - 1) // 2}

        # pass pool_param to the forward pass for the max-pooling layer
        pool_param = {'pool_height': 2, 'pool_width': 2, 'stride': 2}

        scores = None
        ############################################################################
        # TODO: Implement the forward pass for the four-layer convolutional net,    #
        # computing the class scores for X and storing them in the scores          #
        # variable.                                                                #
        #                                                                          #
        # Remember you can use the functions defined in ie590/fast_layers.py and   #
        # ie590/layer_utils.py in your implementation (already imported).          #
        ############################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

        pass

        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        ############################################################################
        #                             END OF YOUR CODE                             #
        ############################################################################

        if y is None:
            return scores

        loss, grads = 0, {}
        ############################################################################
        # TODO: Implement the backward pass for the four-layer convolutional net,  #
        # storing the loss and gradients in the loss and grads variables. Compute  #
        # data loss using softmax, and make sure that grads[k] holds the gradients #
        # for self.params[k]. Don't forget to add L2 regularization!               #
        #                                                                          #
        # NOTE: To ensure that your implementation matches ours and you pass the   #
        # automated tests, make sure that your L2 regularization includes a factor #
        # of 0.5 to simplify the expression for the gradient.                      #
        ############################################################################
        # *****START OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****

        pass

        # *****END OF YOUR CODE (DO NOT DELETE/MODIFY THIS LINE)*****
        ############################################################################
        #                             END OF YOUR CODE                             #
        ############################################################################

        return loss, grads