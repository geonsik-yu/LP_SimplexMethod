import numpy as np

THRESHOLD = 1e-8

def rref(InputMatrix):
	A = InputMatrix.copy()
	rows, cols = A.shape
	r = 0
	pivots_pos = []
	row_exchanges = np.arange(rows)
	for c in range(cols):

		## Find the pivot row:
		pivot = np.argmax (np.abs (A[r:rows,c])) + r
		m = np.abs(A[pivot, c])
		if m <= THRESHOLD:
			## Skip column c, making sure the approximately zero terms are
			## actually zero.
			A[r:rows, c] = np.zeros(rows-r)
		else:
			## keep track of bound variables
			pivots_pos.append((r,c))

			if pivot != r:
				## Swap current row and pivot row
				A[[pivot, r], c:cols] = A[[r, pivot], c:cols]
				row_exchanges[[pivot,r]] = row_exchanges[[r,pivot]]
				

			## Normalize pivot row
			A[r, c:cols] = A[r, c:cols] / A[r, c];

			## Eliminate the current column
			v = A[r, c:cols]
			## Above (before row r):
			if r > 0:
				ridx_above = np.arange(r)
				A[ridx_above, c:cols] = A[ridx_above, c:cols] - np.outer(v, A[ridx_above, c]).T
			## Below (after row r):
			if r < rows-1:
				ridx_below = np.arange(r+1,rows)
				A[ridx_below, c:cols] = A[ridx_below, c:cols] - np.outer(v, A[ridx_below, c]).T
			r += 1
		## Check if done
		if r == rows:
			break;
	return (A, row_exchanges)

def reduceMatrix(InputMatrix):
	rref_results = rref(InputMatrix)
	rowSums = rref_results[0].sum(axis=1, dtype='float')
	remainderIndices = []
	for idx in range(len(rowSums)):
		if rowSums[idx] < THRESHOLD:
			break;
		remainderIndices.append(rref_results[1][idx])

	return InputMatrix[remainderIndices,:]

	print(rowSums)
	print(rref_results)
	print(remainderIndices)
	print(InputMatrix[remainderIndices,:])
matrix = np.array(
		[
				[0, 1 ,0 ,0],
				[0, 0, 1, 0],
				[0, 1, 1, 0],
				[1, 0, 0, 1]
		])
"""
matrix = np.array(
		[
				[0, 1 ,0 ,0, 13, 23,   1, 20, -1],
				[0, 0, 1, 0, 10,  9, 1.1, -3,-12],
				[0, 1, 1, 0, 11, -3,11.8, 13, 12],
				[1, 0, 0, 1,  1, -1,-1.7,  2, 23]
		])
"""
#print(matrix)
#print(matrix.echelon_form())
print(matrix)
reduceMatrix(matrix)

matrix = np.array(
		[
				[0, 1 ,0 ,0],
				[0, 0, 1, 0],
				[0, 0, 1, 0],
				[0, 0, 1, 0],
				[0, 0, 1, 0],
				[0, 1, 1, 0],
				[1, 0, 0, 1],
				[0, 1, 1, 0],
				[0, 1, 1, 0],
				[0, 1, 1, 0],
				[0, 1, 1, 0],
				[0, 1, 1, 0]
		])
print(matrix)
reduceMatrix(matrix)

#print(matrix)
#print(rref(matrix))
