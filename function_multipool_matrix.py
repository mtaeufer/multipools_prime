import numpy as np

# multipool_matrix is a function defining multipool matrices based on the Standard Transversal Design if
# the pool size is a prime number n
# the population size is N = n**2
# the multiplicity k is a number between 1 and n+1

def multipool_matrix(n,k):  
	# check constraint k <= n+1
	if k > n+1: 
		print "Error! k must not be larger than n+1."
		return
	# check constraint n prime if k > 3
	if k > 3:
		for i in range(2,n):
			if (n % i) == 0:
 				print("Error! If k > 3 then n must be a prime number.")
 				return

	A = np.zeros((n*k,n**2))

	for i in range(n):
		for j in range(n):
			A[i][i*n+j] = 1 
	for l in range(1,k):
		for i in range(n):
			for j in range(n):
				A[l*n+i][n*j+((i+(l-1)*j)%n)] = 1
	return A
