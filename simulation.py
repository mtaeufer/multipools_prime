import numpy as np
import matplotlib.pyplot as plt
from function_multipool_matrix import multipool_matrix
# requires the file function_multipool_matrix.py in the same directory 

###
# This script is used for simulations of non-adaptive group testing for the paper
#
# Rapid, Large-Scale, and Effective Detection of COVID-19 Via Non-Adaptive Testing
#
# by Matthias Taeufer.
#
# Visit
# https://doi.org/10.1101/2020.04.06.028431 
# for open-access version of the paper.
###

# maximal_epsilon gives the theoretical false-positive-to-tested-positive ratio in Formula (3.2) of the paper
def maximal_epsilon(rho,n,k):
	if rho == 0:
		return 0
	else:
		return ((1 - rho) * (1 - (1 - rho)**(n-1))**k)/(rho + ((1 - rho) * (1 - (1 - rho)**(n-1))**k))

# n is the pool size
# For the plot in the paper we use n = 31

n = int(raw_input("pool size n:\t"))

# It is assumed that the population size is N = n**2
# See also Section 4 in the paper

N = n**2

# kmin and kmax are the minimal and maximal multiplicity (the number of pools a sample participates in)
# See also Definition 1 in the paper
# We recommend kmin = 3
# We recommend kmax = 3

kmin = int(raw_input("minimal multiplicity k:\t"))
kmax = int(raw_input("maximal multiplicity k:\t"))


# m is the number of populations of N samples simulated
m = 200

# A is a multipool matrix, based on the Shifted Transversal Design
A = multipool_matrix(n,kmax)


# The vector rhovec contains the prevalence values rho considered for the simulation in steps of 0.001
# We need rho <= 1/n , see  Theorem 1 in the paper
rhovec = [i * 0.001 for i in range(1,1000/n)]

print 'Will simulate from rho = 0.001 to rho = ' + str(rhovec[len(rhovec) - 1])

# avvec will store the simulated ratios between false positives and positives
avvec  = np.zeros((kmax - kmin + 1,len(rhovec)))

# epsvec will store the theoretical ratio between false positives and true positives
epsvec = np.zeros((kmax - kmin + 1,len(rhovec)))

# Set a random seed for reproducibility. The next line can be commented out.
np.random.seed(42)


for i in range(len(rhovec)):
	print "rho = " + str(rhovec[i])
	for j in range(kmax - kmin + 1):
		# Number of true positives
		tpos = 0
		# Number of measured positives
		mpos = 0

		epsvec[j][i] = maximal_epsilon(rhovec[i],n,j + kmin)

		for l in range(m):		
			# x is a random {0,1}-vector of samples
			x = np.random.binomial(1,rhovec[i],N)
			# y is the vector of positive pools
			y = [(a != 0) for a in np.dot(A[0:(j + kmin) * n],x)]
			# z is the vector of samples which are tested positive (i.e. all pools of which are tested positive)
			z = [(a == j + kmin) for a in np.dot(A[0:(j + kmin) * n].transpose(),y)]

			tpos  += sum(x)
			mpos  += sum(z)

		avvec[j][i] = float(mpos - tpos)/(float(mpos) + float(mpos == 0)) # correction to avoid dividing by 0



fig = plt.figure()
ax = plt.subplot(111)

# Colors for the plot
colors = ['#ffa500','r','b','g','c'] + ['b','r']*kmax

for j in range(kmax - kmin + 1):
	ax.plot(rhovec, avvec[j], color = colors[j], label = 'Simulated, $k = ' + str(j + kmin) + '$')
	ax.plot(rhovec,epsvec[j], linestyle = '--', color = colors[j], label = 'Theoretical, $k = ' + str(j + kmin) + '$')

box = ax.get_position()

ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 'small')
ax.set_xlabel(r'Incidence $\rho$')
ax.set_ylabel('Ratio of false positives to all positives')
ax.set_title('Ratio of false positives to all positives for $' + str(m*N) + '$ samples in pools of size $' + str(n) + '$', fontsize = 'small', y = 1.02) # increase y coordinate of title to avoid overlap with axis labels

plt.savefig("n=" + str(n) + "_k=" + str(kmin) + "_to_" + str(kmax) + "_m=" + str(m) + ".png")
plt.show()
