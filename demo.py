import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# If running from command line, above should be ran before
# importing normal_corner
import numpy as np
from normal_corner import normal_corner

### EXAMPLE 1: plotting one covariance matrix ###

# Covariance matrix, as a numpy array
covm = np.array([[2.0, 1.9, -0.1], [1.9, 2.0, 0.5], [-0.1,0.1,2.0]])
print(covm)

# Mean matrix, as a numpy array
mean = np.array([2, -2, 0.5])
print(mean)

# Variable labels for plotting, as a list of strings
# in LaTeX format, between $$ symbols
varlabels = ['$var_1$','$var\\frac{2}{3}$','$var3^3$']

# Make a corner plot
fig1 = normal_corner.normal_corner(covm,mean,varlabels)
plt.savefig('example_1.png')
plt.close()

### EXAMPLE 2: plotting a reducer-order covariance matrix on top ###

covm2 = np.array([[1.0, 0.9], [0.9, 1.0]])
mean2 = np.array([-2, -0.5])
fixedvarindex = 1 # Index of value that we skip (1 means 2nd out of 3)
fixedvarvalue = -2 # Value that a skipped variable has

# Make a corner plot
fig2 = normal_corner.normal_corner(covm,mean,varlabels,fixedvarindex=fixedvarindex,fixedvarvalue=fixedvarvalue,covm2=covm2,mean2=mean2)
plt.savefig('example_2.png')
plt.close()
