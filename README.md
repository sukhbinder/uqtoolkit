# uqtoolkit
Simple Uncertainty Quantification Toolkit implemented in pure python from a single file.

Implements polynomial chaos.

## Sample usage

```python

import numpy as np
import matplotlib.pyplot as plt
from uqtoolkit import *

# We will implement a simple model, which squares normally distrubuted random numbers
nsamples=5000
f=np.random.standard_normal(nsamples)
Y=f**2

# Defining Polynomial Chaos (PC) exapansion

polyOrder=2 # Polynomial order
param = polyOrder+1
ndim=1  # No of random variables or dimensions
pc_type='HERMITE'

# Create a PC object
pc=uq_pcset(polyOrder,ndim,pc_type)

# Compute the quadtures
quadrature=uq_quadrature(ndim,param,pc_type)
 
# Evaluate quadrature at the PC
K=uq_getNISP(pc,quadrature)

# Evaluate the model at the quadrature nodes
Ygrid=quadrature['nodes']**2

# Project the response on the PC to compute PC coefs
c=np.dot(K,Ygrid)

# Sample 5000 samples from the pc
U=uq_sample(pc,c,nsamples)

# Mean and Std Dev of response
print 'Mean: ',Y.mean(),U.mean()
print 'St dev: ',Y.std(),U.std()

print 'Mean ',c[0]
print 'stddev ',np.sqrt(c.sum()) 

```
