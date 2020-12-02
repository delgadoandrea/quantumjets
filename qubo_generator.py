from pyqubo import Spin, Array
import numpy as np

def make_qubo(particles, theta):
    # Create an array os spin variables
    n_part = len(particles) 
    s = Array.create('s', shape=n_part, vartype='BINARY')

    # Create an array of zeroes for qubo matrix coefficients
    coeff = [[0] * n_part for _ in range(n_part)] 
    #theta = np.pi/4.

    for i in range(0,n_part):
        for j in range(0,n_part):
            coeff[i][j] = (particles[i]['px']*particles[j]['px'] + particles[i]['py']*particles[j]['py'] + particles[i]['pz']*particles[j]['pz'] - particles[i]['E']*particles[j]['E'] * np.cos(theta))/(1 - np.cos(theta))

    # Construct Hamiltonian
    H = sum([ -1.0*coeff[i][j]*s[i]*s[j] for i in range(0,n_part) for j in range(0,n_part)])

    # Compile model using pyqubo

    model = H.compile()
    qubo, offset = model.to_qubo()

    return qubo, offset