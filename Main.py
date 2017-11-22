import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from OrdDifEq import *
from GillespieAlg import *

# Main file for solving the SIR model equations with the Gillespie algorithm. 
# The code was transformed from Matlab code that was written for my Master thesis project

# total population, N.
N = 7630
# initial number of infected and recovered individuals, I0 and R0
I0, R0 = 10, 0
# Initial number of susceptible individuals
S0 = N - I0 - R0

# contact rate beta
# infection rate gamma
beta, gamma = 0.002, 0.4

# time (in days)
t = np.linspace(0, 14, 14)
tspan = 14

# Initial conditions vector
y0 = S0, I0, R0

# Solve the model equations over the time t.
res = odeint(bbs_ode,y0,t, args=(beta, gamma))
S, I, R = res.T

# solve the equations with Gillespie algorithm
Sg, Ig, Rg, dt, time, alldt = Gillespie(S0, I0, R0, beta, gamma)

# Plot the data for S(t), I(t) and R(t)
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111, axis_bgcolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', lw=1.5, label='ODEInt Susceptible')
ax.plot(t, I, 'r', lw=1.5, label='ODEInt Infected')
ax.plot(t, R, 'g', lw=1.5, label='ODEInt Recovered')
ax.plot(time, Sg, 'y', lw=2, label='Gillespie Susceptible')
ax.plot(time, Ig, 'orange', lw=2, label='Gillespie Infected')
ax.plot(time, Rg, 'm', lw=2, label='Gillespie Recovered')
plt.title('Solution of the SIR model with inner ODE function and Gillespie method')
ax.set_xlabel('Time in days')
ax.set_ylabel('Number of individuals')
ax.set_ylim(0,N+50)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')

legend = ax.legend()
legend.get_frame().set_alpha(0.5)

for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

fig2 = plt.figure(2, facecolor='w')
ax2 = fig2.add_subplot(111, axis_bgcolor='#dddddd', axisbelow=True)
ax2.plot( time, alldt, 'b', alpha=0.5, lw=1, label='$\Delta$t')
ax2.set_xlabel('Time in days')
ax2.set_ylabel('$\Delta$t')
ax2.grid(b=True, which='major', c='w', lw=1.5, ls='-')
plt.title('Distribution of $\Delta$t')
legend = ax2.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax2.spines[spine].set_visible(False)
plt.show()
