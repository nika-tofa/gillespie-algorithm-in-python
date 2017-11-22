from math import log
import random
from numpy import *
def Gillespie(s, i, r, b, g):
    # function that execute Gillespie algorithm for stochastic SIR model
    # contact and infection rates
    beta = b
    gamma = g

    # lists
    rate = [0,0]
    dt = [0,0]
    alldt = [0]

    # initial conditions
    S = [s]
    I = [i]
    R = [r]
    time = 0.0
    timeList = [0]
    tspan = 14

    # counter
    i = 1

    while (time <= tspan):
        # check if I is zero
        if (I[i - 1] == 0):
            S.append(S[i - 1])
            I.append(I[i - 1])
            R.append(R[i - 1])
        else:
            # contact rate
            rate[0] = beta * S[i - 1] * I[i - 1]

            # recovery rate
            rate[1] = gamma * I[i - 1]

            # times to the next event for each event(dt)
            j = 0
            while (j <= 1):
                rnd = random.random()
                dt[j] = -log(rnd) / rate[j]
                j += 1

            # determine the time to the next event
            min_t = min(dt)
            k = dt.index(min(dt))
            alldt.append(min_t)

            # update the system states
            if (k == 0):
                # infection
                S.append(S[i - 1] - 1)
                I.append(I[i - 1] + 1)
                R.append(R[i - 1])
            else:
                # recovery
                S.append(S[i - 1])
                I.append(I[i - 1] - 1)
                R.append(R[i - 1] + 1)

        # increase time to the next time step
        time += dt[k]
        timeList.append(time)
        i += 1

    return S, I, R, dt, timeList, alldt