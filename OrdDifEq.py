def bbs_ode(y, t, beta, gamma):
# This file evaluates the right-hand side of the predator-prey ODE system
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt