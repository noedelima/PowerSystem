'''
Noé de Lima Bezerra
30/12/2016 - 2h53

Funções de Estimação do Estado do Sistema Elétrico de Potência
Derivadas das funções de estado para a matriz jacobiana
'''

import math
from Matrix import *
from PowerSystemTopology import *
from PowerSystemState import *


'''
Estimate the power system state variables using the power system network set
and the measures set of field.
The results can be used to calculate all system parameters.
The equation to solve is:
G = H(x)'(R^-1)H(x)
t = H(x)'(R^-1)(z-h(x))
dx = (G^-1)t
'''


def StateEstimator(PSNetwork, SetMeas, tolerance):
    Top = topology(PSNetwork)
    dim = Top.dimension #The dimension of power system network (number of bus in network)
    V = Voltage(Top) #Create the state vector with Modules and Phase angules of all buses
    Power = StatePowerSystem(Top)
    tolint = tolerance  #set the internal tolerance variable value
    tol = 10**(-tolint)
    inc = 10 #Set the initial residue for incremental iteration
    iteration = 0
    ittot = 0
    H, Rn, z = Jacobian(Top, V, SetMeas)
    HtRn = H.transpose().multiply(Rn) # H'R^-1
    G = HtRn.multiply(H) # H'R^-1H
    try:
        Gn = G.inverse() # G^-1
    except:
        return None
    while inc > tol: #Newton-Raphson Solution Method
        Power.StateActualize(Top, V)
        # The lines below would use the complete formulation of the Newton-Raphson solution method, updating the gain matrix at each iteration
        #H, Rn, z = Jacobian(Top, V, SetMeas)
        #HtRn = H.transpose().multiply(Rn)  # H'R^-1
        #G = HtRn.multiply(H) # H'R^-1H
        #Gn = G.inverse()  # G^-1
        h = StateVector(Top, V, SetMeas, Power)
        zh = z.matrixsum(h.multiplysc(-1)) # zh = z - h(x)
        t = HtRn.multiply(zh)
        inc = Gn.multiply(t) # dx = (G^-1)H'(R^-1)(z-h)
        dV = Voltage(Top)
        for i in range(1, dim):
            dV.Phase.matrix[i][0] = inc.matrix[i-1][0]
        for i in range(dim):
            dV.Magnitude.matrix[i][0] = inc.matrix[dim-1+i][0]
        V.Magnitude = V.Magnitude.matrixsum(dV.Magnitude)
        V.Phase = V.Phase.matrixsum(dV.Phase)
        for i in range(inc.line):
            inc.matrix[i] = math.sqrt(inc.matrix[i][0]**2)
        inc = max(inc.matrix)
        iteration += 1
        ittot += 1
        if iteration > 20:
            if tolint > 1:
                tolint -= 1
                tol = (10**(-tolint))
                iteration = 0
            else:
                print('This system did not converge for', ittot, 'iterations done!')
                return None
    Power.StateActualize(Top, V)
    I = Current(Top)
    I.ActualizeCurrent(Top, V)
    print('System converged for', ittot, 'iterations done!')
    return V, Power, I, tolint


'''
Estimate the power system state variables using the power system network set
and the measures set of field.
The decoupled solution method is used here to Newton-Raphson numerical solution
The results can be used to calculate all system parameters.
The equation to solve is:
Ga = Ha(x)'(Ra^-1)Ha(x)
ta = Ha(x)'(Ra^-1)(za-ha(x))
dx = (Ga^-1)ta
Gr = Hr(x)'(Rr^-1)Hr(x)
tr = Hr(x)'(Rr^-1)(zr-hr(x))
dx = (Gr^-1)tr
with 'a' as Active Power part and 'r' as Reactive Power part
'''


def StateEstimatorDC(PSNetwork, SetMeas, tolerance):
    Top = topology(PSNetwork)
    dim = Top.dimension #The dimension of power system network (number of bus in network)
    V = Voltage(Top) #Create the state vector with Modules and Phase angules of all buses
    Power = StatePowerSystem(Top)
    tolint = tolerance  #set the internal tolerance variable value
    tol = 10**(-tolint)
    inc = 10 #Set the initial residue for incremental iteration
    iteration = 0
    ittot = 0
    # Decoupled formulation for Jacobian and Gain
    Haa, Ran, za = JacobianAA(Top, V, SetMeas)
    HaatRan = Haa.transpose().multiply(Ran) # Haa'(Ra^-1)
    Gaa = HaatRan.multiply(Haa) # Haa'(Ra^-1)Haa
    try:
        Gaan = Gaa.inverse() # Gaa^-1
    except:
        return None
    Hrr, Rrn, zr = JacobianRR(Top, V, SetMeas)
    HrrtRrn = Hrr.transpose().multiply(Rrn) # Hrr'(Rr^-1)
    Grr = HrrtRrn.multiply(Hrr) # Hrr'(Rr^-1)Hrr
    try:
        Grrn = Grr.inverse() # Grr^-1
    except:
        return None
    a, r = Haa.line, Hrr.line
    while inc > tol: #Newton-Raphson Solution Method
        inc = 0
        dV = Voltage(Top)

        # Increment to Voltage Phase Angle
        # The lines below would use the complete formulation of the Newton-Raphson solution method, updating the gain matrix at each iteration
        #Haa, Ran, za = JacobianAA(Top, V, SetMeas)
        #HaatRan = multiply(transpose(Haa), Ran)  # Haa'(Ra^-1)
        #Gaa = multiply(HaatRan, Haa)  # Haa'(Ra^-1)Haa
        #Gaan = inverse(Gaa)  # Gaa^-1
        Power.StateActualize(Top, V)
        h = StateVector(Top, V, SetMeas, Power)
        zha = matrix(a, 1) # zhaa = za - ha(x)
        for i in range(a):
            zha.matrix[i][0] = za.matrix[i][0] - h.matrix[i][0]
        ta = HaatRan.multiply(zha)
        incT = Gaan.multiply(ta) # dTheta = (Gaa^-1)Haa'(Ra^-1)(za-ha)
        for i in range(1, dim):
            dV.Phase[i][0] = incT.matrix[i-1][0]
        iteration += 0.5
        ittot += 0.5


        # Increment to Voltage Magnitude
        # The lines below would use the complete formulation of the Newton-Raphson solution method, updating the gain matrix at each iteration
        #Hrr, Rrn, zr = JacobianRR(Top, V, SetMeas)
        #HrrtRrn = multiply(transpose(Hrr), Rrn)  # Hrr'(Rr^-1)
        #Grr = multiply(HrrtRrn, Hrr)  # Hrr'(Rr^-1)Hrr
        #Grrn = inverse(Grr)  # Grr^-1
        Power.StateActualize(Top, V)
        h = StateVector(Top, V, SetMeas, Power)
        zhr = matrix(r, 1) # zhrr = zr - hr(x)
        for i in range(r):
            zhr.matrix[i][0] = zr.matrix[i][0] - h.matrix[a+i][0]
        tr = HrrtRrn.multiply(zhr)
        incV = Grrn.multiply(tr) # dV = (Grr^-1)Hrr'(Rr^-1)(zr-hr)
        # Put increment values in dV vector
        for i in range(dim):
            dV.Magnitude.matrix[i][0] = incV.matrix[i][0]
        iteration += 0.5
        ittot += 0.5

        V.Magnitude = V.Magnitude.matrixsum(dV.Magnitude)
        V.Phase = V.Phase.matrixsum(dV.Phase)

        # Verify the maximum value of dx
        for i in incT.matrix:
            if math.fabs(i[0]) > inc:
                inc = math.fabs(i[0])
        for i in incV.matrix:
            if math.fabs(i[0]) > inc:
                inc = math.fabs(i[0])

        if iteration > 50:
            if tolint > 1:
                tolint -= 1
                tol = (10**(-tolint))
                iteration = 0
            else:
                print('This system did not converge for', ittot, 'iterations done')
                return None
    Power.StateActualize(Top, V)
    I = Current(Top)
    I.ActualizeCurrent(Top, V)
    print('System converged for', ittot, 'iterations done!')
    return V, Power, I, tolint
