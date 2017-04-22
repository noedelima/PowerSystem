'''
Noé de Lima Bezerra
30/12/2016 - 2h53

Método para obter a topologia da rede elétrica
O objetivo é calcular a matriz Ybus do sistema e armazenar todos os parâmetros
de topologia necessários para calcular as variáveis de estado,
bem como fornecer suporte aos módulos de estimação de estado e análise de
observabilidade de rede.
'''

import math
from Matrix import *
from PowerSystemTopology import *


class Voltage:
    def __init__(self, Topology):
        dim = Topology.dimension
        self.Magnitude = matrix(dim, 1)
        for i in range(dim):
            self.Magnitude.matrix[i][0] = 1
        self.Phase = matrix(dim, 1)

class StatePowerSystem:
    def __init__(self, Top):
        dim = Top.dimension
        self.RealPower = matrix(dim)
        self.ReacPower = matrix(dim)

    def StateActualize(self, Top, V):
        dim = Top.dimension
        for i in range(dim):
            self.RealPower.matrix[i][i] = 0
            self.ReacPower.matrix[i][i] = 0
            for j in range(dim):
                G = Top.Ybus.matrix[i][j].real
                B = Top.Ybus.matrix[i][j].imag
                Vi = V.Magnitude.matrix[i][0]
                Vj = V.Magnitude.matrix[j][0]
                Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
                self.RealPower.matrix[i][i] += Vj*(G*math.cos(Tij)+B*math.sin(Tij))
                self.ReacPower.matrix[i][i] += Vj*(G*math.sin(Tij)-B*math.cos(Tij))
            self.RealPower.matrix[i][i] = self.RealPower.matrix[i][i] * Vi
            self.ReacPower.matrix[i][i] = self.ReacPower.matrix[i][i] * Vi
        for i in range(dim):
            for j in range(dim):
                if not i == j:
                    g = Top.Y.matrix[i][j].real
                    b = Top.Y.matrix[i][j].imag
                    gsh = Top.Ysh.matrix[i][j].real
                    bsh = Top.Ysh.matrix[i][j].imag
                    Vi = V.Magnitude.matrix[i][0]
                    Vj = V.Magnitude.matrix[j][0]
                    Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
                    self.RealPower.matrix[i][j] = (Vi**2)*(g+gsh) - Vi*Vj*(g*math.cos(Tij)+b*math.sin(Tij))
                    self.ReacPower.matrix[i][j] = (-Vi**2)*(b+bsh) - Vi*Vj*(g*math.sin(Tij)-b*math.cos(Tij))


class Current:
    def __init__(self, Top):
        dim = Top.dimension
        self.Current = matrix(dim)

    def ActualizeCurrent(self, Top, V):
        dim = Top.dimension
        for i in range(dim):
            for j in range(dim):
                g = Top.Y.matrix[i][j].real
                b = Top.Y.matrix[i][j].imag
                Vi = V.Magnitude.matrix[i][0]
                Vj = V.Magnitude.matrix[j][0]
                Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
                self.Current.matrix[i][j] = (((g**2)+(b**2))*((Vi**2)+(Vj**2)-2*Vi*Vj*math.cos(Tij)))**0.5


'''
Funções de derivadas para a matriz jacobiana
Base para o estimador de Estado
'''


def dPijdTk(Top, V, i, j, k):
    dim = Top.dimension
    if i == j:
        if i == k: # dPi/dTi
            dPijTk = 0
            for t in range(dim):
                G = Top.Ybus.matrix[i][t].real
                B = Top.Ybus.matrix[i][t].imag
                Tit = V.Phase.matrix[i][0] - V.Phase.matrix[t][0]
                Vi = V.Magnitude.matrix[i][0]
                Vt = V.Magnitude.matrix[t][0]
                dPijTk += Vt*(-G*math.sin(Tit)+B*math.cos(Tit))
            dPijTk = dPijTk*Vi-(Vi**2)*Top.Ybus.matrix[i][i].imag
            #print('\ndP', i, '/dT', i, ' = ', dPijTk)
        else: # dPi/dTj
            G = Top.Ybus.matrix[i][k].real
            B = Top.Ybus.matrix[i][k].imag
            Tik = V.Phase.matrix[i][0] - V.Phase.matrix[k][0]
            Vi = V.Magnitude.matrix[i][0]
            Vk = V.Magnitude.matrix[k][0]
            dPijTk = Vi*Vk*(G*math.sin(Tik)-B*math.cos(Tik))
            #print('\ndP', i, '/dT', j, ' = ', dPijTk)
    else:
        if k == i: # dPij/dTi
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dPijTk = Vi*Vj*(g*math.sin(Tij)-b*math.cos(Tij))
            #print('\ndP', i, j, '/dT', i, ' = ', dPijTk)
        elif k == j: # dPij/dTj
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dPijTk = -Vi*Vj*(g*math.sin(Tij)-b*math.cos(Tij))
            #print('\ndP', i, j, '/dT', j, ' = ', dPijTk)
        else:
            dPijTk = 0
            #print('\ndP', i, j, '/dT', k, ' = ', dPijTk)
    return dPijTk


def dQijdTk(Top, V, i, j, k):
    dim = Top.dimension
    if i == j:
        if i == k: # dQi/dTi
            dQijTk = 0
            for t in range(dim):
                G = Top.Ybus.matrix[i][t].real
                B = Top.Ybus.matrix[i][t].imag
                Tit = V.Phase.matrix[i][0] - V.Phase.matrix[t][0]
                Vi = V.Magnitude.matrix[i][0]
                Vt = V.Magnitude.matrix[t][0]
                dQijTk += Vt*(G*math.cos(Tit)+B*math.sin(Tit))
            dQijTk = dQijTk*Vi-(Vi**2)*Top.Ybus.matrix[i][i].real
            #print('\ndQ', i, '/dT', i, ' = ', dQijTk)
        else: # dQi/dTj
            G = Top.Ybus.matrix[i][k].real
            B = Top.Ybus.matrix[i][k].imag
            Tik = V.Phase.matrix[i][0] - V.Phase.matrix[k][0]
            Vi = V.Magnitude.matrix[i][0]
            Vk = V.Magnitude.matrix[k][0]
            dQijTk = Vi*Vk*(-G*math.cos(Tik)-B*math.sin(Tik))
            #print('\ndQ', i, '/dT', j, ' = ', dQijTk)
    else:
        if k == i: # dQij/dTi
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dQijTk = -Vi*Vj*(g*math.cos(Tij)+b*math.sin(Tij))
            #print('\ndQ', i, j, '/dT', i, ' = ', dQijTk)
        elif k == j: # dQij/dTj
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dQijTk = Vi*Vj*(g*math.cos(Tij)+b*math.sin(Tij))
            #print('\ndQ', i, j, '/dT', j, ' = ', dQijTk)
        else:
            dQijTk = 0
            #print('\ndQ', i, j, '/dT', k, ' = ', dQijTk)
    return dQijTk


def dPijdVk(Top, V, i, j, k):
    dim = Top.dimension
    if i == j:
        if i == k: # dPi/dVi
            dPijVk = 0
            for t in range(dim):
                G = Top.Ybus.matrix[i][t].real
                B = Top.Ybus.matrix[i][t].imag
                Tit = V.Phase.matrix[i][0] - V.Phase.matrix[t][0]
                Vi = V.Magnitude.matrix[i][0]
                Vt = V.Magnitude.matrix[t][0]
                dPijVk += Vt*(G*math.cos(Tit)+B*math.sin(Tit))
            dPijVk += Vi*Top.Ybus.matrix[i][i].real
            #print('\ndP', i, '/dV', i, ' = ', dPijVk)
        else: # dPi/dVj
            G = Top.Ybus.matrix[i][k].real
            B = Top.Ybus.matrix[i][k].imag
            Tik = V.Phase.matrix[i][0] - V.Phase.matrix[k][0]
            Vi = V.Magnitude.matrix[i][0]
            Vk = V.Magnitude.matrix[k][0]
            dPijVk = Vi*(G*math.cos(Tik)+B*math.sin(Tik))
            #print('\ndP', i, '/dV', j, ' = ', dPijVk)
    else:
        if k == i: # dPij/dVi
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dPijVk = -Vj*(g*math.cos(Tij)+b*math.sin(Tij))+2*Vi*(g+gsh)
            #print('\ndP', i, j, '/dV', i, ' = ', dPijVk)
        elif k == j: # dPij/dVj
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dPijVk = -Vi*(g*math.cos(Tij)+b*math.sin(Tij))
            #print('\ndP', i, j, '/dV', j, ' = ', dPijVk)
        else:
            dPijVk = 0
            #print('\ndP', i, j, '/dV', k, ' = ', dPijVk)
    return dPijVk


def dQijdVk(Top, V, i, j, k):
    dim = Top.dimension
    if i == j:
        if i == k: # dQi/dVi
            dQijdVk = 0
            for t in range(dim):
                G = Top.Ybus.matrix[i][t].real
                B = Top.Ybus.matrix[i][t].imag
                Tit = V.Phase.matrix[i][0] - V.Phase.matrix[t][0]
                Vi = V.Magnitude.matrix[i][0]
                Vt = V.Magnitude.matrix[t][0]
                dQijdVk += Vt*(G*math.sin(Tit)-B*math.cos(Tit))
            dQijdVk -= Vi*Top.Ybus.matrix[i][i].imag
            #print('\ndQ', i, '/dV', i, ' = ', dQijdVk)
        else: # dQi/dVj
            G = Top.Ybus.matrix[i][k].real
            B = Top.Ybus.matrix[i][k].imag
            Tik = V.Phase.matrix[i][0] - V.Phase.matrix[k][0]
            Vi = V.Magnitude.matrix[i][0]
            Vk = V.Magnitude.matrix[k][0]
            dQijdVk = Vi*(G*math.sin(Tik)-B*math.cos(Tik))
            #print('\ndQ', i, '/dV', j, ' = ', dQijdVk)
    else:
        if k == i: # dQij/dVi
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dQijdVk = -Vi*(g*math.sin(Tij)-b*math.cos(Tij))-2*Vi*(b+bsh)
            #print('\ndQ', i, j, '/dV', i, ' = ', dQijdVk)
        elif k == j: # dQij/dVj
            g = Top.Y.matrix[i][j].real
            b = Top.Y.matrix[i][j].imag
            gsh = Top.Ysh.matrix[i][j].real
            bsh = Top.Ysh.matrix[i][j].imag
            Tij = V.Phase.matrix[i][0] - V.Phase.matrix[j][0]
            Vi = V.Magnitude.matrix[i][0]
            Vj = V.Magnitude.matrix[j][0]
            dQijdVk = -Vi*(g*math.sin(Tij)-b*math.cos(Tij))
            #print('\ndQ', i, j, '/dV', j, ' = ', dQijdVk)
        else:
            dQijdVk = 0
            #print('\ndQ', i, j, '/dV', k, ' = ', dQijdVk)
    return dQijdVk


def dVidTj(Top, V, i, j):
    return 0


def dVidVj(Top, V, i, j):
    if i == j: # dVi/dVi
        dVidVj = 1
    else: # dVi/dVj
        dVidVj = 0
    return dVidVj


# Jacobian Function to support State Estimation
def Jacobian(Top, V, SetMeas):
    dim = Top.dimension
    H = []
    Sigma = []
    z = []
    line = -1
    for each_meas in SetMeas:
        if each_meas.type == 'Power Flow Measure' and each_meas.Pm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if not i == j:
                line += 1
                H.append([])
                Sigma.append(each_meas.Psigma**(-2))
                z.append([each_meas.Pm])
                for k in range(dim):
                    H[line].append(dPijdTk(Top, V, i, j, k))
                for k in range(dim):
                    H[line].append(dPijdVk(Top, V, i, j, k))
    for each_meas in SetMeas:
        if each_meas.type == 'Power Injection Measure' and each_meas.Pm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if i == j:
                line += 1
                H.append([])
                Sigma.append(each_meas.Psigma**(-2))
                z.append([each_meas.Pm])
                for k in range(dim):
                    H[line].append(dPijdTk(Top, V, i, j, k))
                for k in range(dim):
                    H[line].append(dPijdVk(Top, V, i, j, k))
    for each_meas in SetMeas:
        if each_meas.type == 'Power Flow Measure' and each_meas.Qm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if not i == j:
                line += 1
                H.append([])
                Sigma.append(each_meas.Qsigma**(-2))
                z.append([each_meas.Qm])
                for k in range(dim):
                    H[line].append(dQijdTk(Top, V, i, j, k))
                for k in range(dim):
                    H[line].append(dQijdVk(Top, V, i, j, k))
    for each_meas in SetMeas:
        if each_meas.type == 'Power Injection Measure' and each_meas.Qm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if i == j:
                line += 1
                H.append([])
                Sigma.append(each_meas.Qsigma**(-2))
                z.append([each_meas.Qm])
                for k in range(dim):
                    H[line].append(dQijdTk(Top, V, i, j, k))
                for k in range(dim):
                    H[line].append(dQijdVk(Top, V, i, j, k))
    for each_meas in SetMeas:
        if each_meas.type == 'Voltage Measure' and each_meas.Vm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            line += 1
            H.append([])
            Sigma.append(each_meas.Vsigma**(-2))
            z.append([each_meas.Vm])
            for k in range(dim):
                H[line].append(dVidTj(Top, V, i, k))
            for k in range(dim):
                H[line].append(dVidVj(Top, V, i, k))
    for i in range(line+1):
        H[i].pop(0)
    Rn = matrix(line+1)
    for i in range(line+1):
        Rn.matrix[i][i] = Sigma[i] # Create R^-1 matrix with residual values
    HH = matrix(line+1, 2*dim-1)
    HH.matrix = H
    zz = matrix(line+1, 1)
    for i in range(line+1):
        zz.matrix[i][0] = z[i][0]
    return HH, Rn, zz


# Function to support State Estimation that calculate a vector with a state of system
def StateVector(Top, V, SetMeas, Power):
    dim = Top.dimension
    h = []
    line = -1
    for each_meas in SetMeas:
        if each_meas.type == 'Power Flow Measure' and each_meas.Pm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if not i == j:
                line += 1
                h.append([Power.RealPower.matrix[i][j]])
    for each_meas in SetMeas:
        if each_meas.type == 'Power Injection Measure' and each_meas.Pm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if i == j:
                line += 1
                h.append([Power.RealPower.matrix[i][j]])
    for each_meas in SetMeas:
        if each_meas.type == 'Power Flow Measure' and each_meas.Qm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if not i == j:
                line += 1
                h.append([Power.ReacPower.matrix[i][j]])
    for each_meas in SetMeas:
        if each_meas.type == 'Power Injection Measure' and each_meas.Qm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if i == j:
                line += 1
                h.append([Power.ReacPower.matrix[i][j]])
    for each_meas in SetMeas:
        if each_meas.type == 'Voltage Measure' and each_meas.Vm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            line += 1
            h.append([V.Magnitude.matrix[i][0]])
    hh = matrix(line+1, 1)
    for i in range(line+1):
        hh.matrix[i][0] = h[i][0]
    return hh


'''
Functions to support the decoupled solution method to state estimator
'''
# Jacobian Real Function to support State Estimation
def JacobianAA(Top, V, SetMeas):
    dim = Top.dimension
    Haa = []
    Sigmaa = []
    za = []
    line = -1
    for each_meas in SetMeas:
        if each_meas.type == 'Power Flow Measure' and each_meas.Pm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if not i == j:
                line += 1
                Haa.append([])
                Sigmaa.append(each_meas.Psigma**(-2))
                za.append([each_meas.Pm])
                for k in range(dim):
                    Haa[line].append(dPijdTk(Top, V, i, j, k))
    for each_meas in SetMeas:
        if each_meas.type == 'Power Injection Measure' and each_meas.Pm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if i == j:
                line += 1
                Haa.append([])
                Sigmaa.append(each_meas.Psigma**(-2))
                za.append([each_meas.Pm])
                for k in range(dim):
                    Haa[line].append(dPijdTk(Top, V, i, j, k))
    for i in range(line+1):
        Haa[i].pop(0)
    HHaa = matrix(line+1, dim-1)
    HHaa.matrix = Haa
    Ran = matrix(line+1)
    zza = matrix(line+1, 1)
    for i in range(line+1):
        Ran.matrix[i][i] = Sigmaa[i] # Create Ra^-1 matrix with residual values
        zza.matrix[i][0] = za[i]
    return HHaa, Ran, zza


# Jacobian Reactive Function to support State Estimation
def JacobianRR(Top, V, SetMeas):
    dim = Top.dimension
    Hrr = []
    Sigmar = []
    zr = []
    line = -1
    for each_meas in SetMeas:
        if each_meas.type == 'Power Flow Measure' and each_meas.Qm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if not i == j:
                line += 1
                Hrr.append([])
                Sigmar.append(each_meas.Qsigma**(-2))
                zr.append([each_meas.Qm])
                for k in range(dim):
                    Hrr[line].append(dQijdVk(Top, V, i, j, k))
    for each_meas in SetMeas:
        if each_meas.type == 'Power Injection Measure' and each_meas.Qm:
            i = each_meas.id[0]-1
            j = each_meas.id[1]-1
            if i == j:
                line += 1
                Hrr.append([])
                Sigmar.append(each_meas.Qsigma**(-2))
                zr.append([each_meas.Qm])
                for k in range(dim):
                    Hrr[line].append(dQijdVk(Top, V, i, j, k))
    HHrr = matrix(line+1, dim)
    HHrr.matrix = Hrr
    Rrn = matrix(line+1)
    zzr = matrix(line+1, 1)
    for i in range(line+1):
        Rrn.matrix[i][i] = Sigmar[i] # Create Rr^-1 matrix with residual values
        zzr.matrix[i][0] = zr[i]
    return HHrr, Rrn, zzr
