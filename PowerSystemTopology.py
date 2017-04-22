'''
Noé de Lima Bezerra
30/12/2016 - 2h53

Método para obter a topologia da rede elétrica
O objetivo é calcular a matriz Ybus do sistema e armazenar todos os parâmetros
de topologia necessários para calcular as variáveis de estado,
bem como fornecer suporte aos módulos de estimação de estado e análise de
observabilidade de rede.
'''

from Matrix import *

#Create a network element class
#Each element can be a node or a line
class net_element:
    def __init__(self, fr_line, to_line=None, r_serie=0, x_serie=0, y_shunt=0, y_serie=None):
        if to_line and to_line != fr_line:
            self.type = 'Line'
            self.id = [fr_line, to_line]
        else:
            self.type = 'Node'
            self.id = [fr_line, fr_line]
        self.Rserie = r_serie
        self.Xserie = x_serie
        self.Zserie = r_serie+x_serie*1j
        if y_serie:
            self.Yserie = y_serie
        else:
            try:
                self.Yserie = 1/self.Zserie
            except:
                self.Yserie = 0
        self.Yshunt = y_shunt

        
class topology(net_element):
    def __init__(self, SetNet):
        self.dimension = 0
        for each_element in SetNet:
            if each_element.type == 'Node':
                self.dimension += 1
        self.Y = matrix(self.dimension) #Create self admittance Y matrix with zeros to number of nodes in set of network elements
        self.Ybus = matrix(self.dimension) #Create Ybus matrix with zeros to number of nodes in set of network elements
        self.Ysh = matrix(self.dimension)
        for each_element in SetNet:
            i = each_element.id[0]-1
            j = each_element.id[1]-1
            if each_element.type == 'Node':
                self.Y.matrix[i][j] += each_element.Yshunt
                self.Ysh.matrix[i][j] = each_element.Yshunt
            if each_element.type == 'Line':
                self.Y.matrix[i][j] = each_element.Yserie
                self.Y.matrix[j][i] = each_element.Yserie
                self.Y.matrix[i][i] += each_element.Yshunt
                self.Y.matrix[j][j] += each_element.Yshunt
                self.Ysh.matrix[i][j] = each_element.Yshunt
                self.Ysh.matrix[j][i] = each_element.Yshunt
                
        for i in range(self.dimension):
            for j in range(self.dimension):
                if i == j:
                    for k in range(self.dimension):
                        self.Ybus.matrix[i][j] += self.Y.matrix[i][k]
                else:
                    self.Ybus.matrix[i][j] = -self.Y.matrix[i][j]

class measurement:
    def __init__(self, fr, to=None, p_measure=None, p_sigma=None, q_measure=None, q_sigma=None, v_measure=None, v_sigma=None):
        if to and to != fr:
            self.type = 'Power Flow Measure'
            self.Pm = p_measure
            self.Psigma = p_sigma
            self.Qm = q_measure
            self.Qsigma = q_sigma
            self.id = [fr, to]
        elif v_measure:
            self.type = 'Voltage Measure'
            self.Vm = v_measure
            self.Vsigma = v_sigma
            self.id = [fr, fr]
        else:
            self.type = 'Power Injection Measure'
            self.Pm = p_measure
            self.Psigma = p_sigma
            self.Qm = q_measure
            self.Qsigma = q_sigma
            self.id = [fr, fr]

class PowerSet(measurement):
    def __init__(self, SetMeas, Topology):
        self.PowerNum = 0
        self.VoltageNum = 0
        for each_meas in SetMeas:
            if each_meas.type == 'Power Flow Measure':
                self.PowerNum += 1
            elif each_meas.type == 'Power Injection Measure':
                self.PowerNum += 1
            else:
                self.VoltageNum += 1
        dim = Topology.dimension
        self.Power = matrix(dim)
        self.PowerSigma = matrix(dim)
        self.Voltage = matrix(dim, 1)
        self.VoltageSigma = matrix(dim, 1)
        for measure in SetMeas:
            i = measure.id[0]-1
            j = measure.id[1]-1
            if measure.type == 'Power Flow Measure' or measure.type == 'Power Injection Measure':
                try:
                    self.Power.matrix[i][j] += measure.Pm
                    self.PowerSigma.matrix[i][j] += measure.Psigma
                except:
                    pass
                try:
                    self.Power.matrix[i][j] += measure.Qm*1j
                    self.PowerSigma.matrix[i][j] += measure.Qsigma*1j
                except:
                    pass
            if measure.type == 'Voltage Measure':
                self.Voltage.matrix[i][0] = measure.Vm
                self.VoltageSigma.matrix[i][0] = measure.Vsigma
                
def valid_set_measures(SetMeas, Topology):
    dim = Topology.dimension
    valid = 0
    for each_measure in SetMeas:
        i = each_measure.id[0]-1
        j = each_measure.id[1]-1
        if not (i < dim and j < dim):
            print('Measurement: ' + each_measure.type + str(each_measure.id) + ' is out of Power System set')
            valid += 1
    if not valid == 0:
        return(False)
    else:
        return(True)

