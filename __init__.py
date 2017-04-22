'''
Noé de Lima Bezerra
30/12/2016 - 2h53

Funções de Estimação do Estado do Sistema Elétrico de Potência
Derivadas das funções de estado para a matriz jacobiana
'''

import math
import io
import time

from Matrix import *
from PowerSystemState import *
from PowerSystemStateEstimator import *

# Defining external archives names
# This values shell be obtained in execution time, or as an argument to call this program
# And will supply all parameters of input and output to final operation as this program have to do
Network = input("What is the network topology csv archive's name? ") + '.csv'  # The name of the CSV file with the parameters of bus and branches of the network
Measure = input("What is the measurement set csv archive's name? ") + '.csv'  # The name of the CSV file that contain all measurements values and location id to solve all
Out = 'PSOutput.txt'  # The name of the TXT file to write the results and log of this program; this can be changed later

# Importing System Topology from an external archive
PSNetwork = []
try:
    with open(Network) as PSN:
        print(PSN.readline())
        for branch in PSN:
            element = branch.split(';')
            for i in range(len(element)):
                try:
                    element[i] = float(element[i])
                except:
                    element[i] = 0
            element[0] = int(element[0])
            element[1] = int(element[1])
            if element[5] == 0:
                element[5] = None

            print('Line = ', element)
            if not element[0] == 'From':
                PSNetwork.append(net_element(element[0], element[1], r_serie=element[2], x_serie=element[3], y_shunt=element[4], y_serie=element[5]))
except IOError as err:
    print('File Error in Network Topology: ' + str(err))

# Importing the set of Measurements of another external archive
PSMeasure = []
try:
    with open(Measure) as PSM:
        print(PSM.readline())
        for measure in PSM:
            element = measure.split(';')
            for i in range(len(element)):
                try:
                    element[i] = float(element[i])
                except:
                    element[i] = None
            element[0] = int(element[0])
            if not element[1]:
                element[1] = element[0]
            element[1] = int(element[1])
            if element[5] == 0:
                element[5] = None

            print('Line = ', element)
            if not element[0] == 'From':
                PSMeasure.append(
                    measurement(element[0], element[1], p_measure=element[2], p_sigma=element[3], q_measure=element[4],
                                q_sigma=element[5], v_measure=element[6], v_sigma=element[7]))
except IOError as err:
    print('File Error in Measurement Set: ' + str(err))

key = valid_set_measures(PSMeasure, topology(PSNetwork))

if key:
    print('This System, in principle, is observable, but may not converge.')
else:
    try:
        with open(Out, 'a') as Output:
            print('\n\nThis system is unobservable.\nPleas, try with another measurement set configuration!\n\n', file=Output)
    except IOError as err:
        print('File Error in Results Output: ' + str(err))
    finally:
        exit()

# Asking for the tolerance of convergence
tolerance = 20  #int(input('Set an integer as default tolerance like 10^-'))

check = 0 #Define what to do next
count1 = 0 #count how many times using the same tolerance to increase
count2 = 0 #count how many times using the complet Newton-Raphson method to return to decoupled method

# Here is a State Estimator Call
while check == 0 or check == 1:
    print('\n', time.asctime(time.localtime()), '\n')
    print(count1, 'calls did to actual tolerance', tolerance, ', using this option:', end='')
    if check == 0:
        print('Fast Decoupled Method')
    else:
        print('Fast Newton-Raphson Copled Method')

    if check == 0:
        try:
            V, Power, I, newtolerance = StateEstimatorDC(PSNetwork, PSMeasure, tolerance)
            for i in range(V.Phase.line):
                V.Phase.matrix[i][0] = V.Phase.matrix[i][0] * 180/math.pi
            if tolerance == newtolerance:
                count1 += 1
            else:
                tolerance = newtolerance
                count1 = 0
            if count1 > 10:
                tolerance += 1
                count1 = 0
        except:
            check = 1
            try:
                with open(Out, 'a') as Output:
                    print('\n\n', time.asctime(time.localtime()), '\n\nThe system state was not converge to this '
                        'tolerance:', tolerance, ', using decoupled method!\nTrying complete Newton-Raphson method\n\n',
                        file=Output)
            except IOError as err:
                print('File Error in Results Output: ' + str(err))
    if check == 1:
        try:
            V, Power, I, newtolerance = StateEstimator(PSNetwork, PSMeasure, tolerance)
            for i in range(V.Phase.line):
                V.Phase.matrix[i][0] = V.Phase.matrix[i][0] * 180/math.pi
            if tolerance == newtolerance:
                count1 += 1
            else:
                tolerance = newtolerance
                count1 = 0
            if count1 > 10:
                tolerance += 1
                count1 = 0
            if count2 > 10:
                check = 0
                count2 = 0
            else:
                count2 += 1
        except:
            check = 2
            try:
                with open(Out, 'a') as Output:
                    print('\n\n', time.asctime(time.localtime()), '\n\nThe system state was not converge to this '
                        'tolerance:', tolerance, ', using coupled method!\n\n', file=Output)
            except IOError as err:
                print('File Error in Results Output: ' + str(err))

    # Here is the results
    if check == 0 or check == 1:
        try:
            with open(Out, 'a') as Output:
                print('\n\n', time.asctime(time.localtime()), '\n\n', count1, 'th call and the actual tolerance '
                    'value is', tolerance, '\n\n', file=Output)
                print('The system state variable is:\n', file=Output)
                print('\n\nThe Voltage Magnitude is:\n', file=Output)
                V.Magnitude.printmatrix(fh=Output)
                print('\n\nThe Voltage Phase Angle is:\n', file=Output)
                V.Phase.printmatrix(fh=Output)
                print('\n\nAnd the system state is:\n', file=Output)
                print('\n\nThe Real Power Flow is:\n', file=Output)
                Power.RealPower.printmatrix(fh=Output)
                print('\n\nThe Reactive Power Flow is:\n', file=Output)
                Power.ReacPower.printmatrix(fh=Output)
                print('\n\nThe Current flow is:\n', file=Output)
                I.Current.printmatrix(fh=Output)
        except IOError as err:
            print('File Error in Results Output: ' + str(err))

    print('Choose what to do:\n0: restart estimation using Decoupled Method\n1: restart estimation using Coupled '
        'Method')
    check = input('Press anything to decide what to do next:')
    try:
        check = int(check)
    except:
        check = 2

#StateEstimatorDC(PSNetwork, PSMeasure, tolerance)
#check = input('Press anything to exit process:')
