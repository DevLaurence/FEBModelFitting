import os
import time
import re
import numpy as np
import matplotlib.pyplot as plt
import fileinput
from scipy.optimize import differential_evolution

def follow(filename):
  file = open(filename)
  text = file.read()
  if "N O R M A L   T E R M I N A T I O N" in text:
    return 1
  elif "E R R O R   T E R M I N A T I O N" in text:
    return 0
  else:
    file.close()
    time.sleep(0.05)

def UpdateFEB(filein, fileout, parameters):

  with open(filein, 'r') as file:
    data_original = file.read()

  # data_original = data_original.replace('<c>5.21</c>', '<c>'+str(parameters[0])+'</c>')
  data_original = data_original.replace('<k1>4.05</k1>', '<k1>'+str(parameters[0])+'</k1>')
  data_original = data_original.replace('<k2>10.90</k2>', '<k2>'+str(parameters[1])+'</k2>')
  data_original = data_original.replace('<kappa>0.281</kappa>', '<kappa>'+str(parameters[2])+'</kappa>')
  data_original = data_original.replace('<gamma>5.65</gamma>', '<gamma>'+str(parameters[3])+'</gamma>')
  data_original = data_original.replace('<alpha>1.378</alpha>', '<alpha>'+str(parameters[4])+'</alpha>')
  data_original = data_original.replace('<mu1>0.406</mu1>', '<mu1>'+str(parameters[5])+'</mu1>')
  data_original = data_original.replace('<scale>0.060</scale>', '<scale>'+str(parameters[6])+'</scale>')

  with open(fileout, 'w') as file:
    file.write(data_original)

def RunMonitorFEB(filein, fileout):

  os.system('febio4 -i '+filein+' >/dev/null 2>&1')

  FEB_flag = follow(fileout)

  return FEB_flag

def AnalyzeFEB(filein, data_compare):
  data_out = np.zeros((num_time,5))
  counter_out = 0

  with open(filein, 'r') as file:
    for line in file:
      if "Data Record" in line:
        next(file, '')
        next(file, '')
        temp1 = re.findall("\d*\.?\d+", next(file, ''))
        data_out[counter_out, 0] = float(temp1[0])

        next(file, '')
        temp2 = re.findall("\d+\.?\d*(?:[Ee]\-?\+?\d+)?", next(file, ''))
        data_out[counter_out, 1] = float(temp2[1])
        data_out[counter_out, 2] = float(temp2[2])
        data_out[counter_out, 3] = float(temp2[3])
        data_out[counter_out, 4] = float(temp2[4])

        counter_out += 1

  diffx = (data_out[:, 3] - data_compare[:, 1])**2
  diffy = (data_out[:, 4] - data_compare[:, 2])**2

  residual = np.sqrt(diffx.sum() + diffy.sum())

  residual = residual / num_time

  return residual, data_out

def ResidualFEB(params):

    # parameters = [5.21, 12.548, 48.458, 0.15, 0]

    data_exp = np.loadtxt('jobs/experimentSR2.txt')

    FEB_template = 'jobs/BiaxialSRBiaxial.feb'
    FEB_new = 'jobs/BiaxialSRBiaxialModify.feb'
    FEB_in = 'jobs/BiaxialSRBiaxialModify.feb'
    FEB_out = 'jobs/outputSR.txt'

    UpdateFEB(FEB_template, FEB_new, params)

    flag = RunMonitorFEB(FEB_new, FEB_out)

    if flag == 1:
        [residual, data_out] = AnalyzeFEB(FEB_out, data_exp)
    else:
        residual = 1e6

    return residual

# Convert old MATLAB script to Python
# def DifferentialEvolution():



if __name__ == '__main__':

    num_time = 201

    # bounds = [(1e-6, 10), (1e-6, 50), (1e-6, 50), (1e-6, 0.32)]

    # bounds = [{'name': 'k1', 'type': 'continuous', 'domain': (0, 50)},
    #           {'name': 'k2', 'type': 'continuous', 'domain': (0, 50)},
    #           {'name': 'kappa', 'type': 'continuous', 'domain': (0, 0.33)},
    #           {'name': 'gamma', 'type': 'continuous', 'domain': (0, 90)}]
    #           {'name': 'alpha', 'type': 'continuous', 'domain': (0.3, 2)},
    #           {'name': 'mu1', 'type': 'continuous', 'domain': (0, 0.5)},
    #           {'name': 'scale', 'type': 'continuous', 'domain': (0.01, 0.08)}]

    bounds = [(1e-6, 50), (1e-6, 50), (1e-6, 0.33), (0, 90), (0.2, 2), (1e-6, 0.5), (0.01, 0.08)]

    result = differential_evolution(ResidualFEB, bounds, disp=True, popsize=5, updating='deferred', recombination=0.4,\
                                    mutation=0.3, maxiter=25)
    print(result.x)
    print(result.fun)

    # num_time = 10
    #
    # parameters = [5.21, 12.548, 48.458, 0.15, 0]
    #
    data_exp = np.loadtxt('jobs/experimentSR2.txt')

    FEB_template = 'jobs/BiaxialSRBiaxial.feb'
    FEB_new = 'jobs/BiaxialSRBiaxialModify.feb'
    FEB_in = 'jobs/BiaxialSRBiaxialModify.feb'
    FEB_out = 'jobs/outputSR.txt'

    UpdateFEB(FEB_template, FEB_new, result.x)

    RunMonitorFEB(FEB_new, FEB_out)

    [residual, data_out] = AnalyzeFEB(FEB_out, data_exp)

    data_exp = np.insert(data_exp, 0, [1, 0, 0], axis=0)
    data_out = np.insert(data_out, 0, [0, 1, 1, 0, 0], axis=0)

    print(data_exp)
    print(data_out)

    fig = plt.figure()
    fig.show()
    ax = fig.add_subplot(111)

    ax.plot(data_out[:, 0], data_exp[:, 1], c='k', marker="o", label='Exp. (circ.)')
    ax.plot(data_out[:, 0], data_exp[:, 2], c='k', marker="s", label='Exp. (rad.)')
    ax.plot(data_out[:, 0], data_out[:, 3], c='b', ls='-', label='Model (circ.)')
    ax.plot(data_out[:, 0], data_out[:, 4], c='r', ls='-', label='Model (rad)')

    plt.legend()
    plt.draw()

    plt.savefig('fig/debug.tiff', dpi=75)



