import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def theta(x):
    return np.tanh(1000 * x) if x > 0 else 0

def model(y, t, m, d, r, k, c, a0, l, Ps):
    x1, x2, x3, v1, v2, v3 = y
    
    # Calculo las áreas
    a = {}
    a['1'] = a0['1'] + 2 * l * x1
    a['2'] = a0['2'] + 2 * l * x2
    a['3'] = a0['3'] + 2 * l * x3
    
    # Calculo el área mínima
    a_min = min(a['1'], a['2'], a['3'])
    
    # Calculo las presiones
    P = {}
    P['1'] = Ps * (1 - theta(a_min) * (a_min / a['1']) ** 2) * theta(a['1'])
    P['2'] = Ps * (1 - theta(a_min) * (a_min / a['2']) ** 2) * theta(a['1']) * theta(a['2']) * \
             theta(a['1'] - a['3']) * theta(a['2'] - a['3'])
    P['3'] = 0
    
    # Ecuaciones diferenciales
    dydt = [
        v1,
        v2,
        v3,
        (l * d['1'] * P['1'] - r['1'] * v1 - k['1'] * x1 - theta(-a['1']) * c['1'] * (a['1'] / (2 * l)) - \
         k['1,2'] * (x1 - x2)) / m['1'],
        (l * d['2'] * P['2'] - r['2'] * v2 - k['2'] * x2 - theta(-a['2']) * c['2'] * (a['2'] / (2 * l)) - \
         k['1,2'] * (x2 - x1) - k['2,3'] * (x2 - x3)) / m['2'],
        (l * d['3'] * P['3'] - r['3'] * v3 - k['3'] * x3 - theta(-a['3']) * c['3'] * (a['3'] / (2 * l)) - \
         k['2,3'] * (x3 - x2)) / m['3']
    ]
    
    return dydt