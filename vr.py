def key():
    return ""

def historical(ticker):
    return "https://financialmodelingprep.com/stable/historical-price-eod/full?symbol=" + ticker + "&apikey=" + key()

import numpy as np 
import requests
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ctypes
import time

def TimeIt(f):
    def Handle(*a, **b):
        t0 = int(time.time())
        z = f(*a, **b)
        t1 = int(time.time())
        ds = t1 - t0
        return ds
    return Handle

@TimeIt
def C_GBM(S, drift, vol, t, steps, paths):
    quant = ctypes.CDLL("./vr.so")
    quant.simulation.argtypes = (
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int
    )
    quant.simulation.restype = ctypes.c_double
    steps = int(steps)
    paths = int(paths)
    result = quant.simulation(S, drift, vol, t, steps, paths)
    return result

@TimeIt
def Py_GBM(S, drift, vol, t, steps, paths):
    steps = int(steps)
    paths = int(paths)
    dWT = np.random.randn(steps)
    total = 0
    dt = t / steps
    for i in range(paths):
        S0 = S
        for j in range(steps):
            S0 += drift*S0*dt + vol*S0*dWT[j]
        total += S0
    return total / paths

ticker = 'NVDA'
close = np.array([i['close'] for i in requests.get(historical(ticker)).json()])[::-1]

ror = np.log(close[1:]/close[:-1])

S = close[-1]
t = 1.0/252.0
drift = np.mean(ror)
vol = np.std(ror)

N = 10
steps = np.linspace(100, 10000, N)
paths = np.linspace(50, 1000, N)

x, y = np.meshgrid(paths, steps)
z_py = np.zeros((N, N))
z_c = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        py_result = Py_GBM(S, drift, vol, t, y[i, j], x[i, j])
        c_result = C_GBM(S, drift, vol, t, y[i, j], x[i, j])
        z_py[i, j] = py_result
        z_c[i, j] = c_result
        print(i, j)

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ay = fig.add_subplot(122, projection='3d')

ax.set_title('Seconds Model Ran: Python')
ay.set_title('Seconds Model Ran: C')

ax.set_xlabel('Paths')
ax.set_ylabel('Steps')
ay.set_xlabel('Paths')
ay.set_ylabel('Steps')

the_python = ax.plot_surface(x, y, z_py, cmap='jet_r')
the_c = ay.plot_surface(x, y, z_c, cmap='jet_r')

fig.colorbar(the_python, ax=ax)
fig.colorbar(the_c, ax=ay)

plt.show()

