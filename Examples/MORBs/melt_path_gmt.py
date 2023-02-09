import warnings

import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
import shapely.geometry as shp

import pyMelt as m

Tp=float(sys.argv[1])
Tp_max=float(sys.argv[1])+float(sys.argv[2])
Tp_min=float(sys.argv[1])-float(sys.argv[3])

lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 6., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

path = mantle.adiabaticMelt(Tp, Pstart=6, dP=-0.01)
path_max = mantle.adiabaticMelt(Tp_max, Pstart=6, dP=-0.01)
path_min = mantle.adiabaticMelt(Tp_min, Pstart=6, dP=-0.01)

sol_ar = np.column_stack((P_sol, T_sol))
Tp_ar = np.column_stack((path.P, path.T, path.F))
Tp_ar_max = np.column_stack((path_max.P, path_max.T, path_max.F))
Tp_ar_min = np.column_stack((path_min.P, path_min.T, path_min.F))

np.savetxt('sol.temp', sol_ar)
np.savetxt('Tp_ar.temp', Tp_ar)
np.savetxt('Tp_ar_max.temp', Tp_ar_max)
np.savetxt('Tp_ar_min.temp', Tp_ar_min)
