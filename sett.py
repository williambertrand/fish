''' THIS IS THE FISHE SETTINGS FILE '''
import numpy as np

def init():
    np.set_printoptions(precision = 1, linewidth=220, suppress = True)

    global ROW, COL, V_FLUID, T_DIFFU, RE, U_INF, U_MAX, H, DT, RELAX_FACTOR, TOL
    (ROW, COL)   = (24, 25)
    V_FLUID      = 0.8927
    T_DIFFU      = 0.143e-6
    RE           = 8
    U_INF        = 8
    U_MAX        = 8
    H            = (RE * V_FLUID) / U_INF
    DT           = 0.3 * (H / U_MAX)
    RELAX_FACTOR = 1.5/4.0
    TOL          = 0.001
