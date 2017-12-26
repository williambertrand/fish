###############################################################################
## this is the fishe functions file
###############################################################################

###############################################################################
## imports
###############################################################################
import numpy as np
import sett
import time

###############################################################################
## build array containing node flags
###############################################################################
def buildDomain():
    domain = np.zeros((sett.ROW, sett.COL))

    domain[:         ,          0]    =    1    ## INLET
    domain[:         , sett.COL-1]    =   -1    ## OUTLET
    domain[0         ,          :]    =    2    ## TOP
    domain[sett.ROW-1,          :]    =   -2    ## BOTTOM

    ## design obstructing geometry
    domain[8:15      ,        7:10]   =    7    ## CONSTANT CORE

    domain[8:16      ,          6]    =    3    ## LEFT SIDE
    domain[8:16      ,         10]    =    5    ## RIGHT SIDE
    domain[8         ,       7:10]    =    4    ## TOP
    domain[15        ,       7:10]    =    6    ## BOTTOM

    return domain

###############################################################################
## psi one
###############################################################################
def initPsi(domain):
    psi = np.copy(domain)

    psi[domain ==  1] = sett.U_INF * sett.H * sett.COL
    psi[domain == -1] = 0

    psi[domain ==  2] = 0
    psi[domain == -2] = 0

    psi[domain ==  3] = 0
    psi[domain ==  4] = 0
    psi[domain ==  5] = 0
    psi[domain ==  6] = 0
    psi[domain ==  7] = 0

    return psi

###############################################################################
## solving the laplace equation w/ no vorticity
###############################################################################
def inviscidPsi(psi, domain):
    for row in range(0, sett.ROW):
        for col in range(0, sett.COL):
            if domain[row, col] == 0:
                relaxed = recursive_convergence(psi, row, col)
                psi[row, col] = relaxed

    return psi

def recursive_convergence(psi, row, col):
    now = psi[row, col]

    psi[row, col] += sett.RELAX_FACTOR * (psi[row, col+1] + psi[row, col-1] + psi[row+1, col] + psi[row-1, col] - (4 * psi[row, col]))

    then = psi[row, col+1] + psi[row, col-1] + psi[row+1, col] + psi[row-1, col] - (4 * psi[row, col])

    while(then - now > sett.TOL):
        recursive_convergence(psi, row, col)

    return psi[row, col]

###############################################################################
## enter the Vortex
###############################################################################
def boundaryOmega(psi, domain, omega):
    for row in range(0, sett.ROW):
        for col in range(0, sett.COL):
            if domain[row, col] == 2:
                omega[row, col]  = (-2 * (psi[row+1, col] - psi[row, col])) / (sett.H**2)
            if domain[row, col] == -2:
                omega[row, col]  = (-2 * (psi[row-1, col] - psi[row, col])) / (sett.H**2)

            if domain[row, col] == 3:
                omega[row, col] = (-2 * (psi[row, col-1] - psi[row, col])) / (sett.H**2)
            if domain[row, col] == 4:
                omega[row, col] = (-2 * (psi[row-1, col] - psi[row, col])) / (sett.H**2)
            if domain[row, col] == 5:
                omega[row, col] = (-2 * (psi[row, col+1] - psi[row, col])) / (sett.H**2)
            if domain[row, col] == 6:
                omega[row, col] = (-2 * (psi[row+1, col] - psi[row, col])) / (sett.H**2)

    return omega
###############################################################################
## step u, v
###############################################################################
def updateUV(psi, domain, u, v):
    for col in range(0, sett.COL):
        for row in range(0, sett.ROW):
            if domain[row, col] == 0:
                u[row, col] = (psi[row+1, col] - psi[row-1, col])/ (2*sett.H)
                v[row, col] = (psi[row, col+1] - psi[row, col-1])/ (2*sett.H)

    return u, v

###############################################################################
## step omega ###############################################################################
def updateOmega(psi, omega, domain, u, v):
    for col in range(0, sett.COL):
        for row in range(0, sett.ROW):
            if domain[row, col] == 0:
                if u[row, col] < 0:
                    deluom = ((u[row, col+1] * omega[row, col+1]) -\
                     (u[row, col] * omega[row, col]))/sett.H
                if u[row, col] >= 0:
                    deluom = ((u[row, col] * omega[row, col]) - (u[row, col-1] * omega[row, col-1]))/sett.H

                if v[row, col] < 0:
                    delvom = ((v[row+1, col] * omega[row+1, col]) - (v[row, col] * omega[row, col]))/sett.H
                if v[row, col] >= 0:
                    delvom = ((v[row, col] * omega[row, col]) - (v[row-1, col] * omega[row-1, col]))/sett.H

                delomsq = ((omega[row-1, col] + omega[row+1, col] + omega[row, col-1] + omega[row, col+1]) - (4*omega[row, col])) / (sett.H**2)

                omega[row, col] += sett.DT * ((-deluom) - delvom + (sett.V_FLUID * delomsq))


    return omega

###############################################################################
## step psi ###############################################################################
def nextPsi(psi, omega, domain):
    for row in range(0, sett.ROW):
        for col in range(0, sett.COL):
            if domain[row, col] == 0:

                psi[row, col] += sett.RELAX_FACTOR * (psi[row, col+1] + psi[row, col-1] + psi[row-1, col] + psi[row+1, col] + (4 * (sett.H**2) * omega[row, col]) - (4 * psi[row, col]))

                relaxed2 = recursive_convergence2(psi, omega, row, col)

                psi[row, col] = relaxed2

            if domain[row, col] == -1:
                psi[row, col] = (2*psi[row, col-1]) - psi[row, col-2]

    return psi

def recursive_convergence2(psi, omega, row, col):
    now = sett.RELAX_FACTOR * (psi[row, col+1] + psi[row, col-1] + psi[row+1, col] + psi[row-1, col] + (4 * (sett.H**2) * omega[row, col]) - (4 * psi[row, col]))

    psi[row, col] += sett.RELAX_FACTOR * (psi[row, col+1] + psi[row, col-1] + psi[row+1, col] + psi[row-1, col] + (4 * (sett.H**2) * omega[row, col]) - (4 * psi[row, col]))

    then = sett.RELAX_FACTOR * (psi[row, col+1] + psi[row, col-1] + psi[row+1, col] + psi[row-1, col] + (4 * (sett.H**2) * omega[row, col]) - (4 * psi[row, col]))

    while(then - now > sett.TOL):
        recursive_convergence2(psi, omega, row, col)

    return psi[row, col]

##############################################################################
## init temperature
###############################################################################
def initTemp(domain):
    temp = np.zeros_like(domain)

    for col in range(0, sett.COL):
        for row in range(0, sett.ROW):
            if domain[row, col] == 0:
                temp[row, col] = 10
            if domain[row, col] == 1 or -1:
                temp[row, col] == 10
            if domain[row, col] == 2 or -2:
                temp[row, col] = 10
            if domain[row, col] > 2:
                temp[row, col] = 30

    return temp

##############################################################################
## step temperature
###############################################################################
def stepTemp(temp, domain, u, v):
    for col in range(0, sett.COL):
        for row in range(0, sett.ROW):
            if domain[row, col] == 0:
                if u[row, col] < 0:
                    delute = ((u[row, col+1] * temp[row, col+1]) -\
                     (u[row, col] * temp[row, col]))/sett.H
                if u[row, col] >= 0:
                    delute = ((u[row, col] * temp[row, col]) - (u[row, col-1] * temp[row, col-1]))/sett.H

                if v[row, col] < 0:
                    delvte = ((v[row+1, col] * temp[row+1, col]) - (v[row, col] * temp[row, col]))/sett.H
                if v[row, col] >= 0:
                    delvte = ((v[row, col] * temp[row, col]) - (v[row-1, col] * temp[row-1, col]))/sett.H

                delomsq = ((temp[row-1, col] + temp[row+1, col] + temp[row, col-1] + temp[row, col+1]) - (4*temp[row, col])) / (sett.H**2)

                temp[row, col] += sett.DT * ((-delute) - delvte + (sett.V_FLUID * delomsq))


    return temp
