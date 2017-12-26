###############################################################################
## this is the fishe main file
###############################################################################

###############################################################################
## imports
###############################################################################
import numpy as np
import sett
import fish_func as ff
import time

###############################################################################
## steps
###############################################################################
## 1. create the domain wtih material type flags
## 2. set the boundary and initial conditions for psi, omega
## 3. assume omega = 0 at t = 0, find inviscid solution for PSI. laplace eq.
## 4. invoke no slip condition/vorticity at BOUNDARIES
## 5. bulk computation
## 6. outflow boundary condition, of course
## 7. update temperature in bulk flow from initial conditions of T
## 8. rinse n rerinse


###############################################################################
## constants and printer settings are stored in the sett.py file
###############################################################################
sett.init()

###############################################################################
## build array containing node flags
###############################################################################
domain = ff.buildDomain()

print "\n \t \t THIS IS THE DOMAIN: \n"
print domain
print ('\n')

time.sleep(1)

###############################################################################
## psi one
###############################################################################
psi = ff.initPsi(domain)

print "\n \t \t INITIALIZE PSI: \n"
print psi
print ('\n')

time.sleep(1)

###############################################################################
## solving the laplace equation w/ no vorticity
###############################################################################
print "\n \t \t CONVERGING TO INVISCID SOLUTION . . .\n"
time.sleep(1)

now = sum(sum(i for i in psi))
psi = ff.inviscidPsi(psi, domain)
then = sum(sum(i for i in psi))

while (then - now) > sett.TOL:
    now = sum(sum(i for i in psi))
    psi = ff.inviscidPsi(psi, domain)
    then = sum(sum(i for i in psi))

    print psi

time.sleep(3)

###############################################################################
## enter the Vortex
###############################################################################
omega = np.zeros_like(domain)

omega = ff.boundaryOmega(psi, domain, omega)

print "\n \t \t DETERMINE VORTICITY AT BOUNDARIES: \n"
print omega
time.sleep(3)

###############################################################################
## update u and v from STREAMFUNCTION
###############################################################################
u = np.zeros_like(domain)
v = np.zeros_like(domain)

u, v = ff.updateUV(psi, domain, u, v)

print "\n \t \t UPDATE U AND V FROM STREAMFUNCTION : \n"
print u, "\n"
print v, "\n"
time.sleep(1)

##############################################################################
## step omega ###############################################################################
print "\n \t \t STEP OMEGA : \n"

omega = ff.updateOmega(psi, omega, domain, u, v)

print omega, "\n"
time.sleep(1)


##############################################################################
## step psi ###############################################################################
print "\n \t \t STEP PSI : \n"

prev = sum(sum(i for i in psi))
psi = ff.nextPsi(psi, omega, domain)
later = sum(sum(i for i in psi))

while (later - prev > sett.TOL):
    prev = sum(sum(i for i in psi))
    psi = ff.nextPsi(psi, domain)
    later = sum(sum(i for i in psi))

    print psi

print psi

###############################################################################
## init temperature
###############################################################################
print "\n \t \t INITIALIZED TEMPERATURE : \n"

temp = ff.initTemp(domain)

print temp, "\n"
time.sleep(2)

##############################################################################
## init temperature
###############################################################################
print "\n \t \t STEPPED TEMPERATURE : \n"

temp = ff.stepTemp(temp, domain, u, v)

print temp, "\n"
time.sleep(2)
