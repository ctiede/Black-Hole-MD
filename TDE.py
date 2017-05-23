from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import initialize_BHmd as init

def BHmd(dim, nobj, COM, positions, velocities, M_BH, steps, tau, fileName):
    print "Black Hole MD Simulation"
    print 'Number of Objects:', nobj
    print 'Number of Timesteps:', steps
    print 'Timestep Size:', tau
    print fileName

    pos0 = np.copy(positions)
    pos = np.copy(positions)
    com = np.copy(COM)
    vel = np.copy(velocities)
    acc = np.zeros((nobj,dim))

    r_t = 200*M_BH      #Tidal radius

    print_step = 50
    count = 0
    open(fileName, 'w').close() #clears the file if it exists already
    for step in range(0,steps+1):
        com_dist = np.sqrt(np.sum(com[0,]**2))
        if com_dist > r_t:
            pos, vel, acc = evolve_bulk(dim, nobj, com, vel, acc, M_BH, tau)
        else:
            if count == 0:
                changeStep = step
            pos, vel, acc = evolve(dim, nobj, pos, vel, acc, M_BH, tau)
            count = 1

        if (step+print_step) % print_step == 0:
            print step
            dumpPos(pos, step, fileName)
            dumpVel(vel, step)
            dumpAcc(acc, step)

    print changeStep/print_step

def computeF(dim, nobj, pos0, vel0):
    pos = np.copy(pos0)
    r = np.sqrt(np.sum(pos**2,axis=1)) #take the norm of each row (boject)

    return -pos/r[:,None]**3

def computeF_PW(dim, nobj, pos0, vel0, M_BH):
    r_g = M_BH      #Just to highlight that M_BH is a length (G=c=1)
    r_s = 2*r_g     #Schwarzschild Radius

    pos = np.copy(pos0)
    r = np.sqrt(np.sum(pos**2, axis=1))

    return -pos/( (r[:,None] - r_s)**2 * r[:,None])

def evolve(dim, nobj, pos0, vel0, acc0, M_BH, tau):
    vel = np.copy(vel0)
    pos = np.copy(pos0)
    acc = np.copy(acc0)

    #Velocity Verlet--first half
    vel = vel + 0.5*tau*acc
    pos = pos + tau*vel #+ 0.5*tau*tau*acc

    #Compute accelerations from updated positions
    force = computeF_PW(dim, nobj, pos0, vel0, M_BH)
    acc = force #/ mass vector (mass[:,None])

    #Velocity Verlet--Last Step: Get full step velocities
    vel = vel + 0.5*tau*acc

    return pos, vel, acc

def evolve_bulk(dim, nobj, COM, pos0, vel0, acc0, M_BH, tau):
    vel = np.copy(vel0)
    pos = np.copy(pos0)
    acc = np.copy(acc0)
    com = np.copy(COM)

    #Velocity Verlet--first half
    vel = vel + 0.5*tau*acc
    pos = pos + tau*vel #+ 0.5*tau*tau*acc
    com = com + tau*vel[0]

    #Compute accelerations from updated positions
    force = computeF_PW(dim, nobj, com, vel0, M_BH)
    acc = force #/ mass vector (mass[:,None])

    #Velocity Verlet--Last Step: Get full step velocities
    vel = vel + 0.5*tau*acc

    return pos, vel, acc

def dumpPos(pos, step, fileName):
    f = open(fileName, 'a')
    N = pos.shape[0]
    line = "{0:d} \nAtoms. Timestep: {1:g} \n".format(N+1,step)
    f.write(line)


    type = np.full((N),1, dtype=float)
    type[0] = 1
    line = "{0:g} {1:g} {2:g} {3:g} \n".format(2, 0.0, 0.0, 0.0)
    f.write(line)


    for i in range(N):
        line = "{0:g} {1:g} {2:g} {3:g} \n".format(type[i], pos[i,0], pos[i,1], pos[i,2])
        f.write(line)
    return

def dumpVel(vel, step):
    fileName = 'velTest.txt'
    if step == 0:
        h = open(fileName,'w')
    else:
        h = open(fileName,'a')
    N = vel.shape[0]
    line = "{0:d} \nAtoms. Timestep: {1:g} \n".format(N,step)
    h.write(line)

    type = np.full((N),1, dtype=float)
    type[0] = 2
    for i in range(N):
        line = "{0:g} {1:g} {2:g} {3:g} \n".format(type[i], vel[i,0], vel[i,1], vel[i,2])
        h.write(line)
    return

def dumpAcc(acc, step):
    fileName = 'accTest.txt'
    if step == 0:
        k = open(fileName,'w')
    else:
        k = open(fileName,'a')
    N = acc.shape[0]
    line = "{0:d} \nAtoms. Timestep: {1:g} \n".format(N,step)
    k.write(line)

    type = np.full((N),1, dtype=float)
    type[0] = 2
    for i in range(N):
        line = "{0:g} {1:g} {2:g} {3:g} \n".format(type[i], acc[i,0], acc[i,1], acc[i,2])
        k.write(line)
    return
if __name__ == '__main__':
    dim = 3
    nobj = 1000

    steps = 500000
    dt = 0.01

    R=1
    ecc = 0.9
    r_g = 0.5
    fileName = 'test.xyz'
    pos, vel, com = init.init_star(dim, nobj, R, ecc)


    BHmd(dim, nobj, com, pos, vel, r_g, steps, dt, fileName)
