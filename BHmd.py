from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import initialize_BHmd as init

def cart_to_cyl(vec0):
    #Takes a 2D array of n objects and their cartesian coordinates and
    #converts them to cylindrical coordinates (polar for 2D)
    vec = np.copy(vec0)
    xcoords = vec[:,0]
    ycoords = vec[:,1]

    rs = np.sqrt( xcoords**2 + ycoords**2 )
    phis = np.arccos(xcoords/rs)

    vec[:,0] = rs
    vec[:,1] = phis
    return vec

def cyl_to_cart(vec0):
    #Takes a 2D array of n objects and their cylindrical (polar for 2D)
    #coordinates andconverts them to cartesian coordinates
    vec = np.copy(0)
    rs = vec[:,0]
    phis = vec[:,1]

    xs = rs*np.cos(phis)
    ys = rs*np.sin(phis)

    vec[:,0] = xs
    vec[:,1] = ys
    return vec


def BHmd(dim, nobj, positions, velocities, steps, tau, fileName):
    print "Black Hole MD Simulation"
    print 'Number of Objects:', nobj
    print 'Number of Timesteps:', steps
    print 'Timestep Size:', tau
    print fileName

    pos0 = np.copy(positions)
    pos = np.copy(positions)
    vel = np.copy(velocities)
    acc = np.zeros((nobj,dim))

    print_step = 25
    open(fileName, 'w').close() #clears the file if it exists already
    for step in range(0,steps+1):
        pos, vel, acc = evolve(dim, nobj, pos, vel, acc, tau)

        if (step+print_step) % print_step == 0:
            print step
            dumpPos(pos, step, fileName)
            dumpVel(vel, step)
            dumpAcc(acc, step)

def computeF(dim, nobj, pos0, vel0):
    pos = np.copy(pos0)
    r = np.sqrt(np.sum(pos**2,axis=1)) #take the norm of each row (boject)

    return -pos/r[:,None]**3

def compute_Fbin(dim, nobj, pos0, vel0, sep):
    pos = np.copy(pos0)
    bh1 = np.array([sep/2.0, 0, 0])
    bh2 = np.array([-sep/2.0, 0, 0])

    dist1 = pos - bh1
    dist2 = pos - bh2
    r1 = np.sqrt(np.sum(dist1**2, axis=1))
    r2 = np.sqrt(np.sum(dist2**2, axis=1))

    return -dist1/r1[:,None]**3 - dist2/r2[:,None]**3

def evolve(dim, nobj, pos0, vel0, acc0, tau):
    vel = np.copy(vel0)
    pos = np.copy(pos0)
    acc = np.copy(acc0)

    #Velocity Verlet--first half
    vel = vel + 0.5*tau*acc
    pos = pos + tau*vel #+ 0.5*tau*tau*acc

    #Compute accelerations from updated positions
    force = computeF(dim, nobj, pos0, vel0)
    #print force[100,:]
    #time.sleep(0.5)
    acc = force #/ mass[:,None] ???

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

    steps = 10000
    dt = 0.01

    rmin = 2
    rmax = 20
    theta = np.pi/6.
    fileName = 'test.xyz'
    pos, vel = init.init_disk(dim, nobj, rmin, rmax, theta)


    BHmd(dim, nobj, pos, vel, steps, dt, fileName)
