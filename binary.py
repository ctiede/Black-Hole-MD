from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import initialize_BHmd as init


def BHmd_bin(dim, nobj, positions, velocities, steps, tau, fileName):
    print "Black Hole MD Simulation"
    print 'Number of Objects:', nobj
    print 'Number of Timesteps:', steps
    print 'Timestep Size:', tau

    pos0 = np.copy(positions)
    pos = np.copy(positions)
    vel = np.copy(velocities)
    acc = np.zeros((nobj,dim))
    N = nobj

    #Initialize Black Hole Bindary
    BHpos, BHvel = init.binary(dim)
    BHacc = np.zeros((2, dim))

    print_step = 25
    open(fileName, 'w').close() #clears the file if it exists already
    for step in range(0,steps+1):
        #Update BH positions
        BHpos, BHvel, BHacc = BH_update(dim,BHpos, BHvel, BHacc, tau)
        #Check if particles "too" close to BHs
        N, pos, vel, acc = checkPos(N, BHpos, pos, vel, acc)
        #Update particle positions
        pos, vel, acc = evolve(dim, nobj, BHpos, pos, vel, acc, tau)

        if (step+print_step) % print_step == 0:
            print step
            dumpPos(N, BHpos, pos, step, fileName)
            #dumpVel(vel, step)
            #dumpAcc(acc, step)
            #dump_BH(BHpos, step)

def compute_Fbin(dim, nobj, BHpos, pos0, vel0):
    pos = np.copy(pos0)
    bh1 = BHpos[0]
    bh2 = BHpos[1]

    dist1 = pos - bh1
    dist2 = pos - bh2
    r1 = np.sqrt(np.sum(dist1**2, axis=1))
    r2 = np.sqrt(np.sum(dist2**2, axis=1))

    return -0.5*dist1/r1[:,None]**3 - 0.5*dist2/r2[:,None]**3

def computeBin(dim, pos0, vel0):
    pos = np.copy(pos0)
    r = np.sqrt(np.sum(pos**2,axis=1))
    return -pos/r[:,None]**3

def BH_update(dim, BHpos, BHvel, BHacc, tau):
    pos = np.copy(BHpos)
    vel = np.copy(BHvel)
    acc = np.copy(BHacc)

    #Velocity Verlet--first half
    vel = vel + 0.5*tau*acc
    pos = pos + tau*vel #+ 0.5*tau*tau*acc

    #Compute accelerations from updated positions
    force = computeBin(dim, pos, vel)
    acc = force #/ mass[:,None] ???

    #Velocity Verlet--Last Step: Get full step velocities
    vel = vel + 0.5*tau*acc

    return pos, vel, acc

def evolve(dim, nobj, BHpos, pos0, vel0, acc0, tau):
    vel = np.copy(vel0)
    pos = np.copy(pos0)
    acc = np.copy(acc0)

    #Velocity Verlet--first half
    vel = vel + 0.5*tau*acc
    pos = pos + tau*vel #+ 0.5*tau*tau*acc

    #Compute accelerations from updated positions
    force = compute_Fbin(dim, nobj, BHpos, pos0, vel0)
    #print force[100,:]
    #time.sleep(0.5)
    acc = force #/ mass[:,None] ???

    #Velocity Verlet--Last Step: Get full step velocities
    vel = vel + 0.5*tau*acc

    return pos, vel, acc

def checkPos(N, BHpos, pos0, vel0, acc0):
    vel = np.copy(vel0)
    pos = np.copy(pos0)
    acc = np.copy(acc0)
    bh1 = BHpos[0]
    bh2 = BHpos[1]

    cutoff = 0.1

    dist1 = pos - bh1
    dist2 = pos - bh2
    r1 = np.sqrt(np.sum(dist1**2, axis=1))
    r2 = np.sqrt(np.sum(dist2**2, axis=1))

    mask = (r1<cutoff)|(r2<cutoff)
    pos[mask] = np.NaN
    vel[mask] = np.NaN
    acc[mask] = np.NaN

    Ndel = (r1<cutoff).sum() + (r2<cutoff).sum()
    N = N - Ndel

    return N, pos, vel, acc



def dumpPos(N, BHpos, pos, step, fileName):
    f = open(fileName, 'a')
    #N = pos.shape[0]
    line = "{0:d} \nAtoms. Timestep: {1:g} \n".format(N+2,step)
    f.write(line)


    line = "{0:g} {1:g} {2:g} {3:g} \n".format(2, BHpos[0,0], BHpos[0,1], BHpos[0,2])
    f.write(line)
    line = "{0:g} {1:g} {2:g} {3:g} \n".format(2, BHpos[1,0], BHpos[1,1], BHpos[1,2])
    f.write(line)

    for i in range(N):
        line = "{0:g} {1:g} {2:g} {3:g} \n".format(1, pos[i,0], pos[i,1], pos[i,2])
        f.write(line)
    return


def dump_BH(BHpos, step):
    fileName = 'binaryTest.xyz'
    if step == 0:
        h = open(fileName,'w')
    else:
        h = open(fileName,'a')
    N = BHpos.shape[0]
    line = "{0:d} \nAtoms. Timestep: {1:g} \n".format(N,step)
    h.write(line)

    for i in range(N):
        line = "{0:g} {1:g} {2:g} {3:g} \n".format(1, BHpos[i,0], BHpos[i,1], BHpos[i,2])
        h.write(line)
    return


if __name__ == '__main__':
    dim = 3
    nobj = 2000
    steps = 10000000
    dt = 0.01

    rmin = 2
    rmax = 20
    theta = np.pi/4.
    fileName = 'test_binDisk.xyz'
    pos, vel = init.init_disk(dim, nobj, rmin, rmax, theta)


    BHmd_bin(dim, nobj, pos, vel, steps, dt, fileName)
