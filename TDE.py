from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import initialize_BHmd as init

def BHmd(dim, nobj, COM, positions, velocities, M_BH, r_t, steps, tau, fileName):
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
    N = nobj

    #r_t = Q          #Tidal radius = (M_BH/M_*)^(1/3) * r_*
    print "Tidal radius: ", r_t
    print "Schw Radius: ", M_BH*2
    print "COM Distance: "  ,np.sqrt(np.sum(com[0,]**2))
    #time.sleep(10)

    print_step = 100
    count = 0
    open(fileName, 'w').close() #clears the file if it exists already
    for step in range(0,steps+1):
        com_dist = np.sqrt(np.sum(com[0,]**2))
        #If outside the tidal radius evlolve all particles as COM
        if com_dist > r_t:
            pos, vel, acc, com = evolve_bulk(dim, nobj, com, pos, vel, acc, M_BH, tau)
        #Once inside Tidal radius, evolve each individually
        else:
            if count == 0:
                changeStep = step
                com_pos = com[0]
            #Check if particles "too" close to BHs--if inside r_s get discarded
            N, pos, vel, acc = checkPos(N, pos, vel, acc, M_BH)
            #Update Particle Positions
            pos, vel, acc = evolve(dim, nobj, pos, vel, acc, M_BH, tau)
            count = 1

        if (step+print_step) % print_step == 0:
            print step
            dumpPos(N, pos, step, fileName)
            #dumpVel(vel, step)
            #dumpAcc(acc, step)

    #Create restart files
    np.savetxt('./Restart/restart_pos.dat',pos)
    np.savetxt('./Restart/restart_vel.dat',vel)

    print count
    print changeStep/print_step, com_pos

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
    com = com + tau*vel

    #Compute accelerations from updated positions
    force = computeF_PW(dim, nobj, com, vel0, M_BH)
    acc = force #/ mass vector (mass[:,None])

    #Velocity Verlet--Last Step: Get full step velocities
    vel = vel + 0.5*tau*acc

    return pos, vel, acc, com

def checkPos(N, pos0, vel0, acc0, M_BH):
    vel = np.copy(vel0)
    pos = np.copy(pos0)
    acc = np.copy(acc0)
    r_s = 2*M_BH

    dist = np.sqrt(np.sum(pos**2, axis=1))

    mask = (dist<r_s)
    pos[mask] = np.NaN
    vel[mask] = np.NaN
    acc[mask] = np.NaN

    Ndel = (dist < r_s).sum()
    N = N - Ndel

    return N, pos, vel, acc

def calc_Periapse(Beta, R_T):
    Beta = np.float(Beta)
    R_T = np.float(R_T)
    return R_T/Beta

def calc_Apoapse(Beta, R_T, ecc):
    Beta = np.float(Beta)
    R_T = np.float(R_T)
    ecc = np.float(ecc)

    R_P = calc_Periapse(Beta, R_T)
    fac = R_P/(1-ecc)

    return fac*(1+ecc)

def dumpPos(N, pos, step, fileName):
    f = open(fileName, 'a')
    #N = pos.shape[0]
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

    R=1                         #star radius
    f = .5                     #r_star = f*r_schw
    q = 100                     #(M_BH/M_star)^(1/3)
    beta = 11.6              #Inverse impact parameter: r_t/r_p
    ecc = 0.9                  #Eccentricity

    r_g = 0.5/f                 #r_g = M_BH
    r_t = q                     #tidal radius = (M_BH/M_star)^(1/3)
    r_p = calc_Periapse(beta, r_t)   #Periapse
    r_a = calc_Apoapse(beta, r_t, ecc)    #Apoapse
    #rt_factor = 0.75            #semi-major axis = rt_factor * tidal radius

    fileName = 'test.xyz'
    pos, vel, com = init.init_star_ellip(dim, nobj, R, ecc, beta, r_t)


    BHmd(dim, nobj, com, pos, vel, r_g, r_t, steps, dt, fileName)
