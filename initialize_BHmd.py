from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import h5py
import time


def get_radius(x):
    x = np.asarray(x)
    r = np.sqrt( x[:,0]**2 + x[:,1]**2 + x[:,2]**2 )
    return r

def disk(n, rmin, rmax):
    ga = np.pi*(3-np.sqrt(5))
    th = ga*np.arange(n)
    r = (rmax - rmin)*np.sqrt(np.arange(n)/n) + rmin
    points = np.zeros((n,3))
    points[:,0] = r*np.cos(th)
    points[:,1] = r*np.sin(th)
    return points

def rotate_yax(vec, theta):
    x = np.copy(vec)
    x[:,2] = x[:,0]*np.sin(theta) #give each point a z-component
    x[:,0] = x[:,0]*np.cos(theta) #change x-comp so keep rho (in xz-plane) const
    print "calld"
    return x

def sphere_surf(n):
    ga = np.pi*(3-np.sqrt(5))
    th = ga*np.arange(n)
    z = np.linspace( 1-1.0/n, 1.0/n - 1, n)
    r = np.sqrt(1 - z*z)
    points = np.zeros((n,3))
    points[:,0] = r*np.cos(th)
    points[:,1] = r*np.sin(th)
    points[:,2] = z
    return points

def rand_sphere(n, R):
    #Create array of distances
    U = np.random.rand(n)
    D = R*U**(1./3.)
    #Uniform distribution on surface of sphere
    points = sphere_surf(n)

    return D[:,None]*points


def init_disk(dim, nobj, rmin, rmax, theta):
    #Creates a Keplarian Disk out of the nobj bodies in dim dimensions
    pos = np.zeros((nobj, dim))
    vel = np.zeros((nobj, dim))
    #theta = 0

    #Create positions array of particles in the disk from rmin to rmax
    # at angle theta
    pos = disk(nobj, rmin, rmax)

    #Initialize velocities
    rs = get_radius(pos) #np.sqrt( pos[:,0]**2 + pos[:,1]**2 )
    phis = np.arctan(pos[:,1]/pos[:,0])
    #Array of angular velocities
    vphi = np.sqrt(1/rs)

    vx = -vphi*np.sin(phis)
    vy = vphi*np.cos(phis)
    vx[pos[:,0]>0] = -vx[pos[:,0]>0]
    vy[pos[:,0]>0] = -vy[pos[:,0]>0]

    vel[:,0] = vx
    vel[:,1] = vy

    #Tilt the Disk
    pos = rotate_yax(pos, theta)
    #Rotate Velocities
    vel = rotate_yax(vel, theta)

    print np.max(get_radius(pos))

    print np.max(pos[:,1]), np.min(pos[:,1])

    return pos, vel


def init_star(dim, nobj, R, ecc):
    #Takes dimensions, number of particles, an impact parameter b,
    # the y-distance from the origin, and the y-velocity of the particles
    if dim != 3:
        print "Need to be in 3 dimensions"
    points = rand_sphere(nobj,R)

    fac=50.
    x0 = (1+ecc)*fac
    y0 = 0.
    z0 = 0.
    points[:,0] += x0
    COM = np.array([x0,y0,z0])
    #points[:,1] += y0
    #points[:,2] += z0

    #Init Velocities
    velos = np.zeros((nobj, dim))
    vy = np.sqrt( (1-ecc)/(1+ecc) ) / np.sqrt(fac)
    velos[:,1] = vy

    n_com = np.zeros((nobj, dim))
    n_com[:,] = COM

    np.savetxt('init_pos.dat',points)
    np.savetxt('init_vel.dat',velos)
    return points, velos, n_com

def binary(dim):
    BHpos = np.zeros((2,dim))
    BHvel = np.zeros((2,dim))
    BHpos[:,0] = [0.5, -0.5]
    BHvel[:,1] = [-np.sqrt(2), np.sqrt(2)]

    return BHpos, BHvel



if __name__ == '__main__':
    n = 5000
    b = 0
    rmin = 5
    rmax = 20
    R = 0.1
    ecc = 0.9
    x0 = 0
    y0 = 0
    z0 = 0

    #points = init_disk(3, n, rmin, rmax, np.pi/6.)[0]
    points = init_star(3,n,R,ecc)[0]
    x = points[:,0] + x0
    y = points[:,1] + y0
    z = points[:,2] + z0
    print "Maxs:", np.max(x), np.max(y)   #, np.max(z)
    print "Mins:", np.min(np.abs(x)), np.min(np.abs(y))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z) #, ylim([-1,11]), zlim([-1,11]))
    #ax.set_xlim([-2,12])
    #ax.set_zlim([-2,12])
    #ax.set_ylim([-2,12])
    ax.scatter(0,0,0, c='black')
    plt.show()
