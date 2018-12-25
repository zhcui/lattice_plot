#! /usr/bin/env python
"""
lattice_plot.py
a python module for plot lattice model.
written by Zhihao Cui zcui@caltech.edu
"""

import os, sys
import numpy as np
import scipy.linalg as la
import itertools as it

# import the implementation of the matrix factorization that consider bias
import matplotlib.pyplot as plt
from collections import defaultdict
import random as rd
import operator
from time import time
#from sklearn.metrics.pairwise import euclidean_distances
#from sklearn.random_projection import johnson_lindenstrauss_min_dim
#from sklearn.random_projection import SparseRandomProjection
#from sklearn.datasets import fetch_20newsgroups_vectorized
#from sklearn.datasets import load_digits
#from sklearn.metrics.pairwise import euclidean_distances

#plt.style.use('ggplot')

# TODO
# 2. circle size
# 4. test pure hubbard pairing pattern 

def plot_name(plt, r, name, fontsize=15):
    plt.text(r[0], r[1], name, horizontalalignment='center', verticalalignment='center', fontsize=fontsize)
    return plt

def plot_atom(plt, r, rad, color):
    plt.scatter(r[0], r[1], \
            c = color, s = 1000 * rad, \
            edgecolors='black', linewidths=1)
    return plt

def plot_spin(plt, r, ms, scal=5.0):
    if ms > 0.0:
        dx = 0.0
        dy = np.abs(ms) * 0.5 * scal
        r[1] -= dy*0.5
    else:
        dx = 0.0
        dy = -np.abs(ms) * 0.5 * scal
        r[1] -= dy*0.5
    plt.arrow(r[0], r[1], dx, dy, width=0.05, head_width=0.16, head_length=0.13, length_includes_head=True, color='black')
    return plt

def plot_bond(plt, r0, r1, val, color_list=['C2', 'C4']):
    x, y = zip(r0, r1)
    if val >= 0.0:
        cidx = 0
    else:
        cidx = -1
    plt.plot(x, y, color=color_list[cidx], linestyle='-', linewidth=val*1000, alpha=0.55, zorder=0)
    return plt

def plot_pairing(plt, lattice, rab, idx_list, bond_thr = 2.1):
    s = 0.5**0.5
    neighborDist = lattice.neighborDist
    for i, j in it.combinations_with_replacement(idx_list, 2):
        val = s*(rab[i, j] + rab[j, i])
        r0 = lattice.site_idx2pos(i)
        r1 = lattice.site_idx2pos(j)
        if la.norm(r1 - r0) > bond_thr:
            continue
        plt = plot_bond(plt, r0, r1, val)
    return plt

def plot_lattice(lattice, **kwargs):
    color_list = ['gold', 'C3']
    rad_list = [0.8, 0.4]

    nscsites = lattice.supercell.nsites
    lat_vec = lattice.size
    lat_size = np.diag(lat_vec)
    lat_coords = np.array(lattice.sites).T
    atom_names = lattice.names
    atom_ind, unq_idx = np.unique(atom_names, return_inverse=True)
    color_plot = np.asarray(color_list)[unq_idx]
    rad_plot = np.asarray(rad_list)[unq_idx]
    #color_dict = dict(zip(atom_ind, color_list))
    #color_plot = [color_dict[name] for name in atom_names]

    #lat_coords[0][lat_coords[0] > lat_size[0] // 2] -= lat_size[0]
    #lat_coords[1][lat_coords[1] > lat_size[1] // 2] -= lat_size[1]
   
    circle_rad = 0.5
    fig, ax = plt.subplots()
    #plt.figure(facecolor="white")
    #ax = plt.gca()
    xmin = np.min(lat_coords[0])
    xmax = np.max(lat_coords[0])
    ymin = np.min(lat_coords[1])
    ymax = np.max(lat_coords[1])
    plt.xlim(xmin - 0.15*lat_size[0], xmax + 0.15*lat_size[0])
    plt.ylim(ymin - 0.15*lat_size[0], ymax + 0.15*lat_size[1])
    #plt.xlim(-lat_size[0] / 2 + 1.0, lat_size[0] / 2 + 1.0)
    #plt.ylim(-lat_size[1] / 2 + 1.0, lat_size[1] / 2 + 1.0)
    #plt.xticks(np.arange(-lat_size[0] / 2 , lat_size[0] / 2 + 1.0, 1.0))
    #plt.yticks(np.arange(-lat_size[1] / 2 , lat_size[1] / 2 + 1.0, 1.0))
    
    ax.set_aspect('equal', adjustable = 'box')
    ax.set_facecolor('white') # background color

    ax.axes.get_xaxis().set_visible(False) # do not show tick and labels
    ax.axes.get_yaxis().set_visible(False)

    #ax.spines['top'].set_visible(True)
    #ax.spines['bottom'].set_visible(True)
    #ax.spines['left'].set_visible(True)
    #ax.spines['right'].set_visible(True)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['right'].set_linewidth(1.0)
    ax.spines['top'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
   
    #ax.axhline(linewidth=4, color="g") draw a line along x-axis
    #plt.grid() # add grid
    #plt.set_axis_bgcolor('white')
    plt.scatter(lat_coords[0], lat_coords[1], \
            c = color_plot, s = 1000 * rad_plot, \
            edgecolors='black', linewidths=1)
    return plt
    

if __name__ == '__main__':
    
    import libdmet.utils.logger as log
    import libdmet.dmet.abinitioBCS as dmet
    LatSize = [1, 1]
    ImpSize = [1, 1]
    Lat = dmet.Square3BandSymm(*(LatSize + ImpSize))
    #LatSize = [4, 4]
    #ImpSize = [2, 2]
    #Lat = dmet.SquareLattice(*(LatSize + ImpSize))
    nscsites = Lat.supercell.nsites

    plt = plot_lattice(Lat)

    rab = (np.random.random((nscsites, nscsites)) - 0.5) * 0.01
    idx_list = [0, 3, 6, 9]
    plt = plot_pairing(plt, Lat, rab, idx_list, bond_thr=2.1)
    plt = plot_atom(plt, [1.0, 0.0], rad=0.4, color='C3')
    plt = plot_spin(plt, [1.0, 0.0], -0.3)
    plt = plot_name(plt, [1.0, 1.0], "Cu")
    #plt = pairing_bond(plt, (0.0, 0.0), (1.0, 0.0), 0.005)
    plt.show()
