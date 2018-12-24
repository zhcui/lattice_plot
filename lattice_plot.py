#! /usr/bin/env python
"""
lattice_plot.py
a python module for plot lattice model.
written by Zhihao Cui zcui@caltech.edu
"""

import os, sys
import numpy as np
import scipy.linalg as la

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
    bath_plot = plt.scatter(lat_coords[0], lat_coords[1], \
            c = color_plot, s = 1000 * rad_plot, \
            edgecolors='black', linewidths=1)
    plt.show()
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
    
    plot_lattice(Lat)


