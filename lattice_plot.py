#! /usr/bin/env python

"""
lattice_plot.py
A python module for plot lattice model.
Zhihao Cui zhcui0408@gmail.com
"""

import os, sys
import numpy as np
import scipy.linalg as la
import itertools as it

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

#print(plt.rcParams.get('figure.figsize'))
#plt.style.use('ggplot')

# TODO
# 1. circle size
# 2. test pure hubbard pairing pattern 

def plot_name(plt, r, name, fontsize=15, **kwargs):
    plt.text(r[0], r[1], name, horizontalalignment='center', verticalalignment='center', \
            fontsize=fontsize, **kwargs)
    return plt

def plot_atom(plt, r, rad, color, **kwargs):
    plt.scatter(r[0], r[1], \
            c = color, s = np.sqrt(1000 * rad) * 20, \
            edgecolors='black', linewidths=1, **kwargs)
    return plt

def plot_atom_all(plt, lattice, idx_list, rad_list, color_list, **kwargs):
    for idx, rad, color in zip(idx_list, rad_list, color_list):
        r = lattice.site_idx2pos(idx)
        plt = plot_atom(plt, r, rad, color, **kwargs)
    return plt

def plot_atom_by_pos(plt, pos_list, rad_list, color_list, **kwargs):
    for r, rad, color in zip(pos_list, rad_list, color_list):
        plt = plot_atom(plt, r, rad, color, **kwargs)
    return plt

def plot_p_orb(plt, ax, r, direct='up', phase=["+", "-"]):
    width = 0.4
    height = 0.2
    fontsize = {"+":16, "-":24}
    if direct == 'up':
        r0 = [r[0], r[1] + width * 0.54]
        r1 = [r[0], r[1] - width * 0.54]
        plot_name(plt, r0, phase[0], fontsize=fontsize[phase[0]])
        plot_name(plt, r1, phase[1], fontsize=fontsize[phase[1]])
        angle = 90
    elif direct == 'down':
        r0 = [r[0], r[1] + width * 0.54]
        r1 = [r[0], r[1] - width * 0.54]
        plot_name(plt, r0, phase[1], fontsize=fontsize[phase[1]])
        plot_name(plt, r1, phase[0], fontsize=fontsize[phase[0]])
        angle = 90
    elif direct == 'left':
        r0 = [r[0] + width * 0.54, r[1]]
        r1 = [r[0] - width * 0.54, r[1]]
        plot_name(plt, r1, phase[0], fontsize=fontsize[phase[0]])
        plot_name(plt, r0, phase[1], fontsize=fontsize[phase[1]])
        angle = 0
    elif direct == 'right':
        r0 = [r[0] + width * 0.54, r[1]]
        r1 = [r[0] - width * 0.54, r[1]]
        plot_name(plt, r1, phase[1], fontsize=fontsize[phase[1]])
        plot_name(plt, r0, phase[0], fontsize=fontsize[phase[0]])
        angle = 0
    else:
        raise ValueError

    fc_list = ['C0' if p == "+" else 'C3' for p in phase]
    ellipse1 = Ellipse(xy=r0, width=width, height=height, edgecolor='black', fc=fc_list[0], lw=1.5, angle=angle, alpha=0.8)
    ellipse2 = Ellipse(xy=r1, width=width, height=height, edgecolor='black', fc=fc_list[1], lw=1.5, angle=angle, alpha=0.8)
    ax.add_patch(ellipse1)
    ax.add_patch(ellipse2)

def plot_d_orb(plt, ax, r, direct='up'):
    if direct == 'up':
        plot_p_orb(plt, ax, r, direct='up', phase=["+", "+"])
        plot_p_orb(plt, ax, r, direct='left', phase=["-", "-"])
    else:
        plot_p_orb(plt, ax, r, direct='up', phase=["-", "-"])
        plot_p_orb(plt, ax, r, direct='left', phase=["+", "+"])

def plot_spin(plt, r, ms, scal=4.0, **kwargs):
    if ms > 0.0:
        dx = 0.0
        dy = np.abs(ms) * 0.5 * scal
        r[1] -= dy*0.5
    else:
        dx = 0.0
        dy = -np.abs(ms) * 0.5 * scal
        r[1] -= dy*0.5
    plt.arrow(r[0], r[1], dx, dy, width=0.05, head_width=0.16, head_length=0.13, \
            length_includes_head=False, color='black', **kwargs)
    return plt

def plot_spin_all(plt, lattice, idx_list, ms_list, scal=4.0, **kwargs):
    for idx, ms in zip(idx_list, ms_list):
        r = lattice.site_idx2pos(idx)
        plt = plot_spin(plt, r, ms, scal=scal, **kwargs)
    return plt

def plot_spin_by_pos(plt, pos_list, ms_list, scal=4.0, **kwargs):
    for r, ms in zip(pos_list, ms_list):
        plt = plot_spin(plt, r, ms, scal=scal, **kwargs)
    return plt

def plot_bond(plt, r0, r1, val, color_list=['C2', 'C4'], **kwargs):
    x, y = zip(r0, r1)
    if val >= 0.0:
        cidx = 0
    else:
        cidx = -1
    plt.plot(x, y, color=color_list[cidx], linestyle='-', linewidth=val*1000, \
            alpha=0.65, zorder=0, **kwargs)
    return plt

def plot_pairing(plt, lattice, rab, idx_list, bond_thr = 2.1, bond_min=0.0, **kwargs):
    s = 0.5**0.5
    neighborDist = lattice.neighborDist
    #for i, j in it.combinations_with_replacement(idx_list, 2):
    for i, j in it.combinations(idx_list, 2):
        val = s*(rab[i, j] + rab[j, i])
        r0 = lattice.site_idx2pos(i)
        r1 = lattice.site_idx2pos(j)
        if la.norm(r1 - r0) > bond_thr or la.norm(r1 - r0) < bond_min:
            continue
        plt = plot_bond(plt, r0, r1, val, **kwargs)
    return plt

def plot_pairing_by_pos(plt, lattice, rab, idx_list, bond_thr = 2.1, bond_min=0.0, **kwargs):
    s = 0.5**0.5
    neighborDist = lattice.neighborDist
    #for i, j in it.combinations_with_replacement(idx_list, 2):
    # Oa corresponds to [12, 13, 14, 15]
    idx_dict = {12: 1, 13: 4, 14:7, 15:10}
    Oa_pos = np.asarray([[4.0, 1.0], [1.0, 0.0], [0.0, 3.0], [3.0, 4.0]])
    
    res = []
    #for i_, j_ in it.combinations(idx_list, 2):
    for i_, j_ in it.product(idx_list, idx_list):
        i = idx_dict.get(i_, i_)
        j = idx_dict.get(j_, j_)
        val = s*(rab[i, j] + rab[j, i])
        if i_ >= 12:
            r0 = Oa_pos[i_-12]
        else:
            r0 = lattice.site_idx2pos(i)
        if j_ >= 12:
            r1 = Oa_pos[j_-12]
        else:
            r1 = lattice.site_idx2pos(j)
        if la.norm(r1 - r0) > bond_thr or la.norm(r1 - r0) < bond_min:
            continue
        else:
            res.append((i_, j_, val))
        plt = plot_bond(plt, r0, r1, val, **kwargs)
    
    Oi_list = [2, 5, 8, 11]
    res = [resi for resi in res if (resi[0] < 12 and (not(resi[0] in Oi_list and resi[1]>11)))]
    #print len(res)
    print res
    res.sort()
    dwv = 0.0
    for i in range(0, len(res), 2):
        print
        print res[i][2]
        print res[i+1][2]
        #exit()
        print 0.25 * (2.*res[i][2] - 2.*res[i+1][2])
        dwv += np.abs(0.25 * (2.*res[i][2] - 2.*res[i+1][2]))
    #print res
    print "DWV"
    print dwv
    #exit()
    return plt, dwv

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
   
    #circle_rad = 0.5
    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    #fig, ax = plt.subplots()
    #plt.figure(facecolor="white")
    #ax = plt.gca()
    xmin = np.min(lat_coords[0])
    xmax = np.max(lat_coords[0])
    ymin = np.min(lat_coords[1])
    ymax = np.max(lat_coords[1])
    if "xleft" in kwargs:
        xleft = kwargs["xleft"]
    else:
        xleft = xmin - 0.15*lat_size[0]
    if "xright" in kwargs:
        xright = kwargs["xright"]
    else:
        xright = xmax + 0.15*lat_size[0]
    if "yleft" in kwargs:
        yleft = kwargs["yleft"]
    else:
        yleft = ymin - 0.15*lat_size[1]
    if "yright" in kwargs:
        yright = kwargs["yright"]
    else:
        yright = ymax + 0.15*lat_size[1]

    plt.xlim(xleft, xright)
    plt.ylim(yleft, yright)
    #plt.xlim(xmin - 0.15*lat_size[0], xmax + 0.15*lat_size[0])
    #plt.ylim(ymin - 0.15*lat_size[1], ymax + 0.15*lat_size[1])
    #plt.xlim(-lat_size[0] / 2 + 1.0, lat_size[0] / 2 + 1.0)
    #plt.ylim(-lat_size[1] / 2 + 1.0, lat_size[1] / 2 + 1.0)
    #plt.xticks(np.arange(-lat_size[0] / 2 , lat_size[0] / 2 + 1.0, 1.0))
    #plt.yticks(np.arange(-lat_size[1] / 2 , lat_size[1] / 2 + 1.0, 1.0))
    
    ax.set_aspect('equal', adjustable = 'box')
    ax.set_facecolor('white') # background color

    ax.axes.get_xaxis().set_visible(False) # do not show tick and labels
    ax.axes.get_yaxis().set_visible(False)

    if "noframe" in kwargs and kwargs["noframe"]:
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
    else:
        ax.spines['left'].set_linewidth(1.0)
        ax.spines['right'].set_linewidth(1.0)
        ax.spines['top'].set_linewidth(1.0)
        ax.spines['bottom'].set_linewidth(1.0)

    #ax.axhline(linewidth=4, color="g") draw a line along x-axis
    #plt.grid() # add grid
    #plt.set_axis_bgcolor('white')
    #plt.scatter(lat_coords[0], lat_coords[1], \
    #        c = color_plot, s = 1000 * rad_plot, \
    #        edgecolors='black', linewidths=1)
    return plt
    

if __name__ == '__main__':
    
    import libdmet.utils.logger as log
    import libdmet.dmet.abinitioBCS as dmet
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    LatSize = [1, 1]
    ImpSize = [1, 1]
    Lat = dmet.Square3BandSymm(*(LatSize + ImpSize))
    #LatSize = [4, 4]
    #ImpSize = [2, 2]
    #Lat = dmet.SquareLattice(*(LatSize + ImpSize))
    nscsites = Lat.supercell.nsites

    Cu_list = [0, 3, 6, 9]
    Oi_list = [2, 5, 8, 11]
    Oo_list = [1, 4, 7, 10]
    Oa_list = [12, 13, 14, 15]

    Oa_pos = [[4.0, 1.0], [1.0, 0.0], [0.0, 3.0], [3.0, 4.0]]

    
    #plt = plot_lattice(Lat, noframe=True, xleft=-0.5, xright=2.5, yleft=-0.5, yright=2.5)
    plt = plot_lattice(Lat, noframe=True)
    plt = plot_atom_all(plt, Lat, Cu_list, [0.9]*len(Cu_list), ['gold']*len(Cu_list))
    plt = plot_atom_all(plt, Lat, Oi_list, [0.3]*len(Oi_list), ['C3']*len(Oi_list))
    plt = plot_atom_all(plt, Lat, Oo_list, [0.3]*len(Oo_list), ['C3']*len(Oo_list))
    #plt = plot_atom(plt, [1.0, 0.0], rad=0.3, color='white')
    plt = plot_name(plt, [1.0, 0.4], "$\mathbf{Cu}$", fontsize=25)
    plt = plot_name(plt, [2.0, 0.4], "$\mathbf{O}$", fontsize=25)

#    ax = plt.gca()
#
#    rab = (np.random.random((nscsites, nscsites)) - 0.5) * 0.01
#    idx_list = [0, 3, 6, 9]
#    #plt = plot_pairing(plt, Lat, rab, idx_list, bond_thr=2.1)
#    #plt = plot_spin(plt, [1.0, 0.0], -0.3)
#    plt = plot_name(plt, [0.7, 1.25], "$\mathbf{U_{d}}$", fontsize=25)
#    #plt = pairing_bond(plt, (0.0, 0.0), (1.0, 0.0), 0.005)
#    r_list = [[1.0, 0.0], [0.0, 1.0], [1.0, 2.0], [2.0, 1.0]]
#    phase_list = ['up', 'right', 'up', 'right']
#    for r, ph in zip(r_list, phase_list):
#        plot_p_orb(plt, ax, r, direct=ph)
#    
#    plot_d_orb(plt, ax, [1.0, 1.0], direct=ph)
    #plt.savefig('symmtrized-cluster.png', dpi=600)
    plt.show()
