#! /usr/bin/env python 


import numpy as np
import h5py
from lattice_plot import *

def load_dmet_iter_npy(fname = './dmet_iter.npy'):
    return np.load(fname)

def load_dmet_npy(fname = './dmet.npy'):
    return np.load(fname)

def get_order_param_3band(GRhoImp, AFM_idx = [0, 3, 9, 6]):
    limp = GRhoImp.shape[0] / 2

    AFM_idx_flatten = np.array(AFM_idx).flatten()
    ra = GRhoImp[:limp,:limp]
    rb = np.eye(limp) - GRhoImp[limp:,limp:]
    rab = GRhoImp[:limp,limp:]
    
    rdm_a = GRhoImp[:limp,:limp][np.ix_(AFM_idx_flatten, AFM_idx_flatten)]
    rdm_b = np.eye(len(AFM_idx_flatten)) - GRhoImp[limp:,limp:][np.ix_(AFM_idx_flatten, AFM_idx_flatten)]
    rdm_ab = GRhoImp[:limp, limp:][np.ix_(AFM_idx_flatten, AFM_idx_flatten)]

    m0 = 0.5 * (rdm_a[0,0]-rdm_b[0,0])
    m3 = 0.5 * (rdm_a[3,3]-rdm_b[3,3])
    m1 = 0.5 * (rdm_a[1,1]-rdm_b[1,1])
    m2 = 0.5 * (rdm_a[2,2]-rdm_b[2,2])
    afm = 0.25 * (m0 + m3 - m1 - m2)

    s = 0.5**0.5
    #d01 = s*(rdm_ab[0,1]+rdm_ab[1,0])
    #d23 = s*(rdm_ab[2,3]+rdm_ab[3,2])
    #d02 = s*(rdm_ab[0,2]+rdm_ab[2,0])
    #d13 = s*(rdm_ab[1,3]+rdm_ab[3,1])
    #dwv = 0.25*(d01+d23-d02-d13)


    #
    #
    #
    # 2x2 symmetric supercells
    #        |
    #        4O
    #        |        |
    #     - 3Cu - 5O -6Cu - 7O -
    #        |        |
    #        2O         8O
    #        |        |
    # - 1O - 0Cu - 11O - 9Cu -
    #        |        |
    #                 10O
    #                 |

    def Cu1(rab):
        # Cu-Cu
        d01 = s*(rab[0,3]+rab[3,0])
        d23 = s*(rab[6,9]+rab[9,6])
        d02 = s*(rab[3,6]+rab[6,3])
        d13 = s*(rab[9,0]+rab[0,9])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv

    
    def Cu2(rab):
        # Cu-inner O
        d01 = s*(rab[0,5]+rab[5,0])
        d23 = s*(rab[6,11]+rab[11,6])
        d02 = s*(rab[3,8]+rab[8,3])
        d13 = s*(rab[9,2]+rab[2,9])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv

    def Cu3(rab):
        # Cu-outer O
        d01 = s*(rab[0,4]+rab[4,0])
        d23 = s*(rab[6,10]+rab[10,6])
        d02 = s*(rab[3,7]+rab[7,3])
        d13 = s*(rab[9,1]+rab[1,9])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv
    
    def iO1(rab):
        # iO-Cu
        d01 = s*(rab[2,3]+rab[3,2])
        d23 = s*(rab[8,9]+rab[9,8])
        d02 = s*(rab[5,6]+rab[6,5])
        d13 = s*(rab[11,0]+rab[0,11])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv

    def iO2(rab):
        # iO-iO
        d01 = s*(rab[2,5]+rab[5,2])
        d23 = s*(rab[8,11]+rab[11,8])
        d02 = s*(rab[5,8]+rab[8,5])
        d13 = s*(rab[11,2]+rab[2,11])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv
    
    def iO3(rab):
        # iO-oO
        d01 = s*(rab[2,4]+rab[4,2])
        d23 = s*(rab[8,10]+rab[10,8])
        d02 = s*(rab[5,7]+rab[7,5])
        d13 = s*(rab[11,1]+rab[1,11])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv

    
    def oO1(rab):
        # oO-Cu
        d01 = s*(rab[1,3]+rab[3,1])
        d23 = s*(rab[7,9]+rab[9,7])
        d02 = s*(rab[4,6]+rab[6,4])
        d13 = s*(rab[10,0]+rab[0,10])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv
    
    def oO2(rab):
        # oO-iO
        d01 = s*(rab[1,5]+rab[5,1])
        d23 = s*(rab[7,11]+rab[11,7])
        d02 = s*(rab[4,8]+rab[8,4])
        d13 = s*(rab[10,2]+rab[2,10])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv

    def oO3(rab):
        # oO-oO
        d01 = s*(rab[1,4]+rab[4,1])
        d23 = s*(rab[7,10]+rab[10,7])
        d02 = s*(rab[4,7]+rab[7,4])
        d13 = s*(rab[10,1]+rab[1,10])
        dwv = 0.25*(d01+d23-d02-d13)
        return dwv

    
    f_list = [Cu1, Cu2, Cu3, iO1, iO2, iO3, oO1, oO2, oO3]

    dtot = 0.0
    for f in f_list:
        print f.__name__
        print f.__call__(rab)
        dtot += abs(f.__call__(rab))


    print "total dwv:"
    print dtot
    print "afm parameter:"
    print np.abs(afm)
    rho_Cu_a = np.diag(ra)[[0, 3, 6, 9]]
    rho_Cu_b = np.diag(rb)[[0, 3, 6, 9]]
    rho_Oi_a = np.diag(ra)[[2, 5, 8, 11]]
    rho_Oi_b = np.diag(rb)[[2, 5, 8, 11]]
    rho_Oo_a = np.diag(ra)[[1, 4, 7, 10]]
    rho_Oo_b = np.diag(rb)[[1, 4, 7, 10]]
    
    rho_Cu = rho_Cu_a + rho_Cu_b
    rho_Oi = rho_Oi_a + rho_Oi_b
    rho_Oo = rho_Oo_a + rho_Oo_b

    print "rho_Cu"
    print rho_Cu
    print "rho_Oi"
    print rho_Oi
    print "rho_Oo"
    print rho_Oo

    m_Oi = np.abs(rho_Oi_a - rho_Oi_b) * 0.5
    m_Oo = np.abs(rho_Oo_a - rho_Oo_b) * 0.5


    print "rho_Cu"
    print np.average(rho_Cu)
    print "rho_O"
    print np.average([rho_Oi, rho_Oo])
    print "m_Cu"
    print np.average(np.abs(rho_Cu_a - rho_Cu_b) * 0.5)
    print "m_O"
    print np.average([m_Oi, m_Oo])

    
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
    
    #print rab
    #rab = (np.random.random((nscsites, nscsites)) - 0.5) * 0.01
    idx_list = [0, 3, 6, 9]
    plt = plot_pairing(plt, Lat, rab, idx_list, bond_thr=2.1)
    plt = plot_atom(plt, [1.0, 0.0], rad=0.4, color='C3')
    #plt = pairing_bond(plt, (0.0, 0.0), (1.0, 0.0), 0.005)
    plt.show()



    exit()


    dwv = 0.25*(d01+d23-d02-d13)

    return afm,dwv



if __name__ == '__main__':

    import sys
    import numpy as np
    if len(sys.argv) > 1 :
        fname = sys.argv[1]
    else:
        fname = './dmet.npy'

    #Mu, last_dmu, vcor_param, basis, GRhoEmb, \
    #    GRhoImp, EnergyImp, nelecImp = load_dmet_npy(fname)
    A = load_dmet_npy(fname)
    GRhoImp = A[5]

    np.set_printoptions(3, linewidth =1000, suppress = True)

    print get_order_param_3band(GRhoImp)

