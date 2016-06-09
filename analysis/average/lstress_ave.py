#!/usr/bin/python

import os, sys, glob
import numpy as np
import math

def main(files, fout_name):
    # get average stress
    ave_stress = np.loadtxt(files[0])
    grid_inform = ave_stress[:, 0:3]
    ave_stress = ave_stress[:, 3:]
    for fin_name in files[1:]:
        ave_stress += np.loadtxt(fin_name)[:, 3:]
    ave_stress /= len(files)
        
    # calc stress std
    stdev_stress = np.zeros(ave_stress.shape)
    for fin_name in files:
        lstress       = np.loadtxt(fin_name)[:, 3:]
        stdev_stress += (lstress - ave_stress) * (lstress - ave_stress)
    stdev_stress /= len(files)
    stdev_stress = np.sqrt(stdev_stress)
    stdev_stress /= math.sqrt(len(files))
        
    result = np.hstack((grid_inform, ave_stress, stdev_stress))
    np.savetxt(fout_name, result)
        
if __name__ == "__main__":
    ls_list = set(['all', 'imol', 'bond', 'angle', 'kin'])
    
    argvs = sys.argv
    argc  = len(argvs)
    if argc != 3:
        print "Usage:"
        print "$ python %s dir_name localstress_component" % argvs[0]
        quit()

    ls_comp = argvs[2]
    if not (ls_comp in ls_list) :
        print "Unknown local stress component %s." % ls_comp
        print "The following lists show the supporting local stress components."
        print ls_list
        quit()
    
    cdir   = argvs[1]
    fnames = os.path.join(cdir, "sample*/loc_stress_%s.txt" % ls_comp)
    files  = glob.glob(fnames)
    fout   = os.path.join(cdir, "loc_stress_ave_%s.txt" % ls_comp)
    main(files, fout)