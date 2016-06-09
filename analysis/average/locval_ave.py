#!/usr/bin/python

import os, sys, glob
import numpy as np
import math

def main(fin_name, axis, fout_name):
    locval_data = np.loadtxt(fin_name)
    grid_inform_all = locval_data[:, 0]
    grid_inform = np.unique(grid_inform_all)
    num_grid = len(grid_inform)
    num_frames = len(grid_inform_all) / len(grid_inform)
    
    average = np.zeros(len(grid_inform))
    stdev   = np.zeros(len(grid_inform))
    for i in range(num_frames):
        average += locval_data[i * num_grid : (i + 1) * num_grid, axis]
    average /= num_frames
    
    result = np.transpose(np.vstack((grid_inform, average)))
    np.savetxt(fout_name, result)
    
if __name__ == "__main__":
    argvs = sys.argv
    argc = len(argvs)
    if argc != 2:
        print "Usage:"
        print "$ python %s dir_name" % argvs[0]
        quit()
    cdir = argvs[1]
    fin_name  = os.path.join(cdir, "mlocval.txt")
    fout_name = os.path.join(cdir, "hyphil_mean.txt")
    main(fin_name, 3, fout_name)


