import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path

def istogramma(x,M,colore='g',grandezza='velocity'):
    fig1, ax1 = plt.subplots(1, 1, figsize=(8,4), dpi=100) 
# histogram our data with numpy
    n, bins = np.histogram(x, 100, density=True)
# get the corners of the rectangles for the histogram
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n
# we need a (numrects x numsides x 2) numpy array for the path helper
# function to build a compound path
    XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T
# get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)
# make a patch out of it
    patch = patches.PathPatch(barpath, facecolor=colore, alpha=0.4)
    ax1.add_patch(patch)
# update the view limits
    ax1.set_xlim(left[0], right[-1])
    ax1.set_ylim(bottom.min(), top.max()+.03)
    #ax1.axhline(1.,color='blue', alpha=0.5)
    plt.title(' istogramma %s con %d campioni' % (grandezza,M) )
    plt.show()
#