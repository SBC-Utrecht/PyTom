'''
Created on Dec 5, 2013

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import Particle, ParticleList, SingleTiltWedge
    import os
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Plot angular distribution of a particle list',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-p'], 'Particle list', True, False),
                                   ScriptOption(['-c'], 'Class label to plot', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        pl_filename, class_label, bHelp = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    if bHelp is True:
        print helper
        sys.exit()

    # load the particle list
    pl = ParticleList()
    pl.fromXMLFile(pl_filename)
    
    top_longitude = [] # top sphere
    top_latitude = []
    bottom_longitude = [] # bottom sphere
    bottom_latitude = []
    
    from math import pi
    for p in pl:
        if class_label is not None:
            if str(p.getClass()) != class_label:
                continue
        
        rotation = p.getRotation()
        z2 = rotation.getZ2()
        x = rotation.getX()
        
        # normalize the angles
        if x >= 0 and x <= 180:
            pass
        elif x < 0 and x > -180:
            x = -x
            z2 += 180 # z1 should also plus 180, but in this case it is unused
        elif x > 180 and x < 360:
            x = 360-x
            z2 += 180
        else:
            print 'Ignore particle: ', p.getFilename()
        
        while z2 < 0 or z2 > 360:
            if z2 < 0:
                z2 += 360
            else:
                z2 -= 360
                
        if x <= 90: # top sphere
            top_longitude.append(z2/180.*pi) # transfer to radian
            top_latitude.append(x)
        else: # bottom sphere
            bottom_longitude.append(z2/180.*pi)
            bottom_latitude.append(180-x)
    
    # start plotting
    import numpy as np
    import matplotlib.pyplot as plt
    
    ax = plt.subplot(121, polar=True)
    ax.plot(top_longitude, top_latitude, 'ro')
    ax.set_rmax(90.)
    ax.grid(True)
    ax.set_title("Top View", va='bottom')
    
    ax = plt.subplot(122, polar=True)
    ax.plot(bottom_longitude, bottom_latitude, 'ro')
    ax.set_rmax(90.)
    ax.set_yticks(range(0, 100, 10)) # change the yticks for the bottom sphere
    ax.set_yticklabels(('180', '170', '160', '150', '140', '130', '120', '110', '100', '90'))
    ax.grid(True)
    ax.set_title("Bottom View", va='bottom')
    
    plt.show()
    