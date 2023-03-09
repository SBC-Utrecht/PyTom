#!/usr/bin/env python

'''
Created on Mar 19, 2012

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Cut out the subtomograms from another tomogram given the particle list and the binning factor.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-v'], 'Volume.', True, False),
                                   ScriptOption(['-p'], 'Particle list.', True, False),
                                   ScriptOption(['-b'], 'Binning factor.', True, False),
                                   ScriptOption(['-s'], 'Size in each dimension (in radius).', True, False),
                                   ScriptOption(['-w'], 'Wedge angle.', True, False),
                                   ScriptOption(['-d'], 'Destination directory.', True, True),
                                   ScriptOption(['-r'], 'Resulting particle list for alignment.', True, False),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        vol_filename, pl_filename, binning, radius, w, dest_dir, res_name, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    binning = int(binning)
    radius = int(radius)
    w = int(w)
    if not dest_dir:
        dest_dir = '.'
    
    from pytom.lib.pytom_volume import read, subvolume
    v = read(vol_filename)
    
    from pytom.basic.structures import ParticleList, Particle, WedgeInfo
    pl = ParticleList("./")
    pl.fromXMLFile(pl_filename)
    
    def regulaize(xx, dim):
        if xx*binning-radius < 0:
            if 2*radius > dim:
                raise Exception("Volume size to be cut is too big!")
            return 0
        if xx*binning+radius > dim:
            if dim-2*radius < 0:
                raise Exception("Volume size to be cut is too big!")
            return dim-2*radius
        return xx*binning-radius
    
    res = ParticleList(dest_dir)
    for p in pl:
        x,y,z = p.getPickPosition().toVector()
        x = regulaize(int(x), v.sizeX()); y = regulaize(int(y), v.sizeY()); z = regulaize(int(z), v.sizeZ())
        new_vol = subvolume(v, x, y, z, 2*radius, 2*radius, 2*radius)
        name = dest_dir+'/'+p.getFilename()
        new_vol.write(name) # write the subtomograms to the disk
        res.append(Particle(name, p.getRotation(), None, WedgeInfo(w), 1, p.getPickPosition(), p.getScore())) # append it to the particle list for alignment
    
    res.toXMLFile(dest_dir+'/'+res_name)
