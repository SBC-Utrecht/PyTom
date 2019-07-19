'''
Created on Oct 2, 2012

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Plot angular distribution around a certain axis.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-p'], 'Particle list.', True, False),
                                   ScriptOption(['-a'], 'Euler angle specifying the axis in ZXZ convention', True, False),
                                   ScriptOption(['-c'], 'cut out range', True, True),
                                   ScriptOption(['-o'], 'Output file', True, False),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        pl_filename, axis_angles, cut_out_range, output, bHelp = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
    if bHelp is True:
        print helper
        sys.exit()

    if cut_out_range is not None:
        cut_out_range = [float(i) for i in cut_out_range.split(',')]

    from pytom.basic.structures import Rotation, ParticleList
    pl = ParticleList()
    pl.fromXMLFile(pl_filename)

    angles = []
    for angle in axis_angles.split('/'):
        z1, x, z2 = [float(i) for i in angle.split(',')]
        angles.append(Rotation(z1, z2, x))
    
    from pytom.tools.maths import rotation_distance
    dist = []
    for p in pl:
        rot = p.getRotation()
        d = 180
        for axis in angles:
            dd = rotation_distance(rot, axis)
            if dd < d:
                d = dd
        dist.append(d)
    
    # cut out only the things in the range
    if cut_out_range is not None:
        for i in xrange(len(dist)):
            if dist[i] >= cut_out_range[0] and dist[i] <= cut_out_range[1]:
                dist[i] = 1
            else:
                dist[i] = 0
    
    # write to the disk
    f = open(output, 'w')
    for d in dist:
        f.write(str(d)+'\n')
    f.close()
    