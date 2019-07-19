#!/usr/bin/env pytom

"""
Created on Aug 24, 2011

@author: hrabe
"""



if __name__ == '__main__':
# parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList,Rotation,PointSymmetry
    from pytom_volume import read
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Generate a symmetrized density.',
                          authors='Thomas Hrabe',
                          options= [ScriptOption(['-p','--particleList'], 'Aligned XML particle list.', True, True),
                                    ScriptOption(['-v','--volume'], 'Volume to be symmetrized.', True, True),
                                    ScriptOption(['-r','--result'], 'Resulting symmetry filename.', True, False),
                                    ScriptOption(['-s','--symmetries'], 'How many symmetries (integer)', True, False),
                                    ScriptOption(['--z1Rotation'], 'First rotation around z axis', True, True),
                                    ScriptOption(['--xRotation'], 'Rotation around x axis', True, True),
                                    ScriptOption(['--z2Rotation'], 'Second rotation around z axis', True, True),
                                    ScriptOption(['--help'], 'Print help.', False, False)])
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    
    try:
        particleListName, volumeName, result, numberSymmetries, z1 , x, z2, help = parse_script_options(sys.argv[1:], helper)        
    except Exception as e:
        print e
        sys.exit()
        
    if particleListName:
        pl = ParticleList('.')
        pl.fromXMLFile(particleListName)
    elif volumeName:
        volume = read(volumeName)
    else:        
        raise RuntimeError('You must specify either a particle list or a volume file for symmetrization.')
    
    if not z2:
        z2 = 0
    if not x:
        x = 0
        
    symmetry = PointSymmetry(numberSymmetries,z2,x)
    
    if particleListName:
        newParticleList = symmetry.apply(pl)
        newParticleList.average(resul)
    elif volumeName:
        newVolume = symmetry.applyToParticle(volume,Rotation(0,z2,x))
        newVolume.write(result)