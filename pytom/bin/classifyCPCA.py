#!/usr/bin/env python
'''
Created on Jan 21, 2013

@author: FF
'''

from pytom.classification.CPCAfunctions import *

def doClassification(pl, cpl, ccc, neig, nclass, cName):
    """
    """
    subTomoClust(particleListFilename=pl, classifiedParticleListFilename=cpl,
            cccName=ccc, neig=neig, nclass=nclass)
    averageClasses(particleListFilename=cpl, avName=cName)



if __name__ == '__main__':
    import sys, os
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Classification of particle list using CPCA. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/classification.html',
                          authors='Friedrich Foerster',
                          options=[ScriptOption(['-p','--particleList'], 'ParticleList', True, True),
                                   ScriptOption(['-o','--outputParticleList'], 'classified Particle List.', True, True),
                                   ScriptOption(['-c','--ccc'], 'constrained correlation matrix.', True, True),
                                   ScriptOption(['-e','--neig'],'number of eigenvectors.', True, True),
                                   ScriptOption(['-n','--nclass'],'number of classes.', True, True),
                                   ScriptOption(['-a','--average'], 'name for class averages.', True, True),
                                   ScriptOption(['-t', '--directory'], 'Target Directory.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        pl, cpl, ccc, neig, nclass, cName, directory, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    neig = int(neig)
    nclass = int(nclass)
    if not cName:
        cName = 'class'
    directory = './' if directory is None else directory

    cName = cName if cName[-4:] == '.xml' else cName + '.xml'
    cName = os.path.join(directory, cName)
    cpl = os.path.join(directory, cpl)
    doClassification(pl=pl, cpl=cpl, ccc=ccc, neig=neig, nclass=nclass, cName=cName)
 
