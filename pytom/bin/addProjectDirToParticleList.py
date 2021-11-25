#!/usr/bin/env pytom

def addProjectDir(filename, projectDir, outname):

    from pytom.basic.structures import ParticleList

    pl = ParticleList()
    pl.fromXMLFile(filename)

    for particle in pl:
        s = particle.getInfoGUI()
        s.setProjectDir(projectDir)


    pl.toXMLFile(outname)

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='''Adds project directory to the particles in the particle list''',
        authors='GvdS',
        options=[ScriptOption2(
                     ['-p', '--particleList'], 'Filename of particle list', 'file', 'required'),
                 ScriptOption2(
                     ['-d', '--projectDirectory'], 'Path of project directory.', 'directory', 'optional', ''),
                 ScriptOption2(
                     ['-o', '--outputName'], 'File name of output file', 'string', 'required', 'particleList.xml')])

    filename, directory, outname = parse_script_options2(sys.argv[1:], helper)

    addProjectDir(filename, directory, outname)

