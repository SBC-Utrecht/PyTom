#!/usr/bin/env pytom

def mrc2em(filename,destination):
    #from pytom_volume import read
    from pytom.basic.files import read
    from pytom.tools.files import checkFileExists,checkDirExists
    import os
    if not checkFileExists(filename):
        raise RuntimeError('MRC file not found! ',filename)

    if not checkDirExists(destination):
        raise RuntimeError('Destination directory not found! ', destination)

    emfile = read(filename)

    splitName = filename.split(os.sep)
    filename = splitName[len(splitName)-1]


    newFilename = destination + os.sep + filename[0:len(filename)-4] + '.em'

    emfile.write(newFilename,'em')

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert mrc file to em.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f','--file'], 'Filename', True, True),
                                   ScriptOption(['-d','--directory'], 'A directory of files.', True, True),
                                   ScriptOption(['-t','--targetPath'], 'Path to new file.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        filename, directory, target, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    if filename:
        #convert only one file
        mrc2em(filename,target)
    elif directory:
        import os
        
        fileList = os.listdir(directory)
        for file in fileList:
            if file[-4:] == '.mrc':
                print(directory + os.sep + file , target)
                mrc2em(directory + os.sep + file,target)
