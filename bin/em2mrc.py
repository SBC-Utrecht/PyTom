#!/usr/bin/env pytom

def em2mrc(filename,target):
	from pytom_volume import read
	from pytom.tools.files import checkFileExists,checkDirExists
	import os
	
	if not checkFileExists(filename):
		raise RuntimeError('EM file not found! ',filename)

	if not checkDirExists(target):
		raise RuntimeError('Destination directory not found! ', target)

	emfile = read(filename)
	
	splitName = filename.split(os.sep)
	filename = splitName[len(splitName)-1]
	
	
	newFilename = target + os.sep + filename[0:len(filename)-3] + '.mrc'

	emfile.write(newFilename,'mrc')


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert em file to mrc.',
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
        em2mrc(filename,target)
    elif directory:
        import os
        
        fileList = os.listdir(directory)
        for file in fileList:
            if file[len(file)-3:len(file)] == '.em':
                #print(directory + os.sep + file , target)
                em2mrc(directory + os.sep + file,target)
				
						
		
