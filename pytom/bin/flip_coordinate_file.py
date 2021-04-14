#!/usr/bin/env pytom

import os,sys,glob

def flip(fname,target, size, postfix):
    data = [line for line in open(fname).readlines()]
    a = os.path.basename(fname)
    out = open(os.path.join(target, a[:-4])+(len(postfix)>0)*'_' + '{}.txt'.format(postfix), 'w')

    for line in data:
        if '#' in line or  len(line.split()) != 3:
            out.write(line)
        else:
            x,y,z = map(float, line.split())
            x = size -x
            y = size-y
            z = size - z
            out.write('{:8.0f} {:8.0f} {:8.0f}\n'.format(x,y,z)) 
            
    out.close()

if __name__ == '__main__':

    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert em file to mrc.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f','--file'], 'Filename', True, True),
                                   ScriptOption(['-d','--directory'], 'A directory of files.', True, True),
                                   ScriptOption(['-t','--target'], 'Path to output folder.', True, True),
                                   ScriptOption(['-s','--size'], 'mdoc name.', True, True),
                                   ScriptOption(['-p','--postfix'], 'Prefix to filename.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        filename, directory, target, size, postfix, help = parse_script_options(sys.argv[1:], helper)
    except:#except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()

    try: size =int(size)
    except: size = 464

    if filename:
        #convert only one file
        flip(filename,target, size, postfix)
    elif directory:
        import os

        fileList = os.listdir(directory)
        
        for fname in sorted(fileList):
            if fname.endswith('.txt'):
                print(directory + os.sep + fname , target)
                flip(os.path.join(directory,fname), target, size, postfix)
                
