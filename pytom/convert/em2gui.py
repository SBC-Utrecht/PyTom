#!/usr/bin/env pytom

from pytom.basic.files import read
from pytom.basic.files import read_em_header
import os, sys
from pytom.bin.em2mrc import em2mrc


if __name__ == '__main__':

    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert em file to mrc.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f','--file'], 'Filename', True, True),
                                   ScriptOption(['-d','--directory'], 'A directory of files.', True, True),
                                   ScriptOption(['-t','--target'], 'Path to output folder.', True, True),
                                   ScriptOption(['-m','--mdocname'], 'mdoc name.', True, True),
                                   ScriptOption(['-p','--prefix'], 'Prefix to filename.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    if 1:
        filename, directory, target, mdocname, prefix, help = parse_script_options(sys.argv[1:], helper)
    else:#except Exception as e:
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

        out = open(os.path.join(target, mdocname),'w')

        nr = 0

        d = '''[ZValue = {}]
TiltAngle = {}
SubFramePath = X:\{}

'''        
        fileList = os.listdir(directory)
        print('files: ' , fileList)
        for fname in sorted(fileList):
            if fname.endswith('.em'):
                if prefix and not fname.startswith(prefix): continue
                print(directory + os.sep + fname , target)
                em2mrc(os.path.join(directory,fname), target)
                header = read_em_header(os.path.join(directory,fname))
                out.write(d.format(nr, header.get_tiltangle(), fname.replace('.em','.mrc')))
                nr +=1
        
