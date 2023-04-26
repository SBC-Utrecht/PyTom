#!/usr/bin/env python
'''
Created on Jul 21, 2011

@author: hrabe
'''





if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.filter import bandpassFilter
    from pytom.lib.pytom_volume import read
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Will filter (band/low/high pass) a file. Choose values between 0 and 1.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f','--file'], 'Filename', True, True),
                                   ScriptOption(['-t','--target'], 'Name of filtered file', True, True),
                                   ScriptOption(['-l','--lowestFrequency'], 'The lowest frequency. 0 if not set (in pixels)', True, False),
                                   ScriptOption(['-h','--highestFrequency'], 'The highest frequency. Volume size / 2 if not set (in pixels)', True, True),
                                   ScriptOption(['-s','--smooth'], 'Smoothing of bandpass (in voxels). 0 if not set.', True, True),
                                   ScriptOption(['--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        filename, target, lowestFrequency, highestFrequency, smooth, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    if not filename or not target:
        print(helper)
        sys.exit()

    v = read(filename)
    
    if lowestFrequency:
        lowestFrequency = int(lowestFrequency)
    else:
        lowestFrequency = 0
    
    if highestFrequency:
        highestFrequency = int(highestFrequency)
    else:
        highestFrequency = v.sizeX()/2
    
    if smooth:
        smooth = int(smooth)
    else:
        smooth = 0 
        
    
    r = bandpassFilter(v, lowestFrequency, highestFrequency, None, smooth)
    
    r[0].write(target)
    
    
