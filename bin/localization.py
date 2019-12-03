#!/usr/bin/env pytom

def startLocalizationJob(filename, splitX=0, splitY=0, splitZ=0, doSplitAngles=False, gpuID=None):
    """
    @author: chen
    """
    verbose=True
    from pytom.localization.peak_job import PeakJob
    import os
    job = PeakJob()
    job.fromXMLFile(filename)
    job.check()
    suffix = os.path.basename(job.reference.getFilename()).split('.')[0]
    print(f'suffix: {suffix}')
    if doSplitAngles:
        print('Ignore split volume parameters ...')
        from pytom.localization.parallel_extract_peaks import PeakManager
        manager = PeakManager(suffix=suffix)
        manager.parallelStart_splitAng(job, verbose)
    else:
        from pytom.localization.parallel_extract_peaks import PeakLeader
        leader = PeakLeader(suffix=suffix)
        leader.parallelRun(job, splitX, splitY, splitZ, verbose, gpuID=gpuID)

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run a localization job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/localization.html',
                          authors='Yuxiang Chen, Thomas Hrabe',
                          options=[ScriptOption(['-j','--jobName'], 'Specify job.xml filename', arg=True, optional=False),
                                   ScriptOption(['-x','--splitX'], 'Parts you want to split the volume in X dimension', arg=True, optional=True),
                                   ScriptOption(['-y','--splitY'], 'Parts you want to split the volume in Y dimension', arg=True, optional=True),
                                   ScriptOption(['-z','--splitZ'], 'Parts you want to split the volume in Z dimension', arg=True, optional=True),
                                   ScriptOption(['-g', '--gpuID'], 'gpu index for running job',
                                                arg=True, optional=True),

                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    
    if 1:
        jobName, splitX, splitY, splitZ, gpuID, b_help = parse_script_options(sys.argv[1:], helper)
        
        if b_help is True:
            print(helper)
            sys.exit()
        
        if splitX is None:
            splitX = 0
        else:
            splitX = int(splitX)
        if splitY is None:
            splitY = 0
        else:
            splitY = int(splitY)
        if splitZ is None:
            splitZ = 0
        else:
            splitZ = int(splitZ)
        
        if jobName is None:
            raise RuntimeError()

        if gpuID is None:
            gpuID = None
        else:
            gpuID = list(map(int,gpuID.split(',')))

    else: # backward compatibility
        if len(sys.argv) == 2 or len(sys.argv) == 5:
            pass
        else:
            print(helper)
        
        jobName = sys.argv[1]
        
        if len(sys.argv) == 5:
            splitX = int(sys.argv[2])
            splitY = int(sys.argv[3])
            splitZ = int(sys.argv[4])

    
    from pytom.tools.timing import Timing
    t = Timing(); t.start()
    
    startLocalizationJob(jobName, splitX, splitY, splitZ, doSplitAngles=False, gpuID=gpuID)
    
    time = t.end(); print('The overall execution time: %f' % time)
    
    