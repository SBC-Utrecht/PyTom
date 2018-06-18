#!/usr/bin/env pytom
'''
Created on Nov 3, 2011

@author: hrabe
'''

if __name__ == '__main__':
    # parse command line arguments
    
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options 
    from pytom.alignment.ExMaxAlignment import parallelStart,ExMaxJob 
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run an alignment job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/alignment.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-j','--job'], 'Job', arg=True, optional=False),
                                   ScriptOption(['-v','--verbose'], 'Verbose', arg=False,optional=True),
                                   ScriptOption(['-h', '--help'], 'Help.', arg=False, optional=True)])

    verbose = False
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        jobFile, verbose, helpme = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    if helpme is True:
        print helper
        sys.exit()
    
    # fix the verbose bug
    if verbose is None:
        verbose = False
    
    #from pytom.tools.timing import Timing
    #t = Timing()
    #t.start()

    exj = ExMaxJob()
    exj.fromXMLFile(jobFile)

    
    parallelStart(exj,verbose = verbose)
    
    #print 'Overall execution time: %f s.' % t.end()