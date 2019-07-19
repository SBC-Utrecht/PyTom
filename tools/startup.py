'''
Created on Nov 9, 2010

@author: hrabe
'''


"""
def startup():

    print "Welcome to PyTom v0.1"
    from pytom_volume import *
    from pytom.basic import filter 
    from pytom.basic import fourier
"""


def pytom_help():
    """
    pytom_help: Prints list of commands available after launch
    """
    
    
    
    print 'List of available commands for the start. \nPlease have a look at the API for a complete documentation of PyTom functionality'
    print ''
    print 'This is an overview of key PyTom functions'
    print ''
    print 'You can also always type ?<command> for details about a specific command.'
    print ''
    print '-- Volumes'
    print ''
    print '    vol(x,y,z) : returns volume of size x,y,z      Voxels are NOT initialized'
    print '    vol_compl(x,y,z) : returns volume of size x,y,z Voxels are NOT initialized'
    print '    '
    print '    subvolume'
    print ''
    print '    Explore volume properties with volume.<TAB>'
    print '    '
    print '    for volume arithmetic simply use +,-,*,/'
    print ''
    print '-- IO'
    print ''
    print '    read(str filename) : Read file from disk. Supported formats are EM, MRC, CCP4'
    print ''
    print '    v.write(str filename, str filetype) : Writes a volume object to disk. filetype is optional. Set to \'EM\' \'MRC\' \'CCP4\' if needed.'
    print ''
    print '-- Arithmetic operations'
    print ''
    print '    abs(vol volume (vol_comp volumeComplex) )     : Returns the absolute voxel values of volume. Input can be complex, too.'
    print '    complexDiv(vol_comp volume, vol div)          : '
    print '    limit                                         : '
    print ''
    print '-- Transformations'
    print ''
    print '    rotate(vol volume,vol destination,Z1,Z2,X)        : Rotates volume around Z1,Z2,X angles.'
    print '                                                        Result is stored in the existing!!! destination of same size as volume. ' 
    print '                                                        Uses linear interpolation. Volumes must have same size!'
    print '    rotateCubic(vol volume,vol destination,Z1,Z2,X)   : Rotates volume around Z1,Z2,X angles.'
    print '                                                        Result is stored in the existing!!! destination of same size as volume. '
    print '                                                        Destination of same size as volume.'
    print '                                                        Uses cubic interpolation. Volumes must have same size!'
    print '    rotateSpline(vol volume,vol destination,Z1,Z2,X)  : Rotates volume around Z1,Z2,X angles.'
    print '                                                        Result is stored in the existing!!! destination of same size as volume. '
    print '                                                        Uses spline interpolation. Volumes must have same size!'
    print '    shift(vol volume,vold destination,x,y,z)          : Shifts volume by x,y,z and stores results in destination'
    print '                                                        Uses linear interpolation. Volumes must have same size!'
    print ''
    print '    transform()'
    print '    transformCubic()'
    print '    transformSpline()'
    print ''
    print '    rescale()'
    print '    rescaleCubic()'
    print '    rescaleSpline()'
    print ''
    print '-- Statistics'
    print ''
    print '    mean(vol volume)                                         : Returns mean value in volume'
    print '    variance(vol volume, bool use_sample_standard_deviation) : Returns variance of volume'
    print '                                                               If use_sample_standard_deviation true the denominator is N-1 instead of N (N number of voxels)'
    
    
    
    