DATATYPE_ALIGNMENT_RESULTS_RO = [('TiltAngle',       'f4'),
                              ('Magnification',   'f4'),
                              ('InPlaneRotation', 'f4'),
                              ('AlignmentTransX', 'f4'),
                              ('AlignmentTransY', 'f4'),
                              ('AxisAngle', 'f4'),
                              ('OperationOrder',  'U1000'),
                              ('FileName', 'U1000')]

HEADER_ALIGNMENT_RESULTS_RO = ''
unitsAlignmentResultsRo = ['', 'degrees', 'degrees', 'px', 'px', 'degrees', '', '']
fmtAlignmentResultsRo = FMT_ALIGNMENT_RESULTS_RO ='%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %5s %s'
for n, h in enumerate(DATATYPE_ALIGNMENT_RESULTS_RO):
    HEADER_ALIGNMENT_RESULTS_RO += '{} {}\n'.format(h[0], '({})'.format(unitsAlignmentResultsRo[n])*(unitsAlignmentResultsRo[n]!=''))


DATATYPE_0 = [('DefocusU', 'f4'),
              ('DefocusV', 'f4'),
              ('DefocusAngle', 'f4'),
              ('Voltage', 'i4'),
              ('SphericalAberration', 'f4'),
              ('AmplitudeContrast', 'f4'),
              ('PhaseShift', 'f4'),
              ('PixelSpacing', 'f4'),
              ('MarkerDiameter', 'i4'),
              ('TiltAngle', 'f4'),
              ('RotationTheta', 'f4'),
              ('InPlaneRotation', 'f4'),
              ('TranslationX', 'f4'),
              ('TranslationY', 'f4'),
              ('TranslationZ', 'f4'),
              ('Magnification', 'f4'),
              ('Intensity', 'f4'),
              ('FileName', 'U1000')]

HEADERTEXT_0 = ''
units0 = ['um', 'um', 'deg', 'kV', 'mm', '', 'deg', 'A', 'A', 'deg', 'deg', 'deg', 'px', 'px', 'px', '', '', '' ]
fmt0='%11.6f %11.6f %6.2f %4d %6.2f %4.2f %11.6f %11.6f %4d %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %5.3f %5.3f %s'
for n, h in enumerate(DATATYPE_0):
    HEADERTEXT_0 += '{} {}\n'.format(h[0], '({})'.format(units0[n])*(units0[n]!=''))


DATATYPE_METAFILE = [('DefocusU', 'f4'),
                     ('DefocusV', 'f4'),
                     ('DefocusAngle', 'f4'),
                     ('Voltage', 'i4'),
                     ('SphericalAberration', 'f4'),
                     ('AmplitudeContrast', 'f4'),
                     ('PhaseShift', 'f4'),
                     ('PixelSpacing', 'f4'),
                     ('MarkerDiameter', 'i4'),
                     ('TiltAngle', 'f4'),
                     ('RotationTheta', 'f4'),
                     ('InPlaneRotation', 'f4'),
                     ('TranslationX', 'f4'),
                     ('TranslationY', 'f4'),
                     ('TranslationZ', 'f4'),
                     ('Magnification', 'f4'),
                     ('Intensity', 'f4'),
                     ('ImageSize', 'i4'),
                     ('AcquisitionOrder', 'i4'),
                     ('FileName', 'U1000')]

HEADER_METAFILE = ''
unitsMetaFile = ['um', 'um', 'deg', 'kV', 'mm', '', 'deg', 'A', 'A', 'deg', 'deg', 'deg', 'px', 'px', 'px', '', '','px', '', '' ]
FMT_METAFILE='%11.6f %11.6f %6.2f %4d %6.2f %4.2f %11.6f %11.6f %4d %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %5.3f %5.3f %4d %3d %s'

for n, h in enumerate(DATATYPE_METAFILE):
    HEADER_METAFILE += '{} {}\n'.format(h[0], '({})'.format(unitsMetaFile[n])*(unitsMetaFile[n]!=''))



DATATYPE_MARKER_RESULTS = [('MarkerIndex', 'i4'),
                           ('OffsetX',     'f4'),
                           ('OffsetY',     'f4'),
                           ('OffsetZ',     'f4'),
                           ('PositionX',   'f4'),
                           ('PositionY',   'f4'),
                           ('PositionZ',   'f4')]

HEADER_MARKER_RESULTS = ''
unitsMarkerResults = ['', 'px', 'px', 'px', 'px', 'px', 'px']
fmtMarkerResults='%3d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f'
for n, h in enumerate(DATATYPE_MARKER_RESULTS):
    HEADER_MARKER_RESULTS += '{} {}\n'.format(h[0], '({})'.format(unitsMarkerResults[n])*(unitsMarkerResults[n]!=''))


DATATYPE_ALIGNMENT_RESULTS = [('AlignmentTransX', 'f4'),
                              ('AlignmentTransY', 'f4'),
                              ('TiltAngle',       'f4'),
                              ('InPlaneRotation', 'f4'),
                              ('Magnification',   'f4'),
                              ('FileName', 'U1000')]

HEADER_ALIGNMENT_RESULTS = ''
unitsAlignmentResults = ['px', 'px', 'degrees', 'degrees', '', '']
fmtAlignmentResults = FMT_ALIGNMENT_RESULTS ='%15.10f %15.10f %15.10f %15.10f %15.10f %s'
for n, h in enumerate(DATATYPE_ALIGNMENT_RESULTS):
    HEADER_ALIGNMENT_RESULTS += '{} {}\n'.format(h[0], '({})'.format(unitsAlignmentResults[n])*(unitsAlignmentResults[n]!=''))


DATATYPE_MARKERFILE = [('MarkerIndex', 'i4'),
                       ('TiltAngle', 'f4'),
                       ('PositionX', 'f4'),
                       ('PositionY', 'f4')]

HEADER_MARKERFILE = ''
unitsMarkerfile = ['', 'degrees', 'px', 'px']
FMT_MARKERFILE ='%3d %8.3f %8.2f %8.2f'
for n, h in enumerate(DATATYPE_MARKERFILE):
    HEADER_MARKERFILE += '{} {}\n'.format(h[0], '({})'.format(unitsMarkerfile[n])*(unitsMarkerfile[n]!=''))

LOCAL_ALIGNMENT_RESULTS = [('ParticleIndex', 'i4'),
                           ('AlignmentTransX', 'f4'),
                           ('AlignmentTransY', 'f4'),
                           ('TiltAngle', 'f4'),
                           ('InPlaneRotation', 'f4'),
                           ('Magnification', 'f4'),
                           ('FileName', 'U1000')]

headerLocalAlignmentResults = ''
unitsLAR = ['', 'px', 'px', 'degrees', 'degrees', '', '']
fmtLAR = '%7d %15.3f %15.3f %15.3f %15.3f %15.10f %s'
for n, h in enumerate(LOCAL_ALIGNMENT_RESULTS):
    headerLocalAlignmentResults += '{} {}\n'.format(h[0], '({})'.format(unitsLAR[n]) * (unitsLAR[n] != ''))


RELION31_PICKPOS_STAR = [('CoordinateX', 'f4'),
                         ('CoordinateY', 'f4'),
                         ('CoordinateZ', 'f4'),
                         ('MicrographName', 'U1000'),
                         ('Magnification', 'f4'),
                         ('DetectorPixelSize', 'f4'),
                         ('GroupNumber', 'i4'),
                         ('AngleRot', 'f4'),
                         ('AngleTilt', 'f4'),
                         ('AnglePsi', 'f4')]

headerRelion31Subtomo = '\ndata_\n\nloop_\n\n'
for n, h in enumerate(RELION31_PICKPOS_STAR):
    headerRelion31Subtomo += '_rln{} #{}\n'.format(h[0], n+1)
headerRelion31Subtomo = headerRelion31Subtomo[:-1]


fmtR31S = '%-9.4f %-9.4f %-9.4f %s %9.4f %9.4f %4d %9.4f %9.4f %9.4f'

DATATYPE_PROJECT_ALIGN_RESULTS = [('TomogramName', 'U1000'),
                                  ('AlignmentMethod', 'U1000'),
                                  ('OriginFolder', 'U1000'),
                                  ('AlignmentScore', 'f4'),
                                  ('FirstAngle', 'f4'),
                                  ('LastAngle', 'f4'),
                                  ('ReferenceImageID', 'f4'),
                                  ('ReferenceMarkerID', 'f4'),
                                  ('ExpectedRotationAngle', 'f4'),
                                  ('DeterminedRotationAngle', 'f4')]


DATATYPE_TASOLUTION = [('View', 'i4'),
                        ('Rotation', 'f4'),
                        ('Tilt', 'f4'),
                        ('Deltilt', 'f4'),
                        ('Mag', 'f4'),
                        ('Dmag', 'f4'),
                        ('Skew', 'f4'),
                        ('MeanResid', 'f4')]

HEADER_TASOLUTION = ''
unitsTaSolutionsFile = ['', 'deg', 'deg', 'deg', '', '', '', 'px' ]
FMT_TASOLTUION='%4d %10.1f %10.1f %10.2f %10.4f %10.4f %10.2f %10.2f'

for n, h in enumerate(DATATYPE_TASOLUTION):
    HEADER_TASOLUTION += '{} {}\n'.format(h[0], '({})'.format(unitsTaSolutionsFile[n])*(unitsTaSolutionsFile[n]!=''))


SIMULATED_GROUND_TRUTH = [('ParticleName', 'U1000'),
                          ('x', 'f4'),
                          ('y', 'f4'),
                          ('z', 'f4'),
                          ('ThetaZ', 'f4'),
                          ('PhiX', 'f4'),
                          ('PsiZ', 'f4')]


def generate_default_alignmentresults(folder, metafilename='', prefix='sorted_', filename='alignmentResults.txt'):
    '''This function creates an alignmentResults.txt file and returns a structured array,
    alignmentTransX=alignmentTransY=inPlaneRotation=0, magnification=1
    @param folder: all files found in folder are used.
    @param metafilename: if file is supplied, angles from this file are used as TiltAngle
    @param prefix: only files starting with prefix will be used.. Default sorted_
    @param filename: if filename is supplied, array will be written to file. Default:'alignmentResults.txt'
    @return structured array
    '''
    import os,numpy
    from pytom.basic.files import loadtxt, savetxt

    files = [os.path.join(folder, fname) for fname in os.listdir(folder) if fname.startswith(prefix)]
    data= numpy.zeros((len(files)),dtype=DATATYPE_ALIGNMENT_RESULTS)
    data['Magnification'] = 1
    data['FileName'] = numpy.array(files,dtype=numpy.str)

    if metafilename:
        metadata = loadtxt(metafilename,dtype=DATATYPE_METAFILE)
        m = 0
        for nn, fname in enumerate(data['FileName']):
            if os.path.basename(fname) in files:
                data['TiltAngle'][m] = metadata['TiltAngle'][nn]
                m+=1

    if filename:
        outname = os.path.join(folder,filename)
        savetxt(outname, data, fmt=FMT_ALIGNMENT_RESULTS, header=HEADER_ALIGNMENT_RESULTS)

    return data


# Relion related datatypes

RELION31_DATA_MICROGRAPH_STAR = [ ('MicrographNameNoDW', 'U30'),
                                  ('MicrographName', 'U30'),
                                  ('OpticsGroup', 'i4'),
                                  ('CtfImage', 'U30'),
                                  ('DefocusU', 'f4'),
                                  ('DefocusV', 'f4'),
                                  ('CtfAstigmatism', 'f4'),
                                  ('DefocusAngle', 'f4'),
                                  ('CtfFigureOfMerit', 'f4'),
                                  ('CtfMaxResolution', 'f4')]

RELION31_OPTICS_GROUP = [('OpticsGroupName', 'U200'),
                         ('OpticsGroup', 'i4'),
                         ('MtfFileName', 'U200'),
                         ('MicrographOriginalPixelSize', 'f4'),
                         ('SphericalAberration', 'f4'),
                         ('AmplitudeContrast', 'f4'),
                         ('MicrographPixelSize', 'f4')
                         ]

DATATYPE_RELION31_EXTENDED = [('CoordinateX', 'f4'),
                              ('CoordinateY', 'f4'),
                              ('CoordinateZ', 'f4'),
                              ('MicrographName', 'U1000'),
                              ('Magnification', 'f4'),
                              ('DetectorPixelSize', 'f4'),
                              ('GroupNumber', 'i4'),
                              ('AngleRot', 'f4'),
                              ('AngleTilt', 'f4'),
                              ('AnglePsi', 'f4'),
                              ('CtfMaxResolution', 'f4'),
                              ('ImageName', 'U1000'),
                              ('CtfImage', 'U1000'),
                              ('PixelSize', 'f4'),
                              ('OriginXAngst', 'f4'),
                              ('OriginYAngst', 'f4'),
                              ('OriginZAngst', 'f4'),
                              ('ClassNumber', 'i4'),
                              ('Voltage', 'f4'),
                              ('SphericalAberration', 'f4')]

headerRelion31EXTENDEDSubtomo = '\ndata_\n\nloop_\n\n'
for n, h in enumerate(DATATYPE_RELION31_EXTENDED):
    headerRelion31EXTENDEDSubtomo += '_rln{} #{}\n'.format(h[0], n+1)
headerRelion31EXTENDEDSubtomo = headerRelion31EXTENDEDSubtomo[:-1]

fmtR31EXTENDED = ('%-9.4f %-9.4f %-9.4f %s %.1f %.5f %4d %9.4f %9.4f %9.4f %.1f %s %s %9.4f %9.4f %9.4f %9.4f '
                  '%6d %.3f %.3f')

datatype_list = [DATATYPE_0,DATATYPE_MARKERFILE, DATATYPE_METAFILE, DATATYPE_MARKER_RESULTS, DATATYPE_TASOLUTION,
                 DATATYPE_PROJECT_ALIGN_RESULTS, DATATYPE_ALIGNMENT_RESULTS_RO,
                 RELION31_PICKPOS_STAR, DATATYPE_RELION31_EXTENDED]