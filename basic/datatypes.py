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
