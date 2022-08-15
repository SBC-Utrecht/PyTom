templateXML       = '''<JobDescription Destination="{d[6]}">
  <Volume Filename="{d[0]}" Subregion=" 0,0,{d[7]},{d[8]},{d[9]},{d[10]} "/>
  <Reference Weighting="" File="{d[1]}"/>
  <Mask Filename="{d[2]}" Binning="1" isSphere="True"/>
  <WedgeInfo Angle1="{d[3]}" Angle2="{d[4]}" CutoffRadius="0.0" TiltAxis="custom">
    <Rotation Z1="0.0" Z2="0.0" X="0.0"/>
  </WedgeInfo>
  <Angles Type="FromEMFile" File="{d[5]}"/>
  <Score Type="FLCFScore" Value="-100000000">
    <DistanceFunction Deviation="0.0" Mean="0.0" Filename=""/>
  </Score>
</JobDescription>'''

old = '''
<JobDescription>
  <Volume Filename="{d[0]}">
  </Volume>
  <Reference File="{d[1]}">
  </Reference>
  <Mask Filename="{d[2]}" Binning="1" isSphere="True">
  </Mask>
  <WedgeInfo Angle1="{d[3]}" Angle2="{d[4]}" CutoffRadius="0.0" TiltAxis="custom">
    <Rotation Z1="0.0" Z2="0.0" X="0.0">
    </Rotation>
  </WedgeInfo>
  </Angles>
  <Score Type="FLCFScore" Value="-100000000">
  </Score>
</JobDescription>
 
'''


templateAlignment = '''cd {d[0]}; 

generateAlignedTiltImages.py \\
    --tiltSeriesName {d[2]}  \\
    --firstIndex {d[3]} \\
    --lastIndex {d[4]} \\
    --referenceIndex {d[5]} \\
    --referenceMarkerIndex {d[6]} \\
    --markerFile {d[7]} \\
    --projectionTargets {d[8]} \\
    --projectionBinning {d[9]} \\
    --lowpassFilter 0.9 \\
    --weightingType {d[10]} \\
    --expectedRotationAngle {d[11]} \\
    --numberProcesses 1 '''


# --size {d[9]},{d[13]},{d14}
# switch between sorted and sorted_ctf
# 0 = tomoprocessing dir, 1 = projection dir, 6 = binning, 10, weighting, 7 = output name (tomogram_000),
# 9,13,14 = tomgoram size, 12 = specimen angle
templateWBP       = '''cd {d[0]}

reconstructWB.py \\
    --alignResultFile {d[15]} \\
    --tomogram reconstruction/WBP/{d[7]}_WBP.mrc \\
    --projBinning {d[6]} \\
    --applyWeighting {d[10]}  \\
    --tomogram reconstruction/WBP/{d[7]}_WBP.mrc \\
    --size {d[9]},{d[13]},{d[14]} \\
    {d[12]}


unlink ../../04_Particle_Picking/Tomograms/{d[7]}_WBP.mrc
ln -s {d[0]}/reconstruction/WBP/{d[7]}_WBP.mrc ../../04_Particle_Picking/Tomograms/{d[7]}_WBP.mrc'''


templateWBP_old    = '''cd {d[0]}

reconstructTomogram.py \\
    --tiltSeriesName sorted/sorted \\
    --firstIndex {d[2]} \\
    --lastIndex {d[3]} \\
    --referenceIndex {d[4]} \\
    --markerFile reconstruction/WBP/markerfile.txt \\
    --referenceMarkerIndex {d[5]} \\
    --expectedRotationAngle {d[11]} \\
    --projectionTargets reconstruction/WBP/temp_files_unweighted/sorted_aligned \\
    --projectionBinning {d[6]} \\
    --lowpassFilter 0.9  \\
    --weightingType {d[10]}  \\
    --tomogramFile reconstruction/WBP/{d[7]}_WBP.mrc \\
    --tiltSeriesFormat {d[8]} \\
    --fileType {d[8]}  \\
    --tomogramSizeX {d[9]}  \\
    --tomogramSizeY {d[13]} \\
    --tomogramSizeZ {d[14]} \\
    --reconstructionCenterX 0 \\
    --reconstructionCenterY 0 \\
    --reconstructionCenterZ 0 \\
    {d[12]}


unlink ../../04_Particle_Picking/Tomograms/{d[7]}_WBP.{d[8]}
ln -s {d[0]}/reconstruction/WBP/{d[7]}_WBP.{d[8]} ../../04_Particle_Picking/Tomograms/{d[7]}_WBP.{d[8]}'''


templateINFR      = '''cd {d[0]}

reconstructTomogram.py \\
    --tiltSeriesName reconstruction/INFR/sorted  \\
    --firstIndex {d[2]} \\
    --lastIndex {d[3]} \\
    --referenceIndex {d[4]} \\
    --markerFile reconstruction/INFR/markerfile.txt \\
    --referenceMarkerIndex {d[5]} \\
    --expectedRotationAngle {d[9]} \\
    --projectionTargets reconstruction/INFR/temp_files_unweighted/sorted_aligned \\
    --projectionBinning {d[6]} \\
    --lowpassFilter 0.9 \\
    --weightingType 0 \\
    --tiltSeriesFormat em \\
    --fileType em

unlink ../../04_Particle_Picking/Tomograms/{d[8]}_INFR.em
{d[1]}/bin/pytom {d[7]}/reconstruction/reconstruct_INFR.py -d reconstruction/INFR/temp_files_unweighted/ -o reconstruction/INFR/{d[8]}_INFR.em 
ln -s {d[0]}/reconstruction/INFR/{d[8]}_INFR.em ../../04_Particle_Picking/Tomograms/{d[8]}_INFR.em'''


templateFRMJob    = '''<FRMJob Destination='{d[15]}' BandwidthRange='[{d[0]},{d[1]}]' Frequency='{d[2]}' MaxIterations='{d[3]}' PeakOffset='{d[4]}' RScore='{d[5]}' WeightedAverage='{d[6]}' Binning='1'>
    <Reference PreWedge="" File="{d[7]}" Weighting="{d[8]}"/>
    <Mask Filename="{d[9]}" Binning="{d[10]}" isSphere="{d[11]}"/>
    <SampleInformation PixelSize="{d[12]}" ParticleDiameter="{d[13]}"/>
    <ParticleListLocation Path="{d[14]}"/>
</FRMJob>
'''


templateFRMSlurm  = '''
cd {d[0]}

mpiexec -n {d[3]} FRMAlignment.py -j {d[2]} -v

'''


templateGLocal    = '''cd {d[0]}

mpiexec -n {d[14]} GLocalJob.py \\
--particleList {d[3]} \\
--mask {d[5]} \\
--numberIterations {d[6]} \\
--pixelSize {d[7]} \\
--particleDiameter {d[8]} \\
--binning {d[9]} \\
--destination {d[11]} \\
--SphericalMask \\
--angleShells {d[12]} \\
--angleIncrement {d[13]} \\
--jobName {d[10]} \\
{d[4]}{d[15]}'''


templateCCC       = """cd {d[0]}

mpiexec --tag-output -n {d[7]} calculate_correlation_matrix.py -p {d[2]} -m {d[3]} -f {d[4]} -b {d[5]} -o {d[6]} {d[8]}
"""


templateCPCA      = """cd {d[0]}

classifyCPCA.py -p {d[2]} -o {d[7]}/{d[3]} -c {d[4]} -e {d[5]} -n {d[6]} -a {d[8]} -t {d[7]}
"""


templateAC        = '''cd {d[0]}

mpiexec --tag-output -n {d[13]} auto_focus_classify.py \\
-p {d[2]} \\
{d[3]} {d[4]} \\
-k {d[5]} \\
-f {d[6]} \\
-i {d[7]} \\
-s {d[8]} \\
-b {d[14]} \\
-n {d[9]} --sig {d[10]} -t {d[11]} \\
-o {d[12]}'''


templateTM        = '''cd {d[0]}

mpiexec --tag-output -n {d[1]} localization.py -j {d[3]} -x {d[4]} -y {d[5]} -z {d[6]} {d[7]}
'''


createParticleList = 'coords2PL.py -c {d[0]}  -s {d[1]} -w {d[2]},{d[3]} -p {d[4]} {d[5]}'


extractParticles = '''cd {d[8]}

reconstructWB.py --particleList {d[0]} \\
--projectionDirectory {d[1]} \\
--coordinateBinning {d[2]} \\
--size {d[3]} \\
--applyWeighting {d[9]} \\
--projBinning {d[4]} \\
--recOffset {d[5]},{d[6]},{d[7]} \\
--metafile {d[10]} \\
--numProcesses {d[11]} \\
{d[12]}'''


polishParticles = '''cd {d[0]}

mpiexec -n {d[1]} particlePolishingOrig.py \\
--particleList {d[3]} \\
--projectionDirectory {d[4]} \\
--template {d[5]} \\
--outputDirectory {d[6]} \\
--coordinateBinning {d[7]} \\
--maxParticleShift {d[8]} \\
--recOffset {d[9]},{d[10]},{d[11]} \\
{d[12]} {d[13]} {d[14]}
'''


extractParticlesClosestMarker = '''cd {d[8]}

reconstructWBClosestMarker.py --particleList {d[0]} \\
--projectionDirectory {d[1]} \\
--coordinateBinning {d[2]} \\
--size {d[3]} \\
--applyWeighting {d[9]} \\
--projBinning {d[4]} \\
--recOffset {d[5]},{d[6]},{d[7]} \\
--metafile {d[10]} \\
--logfileReconstruction {d[11]} \\
--numProcesses {d[12]} \\
--prefixFileName {d[13]} \\
'''


multiple_alignment = '''cd {d[0]}

multi_tilt_alignment.py \\
--start {d[2]} \\
--end {d[3]} \\
--tiltSeriesName {d[4]} \\
--projectionTargets {d[5]} \\
--fnames {d[6]}'''


templateExtractCandidates = '''cd {d[0]}

extractCandidates.py \\
--jobFile {d[2]} \\
--result {d[3]} \\
--orientation {d[4]} \\
--particleList {d[5]} \\
--particlePath {d[6]} \\
--size {d[7]} \\
--numberCandidates {d[8]} \\
--minimalScoreValue {d[9]} {d[10]}'''


ParamsFileCTFPlotter = '''#
InputStack           {d[0]}
AngleFile            {d[1]}
DefocusFile          {d[2]}
{d[15]}           {d[3]}
AxisAngle            {d[4]}
PixelSize            {d[5]}
ExpectedDefocus      {d[6]}
AngleRange           {d[7]} {d[8]}
Voltage              {d[9]}
SphericalAberration  {d[10]}
AmplitudeContrast    {d[11]}
DefocusTol           200
PSResolution         101
TileSize             256
LeftDefTol           2000.0
RightDefTol          2000.0
FindAstigPhaseCuton  {d[12]} {d[13]} {d[14]}
Offset 0.02
'''

templateCTFPlotter = '''cd {d[0]}

ctfplotter -pa {d[1]}'''


templateCTFCorrection = '''cd {d[0]}

mpiexec --tag-output -n {d[9]}  ctfCorrection.py \\
-u {d[2]} \\
-c {d[3]} \\
--metafile {d[4]} \\
--rotationAngle {d[5]} \\
--gridSpacing {d[6]} \\
--fieldSize {d[7]} \\
--binningFactor {d[8]}'''


templateFSC = '''cd {d[0]}

fsc.py {d[2]} \\
{d[3]} \\
{d[4]} \\
{d[5]} \\
--outputFolder {d[6]} \\
--fsc {d[7]} \\
--pixelsize {d[8]} \\
--randomizePhases {d[9]} \\
{d[10]}{d[11]} {d[12]}'''


templateFSC2 = '''cd {d[0]}

fsc.py {d[2]} \\
{d[3]} \\
{d[4]} \\
{d[5]} \\
--outputFolder {d[6]} \\
--fsc {d[7]} \\
--pixelsize {d[8]} \\
--randomizePhases {d[9]} \\
{d[10]}{d[11]} {d[12]}'''


templateCTFCorrectionImod = '''cd {d[0]}                                                                                                                                                                                      

ctfphaseflip -inp {d[1]} -o {d[2]} -an {d[3]} -defF {d[4]} \\
-defT {d[5]} -iW {d[6]} -pi {d[7]} -cs {d[8]} \\
-am {d[9]} -vo {d[10]} -AxisAngle {d[11]} {d[15]}

mrcs2mrc.py -f {d[2]} -t {d[12]} -p {d[13]} -o {d[14]}'''


templateMotionCorrection = '''cd {d[0]}

motioncor2 -In{d[1]} {d[2]}/ \\
-InSuffix .{d[8]} \\
-OutMrc {d[3]}/ \\
-Serial 1 \\
{d[7]} \\
{d[4]} {d[5]} {d[6]}
'''


templateAverageParticleList = '''cd {d[0]}

average.py -p {d[2]} -a {d[3]} -c {d[4]} {d[5]}'''


templateConvertData = '''cd {d[0]}

convert.py -t ./ {d[1]}{d[2]}{d[3]}{d[4]} -o {d[5]} {d[6]} --binPyTom {d[7]} --binWarpM {d[8]} --pixelSize {d[9]} \\
{d[10]}{d[11]}{d[12]}'''