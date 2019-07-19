from pytom.basic.structures import ParticleList

alignmentList = 'AlignmentList-1.xml'
particleList = 'alignedParticles.xml'
pl = ParticleList()
pl.fromAlignmentList(alignmentList)
pl.toXMLFile(particleList)
