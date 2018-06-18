'''
Created on Jan 18, 2011

@author: hrabe
'''

def plotFSC(fsc,plotFileName):
    """
    plotFSC: Will plot resolution as it comes in fsc
    @param fsc:
    @param plotFileName 
    """
    
    from pytom.libraries.CairoPlot import dot_line_plot
    
    d = {}
    d[0] = fsc
    
    dot_line_plot(plotFileName,d,500,400,
                  background = None,
                  border = 20,
                  axis = True,
                  grid = True,
                  dots = True,
                  h_labels= 'Resolution',
                  v_labels = 'FSC',
                  h_bounds = None,
                  v_bounds = None)
    
    
def plotScoreDevelopment(particleLists):
    """
    plotScoreDevelopment: Will line plot the sums of scores for the particle lists.
    @param particleLists: particle lists 
    @type particleLists: list of L{pytom.basic.structures.ParticleList}
    """

    if not particleLists.__class__ == list :
        raise TypeError("You must provide a list of ParticleLists to plotScoreDevelopment!")

    scores = []
    
    for particleList in particleLists:
        scores.append(particleList.sumOfScores())
    print scores   
    
    try:
        from pytom.libraries import CairoPlot
        
        data = {}
        data[0] = scores
        
        CairoPlot.dot_line_plot('scorePlot.png',data,400,300,
                  background = None,
                  border = 0,
                  axis = True,
                  grid = True,
                  dots = True,
                  h_labels= 'ParticleList index',
                  v_labels = 'Score',
                  h_bounds = None,
                  v_bounds = None)
        
        
    except Exception:
    
        import matplotlib.pyplot as pyplot
    
        pyplot.figure()    
        pyplot.plot(range(len(scores)),scores,color='b')
        pyplot.show()
    
    
def plotClassSizes(particleList,plotFileName = None):
    """
    plotClassSizes
    @param particleList: 
    """
    
    from pytom.basic.structures import ParticleList
    
    if not particleList.__class__ == ParticleList:
        raise TypeError('Parameter must be a ParticleList!')
    
    if not plotFileName:
        plotFileName = 'ClassSize.png'
    
    classLists = particleList.splitByClass()
    
    sizes = []
    names = []
    for i in xrange(len(classLists)):
        sizes.append(len(classLists[i]))
        names.append(classLists[i][0].getClassName())
    
    try:
        from pytom.libraries import CairoPlot
        
        data = {}
        
        for i in xrange(len(sizes)):
            data[names[i]] = sizes[i]
        
        CairoPlot.donut_plot(plotFileName, data, 500, 500, background = None, gradient = False, shadow = False, colors = None, inner_radius = 0.3)
    except Exception:
        
        print 'CairoPlot failed!'
         
    
    
    
    