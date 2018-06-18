'''
Created on Apr 6, 2011

@author: hrabe
'''

from pytom.basic.structures import PyTomClass
from pytom.parallel.messages import Message


class ReconstructionMessage(Message):
    """
    ReconstructionMessage:
    """

    def __init__(self,sender='',receiver='',particleList='',projectionList='',cubeSize='', binning=1, applyWeighting = False):
        
        super(ReconstructionMessage,self).__init__(sender,receiver)
        
        self._particleList   = particleList
        self._projectionList = projectionList
        self._cubeSize       = cubeSize
        self._binning        = binning
        self._applyWeighting = applyWeighting
                
    def toXML(self):
        from lxml import etree
        
        message_element = etree.Element('ReconstructionMessage',Sender = str(self._sender), Recipient = str(self._recipient), Timestamp = str(self._timestamp))
                                
        message_element.append(self._particleList.toXML())
        message_element.append(self._projectionList.toXML())
        message_element.set('CubeSize',str(self._cubeSize))
        message_element.set('Binning',str(self._binning))
        message_element.set('ApplyWeighting',str(self._applyWeighting))
        
        return message_element
                

    def fromXML(self,xmlObj):
        from lxml.etree import _Element
        from pytom.basic.structures import ParticleList
        from pytom.reconstruction.reconstructionStructures import ProjectionList
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
                
        if xmlObj.tag == 'ReconstructionMessage':
            message_element = xmlObj
        else:
            Exception('Is not a ReconstructionMessage! You must provide a valid ReconstructionMessage object.')
                    
        
        self._sender = message_element.get('Sender')
        self._recipient = message_element.get('Recipient')
        self._timestamp = message_element.get('Timestamp')
        self._cubeSize = int(message_element.get('CubeSize'))
        self._binning = float(message_element.get('Binning'))
        self._applyWeighting = bool(message_element.get('ApplyWeighting'))        
        
        self._particleList = ParticleList('.')
        particleListXML = message_element.xpath('ParticleList')[0]
        self._particleList.fromXML(particleListXML)
               
        self._projectionList = ProjectionList()
        projectionListXML = message_element.xpath('ProjectionList')[0]
        self._projectionList.fromXML(projectionListXML)
        
    
    def getProjectionList(self):
        
        return self._projectionList
    
    def getParticleList(self):
        
        return self._particleList
        
    def getCubeSize(self):
        
        return self._cubeSize
    
    def getBinning(self):
        
        return self._binning
    
    def getApplyWeighting(self):
        
        return self._applyWeighting  
    
        
                       


class ReconstructionWorker(PyTomClass):
    """
    ReconstructionWorker
    """

    def setJob(self,reconstructionJobMessage):
        
        self._projectionList = reconstructionJobMessage.getProjectionList()
        self._particleList   = reconstructionJobMessage.getParticleList()
        self._cubeSize       = reconstructionJobMessage.getCubeSize()
        self._binning        = reconstructionJobMessage.getBinning()
        self._applyWeighting = reconstructionJobMessage.getApplyWeighting() 
        
        
    def run(self,verbose=False):
        
        self._projectionList.generateVolumes(self._particleList,self._cubeSize,self._binning, self._applyWeighting,False,verbose)
        
    def parallelRun(self,verbose=False):
    
        import pytom_mpi
        from pytom.parallel.messages import StatusMessage,MessageError
        from pytom.basic.exceptions import ParameterError
        from pytom.basic.structures import PyTomClassError
        
        end = False
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
        
        mpi_id = pytom_mpi.rank()
        
        while not end:            
                       
            #listen for messages
            
            mpi_msgString = pytom_mpi.receive()
            
            if verbose:
                print mpi_msgString
                
            try:
                
                
                #wait for job and start processing
                msg = ReconstructionMessage()
                msg.fromStr(mpi_msgString)
                
                
                
                self.setJob(msg)
                
                self.run()
                
                
                resultMsg = StatusMessage(mpi_id,'0')
                resultMsg.setStatus('Finished')
                
                pytom_mpi.send(str(resultMsg),0)
               
                               
                
            except (MessageError,PyTomClassError,ParameterError):
                    try:
                        #as StatusMessage and finish
                        msg = StatusMessage('','')
                        msg.fromStr(mpi_msgString)
                        if msg.getStatus() == 'End':
                            end = True
                                               
                        
                    except (MessageError,PyTomClassError,ParameterError):
                        #print mpi_msgString
                        #raise MessageError('Message unknown!')
                        print 'Error parsing message. Message either unknown or invalid.'
                        assert False
            except:
                    print 'wild except'
                    
        
        


def parallelReconstruction(particleList, projectionList, cubeSize, binning, applyWeighting,verbose=False):
    """
    parallelReconstruction
    """
    import pytom_mpi
    from pytom.parallel.messages import StatusMessage
        
    if not pytom_mpi.isInitialised():
        pytom_mpi.init()
        
    mpi_id = pytom_mpi.rank()
    
    if mpi_id == 0:
                
        firstDistribute = False
        numberWorkers = pytom_mpi.size() -1
                        
        #split particleList by number nodes
        
        splitSize = len(particleList) / numberWorkers
                
        pl = []
        
        for i in xrange(0,len(particleList),splitSize):

            pl.append(particleList[i:i+splitSize])
               
        
        for i in xrange(0,numberWorkers):
                                            
            msg = ReconstructionMessage(0,i+1,pl[i],projectionList,cubeSize, binning,applyWeighting)
                                    
            pytom_mpi.send(str(msg),i+1)
            
               
            
           
        finished = False
        msgCounter = 0
        
        while not finished:
            
            mpi_msgString = pytom_mpi.receive()
            msg = StatusMessage(1,'0')
            msg.fromStr(mpi_msgString)
            
            if not msg.getStatus() == 'Finished':
                print 'Worker ' + str(msg.getSender()) + ' sent status: ' + str(msg.getStatus())
                 
            msgCounter += 1
            
            finished = msgCounter == numberWorkers
                        
            
        for i in xrange(0,numberWorkers):
            msg = StatusMessage(mpi_id,'0')
            msg.setStatus('End')
            pytom_mpi.send(str(msg),i+1)
            
        
    else:
        
        worker = ReconstructionWorker()
        worker.parallelRun(verbose)
        
        
    pytom_mpi.finalise()
    
        

            
            
    
    
    

    
