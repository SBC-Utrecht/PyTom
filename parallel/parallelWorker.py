'''
Created on May 13, 2011

@author: luiskuhn
'''


from pytom.basic.structures import PyTomClass
from pytom.parallel.messages import Message

class ParallelWorker(PyTomClass):
    
    def __init__(self):
        
        import pytom_mpi
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
        
        self._mpi_id = pytom_mpi.rank()
        self._numberWorkers = pytom_mpi.size() -1
        
        self._jobList = []
                
    
    def fillJobList(self,jobs):
        """
        fillJobList: Fills current joblist with jobs
        @param jobs: 
        """
        
        if len(jobs) == 0:
            raise RuntimeError('Joblist provided is empty!')
        self._jobList = jobs
    
    def run(self):
        """
        run: Overwrite this function to define specific processing procedures of your jobs 
        """
        print 'You forgot to overwrite run in your worker class'
        assert False
    
    def getMsgObject(self, msgString):
        """
        getMsgObject: Parse message object from string. \
        Must overwrite this function and specify which message to listen to.
        """
        print 'You forgot to overwrite getMsgObject in your worker class'
        assert False
    
    def setJob(self, jobMessage):
        """
        setJob:
        """
        print 'You forgot to overwrite setJob in your worker class'
        assert False
    
      
    def parallelWork(self, verbose=False,doFinalize = True):
        """
        parallelWork: Distribute joblist to workers. Leave as it is.
        @param verbose:
        @param doFinalize:   
        """
        
        import pytom_mpi        
        from pytom.parallel.messages import Message,StatusMessage, MessageError
        from pytom.basic.exceptions import ParameterError
        from pytom.basic.structures import PyTomClassError

        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
        
        if self._mpi_id == 0:
            #if current node == 0, be the master node
            numberJobs = len(self._jobList)
            
            if self._numberWorkers <= numberJobs: 
                numberJobsToSend = self._numberWorkers
            else: 
                numberJobsToSend = numberJobs
            
            for i in xrange(0, numberJobsToSend):
                #send out all first numberJobsToSend jobs
                pytom_mpi.send(str(self._jobList[i]),i+1)
                
                
            numberFinishedJobs = 0
            numberSentJobs = numberJobsToSend
                      
            finished = numberSentJobs == numberJobs and numberFinishedJobs == numberJobs
            
            while not finished:
                #distribute remaining jobs to finished workers
                mpi_msgString = pytom_mpi.receive()
                
                msg = Message('1','0')
                msg.fromStr(mpi_msgString)
                         
                numberFinishedJobs += 1
            
                if numberSentJobs < numberJobs:
                    pytom_mpi.send(str(self._jobList[numberSentJobs]),int(msg.getSender()))
                    numberSentJobs += 1
            
                finished = numberSentJobs == numberJobs and numberFinishedJobs == numberJobs
                            
            if doFinalize:    
                for i in xrange(0,self._numberWorkers):
                    msg = StatusMessage('0',i+1)
                    msg.setStatus('End')
                    pytom_mpi.send(str(msg),i+1)
                    print 'Sending end msg to:', i+1
                
        else:
            #if any other node id, be a worker node           
            end = False              
            
            while not end:            
                           
                #listen for messages
                
                mpi_msgString = pytom_mpi.receive()
                
                if verbose:
                    print mpi_msgString
                    
                try:
                    
                    
                    #wait for job and start processing
                    msg = self.getMsgObject(mpi_msgString)                           
                    
                    self.setJob(msg)
                    
                    self.run()
                    
                    resultMsg = StatusMessage(self._mpi_id,'0')
                    resultMsg.setStatus('Finished')
                    
                    pytom_mpi.send(str(resultMsg),0)
                   
                except (MessageError,PyTomClassError,ParameterError):
                        try:
                            #message is a StatusMessage
                            #if message status is End, finish this worker. 
                            #You can also add other statuses
                            msg = StatusMessage('','')
                            msg.fromStr(mpi_msgString)
                            if msg.getStatus() == 'End':
                                end = True
                                                   
                            
                        except (MessageError,PyTomClassError,ParameterError):
                            #print mpi_msgString
                            #raise MessageError('Message unknown!')
                            raise RuntimeError('Error parsing message. Message either unknown or invalid.')
                            
                except:
                        raise RuntimeError('Something went terribly wrong. Aborting.')

        if doFinalize:
            pytom_mpi.finalise()
            
 