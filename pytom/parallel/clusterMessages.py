'''
Created on Mar 24, 2010

@author: hrabe
'''
from pytom.parallel.messages import Message;

class CorrelationVectorJobMessage(Message):
    """
    CorrelationVectorJobMessage
    @todo: XML UNITTEST , hide attributes
    """
    def setJob(self,job):
        self.job = job
        
    def getJob(self):
        return self.job
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element :
            
            from pytom.basic.exceptions import ParameterError;
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-CorrelationVectorJobMessage object.')
    
        from pytom.cluster.correlationMatrixStructures import CorrelationVectorJob
            
        main = xmlObj.xpath('/CorrelationVectorJobMessage')
    
        message = main[0].xpath('Message')
        message = message[0]
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')
        
        
        jobDescription = main[0].xpath('CorrelationVectorJob')
        
        j = CorrelationVectorJob()
        j.fromXML(jobDescription[0])
        self.setJob(j)
    
    def toXML(self):
        """
        toXML : Compiles a XML file from object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        jobMsg = etree.Element("CorrelationVectorJobMessage")
        
        messageElement = etree.Element("Message",Sender = self._sender.__str__(), Recipient = self._recipient.__str__(), Timestamp = self._timestamp.__str__())
        
        jobMsg.append(messageElement)
        
        jobElement = self.job.toXML()
        jobMsg.append(jobElement)
        
        return jobMsg


class CorrelationVectorMessage(Message):
    """
    CorrelationVectorMessage
    @todo: XML UNITTEST , hide attributes
    """

    def setVector(self,vector):
        self.vector = vector
        
    def getVector(self):
        return self.vector
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element :
            
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-CorrelationVectorMessage object.')
    
        from pytom.cluster.correlationMatrixStructures import CorrelationVector
            
        main = xmlObj.xpath('/CorrelationVectorMessage')
    
        message = main[0].xpath('Message')
        message = message[0]
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')
        
        vectorXML = main[0].xpath('CorrelationVector')
        
        v = CorrelationVector()
        v.fromXML(vectorXML[0])
        self.setVector(v)
    
    def toXML(self):
        """
        toXML : Compiles a XML file from object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        vectorMsg = etree.Element("CorrelationVectorMessage")
        
        messageElement = etree.Element("Message",Sender = self._sender.__str__(), Recipient = self._recipient.__str__(), Timestamp = self._timestamp.__str__())
        
        vectorMsg.append(messageElement)
        vectorMsg.append(self.vector.toXML())
        
        return vectorMsg