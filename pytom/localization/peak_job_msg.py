'''
Created on May 19, 2010

@author: chen
'''

from pytom.parallel.messages import Message

class PeakJobMsg(Message):
    """
    PeakJobMsg : Derived from message, this class constructs the job message that can be used in MPI
    """
    def getJob(self):
        return self.job
    
    def setJob(self,job):
        self.job = job
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: chen
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-PeakJobMsg object.')
    
        from pytom.localization.peak_job import PeakJob
        
        main = xmlObj.xpath('/PeakJobMsg')
        
        if len(main) == 0:
            from pytom.parallel.messages import MessageError
            raise MessageError("This message is not a PeakJobMsg")
        main = main[0]
         
        message = main.xpath('Message')
        message = message[0]
#        jobDescription = main.xpath('JobDescription')  
#        jobDescription = jobDescription[0]
        
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')
        
        jobDescription = main.xpath('/*/JobDescription')
        j = PeakJob()
        j.fromXML(jobDescription[0])
        self.setJob(j)
        
    def toXML(self):        
        """
        toXML : Compiles a XML file from job object
        @author: chen
        """    
        from lxml import etree
        
        jobMsg = etree.Element("PeakJobMsg")
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient, Timestamp = self._timestamp.__str__())
        
        jobMsg.append(messageElement)
        
        jobElement = self.job.toXML()
        jobMsg.append(jobElement)
        
        return jobMsg
    
class PeakResultMsg(Message):

    result = -1

    def setResult(self, result):
        """
        @param result: result
        @type result: L{pytom.localization.peak_job.PeakResult}
        """
        self.result = result
        
    def getResult(self):
        return self.result
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: chen 
        """
        from lxml.etree import _Element;
        
        if xmlObj.__class__ != _Element :
            
            from pytom.basic.exceptions import ParameterError;
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-PeakResultMsg object.');
  
        from pytom.localization.peak_job import PeakResult;

        if xmlObj.tag == "PeakResultMsg":
            main = xmlObj
        else:
            main = xmlObj.xpath('PeakResultMsg');
            if len(main) == 0:
                from pytom.parallel.messages import MessageError
                raise MessageError("This message is not a PeakResultMsg")
            main = main[0]
        
        message = main.xpath('Message');
        message = message[0];
        
        self._sender = message.get('Sender');
        self._recipient = message.get('Recipient');
        self._timestamp = message.get('Timestamp');        
        
        res = main.xpath('Result')[0]
        self.result = PeakResult();
        self.result.fromXML(res);
        
    def toXML(self):       
        from lxml import etree;       
        
        resMsg = etree.Element("PeakResultMsg");
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient, Timestamp = self._timestamp.__str__());
        
        resMsg.append(messageElement);
       
        resultElement = self.result.toXML();
        
        resMsg.append(resultElement);
        
        return resMsg;