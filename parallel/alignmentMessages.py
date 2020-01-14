from pytom.parallel.messages import Message


class GrowingAverageJobMessage(Message):

    def setJob(self,job):
        self.job = job;
        
    def getJob(self):
        return self.job;
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring;
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-MaximisationJobMsg object.')
    
        from pytom.alignment.structures import GrowingAverageJob
            
        main = xmlObj.xpath('/GrowingAverageJobMsg')
    
        message = main[0].xpath('Message')
        message = message[0]
        
        self._sender = int(message.get('Sender'))
        self._recipient = int(message.get('Recipient'))
        self._timestamp = message.get('Timestamp')
        
        jobDescription = main[0].xpath('GrowingAverageJob')
        
        j = GrowingAverageJob()
        j.fromXML(jobDescription[0])
        self.setJob(j)
    
    def toXML(self):
        """
        toXML : Compiles a XML file from object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        jobMsg = etree.Element("GrowingAverageJobMsg")
        
        messageElement = etree.Element("Message",Sender = self._sender.__str__(), Recipient = self._recipient.__str__(), Timestamp = self._timestamp.__str__())
        
        jobMsg.append(messageElement)
        
        jobElement = self.job.toXML()
        jobMsg.append(jobElement)
        
        return jobMsg

class GrowingAverageResultMessage(Message):
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-MaximisationJobMsg object.')
        
        main = xmlObj.xpath('/GrowingAverageResultMsg')
    
        message = main[0].xpath('Message')
        message = message[0]
    
        self._sender = int(message.get('Sender'))
        self._recipient = int(message.get('Recipient'))
        self._timestamp = message.get('Timestamp')

    
    def toXML(self):
        """
        toXML : Compiles a XML file from object
        @author: Thomas Hrabe
        """    
        from lxml import etree;
        
        resultMsg = etree.Element("GrowingAverageResultMsg");
        
        messageElement = etree.Element("Message",Sender = self._sender.__str__(), Recipient = self._recipient.__str__(), Timestamp = self._timestamp.__str__());
        
        resultMsg.append(messageElement);
        
        return resultMsg;
       
class MaximisationJobMsg(Message):
    """
    MaximisationJob : Derived from message, this class stores all infos needed for a maximisation job
    """
    def getJob(self):
        return self.job
    
    def setJob(self,job):
        self.job = job
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring;
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-MaximisationJobMsg object.');
    
        from pytom.alignment.structures import MaximisationJob;
        
        main = xmlObj.xpath('/MaximisationJobMsg');
        
        if len(main) == 0:
            raise TypeError("This message is not an MaximisationJobMsg");
        
        main = main[0];
         
        message = main.xpath('Message');
        message = message[0];
        jobDescription = main.xpath('JobDescription');    
        jobDescription = jobDescription[0];
        
        self._sender = message.get('Sender');
        self._recipient = message.get('Recipient');
        self._timestamp = message.get('Timestamp');
        
        jobDescription = main.xpath('/*/JobDescription');
        j = MaximisationJob();
        j.fromXML(jobDescription[0]);
        self.setJob(j);
        
    def toXML(self):        
        """
        toXML : Compiles a XML file from job object
        @author: Thomas Hrabe
        """    
        from lxml import etree;
        
        jobMsg = etree.Element("MaximisationJobMsg");
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient, Timestamp = self._timestamp.__str__());
        
        jobMsg.append(messageElement);
        
        jobElement = self.job.toXML();
        jobMsg.append(jobElement);
        
        return jobMsg;
    
    def __eq__(self,otherJob):
        
        return self.__str__() == otherJob.__str__();
    
class MaximisationResultMsg(Message):

    result = -1;

    def setResult(self,result):
        self.result = result;
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element;
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-MaximisationResultMsg object.');
  
        from pytom.alignment.structures import MaximisationResult;
  
        message = xmlObj.xpath('Message');
        message = message[0];
        
        self._sender = message.get('Sender');
        self._recipient = message.get('Recipient');
        self._timestamp = message.get('Timestamp');        
        
        self.result = MaximisationResult('','','','','');
        self.result.fromXML(xmlObj);
        
    def toXML(self):
        """
        toXML : Compiles a XML from result object
        @author: Thomas Hrabe
        """        
        from lxml import etree;
        
        resMsg = etree.Element("MaximisationResultMsg");
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient, Timestamp = self._timestamp.__str__());
        
        resMsg.append(messageElement);
       
        resultElement = self.result.toXML();
        
        resMsg.append(resultElement);
        
        return resMsg;

class ExpectationJobMsg(Message):
    
    def getJob(self):
        return self.job;
    
    def setJob(self,job):
        from pytom.alignment.structures import ExpectationJob
        
        if not job.__class__ == ExpectationJob:
            raise TypeError('Object is of wrong class.')
        
        self.job = job
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML ExpectationJobMsg object.')
        
        from pytom.alignment.structures import ExpectationJob
         
        message = xmlObj.xpath('Message')
        message = message[0]
        
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')

        exJobXML = xmlObj.xpath('/ExpectationJobMsg/ExpectationJob')
        
        self.job = ExpectationJob('')
        
        if len(exJobXML) > 0:
            self.job.fromXML(exJobXML[0])
        else:
            raise TypeError('This is not a valid ExpectationJobMsg. ExpectationJob is missing')
        
    def toXML(self):
        """
        toXML : Compiles a XML from job object
        @author: Thomas Hrabe
        """
        from lxml import etree
        
        jobMsg = etree.Element("ExpectationJobMsg")
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient,
                                       Timestamp = str(self._timestamp))
        jobMsg.append(messageElement)
        jobElement = self.job.toXML()
        jobMsg.append(jobElement)
        
        return jobMsg
    
class ExpectationResultMsg(Message):
    
    def setResult(self,result):
        from pytom.alignment.structures import ExpectationResult;
        if not result.__class__ == ExpectationResult:
            raise TypeError('Object is of wrong class.');
        self.result = result;
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-ExpectationResultMsg object.')
        
        from pytom.alignment.structures import ExpectationResult
        
        main = xmlObj.xpath('/ExpectationResultMsg')
        
        if len(main) == 0:
            raise TypeError("This XML is not an ExpectationResultMsg")
        
        main = main[0]
        
         
        message = main.xpath('Message')
        message = message[0]
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')        
        
        self.result = ExpectationResult('')
        self.result.fromXML(xmlObj.xpath('/*/ExpectationResult')[0])
        
        
    def toXML(self):       
        """
        toXML : Compiles a XML from result object
        @author: Thomas Hrabe
        """
        from lxml import etree
        
        resMsg = etree.Element("ExpectationResultMsg")
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient, Timestamp = self._timestamp.__str__())
        
        resMsg.append(messageElement)
       
        resultElement = self.result.toXML()
        
        resMsg.append(resultElement)
        
        return resMsg
         