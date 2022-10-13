class MessageError(Exception):    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        print(self.value)
        

from pytom.basic.structures import PyTomClass


class Message(PyTomClass):
        
    
    def __init__(self,sender='',recipient=''):
        import time
        t = time.localtime()
        self._timestamp = str(t[3])+':'+str(t[4])+':'+str(t[5])+'-'+str(t[2])+'.'+str(t[1])+'.'+str(t[0])
        self._sender = str(sender)
        self._recipient = str(recipient)

    def getSender(self):
        return self._sender
    
    def getRecipient(self):
        return self._recipient
    
    def getTimestamp(self):
        return self._timestamp
    
    def fromXML(self,xmlObj):
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element :
            
            raise RuntimeError('Is not a lxml.etree._Element! You must provide a valid XML-StatusMessage object.')
        
        message = xmlObj.xpath('Message')
        message = message[0]
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')
        
        
    
class StatusMessage(Message):
    """
    StatusMessage : Used for sending simple messages around. 
    """
        
    def setStatus(self,status):
        self._status = status
      
    def getStatus(self):
        return self._status
      
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element :
            
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-StatusMessage object.')
        
        message = xmlObj.xpath('Message')
        message = message[0]
        self._sender = message.get('Sender')
        self._recipient = message.get('Recipient')
        self._timestamp = message.get('Timestamp')
        
        status = xmlObj.xpath('Status')
        status = status[0]
        self._status = status.get('Value')
        
    def toXML(self): 
        """
        toXML : Compiles a XML from result object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        jobMsg = etree.Element("StatusMessage")
        
        messageElement = etree.Element("Message",Sender = self._sender, Recipient = self._recipient, Timestamp = str(self._timestamp))
        
        jobMsg.append(messageElement)
        
        statusElement = etree.Element("Status",Value = self._status)
        
        jobMsg.append(statusElement)
        
        return jobMsg
    