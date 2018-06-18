'''
Created on 27.06.2012

@author: thomas
'''
from pytom.basic.structures import PyTomClass



class PyTomServerMessage(PyTomClass):
    
    
    def __init__(self):
        """
        """


class FileMessage(PyTomServerMessage):
    
    def __init__(self,msgType,path,status):
        
        self._type      = str(msgType)
        self._path      = str(path)
        self._status    = str(status)
        
        
    def toXML(self):
        from lxml import etree
        messageXML = etree.Element("File",type=self._type,path=self._path,status = self._status)
        
        return messageXML


class ErrorMessage(PyTomServerMessage):
    
    def __init__(self,error):
        
        self._error = str(error)
        
        
    def toXML(self):
        from lxml import etree
        messageXML = etree.Element("Error",message=self._error)
        
        return messageXML  
        
class DataMessage(PyTomServerMessage):
    
    def __init__(self,dataType,data,sizeX,sizeY,sizeZ):
        self._type  = dataType
        self._data  = data
        self._sizeX = sizeX        
        self._sizeY = sizeY
        self._sizeZ = sizeZ
         
    def toXML(self):
        from lxml import etree
        messageXML = etree.Element("Data",type = self._dataType, data=self._data,sizeX=str(self._sizeX),sizeY=str(self._sizeY),sizeZ=str(self._sizeZ))
        
        return messageXML
        
         