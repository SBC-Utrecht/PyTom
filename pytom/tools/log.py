class Log:
    
    
    def objectToLogfile(object,logfileName):
        """
        objectToLogfile: Dumps object to logfile
        @param object: Whats to be dumped
        @param logfileName: Destination file name
        @author: Thomas Hrabe  
        """
    
        objString = object.__str__()
    
        f = open(filename,'a')
        
        f.write(objString)
        
        f.close()
        
    
    