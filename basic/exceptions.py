class ParallelError(Exception):    
    """
    PrallelError : Thrown during parallel processing 
    """
    def __init__(self,value):
        self.value = value;
    def __str__(self):
        print self.value;


class ParameterError(Exception):
    """
    ParameterError : Thrown if any input parameters of a arbitrary function are of wrong type, NAN,...
    """    
    def __init__(self,value):
        self.value = value;
    def __str__(self):
        print self.value;
        
        
class FourierSingletonError(Exception):
    """
    FourierSingletonError :  
    """    
    def __init__(self,value):
        self.value = value;
    def __str__(self):
        print self.value;
    