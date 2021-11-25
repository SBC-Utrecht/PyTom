from pytom.angles.localSampling import LocalSampling
 
 
class AV3Sampling(LocalSampling):
    """
    AV3Sampling: Deprecated and for backward compatibility
    @deprecated: 
    """ 
     
    def __init__(self,shells=3,increment=3,z1Start=0.0,z2Start=0.0,xStart=0.0):
        
        super(self.__class__,self).__init__(shells,increment,z1Start,z2Start,xStart)
        

        
                