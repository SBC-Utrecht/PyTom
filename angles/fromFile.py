from pytom.angles.globalSampling import GlobalSampling


class AngleListFromEM(GlobalSampling):
    """
    AngleListFromEM: 
    @deprecated: Use L{pytom.angles.globalSampling.GlobalSampling} instead! 
    """
    
    def __init__(self,filename='',inputIsInRadians=True,refinementParameters=None):
        """
        ___init__:
        @param filename: Filename of em file storing rotations. New files must have the naming (angles_INCREMENT_NUMBERROTATIONS.em) below. Set to either one of these:  
        angles_3_553680.em angles_07_45123.em angles_11_15192.em angles_12.85_7112.em angles_17.86_3040.em angles_18_3040.em angles_19.95_1944.em  angles_25.25_980.em angles_35.76_320.em angles_38.53_256.em angles_50_100.em angles_90_26.em angles_inf_0.em
        @param inputIsInRadians: Are the input values in radiants (True is assumed as default)
        @param refinementParameters: Parameters for next refinement step 
        @author: Thomas Hrabe  
        """
        
        super(AngleListFromEM,self).__init__(filename,inputIsInRadians,refinementParameters)