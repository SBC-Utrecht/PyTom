'''
Created on Nov 27, 2010

@author: chen
'''
import time, sys

class Timing():
    """
    Timing: A really simple class that can be used to determine the runtime of a process.
    """
    def __init__(self):
        if sys.platform == "win32":
            self.timer = time.clock
        else:
            self.timer = time.time
    
    def start(self):
        """
        start: Take first timestamp
        """
        self.t0 = self.timer()
    
    def end(self):
        """
        end: Finish by taking second timestamp
        @return: Time difference between first and second timestamp
        """
        self.t1 = self.timer()
        return self.t1-self.t0