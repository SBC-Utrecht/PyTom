'''
Created on Jul 12, 2010

@author: hrabe
'''

def checkerBoard(sizeX,sizeY,sizeZ):
    """
    checkerBoard: returns volume where voxel values alternate between 0 and 1
    """
    from pytom_volume import vol
    
    board = vol(sizeX,sizeY,sizeZ)
    
    counterX = 0
    counterY = 0
    counterZ = 0
    for x in xrange(sizeX):
        for y in xrange(sizeY):
            for z in xrange(sizeZ):
                board.setV(counterX,x,y,z)
                counterX = (counterX - 1) %2 
            
            counterX = counterX * counterY
            counterY = (counterY -1) %2
            
        counterX = counterX * counterZ
        counterZ = (counterZ -1) %2
            
    return board
    
    

