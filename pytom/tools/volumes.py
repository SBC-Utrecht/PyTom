'''
Created on Jul 12, 2010

@author: hrabe
'''

def checkerBoard(size_x,size_y,size_z):
    """
    checkerBoard: returns volume where voxel values alternate between 0 and 1
    """
    from pytom_volume import vol
    
    board = vol(size_x,size_y,size_z)
    
    counterX = 0
    counterY = 0
    counterZ = 0
    for x in range(size_x):
        for y in range(size_y):
            for z in range(size_z):
                board.setV(counterX,x,y,z)
                counterX = (counterX - 1) %2 
            
            counterX = counterX * counterY
            counterY = (counterY -1) %2
            
        counterX = counterX * counterZ
        counterZ = (counterZ -1) %2
            
    return board
    
    

