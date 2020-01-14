'''
Created on Mar 21, 2011

@author: yuxiangchen
'''

class VolDsp():
    def __init__(self, vol):
        self.vol = vol
    
    def setUpper(self, value=None):
        if value:
            self.upper = value
        else:
            from pytom_volume import max
            self.upper = max(self.vol)
    
    def setLower(self, value=None):
        if value:
            self.lower = value
        else:
            from pytom_volume import min
            self.lower = min(self.vol)
    
    def getSlice(self, n, projectionAxis='z'):
        assert n>= 0
        
        from pytom_numpy import vol2npy
        nvol = vol2npy(self.vol)
        
        from pytom_volume import subvolume
        if projectionAxis == 'x':
            slice = nvol[n,:,:]
        elif projectionAxis == 'y':
            slice = nvol[:,n,:]
        elif projectionAxis == 'z':
            slice = nvol[:,:,n]
        else:
            raise Exception("Projection axis must be either 'x', 'y' or 'z'!")
        
        return slice
    
    def slice2img(self, slice, upper=None, lower=None):
        if upper is None:
            upper = slice.max()
        if lower is None:
            lower = slice.min()
        
        self.setUpper(upper)
        self.setLower(lower)
        
        x, y = slice.shape
        
        # initialize the image
        from PIL import Image
        img = Image.new("RGB", (x,y))
        
        for i in range(x):
            for j in range(y):
                img.putpixel((i,j), self.getRGB(slice[i][j]))
        
        return img
    
    def showSlice(self, n, projectionAxis='z'):
        slice = self.getSlice(n, projectionAxis)
        img = self.slice2img(slice)
        img.show()
    
    def setColorMap(self):
        pass
    
    def getRGB(self, value, cm=None):
        #TODO need to support more color map
        if value < self.lower:
            return (0,0,0)
        if value > self.upper:
            return (255,255,255)
        
        v = 255.*(value-self.lower)/(self.upper-self.lower)
        return (int(v), int(v), int(v))
    
    def drawCircle(self, img, pos, radius, color="#FF0000"):
        from PIL import ImageDraw
        draw = ImageDraw.Draw(img)
        xy = (pos[0]-radius, pos[1]-radius, pos[0]+radius, pos[1]+radius)
        draw.arc(xy, 0, 360, fill=color)
        
        return img