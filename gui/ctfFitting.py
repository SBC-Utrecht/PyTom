from numpy import * 
import os

def CTF(image,U,V,defocus_angle,voltage,Cs,Ac,Mag,pixel_size,phase_shift=0):
    defocus_angle = pi * defocus_angle/180
    wavelength = wavelength_eV2nm(voltage)
    
    X,Y = meshgrid(arange(0,image.shape[0],1.),arange(0,image.shape[1],1.))
    X-=image.shape[0]/2
    Y-=image.shape[1]/2
    
    
    X *= 1/(pixel_size * 1024)
    Y *= 1/(pixel_size * 1024)
    R = sqrt(X**2+Y**2)
    theta = arcsin(Y/R)
    
    z = U*cos(theta-defocus_angle)**2+V*sin(theta-defocus_angle)**2
    z[:,-int(image.shape[1]/2):]= flipud(z[:,-int(image.shape[1]/2):])
    phases = -0.5*pi*Cs*(wavelength**3)*(R**4)+pi*wavelength*z*(R**2)
    
    CTF = -sin(phase_shift+phases)
    return CTF


def wavelength_eV2nm(ev):
    return 6.625 * 10**(-34) / ((2*9.1*10**(-31))* 1.6 *(10**(-19)) * ev *1)**0.5 

if __name__=='__main__':
    from pylab import *
    ctf = CTF(ones((2000,2000)),6908E-9,7115E-9,30.68,200E3,2.7E-3,0.1,80000,3.5E-10)
    imshow(ctf)
    show()
