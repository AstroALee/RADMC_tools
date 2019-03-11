'''
    Set of plotting routines that can be used on RADMC-3D data


    Author: Aaron T. Lee, Summer 2017


    Free to share, use, blah blah.
    I take no responsibility if any errors that exist below screw up your analysis.
    Don't treat code like a black box.

'''

import matplotlib.pyplot as plt
import numpy as np

def PlotSpec(head,wave,spec_array,args):

    logSpec = [ [] for i in range(len(spec_array)) ]
    for i in range(len(spec_array)):
        logSpec[i] = np.log10(spec_array[i])

    nuarr = (3.0e10)  / 1.0e9 / wave

    for i in range(len(spec_array)):
        #if(args.pixel[i]>=0): plt.plot(nuarr,logSpec[i],label=str(args.pixel[i]))
        if(args.pixel[i]>=0): plt.plot(nuarr,spec_array[i],label=str(args.pixel[i]))

    #plt.plot([23.6944955,23.6944955],[0,10],'--',alpha=0.5)
    #plt.plot([23.7226336,23.7226336],[0,10],'--',alpha=0.5)

    plt.legend()
    #plt.axis((min(nuarr), max(nuarr), -17,-14))
    plt.xlabel(r'Frequency   (GHz)')
    #plt.ylabel('Log10 Intensity')
    plt.ylabel('Specific Intensity   (Kelvin)')

    plt.savefig('spectrum.' + args.outfile)
    print("Spectrum file created.")



def PlotIntMap(head,image,args):
    print("Creating integrated intensity map")
    print("For moment " + str(args.moment))

    Cm2Pc = 3.24078e-19

    zoomfac = 2.0

    print(head["Ny"])
    print(head["Nx"])

    dx = head["Pixel_Size_X"]*Cm2Pc
    x_hi = 0.5*dx*head["Nx"]/zoomfac
    x_lo = -x_hi
    dy = head["Pixel_Size_Y"]*Cm2Pc
    y_hi = 0.5*dy*head["Ny"]/zoomfac
    y_lo = -y_hi

    X = np.linspace(x_lo, x_hi, head["Nx"]/zoomfac)
    Y = np.linspace(y_lo, y_hi, head["Ny"]/zoomfac)
    image = image.reshape(head["Nx"],head["Ny"]) # original


    if(args.normalize==0):
        unit_fac = pow(1.0e-5,1+args.moment)
    else:
        unit_fac = pow(1.0e-5,0+args.moment)

    plt.pcolormesh(X, Y, unit_fac*image, cmap = 'viridis')
    cbar = plt.colorbar()
    plt.axis((x_lo, x_hi, y_lo, y_hi))
    plt.xlabel('X   (pc)')
    plt.ylabel('Y   (pc)')
    if(args.normalize==0):
        cbar.set_label(r'[K (km/s)^' + str(1+args.moment) +']',rotation=90)
        plt.title(r'Integrated Intensity $\int I_v v^' + str(args.moment) + ' dv$')
    else:
        cbar.set_label(r'[(km/s)^' + str(args.moment) +']',rotation=90)
        plt.title(r'Normalized Integrated Intensity $\int I_v v^' + str(args.moment) + ' dv$')

    plt.savefig('intmap.' + str(args.moment) + '.' + args.outfile)
    print("Image file created.")
