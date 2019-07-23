import pandas as pd
import numpy as np
import matplotlib.pyplot as plt





''' Global set of constants (cgs)'''
class Csts:
    c = 3.0e10 # speed of light
    kB = 1.3806e-16 # Boltz constant

'''
    Convert specific intensity to Kelvin using Rayleigh-Jeans tail
    inputs must be in cgs
'''
def intToTemp(intensity,freq):
    output = intensity*pow(Csts.c/freq,2)/(2.0*Csts.kB)
    return(output)

def freqToVel(freq,linecenter):
    output = Csts.c*(freq - linecenter)/linecenter
    return(output)

''' Grab the header information from the file '''
def grabRADheader(fname):

    # Variable that will store the header info
    header_data = {}
    # Variable the will store the wavelengths
    wave_data = []
    # Variables for size of header
    header_length = 0
    header_size = 0

    # Open file
    f = open(fname,'r')

    line = f.readline()
    header_data["iformat"] = np.int64(line.strip())
    header_length = header_length + 1
    header_size = header_size + len(line)

    line = f.readline()
    header_data["Nx"],header_data["Ny"] = [np.int64(x) for x in line.strip().split()]
    header_length = header_length + 1
    header_size = header_size + len(line)

    line = f.readline()
    header_data["Nwavelengths"] = np.int64(line.strip())
    header_length = header_length + 1
    header_size = header_size + len(line)

    line = f.readline()
    header_data["Pixel_Size_X"],header_data["Pixel_Size_Y"] = [np.float64(x) for x in line.strip().split()]
    header_length = header_length + 1
    header_size = header_size + len(line)

    # Let's make sure we are in 'observer at infinity' mode
    assert(header_data["iformat"] == 1)

    # The next set of lines will be the camera wavelengths.
    # This will convert to cgs
    for i in range(header_data["Nwavelengths"]):
        line = f.readline()
        wave_data.append( np.float64(line.strip()) )
        header_length = header_length + 1
        header_size = header_size + len(line)
    # converts to cgs, assumed input file was in microns
    wave_data = 1.0e-4*np.array(wave_data)

    # BOOOOO for all those blank lines
    while(True):
        line = f.readline()
        if(len(line.strip()) == 0):
            header_length = header_length + 1
            header_size = header_size + len(line)
        else:
            break

    f.close()
    header_data["header_length"] = header_length
    header_data["header_size"] = header_size

    return header_data, wave_data


''' Load data into a pandas array'''
def LoadData(fname):

    head,wave = grabRADheader(fname)

    # Read in the intensities
    dim = int(head["Nx"]*head["Ny"])
    nwave = int(head["Nwavelengths"])

    # to load into a pandas dataframe, we need two things:
    # list of wavelengths (wave above, done!)
    # list of lists: inner list all intensities of given wavelength, outer list each wavelength
    # Assumes dim lines of data is intensity at each wavelength of the same pixel

    f = open(fname,'r',1048576) # byte number should speed things up a little?
    f.seek(head["header_size"]) # skip past the header

    intData = [ [] for x in range(nwave) ] # empty list of lists
    idx = 0
    for line in f:
        if(idx==dim*nwave): break # should also be end of line
        if(line.isspace()): continue # if for some reason it's blank, skip
        wnum = int(idx*1.0/dim)
        intData[wnum].append( np.float64(line.strip()) ) # assumes one number / line
        idx = idx + 1
    f.close()
    # All information is now read in
    # each list, e.g., intData[5], is a list of length dim and there are nwave of them

    # Convert to a Pandas array
    df = pd.DataFrame(intData,index=wave)

    return head,df

'''
    Function that generates the \Delta\nu and \Delta\lambda array from your
    wavelength array ( list(df.index) will give list of wavelengths )
'''
def Create_DeltaArrays(wavelengths):
    nwave = len(wavelengths)
    dnu = np.zeros(nwave)
    dwa = np.zeros(nwave)

    nu = Csts.c/wavelengths
    for i in range(nwave):
        if(i==0):
            dnu[0] = 2.0*0.5*np.fabs(nu[0]-nu[1])
            dwa[0] = 2.0*0.5*np.fabs(wavelengths[0]-wavelengths[1])
        elif(i==nwave-1):
            dnu[-1] = 2.0*0.5*np.fabs(nu[-1]-nu[-2])
            dwa[-1] = 2.0*0.5*np.fabs(wavelengths[-1]-wavelengths[-2])
        else:
            dnu[i] = 0.5*np.fabs(nu[i]-nu[i-1]) + 0.5*np.fabs(nu[i+1]-nu[i])
            dwa[i] = 0.5*np.fabs(wavelengths[i]-wavelengths[i-1]) + 0.5*np.fabs(wavelengths[i+1]-wavelengths[i])

    return(dnu,dwa)

def Find_MomentAtPixel(df,moment,pixel):
    # get delta nu and lambda arrays for integrals
    # wavelengths
    waves = np.array(df.index)
    nus   = Csts.c/waves
    dnu,dwa = Create_DeltaArrays(waves)
    loc_intensity = np.array(df[pixel]) # gets the specific intensities (per nu)

    val = (loc_intensity * pow(nus,moment) * dnu).sum()
    return(val)

def Find_IntMomentMap(df,moment):
    img = len(df.columns)*[0]
    for i in range(len(df.columns)):
        img[i] = Find_MomentAtPixel(df,moment,i)
    return(img)

def Find_MaxMomentPixel(df,moment):
    img = Find_IntMomentMap(df,moment)
    idx = img.index(max(img))
    return(idx)

# Useful for when plotting in an iPython window (not a notebook)
def pmMap(head,df,moment):
    Nx = head["Nx"]
    Ny = head["Ny"]
    img = np.array(Find_IntMomentMap(df,moment)).reshape(Nx,Ny)
    plt.imshow(img,origin='lower')
    plt.colorbar()
    plt.show()

# Useful for when plotting in an iPython window (not a notebook)
def pSpec(df,pixels):
    if(type(pixels)==int):
        pixels = [pixels]
    for pix in pixels:
        df[pix].plot()
    plt.legend()
    plt.show()
