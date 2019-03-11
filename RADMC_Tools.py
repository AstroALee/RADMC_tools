'''
    Set of analysis and plotting wrappers that can be called for RADMC-3d
    Includes carving out routines for YT


    Author: Aaron T. Lee, Spring/Summer 2018


    Free to share, use, blah blah.
    I take no responsibility if any errors that exist below screw up your analysis.
    Don't treat code like a black box.

'''

# Bringing in relevant modules
import argparse
import numpy as np
import PlotRoutines
import os.path
import sys


''' Global set of constants (cgs)'''
class Csts:
    c = 3.0e10 # speed of light
    kB = 1.3806e-16 # Boltz constant
    NH3lc = 23.7109e9 # Hz



''' Relevant functions are defined here '''

# Code that governs the command line parameters.
parser = argparse.ArgumentParser(
    description='''Let's make some plots!''',
    epilog="Now get to work!")
parser.add_argument("plot_type", help="What do you want to do? 'genintegralmap' generates " + \
                                    "the integrated intensity data and saves it to a file. " + \
                                    "'plotmap' plots the map from the data. If the file cannot " + \
                                    "be found, it generates it first. 'plotspec' plots the 1D spectrum " + \
                                    "for a given pixel, which you must specify with the optional flag.",
    choices=['genintegralmap', 'genmomentmap','plotmap','plotspec'])
parser.add_argument('-if','--infile',help='input file name (default=image.out)',
    default='image.out',type=str)
parser.add_argument('-of','--outfile',help='end of output file name (default=[plottype].outfile.png)',
    default='outfile.png',type=str)
parser.add_argument('-px','--pixel', nargs='*', help='pixel numbers to study (can also provide -1 for max, which is also the default)',
    default=-1,type=int)
parser.add_argument('-mo','--moment', help='Moment map to calculate (0th is the default)',
    default=0,type=int)
parser.add_argument('-nm','--normalize', help='Normalize the moment map by the integrated intensity? (default 0)',
        default='0',choices=[0,1],type=int)


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


''' Function that grabs the header information from your input file '''
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
    wave_data = 1.0e-4*np.array(wave_data) # converts to cgs

    # BOOOOO for all those blank lines
    while(True):
        line = f.readline()
        if(len(line.strip()) == 0):
            header_length = header_length + 1
            header_size = header_size + len(line)
        else:
            break

    f.close()

    return header_data, wave_data, header_length, header_size



'''
    Function that generates the \Delta\nu array from your wavelength array,
    read in from your header
'''
def Create_deltanu(wave):
    dnu = np.zeros(head["Nwavelengths"])
    dwa = np.zeros(head["Nwavelengths"])

    nu = Csts.c/wave
    for i in range(head["Nwavelengths"]):
        if(i==0):
            dnu[0] = 2.0*0.5*np.fabs(nu[0]-nu[1])
            dwa[0] = 2.0*0.5*np.fabs(wave[0]-wave[1])
        elif(i==head["Nwavelengths"]-1):
            dnu[-1] = 2.0*0.5*np.fabs(nu[-1]-nu[-2])
            dwa[-1] = 2.0*0.5*np.fabs(wave[-1]-wave[-2])
        else:
            dnu[i] = 0.5*np.fabs(nu[i]-nu[i-1]) + 0.5*np.fabs(nu[i+1]-nu[i])
            dwa[i] = 0.5*np.fabs(wave[i]-wave[i-1]) + 0.5*np.fabs(wave[i+1]-wave[i])

    return(dnu,dwa)



''' Writes the integrated moment map to a file '''
def SaveMomMap():
    print("Saving integrated moment map to a file.")
    f = open('intmap.' + str(args.moment) + '.file','w')

    # print header
    curline = str(head["Nx"]) + " " + str(head["Ny"]) + " " + \
                str(head["Nwavelengths"]) + " " + str(head["Pixel_Size_X"]) + \
                " " + str(head["Pixel_Size_Y"]) + " \n"
    f.write(curline)

    for i in range(head["Nx"]*head["Nx"]):
        curline = str(intmap_array[i]) + " \n"
        f.write(curline)

    f.close()
    return




''' Writes the integrated intensity map to a file '''
def SaveIntMap():
    print("Saving integrated map to a file.")
    f = open('intmap.file','w')

    # print header
    curline = str(head["Nx"]) + " " + str(head["Ny"]) + " " + \
                str(head["Nwavelengths"]) + " " + str(head["Pixel_Size_X"]) + \
                " " + str(head["Pixel_Size_Y"]) + " \n"
    f.write(curline)

    for i in range(head["Nx"]*head["Nx"]):
        curline = str(intmap_array[i]) + " \n"
        f.write(curline)

    f.close()
    return




''' Reads in the integrated intensity map data from a file '''
def ReadIntMap():

    fname = 'intmap.' + str(args.moment) + '.file'
    print("Reading in " + fname)
    f = open(fname,'r')

    # grab header
    tNx,tNy,tNwave,tPixX,tPixY = f.readline().split()

    # Check to make sure this matches the current header info
    assert(int(tNx)==head["Nx"] and int(tNy)==head["Ny"])
    assert(int(tNwave)==head["Nwavelengths"])
    assert(float(tPixX)==head["Pixel_Size_X"] and float(tPixY)==head["Pixel_Size_X"])

    dim = int(tNx)*int(tNy)
    intmap = np.zeros(dim)

    for idx in range(dim):
        line = f.readline().split()
        intmap[idx] = float(line[0])

    f.close()
    return(intmap)



''' Reads in data to create an 1D spectrum for the given pixel(s) '''
def Create_Spectrum(Sidx):

    dim = int(head["Nx"]*head["Ny"])
    nwave = int(head["Nwavelengths"])

    f = open(args.infile,'r',1048576)
    f.seek(header_size) # Moves file pointer to the line after the header

    #Sidx = pix
    if any(x >= dim for x in Sidx):
        raise ValueError("Index is incorrect for this file (User = " + str(Sidx) + \
                            " vs dim = " + str(dim) )

    for x in Sidx:
        if(x<0):
            if(os.path.isfile('intmap.file')):
                print("Searching for maximum intensity pixel")
                intmap_array = ReadIntMap()
                intmax = max(intmap_array)
                max_idx = list(intmap_array).index(intmax)
                print("Max intensity at pixel " + str(max_idx))
                if any(max_idx==int(y) for y in Sidx):
                    print("Already creating spectrum at max pixel")
                    break # no need to do it twice
                else:
                    Sidx[ Sidx.index(x) ] = max_idx
            else:
                raise ValueError("intmap.file not found to find max!")


    spec = [ [] for i in range(len(Sidx)) ]

    idx=0
    for line in f: # Starts from current location of file pointer
        if(line.isspace()): continue
        for sp in Sidx:
            if(sp == idx%dim):
                sp_idx = Sidx.index(sp)
                #spec[sp_idx].append( np.float64(line.strip()) )
                spec[sp_idx].append( intToTemp(np.float64(line.strip()),Csts.c/wave[sp_idx]) )
        idx=idx+1

    f.close()

    #print(spec)
    return(spec)



''' Creates the integrated intensity moment map from the input file '''
def Create_MomentMap(thismoment):
    print("\nCreating " + str(thismoment) + " moment map. For high-res this could take a while.")
    print("Go refill your coffee.")

    # integrated map array
    dim = int(head["Nx"]*head["Ny"])
    nwave = int(head["Nwavelengths"])

    intmap = np.zeros(dim)

    f = open(args.infile,'r',1048576) # byte number should speed things up a little
    f.seek(header_size) # skip past the header

    NH3lc_wave = Csts.c/Csts.NH3lc

    idx = 0
    for line in f:
        if(idx==dim*nwave): break
        if(line.isspace()): continue
        wnum = int(idx*1.0/dim)
        val = np.float64(line.strip())
        #dnu = delta_nu[wnum]
        dwa = delta_wa[wnum]
        vel = Csts.c*(wave[wnum]-NH3lc_wave)/NH3lc_wave # cgs velocity from line center
        # Converts to Temperature/velocity
        val = intToTemp(val,Csts.c/wave[wnum])
        dvel = dwa*Csts.NH3lc  # Csts.NH3 = freq = c/wave
        #intmap[idx%dim] = intmap[idx%dim] + val*dnu
        if(thismoment==0): # in case 0^0 leads to issues
            intmap[idx%dim] = intmap[idx%dim] + val*dvel
        else:
            intmap[idx%dim] = intmap[idx%dim] + val*pow(vel,thismoment)*dvel
        # Prepare for next cell
        idx = idx + 1

    f.close()

    print("Done creating moment map. Enjoy.")

    if(args.normalize==1 and thismoment>0):
        print("Normalizing")
        intmap0 = np.zeros(dim)
        intmap0 = Create_MomentMap(0)
        for i in range(dim):
            if(intmap0[i]==0):
                intmap[i] = 0 # redundant, how can pos-def int(Idv) = 0 but int(Ivdv) not be?!
            else:
                intmap[i] = intmap[i]/intmap0[i]
        print("Done normalizing")

    return intmap



''' Creates the integrated intensity map from the input file '''
def Create_IntegratedIntMap():
    print("\nCreating integrated map. For high-res this could take a while.")
    print("Go refill your coffee.")

    # integrated map array
    dim = int(head["Nx"]*head["Ny"])
    nwave = int(head["Nwavelengths"])

    intmap = np.zeros(dim)

    f = open(args.infile,'r',1048576) # byte number should speed things up a little
    f.seek(header_size) # skip past the header

    idx = 0
    for line in f:
        if(idx==dim*nwave): break
        if(line.isspace()): continue
        wnum = int(idx*1.0/dim)
        val = np.float64(line.strip())
        dnu = delta_nu[wnum]
        dwa = delta_wa[wnum]
        # Converts to Temperature/velocity
        val = intToTemp(val,Csts.c/wave[wnum])
        dvel = dwa*Csts.NH3lc  # Csts.NH3 = freq = c/wave
        #intmap[idx%dim] = intmap[idx%dim] + val*dnu
        intmap[idx%dim] = intmap[idx%dim] + val*dvel
        # Prepare for next cell
        idx = idx + 1

    f.close()

    print("Done creating integrated map. Enjoy.")
    return intmap




'''
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          Actual 'program' is here
    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
'''

# Parse the command line
args = parser.parse_args()
if(type(args.pixel)==int): args.pixel = [args.pixel]
print(args)


# Parse the header of your output file (only need to do this once)
head , wave , header_length, header_size = grabRADheader(args.infile)

# Create delta nu array (simple and quick; let's just do it now)
delta_nu , delta_wa = Create_deltanu(wave)

# call the relevant routine
if(args.plot_type=='genintegralmap'):
    intmap_array = Create_IntegratedIntMap()
    SaveIntMap()
elif(args.plot_type=='genmomentmap'):
    intmap_array = Create_MomentMap(args.moment)
    SaveMomMap()
elif(args.plot_type=='plotmap'):
    fname = 'intmap.' + str(args.moment) + '.file'
    if(os.path.isfile(fname)):
        intmap_array = ReadIntMap()
    else:
        print("intmap." + str(args.moment) + ".file not found -- generating!")
        intmap_array = Create_MomentMap(args.moment)
        SaveMomMap()
    PlotRoutines.PlotIntMap(head,intmap_array,args)
elif(args.plot_type=='plotspec'):
    spec_array = Create_Spectrum(args.pixel)
    PlotRoutines.PlotSpec(head,wave,spec_array,args)
else:
    print("Something has gone terribly wrong... how did we get here!?")

print("All done.")
