"""
    Varian.py
    Utilities for working with data from the varian systems (now Agilent)
    
    Developed in part by
    Andrew Curtis
    Craig Jones
    Martyn Klassen
    
    2009-2014 (c) Robarts Research Institute, Western University, Canada
    
    Main functions:
    
    readFDF
    readAllIMGs
    readFidDir
    
    also:
        readProcpar
        
    reading and reconstruction fids works like this:
    
    def demoFunction(fidFileName):
        dat, par = readFidDir(fidFileName)

        # in this example we have recievers, echoes, points and phase encodes
        # the ordering will depend on your seqcon and the sequence.
        # The procpars by default don't save information necessary to unambiguously determine
        # the data order, so you need to know the sequence + seqcon interaction.
        tt = dat.reshape( (par['nrcvrs'], par['ne'], par['np']/2.0, par['nv']) )
    
        # reorder dimensions, 
        tt = n.transpose(tt,(0,1,3,2))

        #fft last 2 dims (by default), here we have phaseXread as the last 2
        ft = ft2(tt)

        #combine rcvrs
        sos = sumsq(ft)
    
        #display
        figmri(sos)
    
    also see the associated fdf_io module which is a wrapper of a fast FDF 
    library written by Martyn Klassen.

    See functions at end of file for additional demonstration usage.

"""


import os
import os.path
import re
import numpy as n
import matplotlib.pyplot as plt
import struct

from numpy.fft import fft, ifft, fft2, fftn, fftshift, ifft2, ifftshift, ifftn

# Import python/C FDF I/O library
try:
    import fdf_io.fdflib as fl
except ImportError:
    print ('FDF IO c library not available')


# some utility functions to make running FFTs a little nicer.
# note that this uses fftpack, not FFTW, and is a little slow for large datasizes. 
# Also it is not smart enough to sequentially apply 1d FFT's and transpose for the 
# higher dimensional cases 
def ift(x,axis=[-1]):
    return ifftshift( ifft( ifftshift(x, axes=axis), axis=axis[0] ), axes=axis)

def ft(x,axis=[-1]):
    return fftshift( fft( fftshift(x, axes=axis), axis=axis[0] ), axes=axis)

def ift2(x):
    return ifftshift( ifft2( ifftshift(x, axes=[-2,-1]) , axes=[-2, -1]), axes=[-2,-1])

def ft2(x):
    return fftshift( fft2( fftshift(x, axes=[-2,-1]), axes=[-2, -1] ), axes=[-2,-1])

def ift3(x):
    return ifftshift( ifftn( ifftshift(x, axes=[-3, -2, -1]), axes=[-3, -2, -1] ), axes=[-3, -2,-1])

def ft3(x):
    return fftshift( fftn( fftshift(x, axes=[-3, -2, -1]), axes=[-3, -2, -1] ), axes=[-3, -2,-1])


#square root Sum of squares over axis
def sumsq(x, axis=0):
    return n.sqrt(n.mean(n.abs(x) ** 2, axis=axis))

#normalize
def norml(x):
    return x/n.max(n.abs(x))

# utility fn to display angle using figmri montage code
def da(map, clim=[None, None], cmap=plt.cm.jet, showbar=1, orient='e',interp='nearest'):
    figmri(n.angle(map), clim, cmap, showbar, orient, interp)
    
# utility fn to  display magnitude using figmri montage code
def dm(map, clim=[None, None], cmap=plt.cm.jet, showbar=1, orient='e', interp='nearest'):
    figmri(n.abs(map), clim, cmap, showbar, orient, interp)
   
   
#figmri montage tool for nD data.
# arguments:
# clim: color map limits
# cmap: color map. matplotlib has many to choose from
# showbar: display color bar
# orient: orientation ('v' vertical stack, 'h' horizaontal stack, 't' transpose)
# interp: interpolation mode for imshow(). 
# newfig: generate new figure or no
def figmri(data, clim=[None, None], cmap=plt.cm.gray, showbar=True, orient='e',interp='nearest', newfig=False):
    
    if newfig:
        plt.figure()

    if(n.iscomplexobj(data)):
        data = n.abs(data)
    
    data = data.squeeze()
    ndim = data.ndim
    orient = str.lower(orient)
    N = 1
    
    # find a decent reshaping of the array if it is high dimensional
    if ndim > 3:
        data = n.transpose(data, range(0,ndim-3) + [-2, -3, -1] )
        c = data.shape
        print c
        print ([n.prod(c[:2]), n.prod(c[2:])])
        data = data.reshape([n.prod(c[:2]), n.prod(c[2:])])
    elif ndim == 3:
        c = data.shape
        N=1
        if orient == 'v':
            N = 1
        elif orient == 'h':
            N = c[0]
        else:
            for k in n.arange(int(n.floor(n.sqrt(c[0]))),0,-1):
                if n.remainder(c[0],k) == 0:
                    N=k
                    break
                    
        data = data.reshape([c[0]/N, N] + c[1:] )
        data = n.transpose(data, [0, 2, 1, 3])
        c = data.shape
        data = data.reshape([n.prod(c[:2]), n.prod(c[2:])])
        
    if orient == 't':
        data = n.rot90(data)
    
    plt.imshow(data,cmap=cmap,interpolation=interp)
    plt.clim(clim)
    
    
    
    if showbar:
        plt.colorbar()
    
    


def stripquotes(data):
    "strip whitespace and quotes from vnmrj string params"
    return data.replace("\"", "").strip()


# accessory function to read IMG or FDF
def readImages( filename ):
    if filename.endswith('.fdf'):
        data = readFDF( filename )
    elif filename.endswith('.img'):
        data = readIMG( filename )
    else:
        print "Unknown filename %s " % (filename)

    return data

# Python implementation to read FDFs
def readFDF_py( filename ):
    fp = open( filename, 'rb' )
    xsize = -1
    ysize = -1
    zsize = 1
    bigendian = -1
    done = False
    while not done :
        line = fp.readline()
        if( len( line ) >= 1 and line[0] == chr(12) ):
            break

        if ( len( line) == 0 ):
            break

        if( len( line ) >= 1 and line[0] != chr(12) ):

            if( line.find('bigendian') > 0 ):
                endian = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('echos') > 0 ):
                nechoes = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('echo_no') > 0 ):
                echo_no = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('nslices') > 0 ):
                nslices = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('slice_no') > 0 ):
                sl = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('matrix') > 0 ):
                m = re.findall('(\d+)', line.rstrip())

                if len(m) == 2:
                    xsize, ysize = int(m[0]), int(m[1])
                elif len(m) == 3:
                    xsize, ysize, zsize = int(m[0]), int(m[1]), int(m[2])

    fp.seek(-xsize*ysize*zsize*4,2)

    if bigendian == 1:
        fmt = ">%df" % (xsize*ysize*zsize)
    else:
        fmt = "<%df" % (xsize*ysize*zsize)

    data = struct.unpack(fmt, fp.read(xsize*ysize*zsize*4))
    data = n.array( data ).reshape( [zsize, ysize, xsize ] ).squeeze()

    fp.close()

    return data


# Wrapper for FDF reading. Try to use C based reader (fdflib -- compile shared library separately)
def readFDF( filename ):
    try:
        header = fl.readHeader_c(filename)
        data = fl.readFDF_c(filename, header)
    except:
        print "Something went wrong with c based fdf reader... reverting to python version."
        data = readFDF_py(filename)

    return data


# Load in all *.IMG files in a directory
def readAllIMGs(directory, coil=0, slice=0, image=0):
    (head, tail) = os.path.split(directory)
    dirs = os.listdir(head)

    dirs = [dir for dir in dirs if dir.startswith(tail)]

    if coil > 0:
        coilstr = '{:03}.img'.format(coil)
        dirs = [ dir for dir in dirs if dir.endswith(coilstr) ]

    data = []
    count = 1
    for dir in sorted(dirs):
        print "Reading coil %d/%d..." % (count, len(dirs))
        print "%s" % dir
        count = count + 1
        data.append( readIMG( os.path.join(head, dir), slice, image ) )
       
    #squash the list into an array.
    d2 = n.array(data)
    del data

    return d2



def readIMG(directory, slice=0, image=0, returnNames=False):
    # Get a list of all the FDF files in the directory
    try:
        files = os.listdir(directory)
    except:  
        print "Could not find the directory %s" % directory
        return
        
    files = [ file for file in files if file.endswith('.fdf') ]
    print "Found %d files in %s." % (len(files), directory) 

    if slice > 0:
        slicestr = 'slice{:03}'.format(slice)
        files = [ file for file in files if (file.find(slicestr)>-1) ]

    if image > 0:
        imgstr = 'image{:03}'.format(image)
        files = [ file for file in files if (file.find(imgstr)>-1) ]
        

    count = 1
    data = []
    for file in sorted(files):
        print "Reading file %d/%d..." % (count, len(files))
        print "%s" % file
        count = count + 1
        data.append( readFDF( directory+'/'+file ) )
       
    #squash the list into an array.
    d2 = n.array(data)
    del data

    if returnNames:
        return d2,files
    else:
        return d2



def parseProcpar(inlist):
    """parseProcpar - directly parse vnmrj procpar file
    inputs:
        inlist -- list of strings from file.readlines()
        
    returns:
        paramDict -- dicitonary of key value pairs
                        { 'param name' : 'value' }
                        
    usage:
        pd = parseProcpar(lines)
        
    String parameters are returned as lists of strings.
    Numeric parameters are returned as numpy arrays, which means this works:
        pd['np']/2.0
        
    """
    # Procpar parameter definitions
    # Names and types of paramater info line
    paramLineOne = ['name', 'subtype', 'basictype', 'maxvalue', 'minvalue', 'stepsize', 'dgroup', 'ggroup', 'protection', 'active', 'intptr']
    
    paramTypes = [str, int, int, float, float, float, int, int, int, int, int]
    
    """
    from vnmrj reference, valid subtypes are:
    0 (undefined), 1 (real), 2 (string),
    3 (delay), 4 (flag), 5 (frequency), 6 (pulse), 7 (integer).
    """
    
    #Scalings for param subtypes (not currently used)
    paramScalings = [0, 1, 0, 1, 1, 1, 1e-3, 1]
    paramSubtypeConversions = dict(zip(range(0,len(paramScalings)),paramScalings))
    
    paramDict = {}
    
    #state parameters
    parsedInfo = 0
    parsedData = 0
    parsedDataRem = 0
    
    for item in inlist:
        
        if parsedInfo == 0:
            #print "parsing infos"
            #parse the first line of the parameter entry
            paramInfo = {}
            
            for name, oper, data in zip(paramLineOne, paramTypes, item.strip().split(" ") ):
                paramInfo[name] = oper(data)
                
            parsedInfo = 1
            continue
            
        else:
            #have already parsed info line
            
            if parsedData == 0:
                #print "parsing data first"
                #parse the first data line
                
                splitValues = item.split(" ")
                nvals = int(splitValues[0])
                
                #switch depending on what kind of value is it (stored differently)
                if paramInfo['basictype'] == 1:
                    #convert param values to a numpy array and insert in param dictionary
                    splitValues = stripquotes(item).split(" ")
                    floatData = n.array(map(float, splitValues[1:]))
                    
                    #'pulse' parameters are specified in ms, convert to value in seconds
                    if paramInfo['subtype'] == 6:
                        floatData *= 1e-3
                        
                    paramDict[paramInfo['name']] = floatData
                    parsedDataRem = 0
                    
                elif paramInfo['basictype'] == 2:
                    parsedDataRem = nvals - 1
                    if parsedDataRem>0:
                        #convert param values to a list and insert in param dictionary
                        paramDict[paramInfo['name']] = [stripquotes(splitValues[1]),]
                    else:
                        paramDict[paramInfo['name']] = stripquotes(splitValues[1])
                        
                        
                parsedData = 1
                
                continue
                
            else:
                
                if parsedDataRem > 0:
                    #print "parsing data remaining"
                    
                    #we've started parsing data but haven't finished the multiline string
                    paramDict[paramInfo['name']].append(item.strip())
                    parsedDataRem -= 1
                    continue
                    
                #fallthrough... all other things have been parsed, so we're on the enum line
                # could read it in, but I'm just throwing it out
                
                #last line on entry, reset flags and skip
                parsedInfo = 0
                parsedData = 0
                parsedDataRem = 0
                paramInfo = {}
                
                
    paramDict['nrcvrs'] = paramDict['rcvrs'].count('y')
    return paramDict

# read procpar
def readProcpar(fname = ''):
    
    if fname.endswith('/'):
        fname = fname[:-1]
    
    if os.path.isdir(fname):
        fname = fname + '/procpar'
        
    pfile = file(fname,'r')
    parlines = pfile.readlines()
    pfile.close()
    
    return parseProcpar(parlines)


# readFID
# main utility ofr reading a fid data file
# use readFidDir to handle .fid directories and procpars
def readFID(filename, verbose=False):
    
    dataFileHeaderMap = ["nblocks", "ntraces", "np", "ebytes",
                            "tbytes", "bbytes", "vers_id", "status", "nbheaders" ]

    dataBlockHeaderMap = ["scale", "status", "index", "mode",
                            "ctcount", "lpval", "rpval", "lvl", "tlt" ]

    dataFileHeadString = "6l2hl"
    dataFileHeadString = "6i2hi"

    dataFileHeadSize = struct.calcsize(dataFileHeadString)

    dataBlockHeadString = "4hl4f"
    dataBlockHeadString = "4hi4f"

    dataBlockHeadSize = struct.calcsize(dataBlockHeadString)

    fp = open( filename, 'rb' )

    line = fp.read(dataFileHeadSize)

    endian = ">"

    header = struct.unpack(endian+dataFileHeadString, line)
    fileHeader = dict(zip(dataFileHeaderMap, header))

    if verbose:
        print fileHeader

    if fileHeader['status']<0:
        endian = "<"
        header = struct.unpack(endian+dataFileHeadString, line)
        fileHeader = dict(zip(dataFileHeaderMap, header))
        if verbose:
            print fileHeader

    elementString = "h"
    dtypeString = 'complex'
    
    if fileHeader['ebytes'] == 2:
        elementString = "h"
        dtypeString = 'complex64'
    elif fileHeader['ebytes'] == 4:
        elementString = "f"
        dtypeString = 'complex64'
    elif fileHeader['ebytes'] == 8:
        elementString = "d"
        dtypeString = 'complex128'
    

    traceSize = fileHeader['np']*fileHeader['ebytes']
    blockSize = fileHeader['ntraces']*traceSize

    if verbose:
        print "traceSize: %d  blockSize: %d" % (traceSize, blockSize)

    data = n.zeros( (fileHeader['nblocks'],fileHeader['np']/2.0*fileHeader['ntraces']), dtype=dtypeString )

    for block in xrange(fileHeader['nblocks']):
        if verbose > 1:
            print "Block %d" % block
        line = fp.read(dataBlockHeadSize)
        blockheader = struct.unpack(endian+dataBlockHeadString, line)

        #print dict(zip(dataBlockHeaderMap, blockheader))
        
        #Ignore the block header for now

        #Read all bytes in block
        line = fp.read(blockSize)

        thisBlock = struct.unpack("%s%d%s" % (endian, fileHeader['np']*fileHeader['ntraces'], elementString), line)

        temp = n.array(thisBlock)
        data[block,:] = temp[::2] + temp[1::2]*1j


    fp.close()
    return data



# read .fid directory and get fid data and procpar
def readFidDir(fname):
    """ readFidDir
    
    usage:
        data, par = readFidDir('path/to/my/data.fid')
    
    """
    #check if file exists
    #if not, check if file.fid exists
    if not (fname.endswith('.fid') or fname.endswith('.FID')):
        fname = fname + '.fid'

    if (os.path.exists(fname+'/fid') and os.path.exists(fname+'/procpar')):
        procpar = readProcpar(fname+'/procpar')
        fiddata = readFID(fname+'/fid')
    else:
        return -1

    return (fiddata, procpar)



#Seqcon goes like:
# ne, ns, npe, npe2, npe3 (unused)
def readPETable(fname, shot=0):
    fp = open( fname, 'r' )
    
    line = fp.readlines()
    
    fp.close()
    
    temp = n.empty(0,'int32')
    for l in line[:]:
        temp = n.r_[temp, n.array(map(int,l.strip().split()))]

    permuteOrder = n.argsort(temp)
    
    return permuteOrder, temp
    

# PE table for FSE sequences following Varian's macros
def makePETable(nv, etl, k0=1):
    tabval = n.zeros(nv,dtype='int32')
    
    count = n.arange(nv/2)+1
    tabval[count-1]= -((count*nv/(2*etl)-1) % (nv/2))
    tabval[count-1]= tabval[count-1] + n.floor((count-1)/etl)
    
    count = n.arange(nv/2,nv)+1
    tabval[count-1] = (count-nv/2-1)*nv/(2*etl) % (nv/2)
    tabval[count-1] = tabval[count-1] + n.floor((count-nv/2-1)/etl) +1
    
    #plot(tabval)
    tabval = tabval.reshape((nv/etl,etl))
    
    #cycle columns of table to align k0 traversal
    inds2 = n.r_[n.arange(etl-1,etl-k0-1,-1),n.arange(etl-k0)]
    t = tabval[:, inds2]
    return t




def simpleReconFSE(fname, tabName):
    
    
    dat,par = readFidDir(fname)

    #test reshape
    #permuteOrder,test = readPETable('../data/fse256_8_4')
    permuteOrder,test = readPETable('/Users/andrewcurtis/temp/fse256_8_4')

    #get data sizes from procpar
    shots = par['nv'][0]/par['etl'][0]
    etl = par['etl'][0]
    np = par['np'][0]/2.0
    ns =  par['ns'][0]
    nv = par['nv'][0]
    nrcvrs = par['nrcvrs']*1.0

    sizes = (nrcvrs, shots, ns, etl, np)

    #tt = dat.reshape( sizes ).copy()
    tt = dat.reshape( sizes )
    dat=None

    tt = n.transpose(tt,(0,2,1,3,4))

    tt = tt.reshape((nrcvrs,ns,nv,np))

    #reorder for echoes
    tt =n.take(tt,permuteOrder,axis=2)

    tun = tt[:,10,:,:].squeeze()
    tun[:,1::2,:] = 0
    ft = ft2(tun)

    t2 = ft2(tt)

    sos = sumsq(tt)




def ernst(tr,t1):
    return 180.0/pi*(math.acos(math.exp(-tr/t1)))



if __name__ == '__main__':
    #main()
    pass





