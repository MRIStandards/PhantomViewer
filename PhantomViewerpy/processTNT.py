#!/usr/bin/python

import sys
import io
from collections import OrderedDict
import datetime
from time import gmtime
import re
import numpy as np
from numpy.fft import fftfreq, fftshift
import string

import TNTdtypes
from utils import convert_si, read_pascal_string, make_str as s


class TNTfile:

    def __init__(self, tntfilename):

        self.TNMRComment=''
        self.tnt_sections = OrderedDict()

        with open(tntfilename, 'rb') as tntfile:

            self.tntmagic = np.fromstring(tntfile.read(TNTdtypes.Magic.itemsize),
                                          TNTdtypes.Magic, count=1)[0]

            if not TNTdtypes.Magic_re.match(self.tntmagic):
                raise ValueError("Invalid magic number (is '%s' really a TNMR file?): %s" % (tntfilename, self.tntmagic))

            ##Read in the section headers
            tnthdrbytes = tntfile.read(TNTdtypes.TLV.itemsize)
            while(TNTdtypes.TLV.itemsize == len(tnthdrbytes)):
                #print('unpacking header')
                tntTLV = np.fromstring(tnthdrbytes, TNTdtypes.TLV)[0]
                data_length = tntTLV['length']
                hdrdict = {'offset': tntfile.tell(),
                           'length': data_length,
                           'bool': bool(tntTLV['bool'])}
                if data_length <= 4096:
                    hdrdict['data'] = tntfile.read(data_length)
                    assert(len(hdrdict['data']) == data_length)
                else:
                    tntfile.seek(data_length, io.SEEK_CUR)
                self.tnt_sections[s(tntTLV['tag'])] = hdrdict
                tnthdrbytes = tntfile.read(TNTdtypes.TLV.itemsize)
            # Find and parse the delay tables
            self.DELAY = {}     #Tables with default naming convention
            # This RegExp should match only bytearrays containing
            # things like "deXX:X" or so.
            #_re = re.compile(b'de[0-9]+:[0-9]')
            dealyorloopre=b'((at|lp)[0-9]+:[0-9])|((de|lp)[0-9]+:[0-9])|TI_[1-9]|ti_times|cpmgloop|gr0:2|GxTable|GyTable|gztable|rdArray|seloop|HPMW|TX3:2|teDelay|tiDelay|sliceFrequencies'   #match de0:0 or lp0:0 or TI_0 or ti_times or cpmgloop or gr0:2
            delay_re = re.compile(dealyorloopre)   #SER added to find loop tables
            # seek well past the data section so we read as little
            # of the file into memory as possible
            tntfile.seek(self.tnt_sections["PSEQ"]["offset"])
            search_region = tntfile.read()
            #print (search_region)
            # Do the search, and iterate over the matches
            for match in delay_re.finditer(search_region):
                # Lets go back and read the section properly. Offset back by
                # four to capture the delay table name length
                #print('reading tables and PS params')
                offset = (match.start() - 4)
                # extract the name length, the name, the delay length,
                # and the delay.
                delay_name = read_pascal_string(search_region[offset:])
                offset = match.start() + len(delay_name)
                delay = read_pascal_string(search_region[offset:])
                # Now check for delay tables of length one and discard them
                if len(delay) > 1:
                    delay = delay.split()
                    delay = convert_si(delay)
                    self.DELAY[delay_name] = delay
            try:
                ps=search_region.split(b'Sequence')[1].split(b'INFO')[0]      #kludgy way to pull out T90, T180 and other ps parmaters
                ps=ps[30:]            
                try:
                    T180s=ps.split(b'T180')[1].split(b'u')[0].decode('utf-8')
                    self.T180=[float(st) for st in re.findall(r'-?\d+\.?\d*', T180s)][0]
                except:
                    self.T180=0.0
                try:
                    T90s=ps.split(b'T90')[1].split(b'u')[0].decode('utf-8')
                    self.T90=[float(st) for st in re.findall(r'-?\d+\.?\d*', T90s)][0]
                except:
                    self.T90=0.0
                try:
                    taus=ps.split(b'tau')[1].split(b'm')[0].decode('utf-8')
                    self.tau=[float(st) for st in re.findall(r'-?\d+\.?\d*', taus)][0]
                except:
                    self.tau=0.0
                try:
                    gramp=ps.split(b'GrAmp')[1][:20].decode('utf-8')        #Find gradient amplitude parameter, take 20 bytes after the key 'GrAmp', turn into string
                    self.GrAmp=[float(st) for st in re.findall(r'-?\d+\.?\d*', gramp)][0]       #search for a number
                except:
                    self.GrAmp=0
                try:
                    gs=ps.split(b'Gs')[1][:20].decode('utf-8')        #Find gradient amplitude parameter, take 20 bytes after the key 'GrAmp', turn into string
                    self.Gs=[float(st) for st in re.findall(r'-?\d+\.?\d*', gs)][0]       #search for a number
                    print('gs=',Gs)
                except:
                    self.Gs=0   
                try:
                    graramp=ps.split(b'GradRamp')[1][:20].decode('utf-8')        #Find gradient ramptime parameter, take 20 bytes after the key 'GradRamp', turn into string
                    self.gradPulseRiseTime=[float(st) for st in re.findall(r'-?\d+\.?\d*', graramp)][0]/1000       #search for a number, assume it is in ms, convert to seconds
                except:
                    self.gradPulseRiseTime=-1       #flag to indicate value not found
                    
                try:
                    grad=ps.split(b'Grad')[1][:20].decode('utf-8')        #Find gradient amplitude parameter, take 20 bytes after the key 'GrAmp', turn into string
                    self.Grad=[float(st) for st in re.findall(r'-?\d+\.?\d*', grad)][0]      #search for a number, 
                except:
                    self.Grad=0
                try:
                    delta=ps.split(b'Delta')[1][:20].decode('utf-8')        #Find gradient amplitude parameter, take 20 bytes after the key 'GrAmp', turn into string
                    self.Delta=[float(st) for st in re.findall(r'-?\d+\.?\d*', delta)][0]       #search for a number
                except:
                    self.Delta=0 
                try:
                    NutIncrement=ps.split(b'NutIncrement')[1].split(b'u')[0].decode('utf-8')
                    self.NutIncrement=[float(st) for st in re.findall(r'-?\d+\.?\d*', NutIncrement)][0]
                except:
                    self.NutIncrement=1       
                try:        #pull out comment, get rid of nonascii and line feeds
                    comment=ps.split(b'CMNT')[1]
                    endcomment=comment.find(b'TMG3')
                    self.TNMRComment= comment[7:endcomment].decode('unicode_escape') #('UTF-16LE ')  # ('utf-8') #
                    if self.TNMRComment.find('Slice Offset')>-1:
                        self.TNMRComment=self.TNMRComment.strip(' ')
                    else:
                        self.TNMRComment=self.TNMRComment.replace('\r\n', ';  ')        # replace newlines with semicolons
                        self.TNMRComment=re.sub(r'[^a-zA-Z0-9;=();:.,* ]', '', self.TNMRComment)       #get rid of all nonascii characters
                        #self.TNMRComment=re.sub(r'_{2,20}', '_', self.TNMRComment)
                except:
                    raise
                    self.TNMRComment='No Comment'                  
            except:
                raise                
#        
        assert(self.tnt_sections['TMAG']['length'] == TNTdtypes.TMAG.itemsize)
        self.TMAG = np.fromstring(self.tnt_sections['TMAG']['data'],
                                  TNTdtypes.TMAG, count=1)[0]

        assert(self.tnt_sections['DATA']['length'] ==
               self.TMAG['actual_npts'].prod() * 8)
        ## For some reason we can't set offset and shape together
        #DATA = np.memmap(tntfilename,np.dtype('<c8'), mode='r',
        #                 offset=self.tnt_sections['DATA']['offset'],
        #                 shape=self.TMAG['actual_npts'].tolist(),order='F')
        self.DATA = np.memmap(tntfilename, np.dtype('<c8'), mode='c',
                              offset=self.tnt_sections['DATA']['offset'],
                              shape=self.TMAG['actual_npts'].prod())
        self.DATA = np.reshape(self.DATA,
                               self.TMAG['actual_npts'],
                               order='F')

        assert(self.tnt_sections['TMG2']['length'] == TNTdtypes.TMG2.itemsize)
        self.TMG2 = np.fromstring(self.tnt_sections['TMG2']['data'],
                                  TNTdtypes.TMG2, count=1)[0]                           
                                  
        #print('end of tnt file processing')
        tntfile.closed


#    def writefile(self, outfilename):
#        outfile = open(outfilename, 'wb')
#        outfile.write(self.tntmagic)
#        for tag in self.tnt_sections_order:
#            tlv = np.asarray(self.tnt_sections[tag].items(), dtype=TNTdtypes.TLV)
#

    @property
    def start_time(self):
        """The time when the NMR acquisition was started

        No timezone information is available"""
        time_struct = gmtime(self.TMAG['start_time'])
        return datetime.datetime(*time_struct[:6])

    @property
    def finish_time(self):
        """The time when the NMR acquisition ended

        No timezone information is available"""
        time_struct = gmtime(self.TMAG['finish_time'])
        return datetime.datetime(*time_struct[:6])

    @property
    def date(self):
        """The time when the file was saved

        No timezone information is available"""
        strlen = self.TMAG['date'].index(b'\x00')
        if sys.version_info.major <= 2:
            datestr = str(self.TMAG['date'][:strlen])
        else:
            datestr = str(self.TMAG['date'][:strlen], encoding='ascii')
        return datetime.datetime.strptime(datestr, "%Y/%m/%d %H:%M:%S")
    
    @property
    def gradOrientation(self):
        """gradient oprinetation string XYZ=read, phase, slice"""
        go = str(self.TMAG['grd_orientation'])
        return go   

    def __getattr__(self, name):
        """Expose members of the TMAG and TMG2 structures as attributes"""
        if name in self.TMAG.dtype.names:
            return self.TMAG[name]
        elif name in self.TMG2.dtype.names:
            return self.TMG2[name]
        else:
            raise AttributeError("'%s' is not a member of the TMAG or TMG2 structs" % name)

    def LBfft(self, LB=0, zf=0, phase=None, logfile=None, ph1=0,
              DCoffset=None, altDATA=None):
        if altDATA is None:
            DATA = self.DATA
        else:
            DATA = altDATA
        LBdw = -LB * self.dwell[0] * np.pi  # Multiply by pi to match TNMR
        npts = DATA.shape[0]
        npts_ft = npts * (2 ** zf)

        if DCoffset is None:
            # Taking the last eighth of the points seems to give OK (but not
            # perfect) agreement with the TNMR DC offset correction.
            # This hasn't been tested with enough different values of npts
            # to be sure that this is the right formula.
            DCoffset = np.mean(DATA[int(npts / -8):, :, :, :],
                               axis=0, keepdims=True)
            if logfile is not None:
                logfile.write("average DC offset is %g\n" % np.mean(DCoffset))

        lbweight = np.exp(LBdw * np.arange(npts, dtype=float))
        DATAlb = (DATA - DCoffset) * lbweight[:, np.newaxis, np.newaxis, np.newaxis]

        DATAfft = np.fft(DATAlb, n=npts_ft, axis=0)
        DATAfft = fftshift(DATAfft, axes=[0])
        DATAfft /= np.sqrt(npts_ft)  # To match TNMR behaviour

        if phase is None:  # Phase automatically
            DATAfft *= np.exp(-1j * np.angle(np.sum(DATAfft)))
        else:
            DATAfft *= np.exp(1j * (phase + ph1 * np.linspace(-0.5, 0.5, npts_ft))
                              )[:, np.newaxis, np.newaxis, np.newaxis]

        return DATAfft

    def freq_Hz(self, altDATA=None):
        """Returns the frequency axis (in Hz) for the NMR spectrum"""
        if altDATA is None:
            npts = self.actual_npts[0]
        else:
            npts = altDATA.shape[0]
        dw = self.dwell[0]
        ref_freq = self.ref_freq

        return -(fftshift(fftfreq(npts, dw)) + ref_freq)

    def freq_ppm(self, altDATA=None):
        """Returns the frequency axis (in ppm) for the NMR spectrum"""
        NMR_freq = self.ob_freq[0]
        return self.freq_Hz(altDATA) / NMR_freq

    def fid_times(self, altDATA=None):
        """Returns the time axis (in s) for the FID"""
        if altDATA is None:
            npts = self.actual_npts[0]
        else:
            npts = altDATA.shape[0]
        dw = self.dwell[0]

        return np.arange(npts) * dw

    def ppm_points(self, max_ppm, min_ppm, altDATA=None):
        """Given a maximum and minimum frequency (in ppm), return the indices
        of the points in the spectrum that correspond to the beginning and
        one-past-the-end of that range."""
        ppm = self.freq_ppm(altDATA)
        npts = len(ppm)

        # Account for the situation in which max or min are out of range
        i_max_ppm = 0
        i_min_ppm = npts

        # N.B. the ppm array goes from high to low
        i_max_ppm = npts - np.searchsorted(ppm[::-1], max_ppm, side='right')
        i_min_ppm = npts - np.searchsorted(ppm[::-1], min_ppm, side='left')
        return (i_max_ppm, i_min_ppm)

    def ppm_points_reverse(self, min_ppm, max_ppm, altDATA=None):
        (i_max_ppm, i_min_ppm) = self.ppm_points(max_ppm, min_ppm, altDATA)
        # Does not work if i_max_ppm is 0, because the it returns -1
        return (i_min_ppm - 1, i_max_ppm - 1)

    def spec_acq_time(self):
        """Returns the total time taken to acquire one spectrum

        i.e. number of scans * (acquisition time + delay between scans)"""
        return self.scans * (self.acq_time + self.last_delay)

    def spec_times(self, nspec=None):
        """Return the time at which the acquisition of each spectrum began"""
        if nspec is None:
            nspec = np.prod(self.actual_npts[1:])
        return np.arange(nspec) * self.spec_acq_time()

    def n_complete_spec(self):
        """The number of spectra where all the scans have been completed

        Sometimes acquisition is stopped in the middle of acquiring a
        spectrum. In this case, not all the scans of the last spectrum have
        been acquired, so the summed intensity will be less. It might be
        desirable to omit the last spectrum in this case."""
        assert (self.actual_npts[2:] == 1).all()  # TODO handle general case
        if self.scans == self.actual_scans:
            num_spectra = self.actual_npts[1]
        else:  # The last scan was not finished, so omit it
            num_spectra = self.actual_npts[1] - 1
        return num_spectra
