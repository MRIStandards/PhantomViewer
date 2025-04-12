# -*- coding: utf-8 -*-

import numpy as np

import TNTdtypes


def make_str(b):
    """Convert a bytes object to a str, decoding with latin1 if necessary"""
    if isinstance(b, str):  # Python 2
        return b
    elif isinstance(b, bytes):  # Python 3
        return b.decode('latin1')
    else:
        return b


def unsqueeze(M, new_ndim=4):
    """Add extra dimensions to a matrix so it has the desired dimensionality"""
    newshape = np.ones((new_ndim,))
    newshape[:M.ndim] = M.shape
    return np.reshape(M, newshape, order='A')


def convert_si(si_num_list):
    """takes a list of strings, si_num_lst, of the form xxx.xxx<si suffix>
    and returns an array of floats xxx.xxxe-6 or similar, depending on the si
    suffix
    """
    # A dict of all the SI prefixes
    prefix = {'y': 1e-24,  # yocto
              'z': 1e-21,  # zepto
              'a': 1e-18,  # atto
              'f': 1e-15,  # femto
              'p': 1e-12,  # pico
              'n': 1e-9,   # nano
              'u': 1e-6,   # micro
              'm': 1e-3,   # mili
              'c': 1e-2,   # centi
              'd': 1e-1,   # deci
              's': 1,      # seconds
              'k': 1e3,    # kilo
              'M': 1e6,    # mega
              'G': 1e9,    # giga
              'T': 1e12,   # tera
              'P': 1e15,   # peta
              'E': 1e18,   # exa
              'Z': 1e21,   # zetta
              'Y': 1e24,   # yotta
              }
    # go through the items in the list and try to float them. If they don't
    # float...
    for index, item in enumerate(si_num_list):
        try:
            si_num_list[index] = float(item)
        except ValueError:
            # check if the last character is a valid prefix
            if item[-1] in prefix:
                # if it is, set that list element to the rest of the item
                # multiplied by the value appropriate to the prefix
                si_num_list[index] = prefix[item[-1]] * float(item[:-1])
            else:
                # raise if it doesn't work out
                raise ValueError("Couldn't convert delay table entires\
                                 to float! Make sure your suffixes\
                                 correspond to real SI units.")
    #return it as an array
    return np.array(si_num_list)


def read_pascal_string(data, number_type='<i4', encoding='ascii'):
    number_type = np.dtype(number_type)
    length = np.fromstring(data, dtype=number_type, count=1)[0]
    number_size = number_type.itemsize
    text = data[number_size:number_size + length]
    if not isinstance(text, str):
        text = str(text, encoding)
    return text


def save_gnuplot_matrix(tnt, mat_file, max_ppm=np.Inf, min_ppm=-np.Inf,
                        altDATA=None, times=None, logfile=None):
    """Save a file suitable for use as a gnuplot 'binary matrix'

    Only the real part is saved, and it is converted to 32 bit float.
    The frequency goes in the first row, and the acquisition time goes in
    the first column.

    See http://gnuplot.sourceforge.net/docs_4.2/node330.html for a
    description of the data format."""
    ppm = tnt.freq_ppm(altDATA)
    (i_max_ppm, i_min_ppm) = tnt.ppm_points(max_ppm, min_ppm, altDATA)

    ppm = ppm[i_max_ppm:i_min_ppm]
    if altDATA is None:
        DATAslice = tnt.DATA[i_max_ppm:i_min_ppm, :]
        nspec = tnt.n_complete_spec()
    else:
        DATAslice = altDATA[i_max_ppm:i_min_ppm, :]
        nspec = altDATA.shape[1]

    npts = DATAslice.shape[0]

    gpt_matrix = np.memmap(mat_file, dtype='f4', mode='w+',
                           shape=(npts + 1, nspec + 1), order='F')

    gpt_matrix[0, 0] = npts
    gpt_matrix[1:, 0] = ppm

    if times is None:
        times = tnt.spec_times(nspec)

    for i in range(nspec):
        gpt_matrix[0, i+1] = times[i]
        ## without the 'squeeze', we get some kind of 'output operand
        ## requires a reduction, but reduction is not enabled' error ??
        gpt_matrix[1:, i+1] = DATAslice.real[:, i].squeeze()
        if logfile is not None:
            logfile.write('.')
            logfile.flush()
    if logfile is not None:
        logfile.write('Done\n')
        logfile.flush()

    del(gpt_matrix)  # flush the file to disk


def dump_params_txt(tnt, txtfile):
    """Write a text file with the acquisition and processing parameters"""
    if type(txtfile) == str:
        txtfile = open(txtfile, 'w')

    txtfile.write("TMAG struct (acquisition parameters):\n")
    for fieldname in TNTdtypes.TMAG.names:
        if fieldname.startswith('space'):
            continue
        txtfile.write("{0}:\t{1}\n".format(fieldname, s(tnt.TMAG[fieldname])))

    txtfile.write("\nTMG2 struct (processing parameters):\n")
    for fieldname in TNTdtypes.TMG2.names:
        if fieldname in ['Boolean_space', 'unused', 'space']:
            continue
        txtfile.write("{0}:\t{1}\n".format(fieldname, s(tnt.TMG2[fieldname])))
