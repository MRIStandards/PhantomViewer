# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 10:58:32 2012

@author: Chris Kerr

This file contains the dtype definitions for reading TNT files into numpy
objects. The definitions are based on the documentation in the file
"A1 - TNMR File Format.doc" distributed with the TNMR software.
"""

import re

from numpy import dtype

Magic_re = re.compile(b"^TNT1\.\d\d\d$")

Magic = dtype('a8')
TLV = dtype([('tag', 'a4'), ('bool', '<u4'), ('length', '<u4')])

TMAG = dtype([
    ('npts', '<i4', 4),
    ('actual_npts', '<i4', 4),
    ('acq_points', '<i4'),
    ('npts_start', '<i4', 4),
    ('scans', '<i4'),
    ('actual_scans', '<i4'),
    ('dummy_scans', '<i4'),
    ('repeat_times', '<i4'),
    ('sadimension', '<i4'),
    ('samode', '<i4'),
    # ('space1', 'a0'),

    ## Result of test on 2012-11-29 on meaning of these numbers:
    ## ob_freq = base_freq + offset_freq / 1000
    ## (actually base_freq ~= 0, but base_freq == ob_freq - offset_freq / 1000)
    ## frequency for 0 ppm = ob_freq * 1e6 + ref_freq
    ('magnet_field', '<f8'),
    ('ob_freq', '<f8', 4),
    ('base_freq', '<f8', 4),
    ('offset_freq', '<f8', 4),
    ('ref_freq', '<f8'),
    ('NMR_frequency', '<f8'),
    ('obs_channel', '<i2'),
    ('space2', 'a42'),

    ('sw', '<f8', 4),
    ('dwell', '<f8', 4),
    ('filter', '<f8'),
    ('experiment_time', '<f8'),
    ('acq_time', '<f8'),
    ('last_delay', '<f8'),
    ('spectrum_direction', '<i2'),
    ('hardware_sideband', '<i2'),
    ('Taps', '<i2'),
    ('Type', '<i2'),
    ('bDigRec', '<u4'),
    ('nDigitalCenter', '<i4'),
    ('space3', 'a16'),

    ('transmitter_gain', '<i2'),
    ('receiver_gain', '<i2'),
    ('NumberOfReceivers', '<i2'),
    ('RG2', '<i2'),
    ('receiver_phase', '<f8'),
    ('space4', 'a4'),

    ('set_spin_rate', '<u2'),
    ('actual_spin_rate', '<u2'),

    ('lock_field', '<i2'),
    ('lock_power', '<i2'),
    ('lock_gain', '<i2'),
    ('lock_phase', '<i2'),
    ('lock_freq_mhz', '<f8'),
    ('lock_ppm', '<f8'),
    ('H2O_freq_ref', '<f8'),
    ('space5', 'a16'),

    ('set_temperature', '<f8'),
    ('actual_temperature', '<f8'),

    ('shim_units', '<f8'),
    ('shims', '<i2', 36),
    ('shim_FWHM', '<f8'),

    ('HH_dcpl_attn', '<i2'),
    ('DF_DN', '<i2'),
    ('F1_tran_mode', '<i2', 7),
    ('dec_BW', '<i2'),
    ('grd_orientation', 'a4'),
    ('LatchLP', '<i4'),
    ('grd_Theta', '<f8'),
    ('grd_Phi', '<f8'),
    ('space6', 'a264'),

    ('start_time', '<u4'),
    ('finish_time', '<u4'),
    ('elapsed_time', '<i4'),

    ('date', 'a32'),
    ('nuclei', 'a16', 4),
    ('sequence', 'a32'),
    ('lock_solvent', 'a16'),
    ('lock_nucleus', 'a16')
    ])


GridAndAxis = dtype([
    ('majorTickInc', '<f8', 12),
    ('minorIntNum', '<i2', 12),
    ('labelPrecision', '<i2', 12),
    ('gaussPerCentimeter', '<f8'),
    ('gridLines', '<i2'),
    ('axisUnits', '<i2'),
    ('showGrid', '<u4'),
    ('showGridLabels', '<u4'),
    ('adjustOnZoom', '<u4'),
    ('showDistanceUnits', '<u4'),
    ('axisName', 'a32'),
    ('space', 'a52'),
])


TMG2 = dtype([
    ('real_flag', '<u4'),
    ('imag_flag', '<u4'),
    ('magn_flag', '<u4'),
    ('axis_visible', '<u4'),
    ('auto_scale', '<u4'),
    ('line_display', '<u4'),
    ('show_shim_units', '<u4'),

    ('integral_display', '<u4'),
    ('fit_display', '<u4'),
    ('show_pivot', '<u4'),
    ('label_peaks', '<u4'),
    ('keep_manual_peaks', '<u4'),
    ('label_peaks_in_units', '<u4'),
    ('integral_dc_average', '<u4'),
    ('integral_show_multiplier', '<u4'),
    ('Boolean_space', '<u4', 9),

    ('all_ffts_done', '<u4', 4),
    ('all_phase_done', '<u4', 4),

    ('amp', '<f8'),
    ('ampbits', '<f8'),
    ('ampCtl', '<f8'),
    ('offset', '<i4'),

    ('axis_set', GridAndAxis),

    ('display_units', '<i2', 4),
    ('ref_point', '<i4', 4),
    ('ref_value', '<f8', 4),
    ('z_start', '<i4'),
    ('z_end', '<i4'),
    ('z_select_start', '<i4'),
    ('z_select_end', '<i4'),
    ('last_zoom_start', '<i4'),
    ('last_zoom_end', '<i4'),
    ('index_2D', '<i4'),
    ('index_3D', '<i4'),
    ('index_4D', '<i4'),

    ('apodization_done', '<i4', 4),
    ('linebrd', '<f8', 4),
    ('gaussbrd', '<f8', 4),
    ('dmbrd', '<f8', 4),
    ('sine_bell_shift', '<f8', 4),
    ('sine_bell_width', '<f8', 4),
    ('sine_bell_skew', '<f8', 4),
    ('Trapz_point_1', '<i4', 4),
    ('Trapz_point_2', '<i4', 4),
    ('Trapz_point_3', '<i4', 4),
    ('Trapz_point_4', '<i4', 4),
    ('trafbrd', '<f8', 4),
    ('echo_center', '<i4', 4),

    ('data_shift_points', '<i4'),
    ('fft_flag', '<i2', 4),
    ('unused', '<f8', 8),
    ('pivot_point', '<i4', 4),
    ('cumm_0_phase', '<f8', 4),
    ('cumm_1_phase', '<f8', 4),
    ('manual_0_phase', '<f8'),
    ('manual_1_phase', '<f8'),
    ('phase_0_value', '<f8'),
    ('phase_1_value', '<f8'),
    ('session_phase_0', '<f8'),
    ('session_phase_1', '<f8'),

    ('max_index', '<i4'),
    ('min_index', '<i4'),
    ('peak_threshold', '<f4'),
    ('peak_noise', '<f4'),
    ('integral_dc_points', '<i2'),
    ('integral_label_type', '<i2'),
    ('integral_scale_factor', '<f4'),
    ('auto_integrate_shoulder', '<i4'),
    ('auto_integrate_noise', '<f8'),
    ('auto_integrate_threshold', '<f8'),
    ('s_n_peak', '<i4'),
    ('s_n_noise_start', '<i4'),
    ('s_n_noise_end', '<i4'),
    ('s_n_calculated', '<f4'),

    ('Spline_point', '<i4', 14),
    ('Spline_point_avr', '<i2'),
    ('Poly_point', '<i4', 8),
    ('Poly_point_avr', '<i2'),
    ('Poly_order', '<i2'),

    ('space', 'a610'),

    ('line_simulation_name', 'a32'),
    ('integral_template_name', 'a32'),
    ('baseline_template_name', 'a32'),
    ('layout_name', 'a32'),
    ('relax_information_name', 'a32'),
    ('username', 'a32'),
    ('user_string_1', 'a16'),
    ('user_string_2', 'a16'),
    ('user_string_3', 'a16'),
    ('user_string_4', 'a16')
])
