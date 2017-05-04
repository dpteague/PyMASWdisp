"""
    This script calculates experimental dispersion data for one or more 
    source-offsets. Input file(s) for each source offset should be in seg2
    (.dat) format. 

    Note that the data acquisition paramters (e.g., receiver spacing, source-
    offset, etc) are automatically read from the seg2 file. Thus, if the data
    acquisition parameters were not properly entered into the data acquisition
    system, then these values of the cShotGather class should be manually 
    changed. For example, if the offset was not properly entered into the data 
    acquisition, then you can create a variable called 'offsets' and manually 
    change the offset during each loop by entering the command 
    'cShotGather.offset = offsets[k]' after the cShotGather class is created.
    (See 'shotgathers.py' for more information). 


    This code was developed at the University of Texas at Austin.
    Copyright (C) 2016  David P. Teague, Clinton M. Wood, and Brady R. Cox 
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
"""

# INPUTS************************************************************************
#*******************************************************************************

# Processing inputs.............................................................

# Input file numbers (e.g., [ [1,10], [11, 20], [21,30] ] , where shots
# 1-10 correspond to one offset, 11-20 correspond to another offset, etc.)
infile_numbers = [[11,20],[11,20],[11,20],[11,20],[11,20],[11,20]]

# Directory containing input seg2 (.dat) files. (Enter pwd for current directory)
infile_path = 'C:/DPT_Projects/New Zealand/Preliminary Analysis/Fitzgerald/Active/Raw Data'

# Files to exclude from calculations (enter [] to include all files)
exclude_files=[]

# Name output compressed pickle (.pklz) file
outfile_name = 'FTG_DC_raw'

# Path of output compressed pickle (.pklz) file
outfile_path = infile_path

# Length of recorded trace to use in calculations [s]
timelength = 0.5

# Desired spacing of dispersion data points [Hz]
deltaf = 0.5

# Min and max frequency to consider in calculations
min_freq = 5
max_freq = 100

# Min velocity (FDBF, phase-shift, and tau-p) and max velocity (FDBF and phase-shift) to consider in calculations 
min_vel = 80
max_vel = 400

# Number of velocities or wavenumbers to consider in calculations
n_trial = 2048

# Processing method (can specify single method or one for each source-offset)
#   'fk':           Frequency-wavenumber transformation
#   'fdbf':         Frequency domain beamformer
#   'phase-shift':  Phase-shift transformation
#   'slant-stack':  Slant-stack (linear Radon) transformation
processMethod = [ "fk", "fdbf", "fdbf", "fdbf", "phase-shift", "slant-stack"]

# Weighting technique for FDBF (similar to processing method, can specify multiple)
#   'none' for equal weighting
#   'sqrt' to weight each receiver trace by the square root of the distance from the source
#   'invamp' to weight each receiver by 1/|A(f,x)|
weightType = [ '', 'none', 'sqrt', 'invamp', '', '' ]

# Signal-to-noise ratio inputs
noise_location = 'end'



# Plotting inputs (1 to plot, otherwise 0)......................................
flag_SNR = 1
# Contour plots (1 to plot, otherwise 0)
flag_con_fk = 0                    # Frequency-wavenumber
flag_con_fw = 0                    # Frequency-wavelength
flag_con_fv = 1                    # Frequency-phase velocity
flag_con_fp = 0                    # Frequency-slowness
flag_con_wv = 0                    # Wavelength-phase velocity
# Slices in various domains (1 to plot, otherwise 0)
f_plot_vals = range(8,40,2)        # Frequencies to plot slices in chosen domain
flag_slices_fk = 0                 # Wavenumber-power at select frequencies
flag_slices_fw = 0                 # Wavelength-power at select frequencies
flag_slices_fv = 0                 # Velocity-power at select frequencies
flag_slices_fp = 0                 # Slowness-power at select frequencies

# END OF INPUTS*****************************************************************
#*******************************************************************************


# Load modules
import time
import numpy as np
import pickle
import gzip
import shotgathers
import dcprocessing
import dctypes


# Initialize lists for storing processed dispersion data
offsetRaw = []
frequencyRaw = []
velocityRaw = []
wavelengthRaw = []

# Time the processing
tic = time.clock()

# Compare signal-to-noise ratios for various offsets
if flag_SNR:
    shotgathers.compareSNR( infile_numbers, infile_path, noise_location, timelength, exclude_files, deltaf )

# Loop through all offsets
for k in range( len(infile_numbers) ):

    # Processing***************************************************************
    # Stack records for current offset location
    cShotGather = shotgathers.importAndStackWaveforms( infile_numbers[k][0], infile_numbers[k][1], infile_path, exclude_files )
    # Cut time records at timelength
    cShotGather.cut(timelength)
    # Padd zeros to acheive desired df
    cShotGather.zero_pad(deltaf)
    
    # User-defined processing method for current offset
    if len(processMethod)==1:
        cpMethod = processMethod[0]
    else:
        cpMethod = processMethod[k]

    # Process dispersion data based on user-defined processing method
    if str.lower(cpMethod)=='fdbf':
        # User-defined weighting technique for fdbf
        if len(weightType)==1:
            cWtType = weightType[0]
        else:
            cWtType = weightType[k]
        cDCpower = dcprocessing.fdbf( cShotGather, cWtType, n_trial, min_vel, max_vel, min_freq, max_freq )
    elif str.lower(cpMethod)=='fk':
        cDCpower = dcprocessing.fk( cShotGather, n_trial, min_freq, max_freq )
    elif str.lower(cpMethod)=='phase-shift':
        cDCpower = dcprocessing.phase_shift( cShotGather, n_trial, min_freq, max_freq, min_vel, max_vel )
    elif str.lower(cpMethod)=='slant-stack':
        cDCpower = dcprocessing.tau_p( cShotGather, n_trial, min_freq, max_freq, min_vel, max_vel )
    else:
        raise ValueError('Invalid processing method')

    # Store dispersion data in lists
    offsetRaw.append( cShotGather.offset )
    frequencyRaw.append( cDCpower.freq )
    if str.lower(cDCpower.val_type)=='wavenumber':
        velocityRaw.append( 2*np.pi*cDCpower.freq / cDCpower.peak_vals )
        wavelengthRaw.append( 2*np.pi / cDCpower.peak_vals )
    elif str.lower(cDCpower.val_type)=='velocity': 
        velocityRaw.append( cDCpower.peak_vals ) 
        wavelengthRaw.append( cDCpower.peak_vals / cDCpower.freq )
    else:
        raise ValueError('Invalid val_type. Should be \"wavenumber\" or \"velocity\".')

    # Plotting******************************************************************
    # Contour plots
    if flag_con_fk:
        cDCpower.plotSpect("fk", [0, max_freq, 0, cDCpower.kres ])
    if flag_con_fw:
        cDCpower.plotSpect("fw", [0, max_freq, 0.5*cShotGather.position[1], 2*cShotGather.position[-1] ] )
    if flag_con_fv:
        cDCpower.plotSpect("fv", [0, max_freq, 0, max_vel ])
    if flag_con_fp:
        cDCpower.plotSpect("fp", [0, max_freq, 1.0/max_vel, 1.0/min_vel ])
    if flag_con_wv:
        cDCpower.plotSpect("wv", [0.5*cShotGather.position[1], 2*cShotGather.position[-1], 0, max_vel ])
    # Slices in various domains
    if flag_slices_fk:
        cDCpower.plotSlices("fk", f_plot_vals, (0, cDCpower.kres) )
    if flag_slices_fw:
        cDCpower.plotSlices("fw", f_plot_vals, ( 0.5*cShotGather.position[1], 2*cShotGather.position[-1] ) )
    if flag_slices_fv:
        cDCpower.plotSlices("fv", f_plot_vals, (0, max_vel) )
    if flag_slices_fp:
        cDCpower.plotSlices("fp", f_plot_vals, (1.0/max_vel, 1.0/min_vel) )


# Create class containing all raw dispersion data
cRawDC = dctypes.RawDispersion( frequencyRaw, velocityRaw, offsetRaw )
# Automatically cut data with excessively high velocities
cRawDC.rmvHighVs(max_vel)
# Save results to a compressed pickle file
f = gzip.open(outfile_path+"/"+outfile_name+".pklz", 'wb')
pickle.dump(cRawDC, f)
f.close()


# Print execution time
toc = time.clock()
print "Elapsed time: "+str(toc - tic)+" seconds"


