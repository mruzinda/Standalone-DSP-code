# This script reshapes the one dimensional array from the beamformer and plots the response for analysis
import matplotlib.pyplot as plt
from array import array
import numpy as np

# Open text file containing beamformer output
f = open("output_d.txt", 'r')

# Read file contents
contents = f.read()

# Split elements based on new line i.e. '\n'
contents_tmp = contents.split('\n')

# Convert contents from string to float
contents_float = np.zeros(len(contents_tmp))
for i in range(0,len(contents_tmp)-1):
    contents_float[i] = float(contents_tmp[i])

# Array dimensions
N_beam = 64
N_bin = 10
N_time = 8

# Reshape array to 3D -> Beams X Bins X Time
contents_array = contents_float[0:(N_beam*N_bin*N_time)].reshape(N_beam,N_bin,N_time)

beam_idx = 0 # beam index to plot
time_idx = 0 # time sample index to plot
#contents_array[beam_idx][0:N_bin-1][0:N_time-1]

# Plot intensity map of frequency vs. time
plt.imshow(contents_array[beam_idx][0:N_bin-1][0:N_time-1], extent=[1, N_bin, 1, N_time])
plt.title('Intensity map (Frequency vs. time)')
plt.xlabel('Time samples')
plt.ylabel('Frequency bins')
plt.show()

# Plot of power spectrum
plt.plot(contents_array[beam_idx][0:N_bin-1][time_idx])
plt.title('Power spectrum at a time sample')
plt.xlabel('Frequency bins')
plt.ylabel('Power (arb.)')