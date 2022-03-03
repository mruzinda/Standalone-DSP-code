# This script reshapes the one dimensional array from the beamformer and plots the response for analysis
import matplotlib.pyplot as plt
from array import array
import numpy as np

# Open text file containing beamformer output
#f = open("output_d_c.txt", 'r')
#f = open("/home/mruzinda/beamformer_workspace/src/output_d_cuda.txt", 'r')
#f = open("/home/mruzinda/hpguppi_proc/coherent_beamformer/src/output_d_cuda.txt", 'r')
f = open("/home/mruzinda/hpguppi_proc/upchannelizer/src/output_d_cufft.txt", 'r')
#f = open("output_d_c_simple.txt", 'r')
#f = open("output_d_cuda_simple.txt", 'r')

# Read file contents
contents = f.read()

# Split elements based on new line i.e. '\n'
contents_tmp = contents.split('\n')

# Convert contents from string to float
contents_float = np.zeros(len(contents_tmp))
for i in range(0,len(contents_tmp)-1):
    contents_float[i] = float(contents_tmp[i])

# Array dimensions
N_bin = 64*1024
N_time = 8
N_pol = 2
N_ant = 64 
IQ = 2

# Reshape array to 3D -> Antenna X Polarization X Time X Bins
contents_array = contents_float[0:(N_bin*N_time*N_pol*N_ant*IQ)].reshape(N_bin,N_time,N_pol,N_ant,IQ)

ant_idx = 0 # beam index to plot
pol_idx = 0 # time sample index to plot

# Plot intensity map of frequency vs. time
# "interpolation ='none'" removes interpolation which was there by default. 
# I'm only removing it for the sake of accurate analysis and diagnosis.
#plt.imshow(contents_array[0:N_time,0:N_bin,beam_idx], extent=[1, N_bin, 1, N_time], aspect='auto', interpolation='bicubic')
plt.imshow(contents_array[0:N_bin,0:N_time,pol_idx,ant_idx,1], extent=[1, N_time, 1, N_bin], aspect='auto', interpolation='none')
plt.title('Intensity map (Frequency vs. time)')
plt.ylabel('Time samples')
plt.xlabel('Frequency bins')
plt.show()

time_idx = 0

# Plot of power spectrum
plt.plot(contents_array[0:N_bin,time_idx,pol_idx,ant_idx,1])
plt.title('FFT at a time window')
plt.xlabel('Frequency bins')
plt.ylabel('Raw voltage (arb.)')
plt.show()

fig, axs = plt.subplots(2, 2)
fig.suptitle('FFT at different antennas')
axs[0, 0].plot(contents_array[0:N_bin,time_idx,pol_idx,0,1])
axs[0, 0].set_title('Ant 1')
axs[0, 1].plot(contents_array[0:N_bin,time_idx,pol_idx,1,1], 'tab:orange')
axs[0, 1].set_title('Ant 2')
axs[1, 0].plot(contents_array[0:N_bin,time_idx,pol_idx,2,1], 'tab:green')
axs[1, 0].set_title('Ant 3')
axs[1, 1].plot(contents_array[0:N_bin,time_idx,pol_idx,3,1], 'tab:red')
axs[1, 1].set_title('Ant 4')

# set the spacing between subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.85,
                    wspace=0.4, 
                    hspace=0.6)

for ax in axs.flat:
    ax.set(xlabel='Frequency bins', ylabel='Raw voltage')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()

f.close()

