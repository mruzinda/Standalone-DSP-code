# This script reshapes the one dimensional array from the beamformer and plots the response for analysis
import matplotlib.pyplot as plt
from array import array
import numpy as np

# Open text file containing beamformer output
filename = "/datag/users/mruzinda/o/output_d_cufft.bin"

# Read file contents: np.fromfile(filename, dtype=float, count=- 1, sep='', offset=0)
contents_float = np.fromfile(filename, dtype=np.float32)

# Read file contents
#contents = bytearray(f.read(4*2*131072))

print(len(contents_float))
print(contents_float[0])

# Array dimensions
N_fine = 1024 # 2^17 
N_coarse = 64 # 4
N_win = 8
N_pol = 2
N_ant = 64 
N_iq = 2

# Reshape array to 3D -> Antenna X Polarization X Time X Bins
contents_array = contents_float[0:(N_coarse*N_win*N_pol*N_ant*N_fine*N_iq)].reshape(N_coarse,N_win,N_pol,N_ant,N_fine,N_iq)

contents_restructure = np.zeros(N_win*N_pol*N_ant*N_coarse*N_fine*N_iq).reshape(N_win,N_pol,N_ant,N_coarse*N_fine,N_iq)

# Combine coarse and fine channels
for c in range(0,N_coarse):
    for w in range(0,N_win):
        for p in range(0,N_pol):
            for a in range(0,N_ant):
                for iq in range(0,N_iq):
                    contents_restructure[w,p,a,(0+c*N_fine):(N_fine+c*N_fine),iq] = contents_array[c,w,p,a,0:N_fine,iq]

ant_idx = 0 # beam index to plot
pol_idx = 0 # polarization index to plot
time_idx = 0 # time sample index to plot
coarse_idx = 0 # coarse channel index to plot
iq_idx = 0 # Real or imaginary component
# Plot intensity map of frequency vs. time
# "interpolation ='none'" removes interpolation which was there by default. 
# I'm only removing it for the sake of accurate analysis and diagnosis.
#plt.imshow(contents_array[0:N_win,0:N_coarse,beam_idx], extent=[1, N_coarse, 1, N_win], aspect='auto', interpolation='bicubic')
#plt.imshow(contents_array[coarse_idx,0:N_win,pol_idx,ant_idx,0:N_fine,iq_idx], extent=[1, N_fine, 1, N_win], aspect='auto', interpolation='none')
plt.imshow(contents_restructure[0:N_win,pol_idx,ant_idx,0:N_coarse*N_fine,iq_idx], extent=[1, N_coarse*N_fine, 1, N_win], aspect='auto', interpolation='none')
plt.title('Intensity map (Frequency vs. time)')
plt.ylabel('Time samples')
plt.xlabel('Frequency bins')
plt.show()

# Plot of power spectrum
#plt.plot(contents_array[coarse_idx,time_idx,pol_idx,ant_idx,0:N_fine,iq_idx])
plt.plot(contents_restructure[time_idx,pol_idx,ant_idx,0:N_coarse*N_fine,iq_idx])
plt.title('FFT at a time window')
plt.xlabel('Fine frequency channels')
plt.ylabel('Raw voltage (arb.)')
plt.show()

fig, axs = plt.subplots(2, 2)
fig.suptitle('FFT at different antennas')
axs[0, 0].plot(contents_restructure[time_idx,pol_idx,0,0:N_coarse*N_fine,iq_idx])
axs[0, 0].set_title('Ant 1')
axs[0, 1].plot(contents_restructure[time_idx,pol_idx,1,0:N_coarse*N_fine,iq_idx], 'tab:orange')
axs[0, 1].set_title('Ant 2')
axs[1, 0].plot(contents_restructure[time_idx,pol_idx,2,0:N_coarse*N_fine,iq_idx], 'tab:green')
axs[1, 0].set_title('Ant 3')
axs[1, 1].plot(contents_restructure[time_idx,pol_idx,57,0:N_coarse*N_fine,iq_idx], 'tab:red')
axs[1, 1].set_title('Ant 57')

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

