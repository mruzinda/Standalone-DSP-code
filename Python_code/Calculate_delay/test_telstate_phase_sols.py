import matplotlib.pyplot as plt
from array import array
import numpy as np
import sys
import pickle

with open('test_output.pkl', 'rb') as f:
    p = pickle.load(f)

print(len(p['m037h']))
print(p['m037h'][0])
print(np.angle(p['m037h'][0])) # Phase of the first channel of the 37th antenna in the horizontal polarization

shw = plt.plot(np.angle(p['m037h']))
plt.title('Phase of antenna m037h')
plt.ylabel('Phase')
plt.xlabel('Frequency bins')
plt.show()
