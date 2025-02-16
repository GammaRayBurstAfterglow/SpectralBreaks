import numpy as np
import matplotlib.pyplot as plt

# Loads the outputted file from TimeElapsed.sh
data = np.loadtxt('elapsed.txt')

# Creates a histogram from all of the data in elapsed.txt
plt.hist(data, edgecolor='black')
plt.xlabel("Time (s)")
plt.ylabel("Frequency")

# saves the graph in the current working directory
plt.savefig('histogram.png')
