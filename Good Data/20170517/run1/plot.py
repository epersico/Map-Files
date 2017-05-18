import numpy as np
import matplotlib.pyplot as plt
import sys

args = sys.argv[1:]

for filename in args:
	plt.plotfile(filename, delimiter=' ', cols=(0, 1), 
             names=('y from resonance', 'winding number'), marker='o')
	plt.show()
