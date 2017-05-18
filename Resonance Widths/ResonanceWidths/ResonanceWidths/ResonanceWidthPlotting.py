import numpy as np
import matplotlib.pyplot as plt
import sys
import os


args = sys.argv[1:]

filenames = [ f for f in os.listdir() if f.endswith(".out")]

for filename in filenames:
	try:
		plt.plotfile(filename, delimiter=' ', cols=(0, 1), 
	             names=('y from resonance', 'winding number'), marker='o')
		"""plt.show()"""
		plt.savefig("./plots/"+filename.replace(".out",".png"))
		print("Making plots for",filename)
		plt.close()
	except ValueError:
		print("There was nothing in",filename)
