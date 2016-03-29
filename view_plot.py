import matplotlib.pyplot as plt
import pickle
import sys

filename = sys.argv[1]
pickle.load(file(filename))
plt.show()

