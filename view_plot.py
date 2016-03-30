import matplotlib.pyplot as plt
import pickle
import sys

import json, argparse

# Parse commandline arguments.
parser = argparse.ArgumentParser()
parser.add_argument("plot_file", type=str, help="The file containing the plot (saved in pickled format)")
parser.add_argument("-o", "--output", type=str, default=None, help="Optional output file")

args = parser.parse_args()

filename = args.plot_file
pickle.load(file(filename))

if (args.output == None):
  plt.show()
else:
  plt.savefig(args.output)
