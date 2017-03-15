import numpy as np
from scipy.signal import find_peaks_cwt

def extract_data(dataset):

  # load time-series (TS)
  # **************************************************************************
  if (type(dataset) == str):
    if (dataset.endswith(".txt")):
      fulldata = np.loadtxt(dataset)
    else:
      fulldata = np.load(dataset)
  elif (type(dataset) == numpy.darray):
    fulldata = dataset

  return fulldata

# For now this is really stupid...
def align_signals(dataset1,dataset2,sampleFreq,nSamples):

  nTimeSeries = 2
  
  data1 = extract_data(dataset1)
  data2 = extract_data(dataset2)
  
  data = np.zeros((nSamples, nTimeSeries+1))
  data[:,0] = [ t/sampleFreq for t in range(nSamples)]
  data[:,1] = data1[:,1]
  data[:,2] = data2[:,1]
  
  return data
  
def main():

  import json, argparse

  # Parse commandline arguments.
  parser = argparse.ArgumentParser()
  parser.add_argument("input1", type=str, help="First dataset")
  parser.add_argument("input2", type=str, help="Second dataset")
  parser.add_argument("output", type=str, help="Resulting output files")
  parser.add_argument("-f", "--format", type=str, default="npy", choices=["npy", "txt"], help="Output file format")
  parser.add_argument("-s", "--sample-freq", type=float, default=28000.0, help="Sample frequency")

  args = parser.parse_args()

  output_format = args.format
  outputFile = args.output
  
  data = align_signals(args.input1, args.input2, args.sample_freq, 200000)

  # Save to appropriate file/format.
  if (output_format == "npy"):
    np.save(outputFile, data)
  elif (output_format == "txt"):
    np.savetxt(outputFile, data)
    
if __name__ == "__main__":
  main()
