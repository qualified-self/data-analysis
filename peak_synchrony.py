import numpy as np
from scipy.signal import find_peaks_cwt
import scipy.stats

from matplotlib.pyplot import *

from detect_peaks import detect_peaks

def extract_peaks(dataset,plotFlag=False):

  nTimeSeries = 2
  
  # load time-series (TS)
  # **************************************************************************
  if (type(dataset) == str):
    if (dataset.endswith(".txt")):
      fulldata = np.loadtxt(dataset)
    else:
      fulldata = np.load(dataset)
  elif (type(dataset) == numpy.darray):
    fulldata = dataset

  # Builds a subset by taking only rows firstSample .. lastSample from base dataset
  #if (lastSample != None):
  #  data = fulldata[firstSample:lastSample,0:nTimeSeries]
  #else:
  #  data = fulldata[firstSample:,0:nTimeSeries]

  data = fulldata[:,0:nTimeSeries+1]
  #plot(data[:,0], data[:,1])
  
  nSamples = data.shape[0]
  
  # Find peaks.
  peak_data = np.zeros((nSamples, nTimeSeries+1))
  peak_data[:,0] = data[:,0]
  width = np.array([50, 75, 100, 150])
  for i in range(1,nTimeSeries+1):
    indexes = detect_peaks(data[:,i], mph=600, mpd=200, edge="rising")
    #print indexes
    #indexes = find_peaks_cwt(data[:,i], width, noise_perc=0.1)
    #plot(data[indexes,0], data[indexes,i], "ro")
    for j in indexes:
      peak_data[j,i] = 1.0
  
  #print peak_data

  #plot(peak_data[:,0], peak_data[:,1])
  
  #show()
  
  return peak_data

# From: https://arxiv.org/pdf/1407.5412.pdf
def synchrony_data(peak_data,plotFlag=False):
  nTimeSeries = 2
  
  a0 = 0.5
  n = 100
  var = n/4
  norm = scipy.stats.norm(0, var)
  b = [ norm.cdf(i) for i in range(-n, n+1) ]
  a = [ 0 for i in range(-n, n+1) ]
  for i in range(1,len(b)):
    a[i] = a0 * (b[i]-b[i-1])
  
  nSamples = peak_data.shape[0]
  
  # Compute sync.
  f_data = np.zeros((nSamples, nTimeSeries))
  i_data = np.zeros((nSamples, nTimeSeries))
  sync_data = np.zeros((nSamples, 2))

  for k in range(nTimeSeries):
    for t in range(n,nSamples-n):
      s = 0
      for j in range(-n, n+1):
        s = s + peak_data[t+j, k] * a[j+n]
      f_data[t,k] = s
      i_data[t,k] = (1-peak_data[t,k]) * 0.5
  
  plot(f_data[:,0], f_data[:,1])
  
  sync_data[:,0] = peak_data[:,0]
  for t in range(nSamples):
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for k in range(nTimeSeries):
      fi = f_data[t, k] * i_data[t, k]
      p  = peak_data[t, k]
      sum1 = sum1 + fi
      sum2 = sum2 + p
      sum3 = sum3 + fi * p
    sync_data[t,1] = sum1*sum2-sum3
    if (sync_data[t,1] > 0):
      print sync_data[t,1]
    #print "Sums: {a} {b} {c}".format(a=sum1, b=sum2, c=sum3)
    
  plot(sync_data[:,0], sync_data[:,1])
  show()
  
  return sync_data

def main():

  import json, argparse

  # Parse commandline arguments.
  parser = argparse.ArgumentParser()
  parser.add_argument("input", type=str, help="Dataset (needs to contain two time series)")
  parser.add_argument("output", type=str, help="Resulting output (synchrony measure)")
  parser.add_argument("-f", "--format", type=str, default="npy", choices=["npy", "txt"], help="Output file format")
  parser.add_argument("-s", "--sample-freq", type=float, default=28000, help="Sample frequency")

  args = parser.parse_args()

  output_format = args.format
  outputFile = args.output
  
  # Extract peaks
  peak_data = extract_peaks(args.input, args.sample_freq)
  
  # Perform computation
  sync_data = synchrony_data(peak_data)
  
  # Save to appropriate file/format.
  #if (output_format == "npy"):
  #  np.save(outputFile, data)
  #elif (output_format == "txt"):
  #  np.savetxt(outputFile, data)
    
if __name__ == "__main__":
  main()



