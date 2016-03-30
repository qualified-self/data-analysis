import pyedflib
import datetime
import numpy as np
import pickle

def datetime2timestamp(dt):
  return (dt - datetime.datetime(1970, 1, 1)).total_seconds()

def sampleIndexAtTime(sampleFreq, startTime, timestamp):
  return (timestamp - startTime) * sampleFreq

def startTimestamp(edf, offset):
  return datetime2timestamp(edf.getStartdatetime()) - offset

def combine_edf(filelist, channels, offsets=None):

  nSubjects = len(filelist)
  nChannels = len(channels)

  if (offsets == None):
    offsets = [ 0 for i in range(nSubjects) ]

  # Load all EDF files in memory.
  edf = []
  for i in range(nSubjects):
    edf.append(pyedflib.EdfReader(filelist[i]))

  # Extract common time boundaries.
  maxStartTime = 0
  for i in range(nSubjects):
    startTime = startTimestamp(edf[i], offsets[i])
    maxStartTime = max(maxStartTime, startTime)

  allData = []
  allSampleFreqs = []
  for k in range(nChannels):
    # Load common meta-info.
    channel    = channels[k]
    channelIdx = edf[0].getSignalLabels().index(channel)
    sampleFreq = edf[0].getSampleFrequency(channelIdx)

    # Determine where to start in each file.
    startIdx = []
    signals = []
    nSamples = float("inf")
    for i in range(nSubjects):
      startIdx.append( sampleIndexAtTime(sampleFreq, startTimestamp(edf[i], offsets[i]), maxStartTime) )
      signals.append( edf[i].readSignal(channelIdx)[startIdx[i]:] )
      nSamples = min(nSamples, signals[i].shape[0])

    # Copy data into structure.
    data = np.zeros((nSamples, nSubjects))
    for i in range(nSubjects):
      data[:,i] = signals[i][:nSamples]

    # Append to return structures.
    allData.append( data )
    allSampleFreqs.append( sampleFreq )

  return allData, allSampleFreqs

from matplotlib.pyplot import *

def generate_plot(data, sampleFreq, output):
  nSamples = data.shape[0]
  ax = plot(np.arange(0, nSamples)/sampleFreq, data)
  pickle.dump(ax, file(output, 'w'))

import json
import sys

def main():
  print "Reading JSON file"
  jsonfile = sys.argv[1]
  with open(jsonfile) as f:
    config = json.loads( f.read() )  
  channels = config["channels"]
  
  print "Extracting data"
  data, sampleFreq = combine_edf(config["edf-files"], channels, config["offsets"])

  print "Generating files"
  for i in range(len(channels)):
    basename = "data_{:s}_{:s}".format(config["label"], channels[i])
    generate_plot(data[i], sampleFreq[i], basename + ".plt")
    np.savetxt(basename + ".raw", data[i])
    
  
if __name__ == "__main__":
  main()  
  
  
  
