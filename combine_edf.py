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

def main():

  import json, argparse

  # Parse commandline arguments.
  parser = argparse.ArgumentParser()
  parser.add_argument("config_file", type=str, help="The configuration file (json)")
  parser.add_argument("-i", "--input-dir", type=str, default=".", help="Path to directory containing the input datafiles specified in the configuration file")
  parser.add_argument("-o", "--output-dir", type=str, default=".", help="Path to directory where the output datafiles will be stored")
  parser.add_argument("-f", "--format", type=str, default="npz", choices=["npy", "txt"], help="Output file format")
  parser.add_argument("-p", "--plot", action="store_true", help="Generate plots (as pickled .plt files)")

  args = parser.parse_args()

  # Read config file.
  print "Reading configuration file"
  jsonFile = open(args.config_file)
  config = json.loads( jsonFile.read() )
  channels = config["channels"]

  # Combine EDFs.
  print "Extracting data"
  edfFiles = [ args.input_dir + "/" + f for f in config["edf-files"] ]
  data, sampleFreq = combine_edf(edfFiles, channels, config["offsets"])

  # Generate output files.
  print "Generating files"
  for i in range(len(channels)):

    output_format = args.format

    basename = args.output_dir + "/data_{:s}_{:s}".format(config["label"], channels[i])
    outputFile = basename + "." + output_format

    # Save to appropriate file/format.
    if (output_format == "npy"):
      np.save(outputFile, data[i])
    elif (output_format == "txt"):
      np.savetxt(outputFile, data[i])

    # Generate plot if needed.
    if (args.plot):
      generate_plot(data[i], sampleFreq[i], basename + ".plt")

if __name__ == "__main__":
  main()
