import pyedflib
import datetime
import numpy as np
import pickle

def preprocess_edf(filename, nSubjects, channels):

  nChannels = len(channels)

  print "Trying to open BDF file " + filename
  # Load all EDF files in memory.
  edf = pyedflib.EdfReader(filename)

  allData = []
  allSampleFreqs = []
  signalLabels = edf.getSignalLabels()
  sampleFreq = edf.getSampleFrequency(0)
  for k in range(nChannels):
    # Load common meta-info.
    channel    = channels[k]

    # Determine where to start in each file.
    signals = []
    nSamples = float("inf")
    for i in range(nSubjects):
      # Actual channel label is eg. "Resp.1"
      channelLabel = "{:s}.{:d}".format(channel, i+1)
      channelIdx = signalLabels.index(channelLabel)
      signals.append( edf.readSignal(channelIdx) )
      nSamples = min(nSamples, signals[i].shape[0])

    # Copy data into structure.
    data = np.zeros((nSamples, nSubjects+1))
    # First column: timestamps
    data[:,0] = [ t/sampleFreq for t in range(nSamples)]
    # Other columns: aligned data of the channel
    for i in range(nSubjects):
      data[:,i] = signals[i][:nSamples]

    # Append to return structures.
    allData.append( data )

  return allData

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
  parser.add_argument("-f", "--format", type=str, default="npy", choices=["npy", "txt"], help="Output file format")
  parser.add_argument("-p", "--plot", action="store_true", help="Generate plots (as pickled .plt files)")

  args = parser.parse_args()

  # Read config file.
  print "Reading configuration file"
  jsonFile = open(args.config_file)
  config = json.loads( jsonFile.read() )
  channels = config["channels"]

  # Combine EDFs.
  print "Extracting data"
  edfFile = args.input_dir + "/" + config["edf-file"]
  nSubjects = config["n-subjects"]
  data = preprocess_edf(edfFile, nSubjects, channels)

  # Generate output files.
  print "Generating files"
  for i in range(len(channels)):

    output_format = args.format

    basename = args.output_dir + "/data_{:s}_{:s}".format(config["label"], channels[i])
    outputFile = basename + "." + output_format

    print "Saving file: " + outputFile

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
