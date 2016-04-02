import cluster_phase as cp
import numpy as np
from matplotlib.pyplot import *

# Generate plot based on data compiled by cluster_phase() (used in cluster_phase() with option plotFlag)
def generate_plot(data, sampleRate, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp, plotFlag=True, plotTitle=None, plotMarkers=None):
  fig = figure(1)
  dataLength = indRp.shape[0]
  nTimeSeries = indRp.shape[1]
  t = np.arange(0, dataLength) / sampleRate

  subplot(3,1,1);
  tmpdata = np.zeros((dataLength,nTimeSeries));
  for nts in range(0,nTimeSeries):
    tmpdata[:,nts] = (data[:,nts+1] + (nts*4));
  plot(t, tmpdata, alpha=0.7)
  xlabel('Time');
  ylabel('RAW Data');
  xlim([0, max(t)]);
  #ylim([-185, 185]);
  ymin, ymax = ylim()
  #print ymin, ymax
  vlines(plotMarkers, ymin, ymax, linewidth=2, linestyles='dashed', color='k')

  # plot individ-cluster relative phase
  subplot(3,1,2);
  plot(t, indRp, alpha=0.7);
  xlabel('Time');
  ylabel('IND-Clust Relative Phase');
  xlim([0, max(t)]);
  ylim([-185, 185]);
  ymin, ymax = ylim()
  #print ymin, ymax
  vlines(plotMarkers, ymin, ymax, linewidth=2, linestyles='dashed', color='k')

  # plot group-cluster amplitiude (rho) timeseries
  subplot(3,1,3);
  plot(t, grpRho, alpha=0.7)
  xlabel('Time');
  ylabel('GRP-Clust Amplitude');
  xlim([0, max(t)]);
  ylim([0, 1]);
  ymin, ymax = ylim()
  #print ymin, ymax
  vlines(plotMarkers, ymin, ymax, linewidth=2, linestyles='dashed', color='k')

  # add text
  if (plotTitle != None):
    fig.suptitle(plotTitle)

  text(0, -.4, "Mean GRP Rho: {:.3f}  Mean IND Rhos: {:s}".format(meanGrpRho, np.array_str(meanIndRho,precision=3)))

  # save or display plot
  if plotFlag == True:
    print "Showing graphics"
    show()
  elif plotFlag.endswith(".plt"):
    print "Dumping to " + plotFlag
    pickle.dump(fig, file(plotFlag, 'w'))
  else:
    fig.savefig(plotFlag)

  close(fig)

def save_cluster_phase(filename, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp):
  np.savez(filename, meanGrpRho=meanGrpRho, meanIndRho=meanIndRho, meanIndRp=meanIndRp, grpRho=grpRho, indRp=indRp)

def read_cluster_phase(filename):
  data = np.load(filename)
  return np.asscalar(data["meanGrpRho"]), data["meanIndRho"], data["meanIndRp"], data["grpRho"], data["indRp"]

def read_data(filename):
  if (filename.endswith("txt")):
    data = np.loadtxt(filename)
  else:
    data = np.load(filename)
  return data

def seconds_to_sample(seconds, sampleFreq):
  return int(seconds * sampleFreq)

def sample_to_seconds(sample, sampleFreq):
  return sample / sampleFreq

def get_dataset_n_samples(filename):
  if (filename.endswith(".txt")):
    fulldata = np.loadtxt(filename)
  else:
    fulldata = np.load(filename)
  return fulldata.shape[0]

def main():
  import json, argparse

  # Parse commandline arguments.
  parser = argparse.ArgumentParser()
  parser.add_argument("config_file", type=str, help="The configuration file (json)")
  parser.add_argument("-c", "--channels", type=str, default=None, help="The name(s) of the channel separated by commas")
  parser.add_argument("-i", "--input-dir", type=str, default=".", help="Path to directory containing the input datafiles")
  parser.add_argument("-o", "--output-dir", type=str, default=".", help="Path to directory where the output datafiles will be stored")
  parser.add_argument("-f", "--format", type=str, default="npy", choices=["npy", "txt"], help="File format (for input)")
  parser.add_argument("-P", "--plot-format", type=str, default="view", choices=["pdf", "png", "plt", "view"], help="What to do with the plot (output type)")

  args = parser.parse_args()

  # Read config file.
  print "Reading configuration file"
  jsonFile = open(args.config_file)
  config = json.loads( jsonFile.read() )
  if (args.channels == None):
    channels = config["channels"]
  else:
    channels = args.channels.split(",")
  nSubjects = config["n-subjects"]
  sampleFreq  = config["sample-freq"]
  markers = config["markers"]
  if ("start" in config):
    start = config["start"]
  else:
    start = 0
  firstSample = seconds_to_sample(start, sampleFreq)
  if ("end" in config):
    end = config["end"]
    lastSample  = seconds_to_sample(end, sampleFreq)
  else:
    end = None
    lastSample = None

  for channel in channels:

    print "Processing channel '{:s}'".format(channel)
    basename = "/data_{:s}_{:s}".format(config["label"], channel)
    inputBasename  = args.input_dir  + "/" + basename

    if (end == None):
      lastSample = get_dataset_n_samples(inputBasename + "." + args.format) - 1
      end = sample_to_seconds(lastSample, sampleFreq)

    outputBasename = args.output_dir + "/" + basename + "_cluster_{:d}-{:d}".format(firstSample, lastSample)

    # Extract data.
    fulldata = read_data(inputBasename + "." + args.format)
    meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp = read_cluster_phase(outputBasename + ".npz")
    # Builds a subset by taking only rows firstSample .. lastSample from base dataset
    if (lastSample != None):
      data = fulldata[firstSample:lastSample,]
    else:
      data = fulldata[firstSample:,]

    # Draw plot if needed.
    if (args.plot_format == "view"):
      plotFlag = True
    else:
      plotFlag = outputBasename + "." + args.plot_format
    title = "Experiment: {:s} Channel: {:s} Rate: {:d} Range: {:d}s-{:d}s".format(config["label"], channel, sampleFreq, start, end)
    generate_plot(data, config["sample-freq"], meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp, plotFlag, title, markers)
    print "Done"

if __name__ == "__main__":
  main()
