import cluster_phase as cp
import numpy as np

def save_cluster_phase(filename, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp):
  np.savez(filename, meanGrpRho=meanGrpRho, meanIndRho=meanIndRho, meanIndRp=meanIndRp, grpRho=grpRho, indRp=indRp)

def read_cluster_phase(filename):
  data = np.load(filename)
  return data["meanGrpRho"], data["meanIndRho"], data["meanIndRp"], data["grpRho"], data["indRp"]

def seconds_to_sample(seconds, sampleFreq):
  return int(seconds * sampleFreq)

def sample_to_seconds(sample, sampleFreq):
  return sample / sampleFreq

def main():
  import json, argparse

  # Parse commandline arguments.
  parser = argparse.ArgumentParser()
  parser.add_argument("config_file", type=str, help="The configuration file (json)")
  parser.add_argument("-i", "--input-dir", type=str, default=".", help="Path to directory containing the input datafiles")
  parser.add_argument("-o", "--output-dir", type=str, default=".", help="Path to directory where the output datafiles will be stored")
  parser.add_argument("-f", "--format", type=str, default="npy", choices=["npy", "txt"], help="File format (for both input and output)")
  parser.add_argument("-p", "--plot", action="store_true", help="Generate plots (as pickled .plt files)")

  args = parser.parse_args()

  # Read config file.
  print "Reading configuration file"
  jsonFile = open(args.config_file)
  config = json.loads( jsonFile.read() )
  channels = config["channels"]
  nChannels = len(channels)
  nSubjects = config["n-subjects"]
  sampleFreq = config["sample-freq"]
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

  for i in range(len(channels)):
    print "Processing channel '{:s}'".format(channels[i])
    basename = "/data_{:s}_{:s}".format(config["label"], channels[i])
    inputBasename  = args.input_dir  + "/" + basename

    # Run analysis.
    data, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp = cp.cluster_phase(inputBasename + "." + args.format, nSubjects, sampleFreq, firstSample, lastSample)

    if (end == None):
      lastSample = data.shape[0]-1
      end = int(sample_to_seconds(lastSample, sampleFreq))

    outputBasename = args.output_dir + "/" + basename + "_cluster_{:d}-{:d}".format(firstSample, lastSample)

    # Draw plot if needed.
    if (args.plot):
      plotFlag = outputBasename + ".plt"
      title = "Experiment: {:s} Channel: {:s} Rate: {:d} Range: {:d}-{:d}".format(config["label"], channels[i], sampleFreq, start, end)
      cp.generate_plot(data, config["sample-freq"], meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp, plotFlag, title)

    # Save results.
    save_cluster_phase(outputBasename + ".npz", meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp)

if __name__ == "__main__":
  main()
