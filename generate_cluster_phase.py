import cluster_phase as cp
import numpy as np

def save_cluster_phase(filename, grpRho, indRp):
  data = np.column_stack((grpRho, indRp))
  if (filename.endswith("txt")):
    np.savetxt(filename, data)
  else:
    np.save(filename, data)

def read_cluster_phase(filename):
  if (filename.endswith("txt")):
    data = np.loadtxt(filename)
  else:
    data = np.load(filename)
  return data[:,0], data[:,1:]

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
  nSubjects = len(config["edf-files"])
  rng = config["range"]

  for i in range(len(channels)):
    print "Processing channel '{:s}'".format(channels[i])
    basename = "/data_{:s}_{:s}".format(config["label"], channels[i])
    inputBasename  = args.input_dir  + "/" + basename
    outputBasename = args.output_dir + "/" + basename + "_cluster_{:d}-{:d}".format(rng[0], rng[1])

    # Run analysis.
    data, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp = cp.cluster_phase(inputBasename + "." + args.format, nSubjects, rng[0], rng[1], config["sample_freq"])

    # Draw plot if needed.
    if (args.plot):
      plotFlag = outputBasename + ".plt"
      title = "Experiment: {:s} Channel: {:s} Rate: {:d} Range: {:d}-{:d}".format(config["label"], channels[i], config["sample_freq"], rng[0], rng[1])
      cp.generate_plot(data, config["sample_freq"], meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp, plotFlag, title)

    # Save results.
    save_cluster_phase(outputBasename + "." + args.format, grpRho, indRp)

if __name__ == "__main__":
  main()
