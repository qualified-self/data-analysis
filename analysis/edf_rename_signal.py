## IMPORTANT: This file is still in development
## IT DOES NOT WORK!!! SO DO NOT USE IT!!!

import pyedflib, copy
import numpy as np

def main():

  import json, argparse

  # Parse commandline arguments.
  parser = argparse.ArgumentParser()
  parser.add_argument("signal_operations", type=str, help="Operations to perform over signals")
  parser.add_argument("files", type=str, nargs='+', help="List of EDF files to rename")

  args = parser.parse_args()

  operations = {}
  for op in args.signal_operations.split(','):
    fromSignal, toSignal = op.split(':')
    if (toSignal == ""):
      toSignal = fromSignal
    operations[fromSignal] = toSignal

  for f in args.files:

    # Create input file (for reading).
    reader  = pyedflib.EdfReader(f)
    header = reader.getHeader()
    labels  = reader.getSignalLabels()

    print reader.getFileDuration()

    # Build channel info for output.
    channelInfo = []
    inputChannelIdx = []
    print operations
    for fromSignal, toSignal in operations.iteritems():
      chn = labels.index(fromSignal)
      info = copy.deepcopy( reader.getSignalHeader(chn) )
      info['label'] = toSignal
      print info
      channelInfo.append( info )
      inputChannelIdx.append( chn )

    # Create output file.
    print len(channelInfo)
    outf = "/tmp/" + f + ".out.bdf"
    writer = pyedflib.EdfWriter(outf, len(channelInfo), pyedflib.FILETYPE_BDFPLUS)
    print outf
    writer.setHeader(header)
    writer.setSignalHeaders(channelInfo)
    #print "TEst " + str(pyedflib.set_datarecord_duration(writer.handle, reader.getFileDuration()*1000.0))
    #writer.setDatarecordDuration(1) # pyedflib.EDFLIB_TIME_DIMENSION = 10000000.0
    #print writer.duration

    # Copy signals.
    for i in range(len(channelInfo)):
      duration = reader.getFileDuration()
      for t in range(duration):
        writer.writePhysicalSamples( reader.readSignal(inputChannelIdx[i]) )

    writer.close()

if __name__ == "__main__":
  main()
