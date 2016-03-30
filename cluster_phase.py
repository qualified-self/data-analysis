#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#   cluster_phase.py
#
#   meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp = cluster_phase(filename, nTimeSeries, firstSample, lastSample, sampleRate, plotFlag)
#
#   Input:
#       filename     : data file
#       nTimeSeries  : number of time series
#       sampleRate   : sample rate of the time series
#       firstSample  : first data point in time series used
#       lastSample   : last data point in time series used
#       plotFlag     : do plots (True, False, or filenme.plt to save using pickle)
#
#   Output:
#       meanGrpRho        : mean group rho (0 to 1; 1 = perfect sync)
#       meanIndRho        : mean rho for each TS to group (0 to 1; 1 = perfect sync)
#       meanIndRp         : mean Relative Phase for each TS to group cluster phase
#       grpRho        : group rho time-series
#       indRp         : relative phase time-series for each individual TS to cluster phase
#
#   Example:
#       [meanGrpRho meanIndRho meanIndRp grpRho indRp] = closter_phase('G201EO1.txt', 6, 1, 7200, 120, True);
#
#   ADAPTED TO PYTHON BY (2016):
#   J. S. Senecal (Concordia University)
#
#   BY (2008):
#   Michael J Richardson (Univeristy of Cincinnati) & Till D. Frank (UCONN)
#
#   UPDATED (2011):
#   Michael J Richardson (Univeristy of Cincinnati)
#
#   References:
#   [1]  Frank, T. D., & Richardson, M. J. (2010). On a test statistic for
#        the Kuramoto order parameter of synchronization: with an illustration
#        for group synchronization during rocking chairs.
#
#   [2]  Richardson,M.J., Garcia, R., Frank, T. D., Gregor, M., &
#        Marsh,K. L. (2010). Measuring Group Synchrony: A Cluster-Phase Method
#        for Analyzing Multivariate Movement Time-Series
#
#   Code Contact & References:
#        michael.richardson@uc.edu
#        http://homepages.uc.edu/~richamo/
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

import numpy, scipy, scipy.signal, matplotlib.pyplot, pickle
from numpy import *
from scipy import *
from matplotlib.pyplot import *

def cluster_phase(filename,nTimeSeries,firstSample,lastSample,sampleRate,plotFlag=False):

  filterfreq = 10

  # load time-series (TS)
  # **************************************************************************
  data = loadtxt(filename)

  # Builds a subset by taking only rows firstSample .. lastSample from base dataset
  ts_data = data[firstSample:lastSample,0:nTimeSeries]

  TSlength = ts_data.shape[0]

  delta_t = 1.0/sampleRate
  t = arange(0, TSlength) * delta_t

  # Downsample, Filter and normalize data
  # **************************************************************************
  # linear detrend data to remove drift (chiar moving slightly during trial
  for nts in range(0,nTimeSeries):
    ts_data[:,nts] =      scipy.signal.detrend(ts_data[:,nts])

  # normalize
  for nts in range(0,nTimeSeries):
    ts_data[:,nts] = scipy.stats.mstats.zscore(ts_data[:,nts], ddof=1)

  # filter
  for nts in range(0,nTimeSeries):
    weight_b,weight_a = scipy.signal.butter(2, filterfreq/(sampleRate/2.0));
    irlength = max(weight_b.size-1,weight_a.size-1)
    ts_data[:,nts] = scipy.signal.filtfilt(weight_b, weight_a, ts_data[:,nts])

  # Compute phase for each TS using Hilbert transform
  # **************************************************************************
  phase = zeros((TSlength-1,nTimeSeries))
  for k in range(0,nTimeSeries):
    hrp = scipy.signal.hilbert(ts_data[:,k])
    for n in range(0,TSlength-1):
      phase[n,k] = arctan2( real(hrp[n]), imag(hrp[n]));
    phase[:,k]=unwrap(phase[:,k]);

  # Compute mean running (Cluster) phase
  # **************************************************************************
  clusterphase = zeros(TSlength-1)
  for n in range(0,TSlength-1):
    ztot = complex(0,0);
    for k in range(0,nTimeSeries):
        z = exp(1j * phase[n,k]);
        ztot = ztot+z;
    ztot = ztot/nTimeSeries;
    clusterphase[n] = angle(ztot);
  clusterphase = unwrap(clusterphase);

  # Compute relative phases between phase of TS and cluster phase
  # **************************************************************************
  complexIndRp = zeros((TSlength-1,nTimeSeries),dtype=complex);
  indRp = zeros((TSlength-1,nTimeSeries))
  meanIndRp = zeros(nTimeSeries);
  meanIndRho = zeros(nTimeSeries);
  for k in range(0,nTimeSeries):
    ztot = complex(0,0);
    for n in range(0,TSlength-1):
        z = exp(1j * (phase[n,k] - clusterphase[n]));
        complexIndRp[n,k] = z;
        ztot = ztot+z;
    indRp[:,k] = angle(complexIndRp[:,k]) * 360/(2*numpy.pi); # convert radian to degrees
    ztot = ztot / (TSlength-1);
    meanIndRp[k] = angle(ztot);
    meanIndRho[k] = abs(ztot);
  meanRp = meanIndRp;
  meanIndRp = (meanIndRp / (2*numpy.pi)*360); # convert radian to degrees
  print(' ');
  print('Mean relative phases of individuals to cluster phase')
  print(meanIndRp);
  print('Averaged degree of synchronization of individuals (Rho = 1-circular variance)')
  print(meanIndRho);

  # Compute cluster amplitude rhotot in rotation frame
  # **************************************************************************
  grpRho=zeros(TSlength-1);
  for n in range(0,TSlength-1):
    ztot = complex(0,0);
    for k in range(0,nTimeSeries):
        z = exp(1j * (phase[n,k] - clusterphase[n] - meanRp[k]));
        ztot = ztot+z;
    ztot = ztot / nTimeSeries;
    grpRho[n] = abs(ztot);
  print grpRho

  meanGrpRho = mean(grpRho);
  print('Averaged degree of synchronization of the group')
  print(meanGrpRho);

  # Do Plot
  # **************************************************************************
  # plot data for time-series (separeted on graph for display purposes)
  #scrsz = get(0,'ScreenSize');
  #h = figure('Position',[scrsz(3)/3 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);
  if plotFlag != False:

    fig = figure(1)

    subplot(3,1,1);
    tmpdata = zeros((TSlength,nTimeSeries));
    for nts in range(0,nTimeSeries):
      tmpdata[:,nts] = (ts_data[:,nts] + (nts*4));
    plot(t, tmpdata)
    xlabel('Time');
    ylabel('RAW Data');
    xlim([0, max(t)]);
    #ylim([-185, 185]);

    # plot individ-cluster relative phase
    subplot(3,1,2);
    plot(t[0:TSlength-1], indRp);
    xlabel('Time');
    ylabel('IND-Clust Relative Phase');
    xlim([0, max(t)]);
    ylim([-185, 185]);

    # plot group-cluster amplitiude (rho) timeseries
    subplot(3,1,3);
    plot(t[0:TSlength-1], grpRho)
    xlabel('Time');
    ylabel('GRP-Clust Amplitude');
    xlim([0, max(t)]);
    ylim([0, 1]);

    text(0, -.4, "Mean GRP Rho: {:.3f}  Mean IND Rhos: {:s}".format(meanGrpRho, array_str(meanIndRho,precision=3)))

    if plotFlag == True:
      fig.show()
    else:
      pickle.dump(fig, file(plotFlag, 'w'))

  return meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp

import json

def main():
  print "Reading JSON file"
  jsonfile = sys.argv[1]
  with open(jsonfile) as f:
    config = json.loads( f.read() )
  channels = config["channels"]
  nSubjects = len(config["edf-files"])

  for i in range(len(channels)):
    print "Processing channel '{:s}'".format(channels[i])
    rng = config["range"]
    basename = "data_{:s}_{:s}".format(config["label"], channels[i])
    output = basename + "_cluster_{:d}-{:d}.plt".format(rng[0], rng[1])
    nChannels = len(channels)
    meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp = cluster_phase(basename + ".raw", nSubjects, rng[0], rng[1], config["sample_freq"], False)

if __name__ == "__main__":
  main()
