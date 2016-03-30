#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#   cluster_phase.py
#
#   data, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp = cluster_phase(dataset, nTimeSeries, firstSample, lastSample, sampleRate, plotFlag)
#
#   Input:
#       dataset      : dataset as either a filename OR directly a numpy array
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

# Generate plot based on data compiled by cluster_phase() (used in cluster_phase() with option plotFlag)
def generate_plot(data, sampleRate, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp, plotFlag=True,plotTitle=None):
  fig = figure(1)
  dataLength = data.shape[0]
  nTimeSeries = data.shape[1]
  t = arange(0, dataLength) / sampleRate

  subplot(3,1,1);
  tmpdata = zeros((dataLength,nTimeSeries));
  for nts in range(0,nTimeSeries):
    tmpdata[:,nts] = (data[:,nts] + (nts*4));
  plot(t, tmpdata)
  xlabel('Time');
  ylabel('RAW Data');
  xlim([0, max(t)]);
  #ylim([-185, 185]);

  # plot individ-cluster relative phase
  subplot(3,1,2);
  plot(t[0:dataLength-1], indRp);
  xlabel('Time');
  ylabel('IND-Clust Relative Phase');
  xlim([0, max(t)]);
  ylim([-185, 185]);

  # plot group-cluster amplitiude (rho) timeseries
  subplot(3,1,3);
  plot(t[0:dataLength-1], grpRho)
  xlabel('Time');
  ylabel('GRP-Clust Amplitude');
  xlim([0, max(t)]);
  ylim([0, 1]);

  # add text
  if (plotTitle != None):
    fig.suptitle(plotTitle)

  text(0, -.4, "Mean GRP Rho: {:.3f}  Mean IND Rhos: {:s}".format(meanGrpRho, array_str(meanIndRho,precision=3)))

  # save or display plot
  if plotFlag == True:
    fig.show()
  else:
    print "Dumping to " + plotFlag
    pickle.dump(fig, file(plotFlag, 'w'))
  close(fig)


# Compiles cluster phase.
def cluster_phase(dataset,nTimeSeries,firstSample,lastSample,sampleRate,plotFlag=False):

  filterfreq = 10

  # load time-series (TS)
  # **************************************************************************
  if (type(dataset) == str):
    if (dataset.endswith(".txt")):
      fulldata = loadtxt(dataset)
    else:
      fulldata = load(dataset)
  elif (type(dataset) == numpy.darray):
    fulldata = dataset

  # Builds a subset by taking only rows firstSample .. lastSample from base dataset
  data = fulldata[firstSample:lastSample,0:nTimeSeries]

  dataLength = data.shape[0]

  delta_t = 1.0/sampleRate
  t = arange(0, dataLength) * delta_t

  # Downsample, Filter and normalize data
  # **************************************************************************
  # linear detrend data to remove drift (chiar moving slightly during trial
  for nts in range(0,nTimeSeries):
    data[:,nts] =      scipy.signal.detrend(data[:,nts])

  # normalize
  for nts in range(0,nTimeSeries):
    data[:,nts] = scipy.stats.mstats.zscore(data[:,nts], ddof=1)

  # filter
  for nts in range(0,nTimeSeries):
    weight_b,weight_a = scipy.signal.butter(2, filterfreq/(sampleRate/2.0));
    irlength = max(weight_b.size-1,weight_a.size-1)
    data[:,nts] = scipy.signal.filtfilt(weight_b, weight_a, data[:,nts])

  # Compute phase for each TS using Hilbert transform
  # **************************************************************************
  phase = zeros((dataLength-1,nTimeSeries))
  for k in range(0,nTimeSeries):
    hrp = scipy.signal.hilbert(data[:,k])
    for n in range(0,dataLength-1):
      phase[n,k] = arctan2( real(hrp[n]), imag(hrp[n]));
    phase[:,k]=unwrap(phase[:,k]);

  # Compute mean running (Cluster) phase
  # **************************************************************************
  clusterphase = zeros(dataLength-1)
  for n in range(0,dataLength-1):
    ztot = complex(0,0);
    for k in range(0,nTimeSeries):
        z = exp(1j * phase[n,k]);
        ztot = ztot+z;
    ztot = ztot/nTimeSeries;
    clusterphase[n] = angle(ztot);
  clusterphase = unwrap(clusterphase);

  # Compute relative phases between phase of TS and cluster phase
  # **************************************************************************
  complexIndRp = zeros((dataLength-1,nTimeSeries),dtype=complex);
  indRp = zeros((dataLength-1,nTimeSeries))
  meanIndRp = zeros(nTimeSeries);
  meanIndRho = zeros(nTimeSeries);
  for k in range(0,nTimeSeries):
    ztot = complex(0,0);
    for n in range(0,dataLength-1):
        z = exp(1j * (phase[n,k] - clusterphase[n]));
        complexIndRp[n,k] = z;
        ztot = ztot+z;
    indRp[:,k] = angle(complexIndRp[:,k]) * 360/(2*numpy.pi); # convert radian to degrees
    ztot = ztot / (dataLength-1);
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
  grpRho=zeros(dataLength-1);
  for n in range(0,dataLength-1):
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
  if plotFlag != False:
    generate_plot(data, sampleRate, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp, plotFlag)

  return data, meanGrpRho, meanIndRho, meanIndRp, grpRho, indRp
