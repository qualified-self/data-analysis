#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#   cluster_phase.py
#
#   GRPrhoM, INDrhoM, INDrpM, TSrhoGRP, TSrpIND = cluster_phase(TSfilename, TSnumber, TSfsamp, TSlsamp, TSsamplerate, plotflag)
#
#   Input:
#       TSfilename    : data file
#       TSnumber      : number of time series
#       TSsamplerate  : sample rate of the time series
#       TSfsamp       : first data point in time series used
#       TSlsamp       : last data point in time series used
#       plotflag      : do plots (0=no, 1=yes)
#
#   Output:
#       GRPrhoM        : mean group rho (0 to 1; 1 = perfect sync)
#       INDrhoM        : mean rho for each TS to group (0 to 1; 1 = perfect sync)
#       INDrpM         : mean Relative Phase for each TS to group cluster phase 
#       TSrhoGRP        : group rho time-series
#       TSrpIND         : relative phase time-series for each individual TS to cluster phase
#
#   Example:
#       [GRPrhoM INDrhoM INDrpM TSrhoGRP TSrpIND] = ClusterPhase_do('G201EO1.txt', 6, 1, 7200, 120, 1);
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

def cluster_phase(TSfilename,TSnumber,TSfsamp,TSlsamp,TSsamplerate,plotflag=False):
  
  filterfreq = 10
  
  # load time-series (TS)
  # **************************************************************************
  data = loadtxt(TSfilename)
  
  # Builds a subset by taking only rows TSfsamp .. TSlsamp from base dataset
  ts_data = data[TSfsamp:TSlsamp,0:TSnumber]
  
  TSlength = ts_data.shape[0]
  
  delta_t = 1.0/TSsamplerate
  t = arange(0, TSlength) * delta_t

  # Downsample, Filter and normalize data
  # **************************************************************************
  # linear detrend data to remove drift (chiar moving slightly during trial
  for nts in range(0,TSnumber):
    ts_data[:,nts] =      scipy.signal.detrend(ts_data[:,nts])
 
  # normalize
  for nts in range(0,TSnumber):
    ts_data[:,nts] = scipy.stats.mstats.zscore(ts_data[:,nts], ddof=1)
 
  # filter
  for nts in range(0,TSnumber):
    weight_b,weight_a = scipy.signal.butter(2, filterfreq/(TSsamplerate/2.0));
    irlength = max(weight_b.size-1,weight_a.size-1)
    ts_data[:,nts] = scipy.signal.filtfilt(weight_b, weight_a, ts_data[:,nts])
  
  # Compute phase for each TS using Hilbert transform
  # **************************************************************************
  TSphase = zeros((TSlength-1,TSnumber))
  for k in range(0,TSnumber):
    hrp = scipy.signal.hilbert(ts_data[:,k])
    for n in range(0,TSlength-1):
      TSphase[n,k] = arctan2( real(hrp[n]), imag(hrp[n]));
    TSphase[:,k]=unwrap(TSphase[:,k]);
    
  # Compute mean running (Cluster) phase
  # **************************************************************************
  clusterphase = zeros(TSlength-1)
  for n in range(0,TSlength-1):
    ztot = complex(0,0);
    for k in range(0,TSnumber):
        z = exp(1j * TSphase[n,k]);
        ztot = ztot+z;
    ztot = ztot/TSnumber;
    clusterphase[n] = angle(ztot);
  clusterphase = unwrap(clusterphase);

  # Compute relative phases between phase of TS and cluster phase
  # **************************************************************************
  complexTSrpIND = zeros((TSlength-1,TSnumber),dtype=complex);
  TSrpIND = zeros((TSlength-1,TSnumber))
  INDrpM = zeros(TSnumber);
  INDrhoM = zeros(TSnumber);
  for k in range(0,TSnumber):
    ztot = complex(0,0);
    for n in range(0,TSlength-1):
        z = exp(1j * (TSphase[n,k] - clusterphase[n]));
        complexTSrpIND[n,k] = z;
        ztot = ztot+z;
    TSrpIND[:,k] = angle(complexTSrpIND[:,k]) * 360/(2*numpy.pi); # convert radian to degrees
    ztot = ztot / (TSlength-1);
    INDrpM[k] = angle(ztot);
    INDrhoM[k] = abs(ztot);
  TSRPM = INDrpM;
  INDrpM = (INDrpM / (2*numpy.pi)*360); # convert radian to degrees
  print(' ');
  print('Mean relative phases of individuals to cluster phase')
  print(INDrpM);
  print('Averaged degree of synchronization of individuals (Rho = 1-circular variance)')
  print(INDrhoM);

  print TSrpIND

  # Compute cluster amplitude rhotot in rotation frame
  # **************************************************************************
  TSrhoGRP=zeros(TSlength-1);
  for n in range(0,TSlength-1):
    ztot = complex(0,0);
    for k in range(0,TSnumber):
        z = exp(1j * (TSphase[n,k] - clusterphase[n] - TSRPM[k]));
        ztot = ztot+z;
    ztot = ztot / TSnumber;
    TSrhoGRP[n] = abs(ztot);
  print TSrhoGRP
  
  GRPrhoM = mean(TSrhoGRP);
  print('Averaged degree of synchronization of the group')
  print(GRPrhoM);
  
  # Do Plot
  # **************************************************************************
  # plot data for time-series (separeted on graph for display purposes)
  #scrsz = get(0,'ScreenSize');
  #h = figure('Position',[scrsz(3)/3 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);
  if plotflag != False:

    fig = figure(1)
    
    subplot(3,1,1);
    tmpdata = zeros((TSlength,TSnumber));
    for nts in range(0,TSnumber):
      tmpdata[:,nts] = (ts_data[:,nts] + (nts*4));
    plot(t, tmpdata)
    xlabel('Time');
    ylabel('RAW Data');
    xlim([0, max(t)]);
    #ylim([-185, 185]);
      
    # plot individ-cluster relative phase
    subplot(3,1,2);
    plot(t[0:TSlength-1], TSrpIND);
    xlabel('Time');
    ylabel('IND-Clust Relative Phase');
    xlim([0, max(t)]);
    ylim([-185, 185]);

    # plot group-cluster amplitiude (rho) timeseries
    subplot(3,1,3);
    plot(t[0:TSlength-1], TSrhoGRP)
    xlabel('Time');
    ylabel('GRP-Clust Amplitude');
    xlim([0, max(t)]);
    ylim([0, 1]);

    text(0, -.4, "Mean GRP Rho: {:.3f}  Mean IND Rhos: {:s}".format(GRPrhoM, array_str(INDrhoM,precision=3)))
    
    if plotflag == True:
      fig.show()
    else:
      pickle.dump(fig, file(plotflag, 'w'))

  return GRPrhoM, INDrhoM, INDrpM, TSrhoGRP, TSrpIND

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
    GRPrhoM, INDrhoM, INDrpM, TSrhoGRP, TSrpIND = cluster_phase(basename + ".raw", nSubjects, rng[0], rng[1], config["sample_freq"], False)
 
if __name__ == "__main__":
  main()
