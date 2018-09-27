import tmsanalysis as tms 
import pandas as pd
import time
from ROOT import gROOT

gROOT.ProcessLine('.L fitShaped.C+')

datadir = '~/software/charge-readout-scripts/struck/ngm/'
filename = 'tier1_SIS3316Raw_201807251214228000VCa_3200VX_TopAndBotPMT_ScopeTrigg_62_5MHz_1-ngm.root' 
datafile = datadir + filename

start_time = time.time()
dfout = tms.ProcessFile(datafile,save_waveforms=True,num_events=2)
end_time = time.time()

print('Time elapsed: {}s ({} min)'.format(end_time-start_time,(end_time-start_time)/60.))
