import tmsanalysis as tms 
import pandas as pd
import time

datadir = '/farmshare/user_data/blenardo/struck_data/'
#filename = 'tier1_SIS3316Raw_20180925171528_Cs137_PSDScintTest_TriggThroughDisc_62_5MHz_1-ngm.root'
filename = 'tier1_SIS3316Raw_20181001221826_NaI_1700V_Co_10mVTh_ch1_125MHz__1-ngm.root'
datafile = datadir + filename

start_time = time.time()
dfout = tms.ProcessFile(datafile,save_waveforms=True,num_events=20)
end_time = time.time()

print('Time elapsed: {}s ({} min)'.format(end_time-start_time,(end_time-start_time)/60.))
