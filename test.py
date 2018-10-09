import tmsanalysis as tms 
import pandas as pd
import time

datadir = '/farmshare/user_data/blenardo/struck_data/'
#filename = 'tier1_SIS3316Raw_20180925171528_Cs137_PSDScintTest_TriggThroughDisc_62_5MHz_1-ngm.root'
filename = 'tier1_SIS3316Raw_20181004010243_Co60_CH1NaIOff_CH0PSD1100V_PSDTrigg150mVTh_PSDX10_125MHz__1-ngm.root'
datafile = datadir + filename

start_time = time.time()
dfout = tms.ProcessFile(datafile,save_waveforms=False,num_events=-1)
end_time = time.time()

print('Time elapsed: {}s ({} min)'.format(end_time-start_time,(end_time-start_time)/60.))
