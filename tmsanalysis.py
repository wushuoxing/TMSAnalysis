import os
import sys

import time
import ROOT
from ROOT import TF1
from ROOT import gROOT
from ROOT import *
#ROOT.gROOT.SetBatch(True)

from tms_utilities.digi_filter import digi_hp
from tms_utilities.digi_filter import digi_lp

import numpy as np
import pandas as pd
from scipy.fftpack import fft
from numpy import mean, sqrt, square, arange

import parse_psd_wfm
import parse_nai_wfm
import parse_tms_wfm

import subprocess
root_version = subprocess.check_output(['root-config --version'], shell=True)
isROOT6 = False
if '6.1.0' in root_version or '6.04/06' in root_version:
    print("Found ROOT 6")
    isROOT6 = True


def print_fit_info(fit_result, fit_duration):
        
    chi2 = fit_result.Chi2()
    prob = fit_result.Prob()
    ndf = fit_result.Ndf()
    status = fit_result.Status()
    #print "fit results:"
    #print "\tchi2: %.2f" % chi2
    #print "\tn dof", ndf
    #print "\tchi2 / dof: %.2f" % (chi2/ndf)
    #print "\tprob", prob
    #print "\tstatus", status
    #print "\t%.1f seconds" % fit_duration

########################################################################
def GetChannelTypeMap():
  num_channels = 16
  channel_type_map = ['empty' for i in range(0,num_channels)]

  channel_type_map[0]  = 'Xwire'
  channel_type_map[1]  = 'Xwire'
  channel_type_map[2]  = 'Xwire'
  channel_type_map[3]  = 'Xwire'
  channel_type_map[4]  = 'Xwire'
  channel_type_map[5]  = 'Xwire'
  channel_type_map[6]  = 'Xwire'
  channel_type_map[7]  = 'Xwire'
  channel_type_map[8]  = 'Ywire'
  channel_type_map[9]  = 'Ywire'
  channel_type_map[10] = 'Ywire'
  channel_type_map[11] = 'Ywire'
  channel_type_map[12] = 'Ywire'
  channel_type_map[13] = 'Ywire'
  channel_type_map[14] = 'Ywire'
  channel_type_map[15] = 'Ywire'

  return channel_type_map

########################################################################
def ProcessFile( filename, num_events = -1, save_waveforms = False):

    #print "processing file: ", filename

    channel_type_map = GetChannelTypeMap()
    num_channels = len(channel_type_map)
    #print('Num channels: {}'.format(num_channels))

    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]

    # open the root file and grab the tree
    root_file = ROOT.TFile(filename)
    tree = root_file.Get("HitTree")
    n_entries = tree.GetEntries()
    #print('{} entries in HitTree'.format(n_entries))


    # Figure out how many records to process to get the right
    # number of events
    num_records_to_process = num_events * num_channels
    if num_records_to_process > n_entries or num_events < 0:
       num_records_to_process = n_entries

    output_dataframe = pd.DataFrame()

    i_entry = 0
    i_event = 0
    while i_event < num_events:
    #for i_entry in xrange(n_entries):
        if( i_entry % num_channels == 0 ):

            ch_x_array = []
            energy_x_array = []
            time_x_array = []
            time_err_x_array = []
            ch_y_array = []
            time_y_array = []
            time_err_y_array = []
 
        #tree.GetEntry(i_entry)
        output_series = pd.Series()
 
        for i_channel in xrange(num_channels):
            tree.GetEntry(i_entry)

            slot_num = tree.HitTree.GetSlot()
            ch_num =  tree.HitTree.GetChannel()
            gate_size = tree.HitTree.GetGateCount()

            #print('entry {}, channel {}'.format(i_entry, ch_num))

            # Channel 6 is noisy, so we ignore it for now. 
            # (May have been destroyed by a discharge or something)
            if ch_num==6:
                i_entry = i_entry+1
                continue

            wfm_length = tree.HitTree.GetNSamples()
 
            graph = tree.HitTree.GetGraph()
            channel_wfm = np.array([graph.GetY()[isamp] for isamp in xrange(graph.GetN())])

            if channel_type_map[ch_num] == 'PSD':
               psd_reduced_data = parse_psd_wfm.ParsePSDWfm( channel_wfm, save_waveform=save_waveforms )
               for colname in psd_reduced_data.index:
                  output_series[colname] = psd_reduced_data[colname]

            if channel_type_map[ch_num] == 'Xwire':
               xwire_reduced_data = parse_tms_wfm.ParseTMSXwireWfm( channel_wfm, save_waveform = save_waveforms )
               for colname in xwire_reduced_data.index:
                 output_series[colname] = xwire_reduced_data[colname]

            if channel_type_map[ch_num] == 'Ywire':
               ywire_reduced_data = parse_tms_wfm.ParseTMSYwireWfm( channel_wfm, save_waveform = save_waveforms )

            else:
               # For now, do nothing here.
               i_entry += 1
               continue

            #print(psd_reduced_data)


            i_entry += 1

        output_dataframe = output_dataframe.append( output_series, ignore_index=True )
        i_event += 1   
 
    #basename = os.path.basename(filename)
    #basename = os.path.splitext(basename)[0]
    #basename = basename.split("_")[2]
    #outfilename = "./outfiles/out_new_5kV_%s.root" % basename
    #outfilename = "./outfiles/out_new_%s.root" % basename

    #outfile=ROOT.TFile(outfilename,"RECREATE")

    output_dataframe.to_pickle('psd_output_test.pkl')
    return output_dataframe


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print('arguments: [sis root files]')
        sys.exit(1)

    for filename in sys.argv[1:]:
        process_file(filename)



