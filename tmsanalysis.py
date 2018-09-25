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
from scipy.fftpack import fft
from numpy import mean, sqrt, square, arange

import subprocess
root_version = subprocess.check_output(['root-config --version'], shell=True)
isROOT6 = False
if '6.1.0' in root_version or '6.04/06' in root_version:
    print "Found ROOT 6"
    isROOT6 = True

ROOT.gROOT.SetStyle("Plain")     
ROOT.gStyle.SetOptStat(0)        
ROOT.gStyle.SetPalette(1)        
ROOT.gStyle.SetTitleStyle(0)     
ROOT.gStyle.SetTitleBorderSize(0)       
ROOT.gStyle.SetPadLeftMargin(0.13);
ROOT.gStyle.SetPadBottomMargin(0.1);
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadRightMargin(0.13);
    
gROOT.ProcessLine('.L tms_utilities/fitShaped.C+')

#from ROOT import fitRising
from ROOT import fitShaped

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

  channel_type_map[0]  = 'PSD'
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
def ProcessFile( filename, num_events = -1):

    #print "processing file: ", filename

    channel_type_map = GetChannelTypeMap()
    num_channels = len(channel_type_map)

    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]

    # open the root file and grab the tree
    root_file = ROOT.TFile(filename)
    tree = root_file.Get("HitTree")
    n_entries = tree.GetEntries()
    print "%i entries in HitTree" % n_entries


    # Figure out how many records to process to get the right
    # number of events
    num_records_to_process = num_events * num_channels
    if num_records_to_process > n_entries or num_events < 0:
       num_records_to_process = n_entries

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
          
        for i_channel in xrange(num_channels):
            tree.GetEntry(i_entry)

            ch_num =  tree.HitTree.GetChannel()
            gate_size = tree.HitTree.GetGateCount()

            print "entry %i, channel %i" % (i_entry, channel)

            # Channel 6 is noisy, so we ignore it for now. 
            # (May have been destroyed by a discharge or something)
            if channel==6:
                i_entry = i_entry+1
                continue

            wfm_length = tree.HitTree.GetNSamples()
 
            graph = tree.HitTree.GetGraph()
            channel_wfm = np.array([graph.GetY()[isamp] for isamp in xrange(graph.GetN())])
            if channel_type_map[ch_num] == 'PSD':
               

 
            i_entry += 1

        i_event += 1   
 
    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]
    basename = basename.split("_")[2]
    outfilename = "./outfiles/out_new_5kV_%s.root" % basename
    #outfilename = "./outfiles/out_new_%s.root" % basename

    outfile=ROOT.TFile(outfilename,"RECREATE")
    h_chi2_ind.Write()
    h_chi2_col.Write()
    h_chi2_ind_wfm.Write()
    h_chi2_col_wfm.Write()
    h_amp.Write()
    h_amp_nocorr.Write()
    h_dx.Write()
    p_dQdx.Write()
    p_dQdx_dt.Write()
    h_time_x.Write()
    h_z_x.Write()
    mg_col.Write()
    mg_ind.Write()
    mg_xvsy.Write()
    #time.sleep(40)
    h_tau1.Write()
    h_tau2.Write()
    h_tau3.Write()
    h_tau4.Write()
    h_ratio1.Write()
    h_ratio2.Write()
    h_dqdx_dt.Write()
    h_hit_col.Write()
    h_hit_ind.Write()
    h_theta.Write()
    h_phi.Write()
    outfile.Write()
    '''
    c4.cd()
    h_chi2_ind.Draw()
    c5.cd()
    h_chi2_col.Draw()
    c6.cd()
    h_amp.Draw()
    c2.cd()
    h_dx.Draw()
    c10.cd()
    p_dQdx.Draw()
 
    c2.Print("dx.pdf")
    c4.Print("chi2_ind.pdf")
    c5.Print("chi2_col.pdf")
    c6.Print("amplitude.pdf")
    c10.Print("dQdx_z.pdf")
    '''
    #c4.Update()
        #if i_entry > 2: break # debugging

    #legend.Draw()
    #canvas.Update()
    #canvas.Print("%s_spectrum.png" % basename)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "arguments: [sis root files]"
        sys.exit(1)

    for filename in sys.argv[1:]:
        process_file(filename)



