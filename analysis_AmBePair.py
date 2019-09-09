"""
This script draws a spectrum from an NGM root tree. 

arguments [NGM sis tier1 root files]
"""

import os
import sys

import time
import ROOT
from ROOT import TF1
from ROOT import gROOT
from ROOT import *
#ROOT.gROOT.SetBatch(True)

from digi_filter import digi_hp
from digi_filter import digi_lp

from graph_trans import graph_fft
from graph_trans import graph_fft_hp_lp

import numpy as np
from scipy.fftpack import fft
from numpy import mean, sqrt, square, arange

import subprocess
root_version = subprocess.check_output(['root-config --version'], shell=True)
isROOT6 = False
if '6.1.0' in root_version or '6.04/06' in root_version:
    print "Found ROOT 6"
    isROOT6 = True

if os.getenv("EXOLIB") is not None and not isROOT6:
    try:
        gSystem.Load("$EXOLIB/lib/libEXOROOT")
        #from ROOT import CLHEP
        #microsecond = CLHEP.microsecond
        #second = CLHEP.second
        #print "imported CLHEP/ROOT"
    except (ImportError, AttributeError) as e:
        print "couldn't import CLHEP/ROOT"

ROOT.gROOT.SetStyle("Plain")     
ROOT.gStyle.SetOptStat(0)        
ROOT.gStyle.SetPalette(1)        
ROOT.gStyle.SetTitleStyle(0)     
ROOT.gStyle.SetTitleBorderSize(0)       
ROOT.gStyle.SetPadLeftMargin(0.13);
ROOT.gStyle.SetPadBottomMargin(0.1);
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadRightMargin(0.13);
ROOT.gStyle.SetEndErrorSize(0) 
  
#gROOT.ProcessLine('.L fitRising.C+')
gROOT.ProcessLine('.L fitShaped.C+')

#from ROOT import fitRising
from ROOT import fitShaped
#from ROOT import EXODoubleWaveform
from ROOT import EXOBaselineRemover
from ROOT import EXORisetimeCalculation
from ROOT import EXOSmoother
from ROOT import EXOPoleZeroCorrection
from ROOT import EXOExtremumFinder
from ROOT import EXOTrapezoidalFilter
from ROOT import EXOMatchedFilter

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

def process_file(filename):

    # options
    #channel = 0

    # maw options
    #gap_time = 30 # delay
    #peaking_time = 10 # length of moving average unit
 
    rise_time_calculator = EXORisetimeCalculation()
    # gt 30, pt 10 seems ok for 250 MS/s
    # gt 30, pt 2 to 4 seems ok for 25 MS/s
    gap_time = 30 # delay
    peaking_time = 2 # length of moving average unit

    #print "processing file: ", filename

    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]

    # open the root file and grab the tree
    root_file = ROOT.TFile(filename)
    tree = root_file.Get("HitTree")
    n_entries = tree.GetEntries()
    #n_entries = 400
    #if n_entries > 9000:
    #    n_entries = 9000
    #n_entries = 3000
    print "%i entries" % n_entries
    
    tree.SetLineWidth(2)
    tree.SetLineColor(ROOT.kBlue+1)
    tree.SetMarkerStyle(8)
    tree.SetMarkerSize(0.3)

    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]
    #basename = basename.split("_")[-1]
    #outfilename = "./outfiles/out_new_5kV_%s.root" % basename
    #outfilename = "/p/lscratchh/wu41/software/shuoxing_charge-readout-scripts/struck/ngm/tms_system/co60_calib/outfiles/new_out_8kV_%s.root" % basename
    outfilename = "./outfiles/new_out_8kV_%s.root" % basename
    outfile=ROOT.TFile(outfilename,"RECREATE")

    n_candidate = 0
    # set up a canvas
    #c1 = ROOT.TCanvas("canvas1","canvas1",600,400)
    #c2 = ROOT.TCanvas("canvas2","canvas2",600,400)
    #c3 = ROOT.TCanvas("canvas3","canvas3",600,400)
    #c1.SetGrid()

    h_chi2_ind = ROOT.TH1F("chi2_ind","chi2_ind",30,0,80)
    h_chi2_col = ROOT.TH1F("chi2_col","chi2_col",30,0,80)

    h_chi2_ind_wfm = ROOT.TH1F("chi2_ind_wfm","chi2_ind_wfm",50,0,50)
    h_chi2_col_wfm = ROOT.TH1F("chi2_col_wfm","chi2_col_wfm",50,0,50)

    h_time = ROOT.TH1F("Signal arrival time","Signal arrival time",1000,950,1050)

    #h_amp = ROOT.TH1F("Signal amplitude","Signal amplitude",35,0,100)
    #h_amp_nocorr=ROOT.TH1F("Signal amplitude no corr","Signal amplitude no corr",35,0,100)

    h_amp = ROOT.TH1F("Signal amplitude","Signal amplitude",1500,0,150)
    h_amp_nocorr=ROOT.TH1F("Signal amplitude no corr","Signal amplitude no corr",100,0,50)

    h_amp0 = ROOT.TH1F("Signal amplitude ch0","Signal amplitude ch0",400,0,50)
    h_amp1 = ROOT.TH1F("Signal amplitude ch1","Signal amplitude ch1",400,0,50)
    h_amp2 = ROOT.TH1F("Signal amplitude ch2","Signal amplitude ch2",400,0,50)
    h_amp3 = ROOT.TH1F("Signal amplitude ch3","Signal amplitude ch3",400,0,50)
    h_amp4 = ROOT.TH1F("Signal amplitude ch4","Signal amplitude ch4",400,0,50)
    h_amp5 = ROOT.TH1F("Signal amplitude ch5","Signal amplitude ch5",400,0,50)
    h_amp6 = ROOT.TH1F("Signal amplitude ch6","Signal amplitude ch6",400,0,50)
    h_amp7 = ROOT.TH1F("Signal amplitude ch7","Signal amplitude ch7",400,0,50)

    h_amp_NaI = ROOT.TH1F("Signal amplitude NaI","Signal amplitude NaI",10000,0,5000)
    h_amp_Plastic = ROOT.TH1F("Signal amplitude Plastic","Signal amplitude Plastic",500,0,500)
    h_amp_Che = ROOT.TH1F("Signal amplitude CK","Signal amplitude CK",10000,0,10000)

    h_amp_Che_NaI = ROOT.TH2F("Signal amp CK vs NaI","Signal amp CK vs NaI",10000,0,5000,10000,0,10000)

    h_amp_NaI_nY = ROOT.TH2F("Signal amp NaI vs nY","Signal amp NaI vs nY",10,0,10,10000,0,5000)

    h_dx = ROOT.TH1F("DeltaX","DeltaX",100,0,20)

    p_dQdx = ROOT.TProfile("p_dQ_dx_vs_Z","p_dQ_dx_vs_Z",20,-5,35,0,100,"g")
    #p_dQdx_dt = ROOT.TProfile("p_dQ_dx_vs_dt","p_dQ_dx_vs_dt",30,-1,5,0,100,"g")
    p_dQdx_dt = ROOT.TProfile("p_dQ_dx_vs_dt","p_dQ_dx_vs_dt",60,-1,5,0,100,"g")

    #h_dqdx_dt = ROOT.TH2F("h_dQ_dx_vs_dt","h_dQ_dx_vs_dt",15,-1,5,50,0,100)
    h_dqdx_dt = ROOT.TH2F("h_dQ_dx_vs_dt","h_dQ_dx_vs_dt",60,-1,5,50,0,100)

    h_time_x = ROOT.TH1F("arrival_time_distribution_X","arrival_time_distribution_X",60,-1,5)
    h_z_x = ROOT.TH1F("arrival_z_distribution_X","z_time_distribution_X",15,0,15)

    h_time_y = ROOT.TH1F("arrival_time_distribution_Y","arrival_time_distribution_Y",60,-1,5)

    h_tau1 = ROOT.TH1F("h_tau1","h_tau1",100,0,300)
    h_tau2 = ROOT.TH1F("h_tau2","h_tau2",100,0,300)
    h_tau3 = ROOT.TH1F("h_tau3","h_tau3",100,0,300)
    h_tau4 = ROOT.TH1F("h_tau4","h_tau4",100,0,300)
     
    h_hit_col = ROOT.TH1F("h_hit_col","h_hit_col",10,0,10)
    h_hit_ind = ROOT.TH1F("h_hit_ind","h_hit_ind",10,0,10)

    h_theta = ROOT.TH1F("h_theta","h_hit_col",90,-90,90)
    h_phi = ROOT.TH1F("h_phi","h_phi",90,-90,90)

    h_ratio1 = ROOT.TH1F("h_tau_ratio1","h_tau_ratio1",60,0,3)
    h_ratio2 = ROOT.TH1F("h_tau_ratio2","h_tau_ratio2",60,0,3)

    
    i_entry = 0
    n_track = 0
    while i_entry < n_entries:
    #for i_entry in xrange(n_entries):
        if(i_entry%32==0):
	    sum_x_energy = 0
            ch_x_array = []
            energy_x_array = []
            time_x_array = []
            time_err_x_array = []
            ch_y_array = []
            time_y_array = []
            time_err_y_array = []
            energy_NaI = 0
	    energy_Plastic = 0
	    energy_CK = 0
	    hist_x_collection = []
	    hist_x_collection_fft = []

        #tree.GetEntry(i_entry)

        n_channels=32
          
        for i_channel in xrange(n_channels):
            tree.GetEntry(i_entry)

            channel =  tree.HitTree.GetChannel()
            gate_size = tree.HitTree.GetGateCount()

            slot = tree.HitTree.GetSlot()

            print "entry %i, channel %i, slot %i" % (i_entry, channel, slot)

	    if slot ==1:
	        i_entry=i_entry+1

                graph = tree.HitTree.GetGraph()

		rms_wfm_pmt = np.array([graph.GetY()[isamp] for isamp in xrange(1850,1950,1)])

                rms_wfm_pmt_mean = np.mean(rms_wfm_pmt)

		rms_wfm_pmt = rms_wfm_pmt-rms_wfm_pmt_mean

		signal_wfm_pmt = np.array([graph.GetY()[isamp] for isamp in xrange(1950,2100,1)])

                signal_wfm_pmt = signal_wfm_pmt - rms_wfm_pmt_mean

		noise_rms = np.sqrt(rms_wfm_pmt.dot(rms_wfm_pmt)/rms_wfm_pmt.size)

		signal_amp_pmt = np.abs(np.min(signal_wfm_pmt))
		#if signal_amp_pmt>3*noise_rms:
		#    signal_amp_pmt = np.abs(np.min(signal_wfm_pmt))
		#else:
		#    signal_amp_pmt = rms_wfm_pmt[0]

		#signal_amp_pmt = np.abs(np.min(signal_wfm_pmt))

                if(channel==0):
		    h_amp_Plastic.Fill(signal_amp_pmt)
		    energy_Plastic = signal_amp_pmt
		elif(channel==1):
		    h_amp_NaI.Fill(signal_amp_pmt)
		    energy_NaI = signal_amp_pmt
		elif(channel==2):
		    h_amp_Che.Fill(signal_amp_pmt)
		    energy_CK = signal_amp_pmt

	        continue

            wfm_length = tree.HitTree.GetNSamples()
 
            #print "sample length: ", wfm_length

            graph = tree.HitTree.GetGraph()

            graph.SetTitle("")
            graph.GetXaxis().SetTitle("Samples [16 ns]")
            graph.GetYaxis().SetTitle("ADC count")
            graph.GetXaxis().CenterTitle(1)
            graph.GetYaxis().CenterTitle(1)
            graph.GetYaxis().SetTitleOffset(1.4)

            legend = ROOT.TLegend(0.7, 0.92, 0.9, 0.99)
            legend.AddEntry(graph, "event %i, channel %i" % (i_entry/32+1,channel), "pl")
            
            graph_fftcut=ROOT.TGraph()

	    graph_fftcut = graph_fft(graph)

            graph_fftcut_hp_lp = ROOT.TGraph()
            graph_fftcut_hp_lp = graph_fft_hp_lp(graph)

            graph_fftcut_hp_lp.SetName("fft hp lp event %i channel %i" % (i_entry/32+1, channel))
            graph_fftcut.SetName("fft event %i channel %i" % (i_entry/32+1, channel))

            if(channel<8):
		fitFunc = TF1("fitLine","[0]*x[0]+[1]", 0, 1800, 2)
		graph_fftcut.Fit(fitFunc,"QS","",0,1800)
                par0=fitFunc.GetParameter(0)
                par1=fitFunc.GetParameter(1)

                print par0, par1

		for i_point in xrange(graph_fftcut.GetN()): 
	            x = graph_fftcut.GetX()[i_point]
		    y = graph_fftcut.GetY()[i_point]
		    graph_fftcut.SetPoint(i_point,x,y-(par0*x+par1)) 

            noise_wfm_fft_hp_lp = np.array([graph_fftcut_hp_lp.GetY()[isamp] for isamp in xrange(1750,1950,1)])

            signal_wfm_fft_hp_lp = np.array([graph_fftcut_hp_lp.GetY()[isamp] for isamp in xrange(1950,2400,1)])

            noise_baseline_fft_hp_lp = np.mean(noise_wfm_fft_hp_lp)
	    noise_wfm_fft_hp_lp = noise_wfm_fft_hp_lp - noise_baseline_fft_hp_lp
            noise_rms_fft_hp_lp=np.sqrt(noise_wfm_fft_hp_lp.dot(noise_wfm_fft_hp_lp)/noise_wfm_fft_hp_lp.size)
                     
	    signal_wfm_fft_hp_lp = signal_wfm_fft_hp_lp - noise_baseline_fft_hp_lp
  
            noise_wfm_fft = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(1600,1800,1)])
            #signal_wfm_fft = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(2600,2800,1)])
            signal_wfm_fft = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(3000,3500,1)])
            
            noise_rms_fft=np.sqrt(noise_wfm_fft.dot(noise_wfm_fft)/noise_wfm_fft.size)

            if channel<8 and slot==0:
                signal_amp_fft_hp_lp = np.max(signal_wfm_fft_hp_lp)  
		signal_amp_fft = np.mean(signal_wfm_fft)-np.mean(noise_wfm_fft)

	    elif channel>7 and slot==0 :
                signal_amp_fft_hp_lp = np.min(signal_wfm_fft_hp_lp)
		signal_amp_pos_fft_hp_lp = np.max(signal_wfm_fft_hp_lp)

            #if (((channel<8) and (signal_amp)>5*noise_rms) or ((channel>7) and ((-signal_amp) > 3.5*noise_rms) and (signal_amp_pos > 3.5*noise_rms))) :
            if (((channel<8) and (signal_amp_fft)>4*noise_rms_fft) or ((channel>7) and ((-signal_amp_fft_hp_lp) > 3.5*noise_rms_fft_hp_lp))):
                if(channel<8 and slot==0): 

		    energy_x_array.append(signal_amp_fft)
		    ch_x_array.append(channel)
		    hist_x_collection.append(graph_fftcut_hp_lp)
		    hist_x_collection_fft.append(graph_fftcut)

		    sum_x_energy = sum_x_energy+signal_amp_fft
    
		elif(channel>7 and slot==0):
		    ch_y_array.append(channel-8)
            i_entry += 1
        
        #if(i_entry%16==0):   
        if((i_channel+1)%32==0):   
            print "event %i" % (i_entry/32+1)
            print "length of X wire signals %i: " % len(ch_x_array) 
            print "length of Y wire signals %i: " % len(ch_y_array) 
	    h_amp_Che_NaI.Fill(energy_NaI,energy_CK)
	    h_amp_NaI_nY.Fill(len(ch_y_array),energy_NaI)
            #if len(ch_x_array)>3 :
            #if (len(ch_x_array)>1) and (len(ch_y_array)>2) :
            #if (len(ch_x_array)>2) and (len(ch_y_array)>4) and ch_y_array[0]<10 and ch_y_array[len(ch_y_array)-1]>50 :
            #if (len(ch_x_array)<3 and len(ch_x_array)>0 and len(ch_y_array)<3 and len(ch_y_array)>0):
            #if (len(ch_x_array)<3 and len(ch_x_array)>0 and len(ch_y_array)<2 and energy_NaI<50 and energy_NaI>16):
            if (len(ch_x_array)<3 and len(ch_x_array)>0 and len(ch_y_array)<2 and energy_NaI<500):
	        if(ch_x_array[0]==0):
		    h_amp0.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==1):
		    h_amp1.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==2):
		    h_amp2.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==3):
		    h_amp3.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==4):
		    h_amp4.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==5):
		    h_amp5.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==6):
		    h_amp6.Fill(energy_x_array[0])
	        elif(ch_x_array[0]==7):
		    h_amp7.Fill(energy_x_array[0])

                #h_amp.Fill(energy_x_array[0])
		#if(len(ch_x_array)==1 or (len(ch_x_array)==2 and np.abs(ch_x_array[0]-ch_x_array[1])==1)):
		if(len(ch_x_array)==1):
                    h_amp.Fill(sum_x_energy)
                #for i in xrange(len(ch_x_array)):
		#    h_time_x.Fill(time_x_array[i]/5.2)
                    for id in xrange(len(ch_x_array)):
	                hist_x_collection[id].SetLineColor(n_candidate)
	                hist_x_collection[id].SetMarkerColor(n_candidate)
                        hist_x_collection[id].Write()

	                hist_x_collection_fft[id].SetLineColor(n_candidate)
	                hist_x_collection_fft[id].SetMarkerColor(n_candidate)
                        hist_x_collection_fft[id].Write()

                    n_candidate = n_candidate+1

    h_amp.Write()
    h_amp0.Write()
    h_amp1.Write()
    h_amp2.Write()
    h_amp3.Write()
    h_amp4.Write()
    h_amp5.Write()
    h_amp6.Write()
    h_amp7.Write()
    h_amp_Plastic.Write()
    h_amp_NaI.Write()
    h_amp_Che.Write()
    h_amp_Che_NaI.Write()
    h_amp_NaI_nY.Write()
    #for id in xrange(len(ch_x_array)):
    #	hist_x_collection[id].SetLineColor(n_candidate)
    #    hist_x_collection[id].Write()
    
    #mg_col.Write()
    #mg_ind.Write()
    #mg_xvsy.Write()
    #time.sleep(40)
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



