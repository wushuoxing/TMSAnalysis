#!/usr/bin/env python

"""
This script draws events from NGM root file(s). 

arguments [NGM root files of events: tier1_SIS3316Raw_*.root]
"""

import os
import sys
import math
import commands
#from scipy.fftpack import fft

import ROOT
ROOT.gROOT.SetBatch(True) # uncomment to draw multi-page PDF

#from struck import struck_analysis_parameters
import struck_analysis_parameters
from graph_trans import graph_fft_hp_lp
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

# set the ROOT style
ROOT.gROOT.SetStyle("Plain")     
ROOT.gStyle.SetOptStat(0)        
ROOT.gStyle.SetPalette(1)        
ROOT.gStyle.SetTitleStyle(0)     
ROOT.gStyle.SetTitleBorderSize(0)       
print "title align:", ROOT.gStyle.GetTitleAlign() 
ROOT.gStyle.SetTitleX(.5)
# default: 13
ROOT.gStyle.SetTitleAlign(23) 
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadTickX(1)
#ROOT.gStyle.SetGridStyle(0);
#h2->GetXaxis()->SetNdivisions(-205);

# set up a canvas
canvas = ROOT.TCanvas("canvas","", 1000, 800)
canvas.SetGrid(1,1)
canvas.SetLeftMargin(0.12)
canvas.SetTopMargin(0.12)
canvas.SetBottomMargin(0.12)

# for RHS legend:
canvas.SetTopMargin(0.05)
canvas.SetBottomMargin(0.1)
canvas.SetRightMargin(0.11)

ROOT.gStyle.SetTitleFontSize(0.04)

nchannels = len(struck_analysis_parameters.channel_map)
#nchannels = 16 # FIXME -- for DT unit!!

print "%i channels" % nchannels

def process_file(filename=None, n_plots_total=0):

    # options ------------------------------------------
    is_for_paper = True # change some formatting
    threshold = 550 # keV
    energy_offset = 700 # keV, space between traces

    units_to_use = 1 # 0=keV, 1=ADC units, 2=mV

    # if this is 1.0 there is no effect:
    pmt_shrink_scale = 0.5 # shrink PMT signal by an additional factor so it doesn't set the graphical scale

    #------------------------------------------------------

    # figure out which units we are using, to construct the file name
    units = "keV"
    if units_to_use == 1:
        units = "AdcUnits"
    elif units_to_use == 2:
        units = "mV"

    y_min = -200 # adc

    y_max = 17*energy_offset # +250 # keV

    #If input filename is null assume we want to examine the most recent file
    if(filename == None):
        # example name: SIS3316Raw_20160712204526_1.bin
        output  = commands.getstatusoutput("ls -rt tier1_SIS3316Raw_*.root | tail -n1")
        filename = output[1]
        print "--> using most recent NGM file, ", filename

    print "--> processing file", filename

    channel_map = struck_analysis_parameters.channel_map
    top_pmt_channel = struck_analysis_parameters.top_pmt_channel
    print "top_pmt_channel:", top_pmt_channel
    bot_pmt_channel = struck_analysis_parameters.bot_pmt_channel
    print "bot_pmt_channel:", bot_pmt_channel
    che_pmt_channel = struck_analysis_parameters.che_pmt_channel
    print "che_pmt_channel:", che_pmt_channel

    charge_channels_to_use = struck_analysis_parameters.charge_channels_to_use
    channels = []
    for (channel, value) in enumerate(charge_channels_to_use):
        #print channel, value
        if value is not 0: channels.append(channel)
    print "there are %i charge channels" % len(channels)
    channels.append(top_pmt_channel)
    channels.append(bot_pmt_channel)
    channels.append(che_pmt_channel)
    colors = struck_analysis_parameters.get_colors()

    basename = os.path.basename(filename)
    basename = os.path.splitext(basename)[0]
    #print "basename:", basename

    # open the root file and grab the tree
    root_file = ROOT.TFile(filename)
    tree = root_file.Get("HitTree")
    n_entries = tree.GetEntries()
    print "%i entries in tree" % n_entries

    # get NGM system config
    sys_config = root_file.Get("NGMSystemConfiguration")
    card0 = sys_config.GetSlotParameters().GetParValueO("card",0)
    card1 = sys_config.GetSlotParameters().GetParValueO("card",1)

    #sampling_freq_Hz = struck_analysis_parameters.get_clock_frequency_Hz_ngm(card0.clock_source_choice)
    sampling_freq_Hz = 62.5e6
    print "sampling_freq_Hz: %.1f MHz" % (sampling_freq_Hz/1e6)

    channel_map = struck_analysis_parameters.channel_map

    tree.GetEntry(0) # use 0th entry to get some info that we don't expect to change... 
    trace_length_us = tree.HitTree.GetNSamples()/sampling_freq_Hz*1e6
    print "length in microseconds:", trace_length_us

    trigger_time = card0.pretriggerdelay_block[0]/sampling_freq_Hz*1e6
    print "trigger time: [microseconds]", trigger_time

    frame_hist = ROOT.TH1D("hist", "", tree.HitTree.GetNSamples(), 0, trace_length_us+0.5)
    #frame_hist = ROOT.TH1D("hist", "", int(tree.HitTree.GetNSamples()*33.0/32.0), 0, 25.0)
    frame_hist.SetXTitle("Time [#mus]")
    sum_offset = 400
    frame_hist.SetYTitle("Energy (with arbitrary offsets) [keV]")
    if units_to_use == 1: # ADC units
        frame_hist.SetYTitle("ADC units")

    frame_hist.GetYaxis().SetTitleOffset(1.6)
    frame_hist.GetYaxis().SetLimits(-100,13000)
    frame_hist.GetYaxis().SetRangeUser(-100,13000)
    frame_hist.GetXaxis().SetRangeUser(20,70)
    #frame_hist.GetXaxis().SetNdivisions(-205);

    #legend = ROOT.TLegend(0.15, 0.86, 0.9, 0.99)
    legend = ROOT.TLegend(
        #1.005 - canvas.GetRightMargin(), 
        1.005 - canvas.GetRightMargin(), 
        canvas.GetBottomMargin(), 
        0.985, 
        1.0-canvas.GetTopMargin()
    )
    #legend.SetNColumns(7)
    #legend.SetBorderSize(0)

    # set up some placeholder hists for the legend
    hists = []
    for (i, channel) in enumerate(channels):
        #print "i=%i, ch=%i (%s)" % (i, channel, channel_map[channel])
	#histname = "hist %i" % (channel)
        #hist = ROOT.TH1D(histname,"",10,0,10)
        #hist = ROOT.TH1D("hist %i" % channel,"",10,0,10)
        hist = ROOT.TH1D("hist %i" % i,"",10,0,10)
        try:
            color = colors[i]
        except IndexError:
            color = ROOT.kBlack
        if(channel==bot_pmt_channel): color=1

        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hists.append(hist)
    
    # loop over all events in file
    i_entry = 0
    n_plots = 0
    y_max_old = y_max
    y_min_old = y_min
    # use while loop instead of for loop so we can modify i_entry if needed
    while i_entry < n_entries:

        canvas.SetLogy(0)
        canvas.SetLogx(0)

        frame_hist.Draw()
        wfm_length = tree.HitTree.GetNSamples()

        #print "==> entry %i of %i | charge energy: %i" % ( i_entry, n_entries, chargeEnergy,)
      
        legend.Clear()
        legend_entries = [""]*nchannels
	#legend_entries = [""]*17

        graph_total = ROOT.TMultiGraph()

        #for i in xrange(nchannels):
        for i in xrange(32):
            
            tree.GetEntry(i_entry)
            i_entry += 1
            slot = tree.HitTree.GetSlot()
            card_channel = tree.HitTree.GetChannel() # 0 to 16 for each card

            channel = card_channel + 16*slot # 0 to 31
            card = sys_config.GetSlotParameters().GetParValueO("card",slot)

            if(channel>18): continue

            gain = card.gain[card_channel]
            # gain: 1 = 2V; 0 = 5V
            voltage_range_mV = struck_analysis_parameters.get_voltage_range_mV_ngm(gain)

            try:
                color = colors[channel]
            except IndexError:
                color = ROOT.kBlack

            multiplier = 1
            if units_to_use == 1: # ADC units
                multiplier = 1.0
            elif units_to_use == 2: # mV
                multiplier = voltage_range_mV/pow(2,16)
	    if channel>7 and channel != top_pmt_channel and channel != bot_pmt_channel and channel != che_pmt_channel:
                multiplier = 50
	    elif channel<8:
	        multiplier = 25

            graph_orig = tree.HitTree.GetGraph()

            graph = ROOT.TGraph()
	   
	    if(channel!=top_pmt_channel and channel!=bot_pmt_channel and channel!=che_pmt_channel):
	        graph = graph_fft_hp_lp(graph_orig)
            else:
		graph = graph_orig

            # as in http://exo-data.slac.stanford.edu/exodoc/src/EXOBaselineRemover.cxx.html#48
            baseline = 0.0
            energy = 0.0
            baseline_avg_sq = 0.0
            energy_avg_sq = 0.0
            #for i_sample in xrange(n_samples_to_avg):
	    for i_sample in xrange(800,1800):
                y = graph.GetY()[i_sample]
                y2 = graph.GetY()[i_sample+2000]
                baseline += y / 1000.0
                energy += y2 / 1000.0
                baseline_avg_sq += y*y/1000.0
                energy_avg_sq += y2*y2/1000.0

            # add an offset so the channels are drawn at different levels
	    if(channel!=top_pmt_channel and channel!=bot_pmt_channel and channel!=che_pmt_channel):
                offset = channel*energy_offset+200
            else:
                offset = channel*energy_offset+200

            if(channel==top_pmt_channel or channel==bot_pmt_channel):
                signal_wfm = np.array([graph.GetY()[isamp] for isamp in xrange(1800,2100,1)])
		energy=np.min(signal_wfm)
                multiplier = 12000.0/np.abs(energy-baseline)
                multiplier = 3
            if(channel==che_pmt_channel):
                multiplier = 5

            print channel, baseline, offset, multiplier, channel_map[channel]

            #graph.SetLineColor(color)
            fcn_string = "(y - %s)*%s + %s" % ( baseline, multiplier, offset)
            #print "fcn_string:", fcn_string
            fcn = ROOT.TF2("fcn",fcn_string)
            #graph.Apply(fcn)
            
	    graph_new=ROOT.TGraph()
            # convert x axis to microseconds
	    for i_point in xrange(graph.GetN()):
                x = graph.GetX()[i_point]
                y = graph.GetY()[i_point]
                graph_new.SetPoint(i_point, x/sampling_freq_Hz*1e6, ((y-baseline)*multiplier+offset))
                #graph.SetPoint(i_point, x/sampling_freq_Hz*1e6, y)

            graph_new.SetLineWidth(1)
            graph_new.SetLineColor(color)
	    if(channel==bot_pmt_channel): graph_new.SetLineColor(1)
            
	    #graph_new.Draw("APL")
            if(channel<19):
                graph_new.GetXaxis().SetRangeUser(30,40)
		graph_total.Add(graph_new)
	        #graph_new.Draw("L,SAME")

            leg_entry= "%s" % channel_map[channel]
	    	
            i_legend_entry = channel
	    if(i_legend_entry<19):
                legend_entries[i_legend_entry] = leg_entry

            # end loop over channels

        print "---------> entry %i of %i (%.1f percent), %i plots so far" % (
            i_entry/nchannels,  # current event
            n_entries/nchannels,  # n events in file                       
            100.0*i_entry/n_entries, # percent done
            n_plots,
        ) 

        for i in xrange(nchannels):
        #for i in xrange(17):
            index = (nchannels-1) - i # fill from top down
            if legend_entries[index] != "":
                legend.AddEntry( hists[index], legend_entries[index], "f")

        graph_total.GetXaxis().SetRangeUser(30,40)
        graph_total.Draw()

        # line to show trigger time
        line = ROOT.TLine(trigger_time, y_min, trigger_time,y_max)
        line.SetLineWidth(2)
        line.SetLineStyle(7)
        #line.Draw()

        max_drift_time = struck_analysis_parameters.max_drift_time
        line1 = ROOT.TLine(trigger_time+max_drift_time, y_min, trigger_time+max_drift_time,y_max)
        line1.SetLineWidth(2)
        line1.SetLineStyle(7)
        #line1.Draw()

        legend.Draw()
        canvas.Update()
        n_plots += 1
        print "--> Event %i, %i of %i plots so far" % (i_entry/nchannels, n_plots, n_plots_total)

        if not ROOT.gROOT.IsBatch(): 

            val = raw_input("--> entry %i | enter to continue (q to quit, p to print, or entry number) " % i_entry)

            if val == 'q': sys.exit()
            elif val == 'p':
                canvas.Update()
                canvas.Print("%s_event_%i.pdf" % (basename, n_plots-1))
            try:
                i_entry = int(val)
                print "getting entry %i" % i_entry
                continue
            except: 
                pass

        else: # batch mode
            # if we run in batch mode, print a multi-page canvas
            #plot_name = "EventsWithChargeAbove%ikeV_6thLXe.pdf" % threshold
            
            plot_name = "%s_%s.pdf" % (basename, units) # the pdf file name

            if n_plots == 1: # start the file
                canvas.Print("%s[" % plot_name)
            canvas.Print(plot_name)


        if n_plots >= n_plots_total:
            print "quitting..."
            break

        # end loop over entries
    
    if ROOT.gROOT.IsBatch(): # end multi-page pdf file
        canvas.Print("%s]" % plot_name)
    return n_plots


if __name__ == "__main__":

    n_plots_total = 1000
    #n_plots_total = 3
    n_plots_so_far = 0
    if len(sys.argv) > 1:
        for filename in sys.argv[1:]:
            n_plots_so_far += process_file(filename, n_plots_total)
    else:
        process_file(n_plots_total=n_plots_total)

