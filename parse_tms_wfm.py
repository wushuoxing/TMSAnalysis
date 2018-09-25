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
    
#gROOT.ProcessLine('.L fitRising.C+')
gROOT.ProcessLine('.L fitShaped.C+')

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

def process_file(filename):

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
    #n_entries = 50000
    print "%i entries" % n_entries
    
    i_entry = 0
    while i_entry < n_entries:
    #for i_entry in xrange(n_entries):
        if(i_entry%16==0):
            ch_x_array = []
            energy_x_array = []
            time_x_array = []
            time_err_x_array = []
            ch_y_array = []
            time_y_array = []
            time_err_y_array = []
 
        #tree.GetEntry(i_entry)

        n_channels=16
          
        for i_channel in xrange(n_channels):
            tree.GetEntry(i_entry)

            channel =  tree.HitTree.GetChannel()
            gate_size = tree.HitTree.GetGateCount()

            print "entry %i, channel %i" % (i_entry, channel)

            # Channel 6 is noisy, so we ignore it for now. 
            # (May have been destroyed by a discharge or something)
            if channel==6:
                i_entry = i_entry+1
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
            legend.AddEntry(graph, "event %i, channel %i" % (i_entry/16+1,channel), "pl")
            
            channel_wfm = np.array([graph.GetY()[isamp] for isamp in xrange(graph.GetN())])
 
            channel_fft = np.fft.rfft(channel_wfm)
             
            graph_fft=ROOT.TGraph() 
            
            for isamp in xrange(len(channel_fft)):
                graph_fft.SetPoint(isamp,2.32*isamp,np.abs(channel_fft[isamp]))
          
            #c1.Clear()
            #c2.Clear()
            #c3.Clear()

            #c1.cd()
            #graph.Draw()

            #legend.Draw()

            #c2.cd()
 
            #c2.SetLogy(1)
            #graph_fft.GetYaxis().SetRangeUser(-2000000,3000000)

            graph_fft.GetYaxis().SetRangeUser(50,30000000)
            graph_fft.GetYaxis().SetTitle("Magnitude (arbtuary unit)")
            graph_fft.GetXaxis().SetTitle("Frequency (kHz)")
            graph_fft.GetXaxis().CenterTitle(1)
            graph_fft.GetYaxis().CenterTitle(1)
            #graph_fft.Draw()

            #c3.cd()

            #channel_fft[1]=0

            channel_fft[100:]=0

            channel_wfm_fft=np.fft.irfft(channel_fft)
                    
            channel_wfm_fft_hp=digi_hp(channel_wfm_fft)

            channel_wfm_fft_hp_lp=digi_lp(channel_wfm_fft_hp)

            graph_fftcut=ROOT.TGraph()

            graph_fftcut.SetTitle("")
            graph_fftcut.GetXaxis().SetTitle("Samples [16 ns]")
            graph_fftcut.GetYaxis().SetTitle("ADC count")
            graph_fftcut.GetXaxis().CenterTitle(1)
            graph_fftcut.GetYaxis().CenterTitle(1)
            graph_fftcut.GetYaxis().SetTitleOffset(1.4)
            
            #channel_wfm_fft2 = single_pole_hp( channel_wfm_fft, 0.99 )   # Brian's recursive filter code
            #channel_wfm_lp = single_pole_lp( channel_wfm, 0.985 )  # Brian's recursive filter code
            #channel_wfm_lp = maw( channel_wfm_lp, 50 )  # Brian's recursive filter code
            #channel_wfm_fft2 = single_pole_hp( channel_wfm_lp, 0.985 ) # Brian's recursive filter code

            #graph_fftcut.GetYaxis().SetLimits(graph.GetHistogram().GetMinimum(),graph.GetHistogram().GetMaximum())
            #graph_fftcut.GetYaxis().SetRangeUser(graph.GetHistogram().GetMinimum(),graph.GetHistogram().GetMaximum())
 
            for isamp in xrange(graph.GetN()):
                x = graph.GetX()[isamp]
                #graph_fftcut.SetPoint(isamp,x,channel_wfm_lp[isamp])
                graph_fftcut.SetPoint(isamp,x,channel_wfm_fft_hp_lp[isamp])

            #graph_fftcut.GetYaxis().SetLimits(graph.GetHistogram().GetMinimum(),graph.GetHistogram().GetMaximum())
            #graph_fftcut.GetYaxis().SetRangeUser(graph.GetHistogram().GetMinimum(),graph.GetHistogram().GetMaximum())

            graph_fftcut.SetTitle("")
            graph_fftcut.GetXaxis().SetTitle("Samples [40 ns]")
            graph_fftcut.GetYaxis().SetTitle("ADC count")
            graph_fftcut.GetXaxis().CenterTitle(1)
            graph_fftcut.GetYaxis().CenterTitle(1)
            graph_fftcut.GetYaxis().SetTitleOffset(1.4)
            graph_fftcut.GetXaxis().SetRangeUser(1500,3000)

            #graph_fftcut.Draw("APL")

            rms_wfm = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(8500,9000,1)])
 
            rms_wfm=rms_wfm-np.mean(rms_wfm)
            
            rms_value=np.sqrt(rms_wfm.dot(rms_wfm)/rms_wfm.size)

            legend1=ROOT.TLegend(0.4, 0.75, 0.9, 0.88)
            legend1.SetTextSize(0.04)
            legend1.AddEntry(graph_fftcut, "RMS value (Ne-) is %f" % (rms_value*240), "pl")
            #legend1.Draw()

            legend2=ROOT.TLegend(0.4, 0.88, 0.9, 0.99)
            legend2.SetTextSize(0.04)
            legend2.AddEntry(graph_fftcut, "Event %i, Channel %i" % (i_entry/16+1,channel), "pl")
            #legend2.Draw()

            #graph_dummy=ROOT.TGraph()

            #for isamp in xrange(500,9500,1):
            #    graph_dummy.SetPoint(isamp,graph_fftcut.GetX()[isamp],graph_fftcut.GetY()[isamp])

            #print "RMS value (Ne-) is: ", rms_value*240

            noise_wfm = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(700,1700,1)])

            signal_wfm = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(2000,2400,1)])

            noise_baseline = np.mean(noise_wfm)
            noise_rms=np.sqrt(noise_wfm.dot(noise_wfm)/noise_wfm.size)
           
            # Channels 0-7 are the collection wires (positive pulses), 8-15 are induction wires (negative) 
            if channel<8 :
                signal_amp = np.max(signal_wfm)
            else :
                signal_amp = np.min(signal_wfm)

            fitFunc = TF1("fitShapedTT",fitShaped, 1500, 3500, 5)

            fitFunc.SetParName(0, "integral")
            fitFunc.SetParName(1, "starting time")
            fitFunc.SetParName(2, "differential time")
            fitFunc.SetParName(3, "integration time")
            fitFunc.SetParName(4, "noise floor")

            fitFunc.SetParameter(1, float(2000.0))
            fitFunc.SetParameter(2, float(125.0*1))
            fitFunc.SetParameter(3, float(31.25*2))
            fitFunc.SetParameter(4, noise_baseline)

            if ((channel<8) and (signal_amp-noise_baseline)>4*noise_rms) or ((channel>7) and (noise_baseline-signal_amp) > 3*noise_rms) :
                if channel==6 :    
                    continue
                elif channel<8 and (np.argmax(signal_wfm)+2000)>2000 : 
                    fitFunc.SetParameter(0,np.max(signal_wfm)*100)
                    fit_result = graph_fftcut.Fit(fitFunc,"QS","",np.argmax(signal_wfm)+2000-150,np.argmax(signal_wfm)+2000+400)
                    #fit_result = graph_fftcut.Fit(fitFunc,"QS","",np.argmax(signal_wfm)+2000-250,np.argmax(signal_wfm)+2000+500)
                    deltaT=0.1*(fitFunc.GetParameter(4)-fitFunc.GetMinimum())*60
                    startT=fitFunc.GetParameter(1)
                    h_tau1.Fill(fitFunc.GetParameter(2))
                    h_tau2.Fill(fitFunc.GetParameter(3))
		    if(fitFunc.GetParameter(3)!=0):
		        h_ratio1.Fill(fitFunc.GetParameter(3)/fitFunc.GetParameter(2))
                    #h_tau1.Fill(fitFunc.GetParameter(2))
                    #h_tau2.Fill(fitFunc.GetParameter(3))
                    if(fit_result.Ndf()!=0) :
                        goodness_fit_col=fit_result.Chi2()/fit_result.Ndf()
                        h_chi2_col_wfm.Fill(goodness_fit_col)
                    #if((deltaT+startT)>2000) and fit_result.Chi2()/fit_result.Ndf() <20:
                    #if(startT>2000) and (startT-2000)*0.0832<20 and fit_result.Chi2()/fit_result.Ndf() <5 and fitFunc.GetParameter(3)/fitFunc.GetParameter(2) < 0.9 and fitFunc.GetParameter(3)/fitFunc.GetParameter(2) > 0.8 :
                    if(startT>2000) and (startT-2000)*0.052<20 and fit_result.Chi2()/fit_result.Ndf() <5 :
                        ch_x_array.append(channel*9)
                        #energy_x_array.append(fitFunc.GetMaximum())
                        energy_x_array.append(signal_amp)
                        #time_x_array.append((deltaT+startT-2000)*0.0832)
                        time_x_array.append((startT-2000)*0.052) # 0.052 is the drift velocity times sampling period
                        time_err_x_array.append(fitFunc.GetParError(1))
                        #h_tau1.Fill(fitFunc.GetParameter(2))
                        #h_tau2.Fill(fitFunc.GetParameter(3))
			#if(fitFunc.GetParameter(3)!=0):
			#    h_ratio1.Fill(fitFunc.GetParameter(3)/fitFunc.GetParameter(2))
                    #ch_x_array.append(channel*9)
                    '''
                    risesample=0
                    for i in xrange(0,np.argmax(signal_wfm)):
                        if((signal_wfm[np.argmax(signal_wfm)-i]-noise_baseline)<5*noise_rms):
                            risesample=i
                            break
                    if (np.argmax(signal_wfm)-risesample)*0.0832<15.5:
                        ch_x_array.append(channel*9)
                        energy_x_array.append(signal_amp)
                        #time_x_array.append((np.argmax(signal_wfm)+1800-2000)*0.0832)
                        time_x_array.append((np.argmax(signal_wfm)-risesample)*0.0832)
                        time_err_x_array.append(0)
                     '''

                elif channel>7 and (np.argmin(signal_wfm)+2000)>2000 :
                    fitFunc.SetParameter(0,np.min(signal_wfm)*100)
                    fit_result = graph_fftcut.Fit(fitFunc,"QS","",np.argmin(signal_wfm)+2000-150,np.argmin(signal_wfm)+2000+400)
                    #fit_result = graph_fftcut.Fit(fitFunc,"QS","",np.argmin(signal_wfm)+2000-250,np.argmin(signal_wfm)+2000+500)
                    deltaT=0.1*(fitFunc.GetMaximum()-fitFunc.GetParameter(4))*60
                    startT=fitFunc.GetParameter(1)
                    h_tau3.Fill(fitFunc.GetParameter(2))
                    h_tau4.Fill(fitFunc.GetParameter(3))
		    if(fitFunc.GetParameter(3)!=0):
		        h_ratio2.Fill(fitFunc.GetParameter(3)/fitFunc.GetParameter(2))
                    if(fit_result.Ndf()!=0) :
                        goodness_fit_ind=fit_result.Chi2()/fit_result.Ndf()
                        h_chi2_ind_wfm.Fill(goodness_fit_ind)
                    #if((deltaT+startT)>2000) and fit_result.Chi2()/fit_result.Ndf() <20:
                    #if(startT>2000) and (startT-2000)*0.0832<20 and fit_result.Chi2()/fit_result.Ndf() <5 and fitFunc.GetParameter(3)/fitFunc.GetParameter(2) < 0.9 and fitFunc.GetParameter(3)/fitFunc.GetParameter(2) > 0.8:
                    if(startT>2000) and (startT-2000)*0.052<20 and fit_result.Chi2()/fit_result.Ndf() <5 :
                        ch_y_array.append((channel-8)*9)
                        #time_y_array.append((deltaT+startT-2000)*0.0832)
                        time_y_array.append((startT-2000)*0.052)
                        time_err_y_array.append(fitFunc.GetParError(1))
                        #h_tau3.Fill(fitFunc.GetParameter(2))
                        #h_tau4.Fill(fitFunc.GetParameter(3))
			#if(fitFunc.GetParameter(3)!=0):
			#    h_ratio2.Fill(fitFunc.GetParameter(3)/fitFunc.GetParameter(2))
                    #ch_y_array.append((channel-8)*9)
                    '''
                    risesample=0
                    for i in xrange(0,np.argmin(signal_wfm)):
                        if((noise_baseline-signal_wfm[np.argmin(signal_wfm)-i])<3*noise_rms):
                            risesample=i
                            break
                    if (np.argmin(signal_wfm)-risesample)*0.0832<15.5 :
                        ch_y_array.append((channel-8)*9)
                        #time_y_array.append((np.argmin(signal_wfm)+1800-2000)*0.0832)
                        time_y_array.append((np.argmin(signal_wfm)-risesample+2000-2000)*0.0832)
                        time_err_y_array.append(0)
                    ''' 
            #c3.Update()
            i_entry += 1
        
###############################################################################
# Below this point is the track-fitting code


        #if(i_entry%16==0):   
        if((i_channel+1)%16==0):   
            print "event %i" % (i_entry/16+1)
            print "length of X wire signals %i: " % len(ch_x_array) 
            print "length of Y wire signals %i: " % len(ch_y_array) 
            #if len(ch_x_array)>3 :
            if (len(ch_x_array)>1) and (len(ch_y_array)>2) :
                g_track_col = ROOT.TGraphErrors()
                g_track_ind = ROOT.TGraphErrors()
                for i in xrange(len(ch_x_array)):
                    g_track_col.SetPoint(i,ch_x_array[i],time_x_array[i])
                    g_track_col.SetPointError(i,0,time_err_x_array[i])
                for i in xrange(len(ch_y_array)):
                    g_track_ind.SetPoint(i,ch_y_array[i],time_y_array[i])
                    g_track_ind.SetPointError(i,0,time_err_y_array[i])
                #c7.cd()    
                g_track_col.SetTitle("")
                g_track_col.GetYaxis().SetTitle("Z [mm]")
                g_track_col.GetXaxis().SetTitle("X [mm]")
                g_track_col.GetXaxis().CenterTitle(1)
                g_track_col.GetYaxis().CenterTitle(1)
                g_track_col.GetYaxis().SetTitleOffset(1.4)
                g_track_col.SetMarkerStyle(20)
                g_track_col.SetMarkerSize(1)
                g_track_col.SetMarkerColor(2)
                g_track_col.SetLineColor(2)
                #g_track_col.Draw("AP")
                #c8.cd()
                g_track_ind.SetTitle("")
                #g_track_ind.GetYaxis().SetTitle("Samples [16 ns]")
                #g_track_ind.GetXaxis().SetTitle("Channel ID (Y)")
                g_track_ind.GetYaxis().SetTitle("Z [mm]")
                g_track_ind.GetXaxis().SetTitle("Y [mm]")
                g_track_ind.GetXaxis().CenterTitle(1)
                g_track_ind.GetYaxis().CenterTitle(1)
                g_track_ind.GetYaxis().SetTitleOffset(1.4)
                g_track_ind.SetMarkerStyle(20)
                g_track_ind.SetMarkerSize(1)
                g_track_ind.SetMarkerColor(2)
                g_track_ind.SetLineColor(2)
                #g_track_ind.Draw("AP")
                #c7.Update()
                #c8.Update()
                #time.sleep(2)
                fitLineInd = TF1("FLInd","[0]*x+[1]",0,72)      
                fitLineCol = TF1("FLCol","[0]*x+[1]",0,72)      
                fitLineCol.SetParameter(0,((time_x_array[len(ch_x_array)-1]-time_x_array[0])/len(ch_x_array)))
                fitLineCol.SetParameter(1,np.mean(time_x_array))
                fitLineInd.SetParameter(0,((time_y_array[len(ch_y_array)-1]-time_y_array[0])/len(ch_y_array)))
                fitLineInd.SetParameter(1,np.mean(time_y_array))

                fitLineInd.SetLineColor(2)
                fitLineInd.SetLineWidth(2)
                fitLineCol.SetLineColor(2)
                fitLineCol.SetLineWidth(2)
                
                fitresult_ind=g_track_ind.Fit(fitLineInd,"QSR")
                fitresult_col=g_track_col.Fit(fitLineCol,"QSR")
                time.sleep(0)
                '''
                g_track_col_noerr = ROOT.TGraph()
                g_track_ind_noerr = ROOT.TGraph()
                for i in xrange(len(ch_x_array)):
                    g_track_col_noerr.SetPoint(i,ch_x_array[i],time_x_array[i])
                for i in xrange(len(ch_y_array)):
                    g_track_ind_noerr.SetPoint(i,ch_y_array[i],time_y_array[i])
                '''
                #mg_col.Add(g_track_col_noerr)
                #mg_ind.Add(g_track_ind_noerr)
                
                #c7.Update()
                #c8.Update()

                goodness_fit_ind=100
                goodness_fit_col=100
                if(fitresult_ind.Ndf()!=0) :
                    goodness_fit_ind=fitresult_ind.Chi2()/fitresult_ind.Ndf()
                    h_chi2_ind.Fill(goodness_fit_ind)
                    print "Ind chi2/ndf: %f " % (fitresult_ind.Chi2()/fitresult_ind.Ndf())
                if(fitresult_col.Ndf()!=0) :
                    goodness_fit_col=fitresult_col.Chi2()/fitresult_col.Ndf()
                    h_chi2_col.Fill(goodness_fit_col)
                    print "Col chi2/ndf: %f " % (fitresult_col.Chi2()/fitresult_col.Ndf())
                
                if fitLineInd.GetParameter(0)!=0 and (fitresult_ind.Ndf()!=0 and fitresult_ind.Chi2()/fitresult_ind.Ndf() < 8) and ((fitresult_col.Ndf()!=0 and fitresult_col.Chi2()/fitresult_col.Ndf() < 8) or fitresult_col.Ndf()==0): 
                    slopexvsy=fitLineCol.GetParameter(0)/fitLineInd.GetParameter(0)
                    constxvsy=(fitLineCol.GetParameter(1)-fitLineInd.GetParameter(1))/fitLineInd.GetParameter(0)
                    
                    print "slopexvsy: %f; constxvsy: %f" % (slopexvsy,constxvsy)

                    #xvsy = TF1("xvsy","slopexvsy*x+constxvsy",0,8)      
                    xvsy = ROOT.TGraph()
                    for i in xrange(0,720,1):
                        xvsy.SetPoint(int(i),float(i)/10.0,slopexvsy*float(i)/10.0+constxvsy)

                    #TF1("xvsy","[0]*x+[1]",0,8)
                    xvsy.GetYaxis().SetLimits(0,72)
                    xvsy.GetYaxis().SetRangeUser(0,72)
                    xvsy.GetXaxis().SetLimits(0,72)
                    xvsy.GetXaxis().SetRangeUser(0,72)
                    #xvsy.SetParameter(0,slopexvsy)
                    #xvsy.SetParameter(1,constxvsy)
                    xvsy.GetYaxis().SetTitle("Induction Ch [mm]")
                    xvsy.GetXaxis().SetTitle("Collection Ch [mm]")
                    xvsy.GetXaxis().CenterTitle(1)
                    xvsy.GetYaxis().CenterTitle(1)
                    xvsy.SetLineColor(2)
                    xvsy.SetLineWidth(2)
                    
                    #c9.cd()
                    #xvsy.Draw("AL")
                    #c9.Update()
 
                    if fitLineInd(36)>0 and fitLineInd(36)<16 and (((slopexvsy*0+constxvsy)>0 and (slopexvsy*0+constxvsy)<72) or ((slopexvsy*72+constxvsy)>0 and (slopexvsy*72+constxvsy)<72) or ((slopexvsy*0+constxvsy)<0 and (slopexvsy*72+constxvsy)>72) or ((slopexvsy*0+constxvsy)>72 and (slopexvsy*72+constxvsy)<0)) :
                        mg_xvsy.Add(xvsy)
                        for i in xrange(len(ch_x_array)):
                            term1=ch_x_array[i]/np.sqrt(ch_x_array[i]*ch_x_array[i]+(slopexvsy*ch_x_array[i]+constxvsy)*(slopexvsy*ch_x_array[i]+constxvsy)+fitLineCol(ch_x_array[i])*fitLineCol(ch_x_array[i]))
                            term2=np.sqrt(0.9*0.9+slopexvsy*slopexvsy*0.9*0.9+fitLineCol.GetParameter(0)*fitLineCol.GetParameter(0)*0.9*0.9)
                            #term2=np.sqrt((1/slopexvsy)*(1/slopexvsy)*0.9*0.9+0.9*0.9+fitLineInd.GetParameter(0)*fitLineInd.GetParameter(0)*0.9*0.9)
                            deltaX=term1*term2
                            h_dx.Fill(deltaX)
                            #energy_corrected=energy_x_array[i]*np.exp(time_x_array[i]/29.7)
                            energy_corrected=energy_x_array[i]*np.exp(time_x_array[i]/27)
                            #h_amp.Fill(energy_x_array[i]/deltaX)
                            if time_x_array[i]>0 and time_x_array[i] < 7 :
                                h_amp.Fill(energy_corrected/deltaX)
                                h_amp_nocorr.Fill(energy_x_array[i]/deltaX)
                                g_track_col_noerr = ROOT.TGraph()
                                g_track_ind_noerr = ROOT.TGraph()
                                for m in xrange(len(ch_x_array)):
                                    g_track_col_noerr.SetPoint(m,ch_x_array[m],time_x_array[m])
                                for m in xrange(len(ch_y_array)):
                                    g_track_ind_noerr.SetPoint(m,ch_y_array[m],time_y_array[m])
                                g_track_col_noerr.SetLineColor((i_channel+1)/16)
                                g_track_ind_noerr.SetLineColor((i_channel+1)/16)
                                mg_col.Add(g_track_col_noerr)
                                mg_ind.Add(g_track_ind_noerr)
				h_hit_col.Fill(len(ch_x_array))
				h_hit_ind.Fill(len(ch_y_array))

				thetaangle=57.3*np.arctan(slopexvsy)
				if thetaangle>0:
				    thetaangle=thetaangle-90
				else:
				    thetaangle=thetaangle+90
                                h_theta.Fill(thetaangle)
                                h_phi.Fill(57.3*np.arctan(fitLineCol.GetParameter(0)))
                                #mg_xvsy.Add(xvsy)
                                #p_dQdx.Fill(energy_x_array[i]/deltaX,time_x_array[i])

                                p_dQdx.Fill(time_x_array[i],energy_x_array[i]/deltaX)
                                p_dQdx_dt.Fill(time_x_array[i]/3.25,energy_x_array[i]/deltaX)
				h_dqdx_dt.Fill(time_x_array[i]/3.25,energy_x_array[i]/deltaX)
				h_time_x.Fill(time_x_array[i]/3.25)
				h_z_x.Fill(time_x_array[i])
 
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



