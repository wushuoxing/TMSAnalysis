import numpy as np
from scipy.fftpack import fft
import ROOT
from ROOT import TF1
from ROOT import gROOT
from ROOT import *

from digi_filter import digi_hp
from digi_filter import digi_lp

def graph_fft( graph_input ):
    
    channel_wfm = np.array([graph_input.GetY()[isamp] for isamp in xrange(graph_input.GetN())]) 
    channel_fft = np.fft.rfft(channel_wfm)
    channel_fft_filter_array = np.zeros_like(channel_fft)   

    sampling_freq_Hz=62.5*1.e6
    fft_freq = np.fft.rfftfreq(len(channel_wfm),d=1./sampling_freq_Hz)                                                       
    ch_low_pass = 8.0e6                                                                                                       
    fft_freq_pass = np.logical_and(fft_freq > 5e3, fft_freq < ch_low_pass)

    channel_fft_filter_array[fft_freq_pass] = channel_fft[fft_freq_pass] 
    channel_wfm_fft = np.fft.irfft(channel_fft_filter_array)

    graph_fftcut=ROOT.TGraph()  
                             
    for isamp in xrange(graph_input.GetN()):                                                                                        
	 x = graph_input.GetX()[isamp]                                                                                               
	 graph_fftcut.SetPoint(isamp,x,channel_wfm_fft[isamp])      
 
    return graph_fftcut


def graph_fft_hp_lp(graph_input):

    
    channel_wfm = np.array([graph_input.GetY()[isamp] for isamp in xrange(graph_input.GetN())]) 
    channel_fft = np.fft.rfft(channel_wfm)
    channel_fft_filter_array = np.zeros_like(channel_fft)   

    sampling_freq_Hz=62.5*1.e6
    fft_freq = np.fft.rfftfreq(len(channel_wfm),d=1./sampling_freq_Hz)                                                       
    ch_low_pass = 8.0e6                                                                                                       
    fft_freq_pass = np.logical_and(fft_freq > 5e3, fft_freq < ch_low_pass)

    channel_fft_filter_array[fft_freq_pass] = channel_fft[fft_freq_pass] 
    channel_wfm_fft = np.fft.irfft(channel_fft_filter_array)

    channel_wfm_fft_hp=digi_hp(channel_wfm_fft)
    channel_wfm_fft_hp_lp=digi_lp(channel_wfm_fft_hp)

    graph_fftcut_hp_lp=ROOT.TGraph()  
                             
    for isamp in xrange(graph_input.GetN()):                                                                                        
	 x = graph_input.GetX()[isamp]                                                                                               
	 graph_fftcut_hp_lp.SetPoint(isamp,x,channel_wfm_fft_hp_lp[isamp])      
 
    return graph_fftcut_hp_lp

