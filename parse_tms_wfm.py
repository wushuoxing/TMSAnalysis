import os
import sys
import time
import ROOT
from ROOT import TF1
from ROOT import gROOT
from ROOT import *

from tms_utilities.digi_filter import digi_hp
from tms_utilities.digi_filter import digi_lp
from tms_utilities.useful_function_shapes import fitShaped

import numpy as np
from scipy.fftpack import fft
import scipy.optimize as opt
from numpy import mean, sqrt, square, arange
import pandas as pd
from iminuit import Minuit


#gROOT.ProcessLine('.L fitShaped.C+')

#from ROOT import fitShaped

#def print_fit_info(fit_result, fit_duration):
#        
#    chi2 = fit_result.Chi2()
#    prob = fit_result.Prob()
#    ndf = fit_result.Ndf()
#    status = fit_result.Status()
#    #print "fit results:"
#    #print "\tchi2: %.2f" % chi2
#    #print "\tn dof", ndf
#    #print "\tchi2 / dof: %.2f" % (chi2/ndf)
#    #print "\tprob", prob
#    #print "\tstatus", status
#    #print "\t%.1f seconds" % fit_duration

def ParseTMSXwireWfm( raw_waveform, ch_num, save_waveform = False ):

    output_series = pd.Series()

    # Denoising/filtering happens here.
    channel_fft = np.fft.rfft( raw_waveform )
    channel_fft[100:]=0
    channel_wfm_fft=np.fft.irfft(channel_fft)
    channel_wfm_fft_hp=digi_hp(channel_wfm_fft)
    channel_wfm_fft_hp_lp=digi_lp(channel_wfm_fft_hp)

    graph_fftcut=ROOT.TGraph()
    for isamp in xrange( len(raw_waveform) ):
        x = isamp
        #graph_fftcut.SetPoint(isamp,x,channel_wfm_lp[isamp])
        graph_fftcut.SetPoint(isamp,x,channel_wfm_fft_hp_lp[isamp])

    noise_wfm = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(700,1700,1)])
    signal_wfm = np.array([graph_fftcut.GetY()[isamp] for isamp in xrange(2000,2400,1)])

    noise_baseline = np.mean(noise_wfm)
    noise_rms=np.sqrt(noise_wfm.dot(noise_wfm)/noise_wfm.size)
           
    signal_amp = np.max(signal_wfm)

    # fit_init is a tuple 
    fit_init_list = list()
    fit_init_list.append( 1. )             # integral
    fit_init_list.append( float(2000.0) )  # starting time
    fit_init_list.append( float(125.*1) )  # differential time
    fit_init_list.append( float(31.25*2))  # integration time
    fit_init_list.append( noise_baseline)  # noise floor
    fit_init = tuple( fit_init_list )


#    fitFunc = TF1("fitShapedTT",fitShaped, 1500, 3500, 5)

#    fitFunc.SetParName(0, "integral")
#    fitFunc.SetParName(1, "starting time")
#    fitFunc.SetParName(2, "differential time")
#    fitFunc.SetParName(3, "integration time")
#    fitFunc.SetParName(4, "noise floor")

#    fitFunc.SetParameter(1, float(2000.0))
#    fitFunc.SetParameter(2, float(125.0*1))
#    fitFunc.SetParameter(3, float(31.25*2))
#    fitFunc.SetParameter(4, noise_baseline)

    if (signal_amp-noise_baseline)>4*noise_rms and (np.argmax(signal_wfm)+2000)>2000 :

         fitFunc.SetParameter(0,np.max(signal_wfm)*100)

         fit_start = np.argmax(signal_wfm)+2000-150
         fit_end   = np.argmax(signal_wfm)+2000+400

         fit_mask = range(fit_start,fit_end+1)
         
         fit_x = np.arange(float(fit_start),float(fit_end+1))
         fit_y = channel_wfm_fft_hp_lp[fit_mask]  

        
         popt, pcov = opt.curve_fit( fitShaped, fit_x, fit_y, p0=fit_init )
         perr = np.sqrt( np.diag( pcov ) )

         y_from_fit = fitShaped( fit_x, popt[0], popt[1], popt[2], popt[3], popt[4] )
         redChi2 = np.sum( (y_from_fit - fit_y)**2 ) / ( len(fit_x) - len(popt) )

         deltaT=0.1*(popt[4]-fitFunc.GetMinimum())*60
         startT= popt[1]

         if(startT>2000) and (startT-2000)*0.052<20 and fit_result.Chi2()/fit_result.Ndf() <5 :
                output_series['Drift'] = (startT-2000)*0.052
                output_series['Energy'] = signal_amp
                output_series['Drift Err'] = perr[1]*0.052
                output_series['Fit RedChi2'] = redChi2
         else: 
                output_series['Drift'] = -1.
                output_series['Energy'] = -1.
                output_series['Drift Err'] = -1.
                output_series['Fit RedChi2'] = redChi2

    else:
         output_series['Drift'] = -1.
         output_series['Drift Err'] = -1.
         output_series['Energy'] = -1.
         output_series['Fit RedChi2'] = -1.

    if save_waveform:
       output_series['Data'] = raw_waveform
       output_series['Denoised Data'] = channel_wfm_fft_hp_lp

    return output_series
