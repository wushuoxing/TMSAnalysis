import os
import sys
import time
import ROOT
from ROOT import TF1
from ROOT import gROOT
from ROOT import *

from tms_utilities.digi_filter import digi_hp
from tms_utilities.digi_filter import digi_lp
from tms_utilities.useful_function_shapes import fitShapedConst

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

    print('Parsing channel {}'.format(ch_num))

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
#    fit_init_list.append( float(125.*1) )  # differential time
#    fit_init_list.append( float(31.25*2))  # integration time
    fit_init_list.append( noise_baseline)  # noise floor


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

         #fitFunc.SetParameter(0,np.max(signal_wfm)*100)

         fit_init_list[0] = np.max(signal_wfm)*8.
         fit_init = tuple( fit_init_list )
         print('Init: [{},{},{}]'.format(fit_init_list[0],fit_init_list[1],fit_init_list[2]))

         fit_start = 1600#np.argmax(signal_wfm)+2000-150
         fit_end   = 4000#np.argmax(signal_wfm)+2000+400

         fit_mask = range(fit_start,fit_end+1)
         
         fit_x = np.arange(float(fit_start),float(fit_end+1))
         fit_y = channel_wfm_fft_hp_lp[fit_mask]  
#         print('fit_x: {}'.format(fit_x))
#         print('fit_y: {}'.format(fit_y))

         try:
           popt, pcov = opt.curve_fit( fitShapedConst, fit_x, fit_y, p0=fit_init_list )
         except RuntimeError as e:
           print(e)
           popt = np.zeros(len(fit_init))
           pcov = np.zeros(len(fit_init))
         print('Fit: [{},{},{}]'.format(popt[0],popt[1],popt[2]))
#         print('Popt {}'.format(popt))
#         print('Pcov {}'.format(pcov))
#         print('pcov.shape {}'.format(pcov.shape))
#         if pcov.shape
         try:
           perr = np.sqrt( np.diag( pcov ) )
         except ValueError:
           perr = np.zeros(len(popt))


         y_from_fit = fitShapedConst( fit_x, popt[0], popt[1], popt[2] )
         redChi2 = np.sum( (y_from_fit - fit_y)**2 ) / ( len(fit_x) - len(popt) )

         deltaT=0.1*(popt[2]-np.min(y_from_fit))*60
         startT= popt[1]

         if(startT>2000) and (startT-2000)*0.052<20 and redChi2 <5 :
                output_series['Drift {}'.format(ch_num)] = (startT-2000)*0.052
                output_series['Energy {}'.format(ch_num)] = signal_amp
                output_series['Drift Err {}'.format(ch_num)] = perr[1]*0.052
                output_series['Fit RedChi2 {}'.format(ch_num)] = redChi2
         else: 
                output_series['Drift {}'.format(ch_num)] = -1.
                output_series['Energy {}'.format(ch_num)] = -1.
                output_series['Drift Err {}'.format(ch_num)] = -1.
                output_series['Fit RedChi2 {}'.format(ch_num)] = redChi2

    else:
         output_series['Drift {}'.format(ch_num)] = -1.
         output_series['Drift Err {}'.format(ch_num)] = -1.
         output_series['Energy {}'.format(ch_num)] = -1.
         output_series['Fit RedChi2 {}'.format(ch_num)] = -1.

    if save_waveform:
       output_series['TMS Data {}'.format(ch_num)] = raw_waveform
       output_series['TMS Denoised Data {}'.format(ch_num)] = channel_wfm_fft_hp_lp

    return output_series
