import pandas as pd
import numpy as np
from scipy import optimize as opt
from tms_utilities.useful_function_shapes import TwoExpConv
from tms_utilities.useful_function_shapes import Gaussian

########################################################################################
def ParseNaIWfm( raw_waveform, save_waveform=False ):
 
  colnames = [\
             'Baseline',\
             'Baseline RMS',\
#             'Pulse Area',\
             'Pulse Height',\
             'Pulse Time']
             
  if save_waveform:
     colnames.append('Data')

  reduced_data = pd.Series()

  if save_waveform:
     reduced_data['Data'] = raw_waveform

  reduced_data['Baseline'], reduced_data['Baseline RMS'] = GetBaseline( raw_waveform )    
  reduced_data.name = 'NaI Event'

  try:
    reduced_data['Pulse Height'],\
    reduced_data['Pulse Time'] = NaICrossCorrelationPulseFinder( raw_waveform - reduced_data['Baseline'] )
  except ValueError:
      dummy=0

  reduced_data['Data (BSS)'] = raw_waveform - reduced_data['Baseline']

  return reduced_data
    
########################################################################################
def GetBaseline( data ):

   # Baseline is computed from first 50 samples. This may need to be fine-tuned.  

   fullmean = np.mean(data[0:50])
   fullrms = np.std(data[0:50])
   mask = (data - fullmean)**2 < (6 * fullrms)**2
 
   baseline = np.mean(data[mask])
   baselineRMS = np.std(data[mask])

   #print('{} {} {} {}'.format(tempmean,tempRMS,baseline,baselineRMS))
      
   return baseline, baselineRMS


#######################################################################################
def DoubleExpConvConstNaI(x,A,mu):                                                                                        

    # These constants will need to be fine-tuned for the new digitizer

    t1 = 3.0
    t2 = 30.0
    mask = (x-mu)<0
    y = A/(t1-t2)*( np.exp(-(x-mu)/t1) - np.exp(-(x-mu)/t2) )
    y[mask] = np.zeros(len(y))[mask]
    return y

#######################################################################################
def NaICrossCorrelation(x,y):
    cross_cor = np.zeros(len(y))
    for i in range(0,len(y)):
        cross_cor[i] = np.sum( DoubleExpConvConstNaI(x,-1.,float(i))*y )

    return cross_cor*1.654311 # Empirical scaling factor so that the cross correlation 
                              # height equals the true pulse height in ADC counts.

#######################################################################################
def NaICrossCorrelationPulseFinder( data ):

    x = np.linspace(0.,len(data)-1.,len(data))
    y = data               

    cross_cor = NaICrossCorrelation(x,y)
    pulse_heights = np.array([])                
    pulse_times = np.array([])
    pulse_idxs = np.array([])
    
    pulse_idxs, pulse_heights = FindLocalMaxima( cross_cor )
    
    if len(pulse_idxs) > 1:
          #print('Pulse indexes: {}'.format(pulse_idxs))
          raise ValueError('Multiple pulses found in event. Skipping...')
          return  0,0
    if len(pulse_idxs) == 0:
        raise ValueError('No pulses found in event')
        return 0,0
        
    start = int(pulse_idxs[0]) - 1
    end = int(pulse_idxs[0]) + 4
    
    (Afit, mufit, sigfit),_ = opt.curve_fit( Gaussian,\
                                            x[start:end],\
                                            cross_cor[start:end],\
                                            p0=(pulse_heights[0],pulse_idxs[0],10.) )
    return Afit, mufit, cross_cor
    
    
#######################################################################################
def FindLocalMaxima(data):
    mask = data>3.5
    local_maxima = np.array([])
    values = np.array([])
    for i in range(1,len(data)-1):
        if data[i] > data[i-1] and data[i] > data[i+1] and mask[i]:
            local_maxima = np.append(local_maxima,i)
            values = np.append(values,data[i])
    return local_maxima, values     
                   
