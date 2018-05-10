import numpy as np
import math  as mt

def statistics(_nbin,_x,_y,_min,_max):
    x=_x
    y=_y
    nbin = _nbin
    
    min_x               = _min             # min value for binning range 
    max_x               = _max             # max value for binning range
    ymatrix             = []               # values of 'y' for each bin. This matrix is use to obtain the median value per bin
    sum_mean            = np.zeros(nbin)   # sum of values inside a bin
    binx                = np.array([])     # values for the x-axis, mean value of the binning interval
    medianx             = np.array([])     # median for each bin
    elements_per_bin    = np.zeros(nbin)   # number of elements per bin
    
    # generation of ymatrix
    for i in range(nbin):
        ymatrix.append([])

    # bin length    
    deltax = (max_x-min_x)/float(nbin)
    print deltax   
 
    for i in range(len(x)):
        bin_index = int(mt.floor((x[i] - min_x)/deltax))
        if bin_index>=nbin:
            continue
        elements_per_bin[bin_index] += 1
        ymatrix[bin_index].append(y[i])
        sum_mean[bin_index]=sum_mean[bin_index]+y[i]
    
    for i in range(nbin):
        binx        = np.append(binx, min_x+(2*i+1)*deltax/2.)
        medianx     = np.append(medianx, np.median(ymatrix[i]))
        sum_mean[i] = sum_mean[i]/elements_per_bin[i]

    sigma=np.array([])
    for b in range(nbin):
        sum_sigma=0.
        for val in range(len(ymatrix[b])):
            sum_sigma= sum_sigma+(ymatrix[b][val]-sum_mean[b])**2
        if elements_per_bin[i] < 30:
            sum_sigma = (sum_sigma/elements_per_bin[b])**0.5
        else:
            sum_sigma = (sum_sigma/(elements_per_bin[b]-1))**0.5
        sigma = np.append(sigma,sum_sigma)

    return binx, medianx, sum_mean, sigma 


def internal_stat(_var):
    return np.mean(_var),np.median(_var),np.std(_var),np.percentile(_var,25),np.percentile(_var,75)

def internal_stat_3(_var):
    return np.mean(_var),np.median(_var),np.std(_var)




