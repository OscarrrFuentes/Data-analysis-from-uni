'''
plot(galaxies)

    Parameters
    ----------
    galaxies : A list of strings which are the NGC number of the galaxies
    
    colours (optional) : A list of strings which are the letters denoting 
    what colour marker you'd like for each galaxy

    Returns
    -------
    Two plots of: 
     - log(IR flux density) against log(UV flux density)
     - log(IR flux density/UVflux density) against galactocentric radius 
     - The Pearson Correlation Coefficient either in the graph if no colour 
    specified, or printed in the console

Must have files saved in the format:
    
    3938_emission_npix.txt
    3938_Stdev_mean.txt
    
which will both have the columns:
    
    Right Ascension, Declination, Sum UV, Npix UV, Sum IR, Npix IR, Galactocentric Radius
    Mean UV, StDev UV, Mean IR, StDev IR
    
respectively, where emission_npix are the data for the bright spots 
and Stdev_mean are the data for the background
'''


import numpy as np
import matplotlib.pyplot as plt

def corr_coeff(xarr, yarr):
    '''
    

    Parameters
    ----------
    xarr : A numpy array withall the x-values
    yarr : A numpy array with the corresponding y-values

    Returns
    -------
    The Pearson Correlation Coefficient
    '''
    
    
    xmean = sum(xarr)/len(xarr)
    ymean = sum(yarr)/len(yarr)
    top = sum((xarr-xmean)*(yarr-ymean))
    bottom = np.sqrt(sum((xarr-xmean)**2)*sum((yarr-ymean)**2))
    return top/bottom

def data(galaxy):

    data = np.genfromtxt(galaxy+'_emission_npix.txt', dtype = 'float', delimiter = ',', skip_header = 1)
    stdevs = np.genfromtxt(galaxy+'_Stdev_mean.txt', dtype = 'float', delimiter = ',', skip_header = 1)

    UV_bckg_mean = sum(stdevs[:,0])/len(stdevs[:,0])
    IR_bckg_mean = sum(stdevs[:,2])/len(stdevs[:,2])
    data[:,2] -= UV_bckg_mean*data[:,3]
    data[:,4] -= IR_bckg_mean*data[:,5]
    
    UV_unc_perpix = (sum(stdevs[:,1])/len(stdevs[:,1]))/UV_bckg_mean
    IR_unc_perpix = (sum(stdevs[:,3])/len(stdevs[:,3]))/IR_bckg_mean
    
    
    UV_perc_unc = (np.sqrt(data[:,3]))*UV_unc_perpix
    IR_perc_unc = (np.sqrt(data[:,5]))*IR_unc_perpix
    
    data[:,2] *= (3.365*(10**(-5)))/(np.pi*5**2)
    UV_uncertainty = np.sqrt((UV_perc_unc**2)+(0.03**2))
    
    data[:,4] *= (10**6)/(4.25*10**10)
    IR_uncertainty = np.sqrt((IR_perc_unc**2)+(0.04**2))
    
    log_UV_unc = (1/np.log(10))*UV_uncertainty
    log_IR_unc = (1/np.log(10))*IR_uncertainty
    
    return [data[:,2], data[:,4], data[:,6], log_UV_unc, log_IR_unc]
                                    

def plot(galaxies, colours = None):
    '''

    Parameters
    ----------
    galaxies : A list of strings which are the NGC number of the galaxies
    
    colours (optional) : A list of strings which are the letters denoting 
    what colour marker you'd like for each galaxy

    Returns
    -------
    Two plots of: 
     - log(IR flux density) against log(UV flux density)
     - log(IR flux density/UVflux density) against galactocentric radius 
     - The Pearson Correlation Coefficient either in the graph if no colour 
    specified, or printed in the console

    '''
    
    x, y, r = np.array([]), np.array([]), np.array([])
    title = ''
    for i in galaxies:
        x, y, r = np.append(x, data(i)[0]), np.append(y, data(i)[1]), np.append(r, data(i)[2])
        title += 'NGC '+i+'  '
    r1 = 'r = '+str(corr_coeff(np.log10(x), np.log10(y)))
    r2 = 'r = '+str(corr_coeff(r, np.log10(y/x)))
    
    count = 0
    for i in galaxies:
        x = data(i)[0]
        y = data(i)[1]
        if colours:
            plt.plot(np.log10(x), np.log10(y), 'x'+colours[count], label = 'NGC '+i)
        else:
            plt.plot(np.log10(x), np.log10(y), 'xk')
        count += 1
    plt.title(title)
    plt.xlabel('log(UV flux density)')
    plt.ylabel('log(IR flux density)')
    if colours:
        plt.legend()
        print('First graph: '+r1)
    else:
        plt.legend([r1])
    plt.show()  
    
    
    count = 0
    for i in galaxies:
        x = data(i)[2]
        y = data(i)[1]/data(i)[0]
        if colours:
            plt.plot(x, np.log10(y), 'x'+colours[count], label = 'NGC '+i)
        else:
            plt.plot(x, np.log10(y), 'xk')
        count += 1
    plt.title(title)
    plt.xlabel('galactocentric radius (arcsec)')
    plt.ylabel('log(IR flux density/UV flux density)')
    if colours:
        plt.legend()
        print('Second graph: '+r2)
    else:
        plt.legend([r1])
    plt.show()

    

# plt.errorbar(np.log10(data[:,2]), np.log10(data[:,4]), log_IR_unc, log_UV_unc, fmt = 'xk')
# plt.xlabel('log(UV flux density)')
# plt.ylabel('log(IR flux density)')
# plt.show()
