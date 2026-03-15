import numpy as np
import matplotlib.pyplot as plt 
from scipy.io import loadmat

def ease_plot(XLag, YLag, start, end, title, figname) : 
    
    """Function that plots the Sea Ice concentrations on EASEgrid
    for 1 specific year. 
    
    Inputs: 
        XLag, YLag: position vectors on ease_grid
        Start, end: start and end week
        title : title of each Figure
        figname: figname name when saved
    
    Returns:
        Return n figures with the concentration for the specific week/date

    """
    
    bathyEASE = loadmat('./Grid/bathyEASE.mat')
    Land = np.zeros((361, 361)) * np.nan 
    bathyEASE = bathyEASE['bathyEASE']
    Land[np.where(bathyEASE > 0)] = 1
    Land[np.where(bathyEASE < 0)] = 0

    
    for i in range(start,end) : 
        plt.figure()
        plt.pcolor(Land, cmap = 'winter')
        plt.gca().invert_yaxis()
        plt.scatter(YLag[i], XLag[i], color = 'black', marker = '.', alpha = 0.2, s = 1)
        plt.title(title.format(i+1))
        plt.savefig(figname.format(i))
        # plt.show()
        
        
def add_land(EASE_path) : 
    
    bathyEASE = loadmat(EASE_path)
    Land = np.zeros((361, 361)) * np.nan

    bathyEASE = bathyEASE['bathyEASE']

    Land[np.where(bathyEASE > 0)] = 1
    # Land[np.where(bathyEASE < 0)] = 0
    
    return Land