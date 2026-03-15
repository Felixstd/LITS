import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matlab.engine
import time 
from maps import *
from LITS_age import *
from scipy.io import savemat, loadmat

X, Y, i, j, ActiveFlag, Melt, Age, Ageround, ActiveGeo, nWeek, TracerNum \
     = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
     
def saveLITS(FigPath, SeaIce_Tracers, name) : 
    
    """
    This is a function to save the lits ouput, it just saves many files depending on the weeks. 
    """
    
    idx_1fourth_seaIce_Tracers = np.where(SeaIce_Tracers[:, nWeek] < 520)
    SeaIce_Tracers_1fourth = SeaIce_Tracers[idx_1fourth_seaIce_Tracers[0], :]


    idx_2_seaIce_Tracers = np.where((SeaIce_Tracers[:, nWeek] >= 520) & (SeaIce_Tracers[:, nWeek] < 820))
    SeaIce_Tracers_2 = SeaIce_Tracers[idx_2_seaIce_Tracers[0], :]

    idx_3_seaIce_Tracers = np.where((SeaIce_Tracers[:, nWeek] >= 820) & (SeaIce_Tracers[:, nWeek] < 1020))
    SeaIce_Tracers_3 = SeaIce_Tracers[idx_3_seaIce_Tracers[0], :]

    idx_4_seaIce_Tracers = np.where((SeaIce_Tracers[:, nWeek] >= 1020) & (SeaIce_Tracers[:, nWeek] < 1220))
    SeaIce_Tracers_4= SeaIce_Tracers[idx_4_seaIce_Tracers[0], :]

    idx_5_seaIce_Tracers = np.where((SeaIce_Tracers[:, nWeek] >= 1220) & (SeaIce_Tracers[:, nWeek] < 1420))
    SeaIce_Tracers_5= SeaIce_Tracers[idx_5_seaIce_Tracers[0], :]

    idx_6_seaIce_Tracers = np.where((SeaIce_Tracers[:, nWeek] >= 1420) & (SeaIce_Tracers[:, nWeek] < 1620))
    SeaIce_Tracers_6 = SeaIce_Tracers[idx_6_seaIce_Tracers[0], :]

    idx_7_seaIce_Tracers = np.where((SeaIce_Tracers[:, nWeek] >= 1620) & (SeaIce_Tracers[:, nWeek] < 1820))
    SeaIce_Tracers_7 = SeaIce_Tracers[idx_7_seaIce_Tracers[0], :]

    idx_8_seaIce_Tracers = np.where(SeaIce_Tracers[:, nWeek] >= 1820)
    SeaIce_Tracers_8 = SeaIce_Tracers[idx_8_seaIce_Tracers[0], :]

    dic_8 = {'Tracers' : SeaIce_Tracers_8}
    dic_7 = {'Tracers' : SeaIce_Tracers_7}
    dic_6 = {'Tracers' : SeaIce_Tracers_6}
    dic_5 = {'Tracers' : SeaIce_Tracers_5}
    dic_4 = {'Tracers' : SeaIce_Tracers_4}
    dic_3 = {'Tracers' : SeaIce_Tracers_3}
    dic_2 = {'Tracers' : SeaIce_Tracers_2}
    dic_1fourth = {'Tracers' : SeaIce_Tracers_1fourth}

    savemat(FigPath+'SeaIce_Tracers1fourth_elim_{}.mat'.format(name), dic_1fourth)
    savemat(FigPath+'SeaIce_Tracers2_elim_{}.mat'.format(name), dic_2)
    savemat(FigPath+'SeaIce_Tracers3_elim_{}.mat'.format(name), dic_3)
    savemat(FigPath+'SeaIce_Tracers4_elim_{}.mat'.format(name), dic_4)
    savemat(FigPath+'SeaIce_Tracers5_elim_{}.mat'.format(name), dic_5)
    savemat(FigPath+'SeaIce_Tracers6_elim_{}.mat'.format(name), dic_6)
    savemat(FigPath+'SeaIce_Tracers7_elim_{}.mat'.format(name), dic_7)
    savemat(FigPath+'SeaIce_Tracers8_elim_{}.mat'.format(name), dic_8)
    
def saveAge_Grid(FigPath, Age_Grid_list) : 
    
    Age_Grid_1 = {'Grid' : Age_Grid_list[:500]}
    Age_Grid_2 = {'Grid' : Age_Grid_list[500:1000]}
    Age_Grid_3 = {'Grid' : Age_Grid_list[1000:1500]}
    Age_Grid_4 = {'Grid' : Age_Grid_list[1500:2000]}
    Age_Grid_5 = {'Grid' : Age_Grid_list[2000:]}
    
    savemat(FigPath+'Age_grid_list_1.mat', Age_Grid_1)
    savemat(FigPath+'Age_grid_list_2.mat', Age_Grid_2)
    savemat(FigPath+'Age_grid_list_3.mat', Age_Grid_3)
    savemat(FigPath+'Age_grid_list_4.mat', Age_Grid_4)
    savemat(FigPath+'Age_grid_list_5.mat', Age_Grid_5)
    
    
    
eng = matlab.engine.start_matlab()


dataDir = '/aos/home/fstdenis/LITS/'
SIC_dir = '/aos/home/fstdenis/LITS/SIC_V4/*1979_0*'
SIV_dir = '/aos/home/fstdenis/LITS/SIV/'
GridDir = '/aos/home/fstdenis/LITS/Grid/Poly_Masks/'
FigPath1 = '/storage/fstdenis/1979-2022/Bilans17_AgeGrid_september_WithoutElim/Growth/'
FigPath2 = './'


ArcOceanMask = scipy.io.loadmat('./ArcticOceanMask4IceEdge.mat')
ArcOceanMask = ArcOceanMask['ArcticOceanMask4IceEdge'] 

SIC_data, SIC_mask, SIC_file = get_files_and_data(SIC_dir)

SIC_indice, SIC_indice_x, SIC_indice_y = idx_SI_x_y(SIC_data)
SIC_indice_mask, SIC_indice_x_mask, SIC_indice_y_mask = idx_SI_x_y(SIC_mask, 0.15)


year = list(range(1979, 2022))
weeks = np.arange(1,53)


##################################### <LITS> #######################################
 
xstart = SIC_indice_x_mask[0][:]
ystart = SIC_indice_y_mask[0][:]
agestart = np.zeros(len(xstart))
ActiveFlagstart = np.zeros(len(xstart))*np.nan


startyear = 1979
startweek = 1
endyear = 1985
endweek = 1

#Main calculation with LITS
SeaIce_Tracers, Age_grid_list, SeaIce_Tracers_melted, SeaIce_melt_Beaufort, SeaIce_MYI_aniv_Beaufort, SeaIce_MYI_Melt_Arctic, SeaIce_MYI_aniversary,\
 SeaIce_MYI_import_Beaufort, SeaIce_MYI_export_Beaufort, SeaIce_Promotion_Arctic, SeaIce_Promotion_Beaufort, MYI_Extent_Beaufort, OpenWater2FYI_prom_Beaufort\
     = LITS_Age2(xstart, ystart, startyear, startweek, endyear, endweek, 0.15)


#saving data
saveLITS(FigPath1, SeaIce_Tracers, 'active')
saveLITS(FigPath1, SeaIce_Tracers_melted, 'melted')
saveAge_Grid(FigPath1, Age_grid_list)


try: 
    dic_Beaufort_Growth = {'Melt' : SeaIce_melt_Beaufort, 'Aniv' : SeaIce_MYI_aniv_Beaufort, 'Prom' : SeaIce_Promotion_Beaufort, 'OW_Prom' : OpenWater2FYI_prom_Beaufort, 'Extent' : MYI_Extent_Beaufort}
    dic_Beaufort_ImportExport = {'Import' : SeaIce_MYI_import_Beaufort, 'Export' : SeaIce_MYI_export_Beaufort}
    dic_Arctic_Growth = {'Melt' : SeaIce_MYI_Melt_Arctic, 'Aniv' : SeaIce_MYI_aniversary, 'Prom' : SeaIce_Promotion_Arctic}


    savemat(FigPath1+'Beaufort_Growth.mat', dic_Beaufort_Growth)
    savemat(FigPath1+'Beaufort_ImpExp.mat', dic_Beaufort_ImportExport)
    savemat(FigPath1+'Arctic_Beaufort.mat', dic_Arctic_Growth)

except:
    print('Saving did not work')


MaxWeek = SeaIce_Tracers[:,nWeek].max() 
Land = add_land('./Grid/bathyEASE.mat')

bathyEASE = loadmat('./Grid/bathyEASE.mat')


cmap_1 = mpl.colors.ListedColormap(['peru'])
cmap_2 = mpl.colors.ListedColormap(['midnightblue', 'cyan', 'limegreen', 'darkorange', 'red'])
cmap_3 = mpl.colors.ListedColormap(['slateblue', 'midnightblue', 'cyan', 'limegreen', 'darkorange', 'red'])
bounds_2 = [1, 2, 3, 4, 5, 6]
norm_2 = mpl.colors.BoundaryNorm(bounds_2, cmap_2.N)

        
    
for m in range(250,int(MaxWeek),1) :
    year, week = JWeek2YearWeek(m+1)
    
    plt.figure()
    plt.rcParams['axes.facecolor'] = 'slateblue'
    plt.gca().invert_yaxis()
    # age = np.ma.masked_where(np.isnan(Age_grid_list[m]),Age_grid_list[m])
    plt.pcolor(Land, cmap = cmap_1)
    plt.pcolor(Age_grid_list[m], cmap = cmap_2, norm = norm_2)
    plt.axis('equal')
    plt.xlim((75,275))
    plt.ylim((75,275))
    bounds = [0, 1, 2, 3, 4, 5, 6]
    norm = mpl.colors.BoundaryNorm(bounds, cmap_3.N, clip = True)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap = cmap_3, norm = norm), orientation = 'vertical', label = 'Ice Age (years)')
    plt.gca().invert_yaxis()
    cbar.set_ticks(np.arange(0,6) + 0.5)
    cbar.ax.set_yticklabels(['OW','1', '2', '3', '4', '5+'])
    plt.title('{} week, {}'.format(week, year))
    plt.savefig('{}.png'.format(m), dpi = 500)
    plt.close()
    