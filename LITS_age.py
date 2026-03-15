import scipy.io
import pickle
import matlab.engine 
import numpy as np
from glob import glob
from scipy.io import loadmat
import math
import scipy.interpolate as spint
import matplotlib.pyplot as plt
import matplotlib as mpl
from maps import add_land
from scipy.interpolate import NearestNDInterpolator


eng = matlab.engine.start_matlab()

X, Y, i, j, ActiveFlag, Melt, Age, Ageround, ActiveGeo, nWeek, TracerNum \
     = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

Land = add_land('./Grid/bathyEASE.mat')
    
def get_files_and_data(SIC_dir) : 
    
    """Function that get the files and the data for the concentrations of sea ice. 
    
    Inputs: 
        The directory where the concentrations files are located. 

    Returns:
        Two numpy arrays one with all the data that it got and all the files that it read. 
    """
    
    SIC_file = []
    SIC_data = []
    
    files_SIC = sorted(glob(SIC_dir))
    
    for file in files_SIC : 
        SIC_mat  = loadmat(file)
        SIC_file.append(SIC_mat)
    
    for i in range(len(SIC_file)):
        SIC_data.append(SIC_file[i]['sicEASE'])
        
    ArcOceanMask = scipy.io.loadmat('./ArcticOceanMask4IceEdge.mat')
    ArcOceanMask = ArcOceanMask['ArcticOceanMask4IceEdge'] 
    
    SIC_data = np.asarray(SIC_data)
    SIC_mask = SIC_data * ArcOceanMask
    
    
    return SIC_data, SIC_mask, SIC_file

def idx_SI_x_y(SIC_data, eps = 0.15) : 
    
    """Function that get the indices to where there is sea ice. 
    The indices correspond to the sea ice positon on the EASE grid too. 

    Inputs: 
    SIC_data : the array corresponding to all the concentrations data for each week/year. 

    Returns:
        SIC_indice: all of the indices
        SIC_indice_x: where there's sea ice on the x axis
        SIC_indice_y: where there's sea ice on the y axis
    """
    
    SIC_indice = []
    for i in range(len(SIC_data)):  #1979 = 0, 2018 = 39
    
        SIC_ind = (np.where(SIC_data[i,:,:] >  eps)) #finding SIC>15% indices for each year
        SIC_indice.append(SIC_ind) #appending indices for each year
    
    SIC_indice = np.asarray(SIC_indice)
    SIC_indice_x = SIC_indice[:, 0]
    SIC_indice_y = SIC_indice[:, 1]
    
    return SIC_indice, SIC_indice_x, SIC_indice_y

def YearWeek2JWeek( Year, Week ):
    
    """Function to create "Julian" weeks from year, week
    with 1/1/1979 as the starting point.  
    
    Needed to work with weekly SSMI velocity data,
    where we need to create lagrangian tracks across
    year boundaries.  

    Returns:
        Jweek
    """
    
    refyear = 1979 
    jyear = Year - refyear
    return 52* jyear + Week 

def JWeek2YearWeek( JWeek ):
    
    
    """function to create Year and week from "Julian" week, 
    with 1/1/1979 as the starting point.  
    
    Needed to work with weekly SSMI velocity data,
    where we need to create lagrangian tracks across
    year boundaries. 

    Returns:
        Year
        Week
    """

    Week = np.mod(JWeek,52) #divide by 52
    Year = math.floor(JWeek/52) + 1979 
    
    if Week == 0: # Avoid errors when the modulus is 0 for week 52
        Week = 52 
        Year = Year - 1 
    return Year, Week

def Prom_ExpImp(SeaIce_Tracers, SeaIce_MYI_aniversary, SeaIce_MYI_aniv_Beaufort, SeaIce_Promotion_Arctic, SeaIce_Promotion_Beaufort, SeaIce_MYI_import_Beaufort, SeaIce_MYI_export_Beaufort, n) : 
    if n == 0 : 
        print('n = 0 no import/export, no prom/melt')
        
    else : 

        #*promotion 1->2 2->3 3->4 4->5...
        #finding the sea ice of this week and last week
        idx_week3, idx_lastweek = np.where(SeaIce_Tracers[:, nWeek] == n), np.where(SeaIce_Tracers[:, nWeek] == n-1)
        SeaIce_Week1, SeaIce_lastweek1 = SeaIce_Tracers[idx_week3[0], :], SeaIce_Tracers[idx_lastweek[0], :]
        
        Domain_SeaIceWeek1 = np.where(SeaIce_Week1[:, ActiveGeo] > 0)
        Domain_SeaIcelastWeek1 = np.where(SeaIce_lastweek1[:, ActiveGeo] > 0)
        
        SeaIce_Week = SeaIce_Week1[Domain_SeaIceWeek1[0], :]
        SeaIce_lastweek = SeaIce_lastweek1[Domain_SeaIcelastWeek1[0], :]
        
        
        #finding MYI
        idx_myi = np.where(SeaIce_Week[:,Ageround] >= 2)
        SeaIce_MYI_Week = SeaIce_Week[idx_myi[0], :]
        
        #finding SYI
        idx_syi = np.where(SeaIce_Week[:, Ageround] == 2)
        SeaIce_SYI_Week = SeaIce_Week[idx_syi[0], :]
        
        #finding FYI 
        idx_fyi = np.where(SeaIce_lastweek[:, Ageround] == 1)
        FYI_lastweek = SeaIce_lastweek[idx_fyi[0], :]

        #finding the ice floes that were advected and the did not melt
        MYI_lastweek_find = np.isin(SeaIce_lastweek[:, TracerNum], SeaIce_MYI_Week[:,TracerNum])
        MYI_week_find = np.isin(SeaIce_MYI_Week[:, TracerNum], SeaIce_lastweek[:,TracerNum])

        SYI_lastweek_find = np.isin(SeaIce_lastweek[:, TracerNum], SeaIce_SYI_Week[:,TracerNum])
        SYI_week_find = np.isin(SeaIce_SYI_Week[:, TracerNum], SeaIce_lastweek[:,TracerNum])
        
        Prom_icefloe_find = np.isin(SeaIce_SYI_Week[:, TracerNum], FYI_lastweek[:, TracerNum])
        
        idx_prom_week = np.where(Prom_icefloe_find == True)
        
        idx_icefloe_lastweek = np.where(MYI_lastweek_find == True)
        idx_icefloe_week = np.where(MYI_week_find == True)

        idx_SYI_lastweek = np.where(SYI_lastweek_find == True)
        idx_SYI_week = np.where(SYI_week_find == True)

        MYI_lastweek = SeaIce_lastweek[idx_icefloe_lastweek[0], :]
        MYI_week = SeaIce_Week[idx_icefloe_week[0], :]

        SYI_lastweek = SeaIce_lastweek[idx_SYI_lastweek[0], :]
        SYI_week = SeaIce_Week[idx_SYI_week[0], :]
        
        Prom_Week = SeaIce_SYI_Week[idx_prom_week[0], :]

        #*calculating the import and export 

        #MYI import Beaufort this week 
        idx_beaufort_myi_week = np.where(MYI_week[:, ActiveGeo] == 1)
        MYI_Beaufort_week = MYI_week[idx_beaufort_myi_week[0], :]
        
        idx_beaufort_last_week = np.where(MYI_lastweek[:, ActiveGeo] == 1)
        MYI_lastWeek_Beaufort = MYI_lastweek[idx_beaufort_last_week[0], :]

        Beaufortweek_last = np.isin(MYI_lastweek[:, TracerNum], MYI_Beaufort_week[:, TracerNum])
        idx_beaufortweek_last = np.where(Beaufortweek_last == True)

        idx_import_beaufort_MYI = np.where(MYI_lastweek[idx_beaufortweek_last[0], ActiveGeo] != 1)[0]


        #MYI export Beaufort this week
        idx_not_beaufort_myi_week = np.where(MYI_week[:, ActiveGeo] != 1)
        MYI_not_Beaufort_week = MYI_week[idx_not_beaufort_myi_week[0], :]

        Beaufortweek_last = np.isin(MYI_lastWeek_Beaufort[:, TracerNum], MYI_not_Beaufort_week[:, TracerNum])
        idx_not_beaufortweek_last = np.where(Beaufortweek_last == True)

        idx_export_beaufort_MYI = np.where(MYI_lastweek[idx_not_beaufortweek_last[0], ActiveGeo] == 1)[0]


        #* MYI promotion
        
        #TODO need to add here the promotion of only 1 - > 2 year ice,
        #TODO need to change the anniversaries for only promotion 
        #if the Difference is equal to one, the ice gained one year
        Difference_Age = MYI_week[:, Ageround] - MYI_lastweek[:, Ageround]

        idx_aniv_arctic = np.where(Difference_Age == 1)[0]
        
        #finding out of SYI, the one that were FYI last week
        idx_prom_arctic = np.where(SYI_lastweek[:, Ageround] == 1)[0]
        
        idx_prom_beaufort = np.where(Prom_Week[:, ActiveGeo] == 1)[0]

        #promotion in Beaufort Sea
        MYI_aniv_week = MYI_week[idx_aniv_arctic, :]

        idx_aniv_beaufort = np.where(MYI_aniv_week[:, ActiveGeo] == 1)[0]

        SeaIce_MYI_aniversary[n] = len(idx_aniv_arctic)
        SeaIce_MYI_aniv_Beaufort[n] = len(idx_aniv_beaufort)
        SeaIce_Promotion_Arctic[n] = len(idx_prom_arctic)
        SeaIce_Promotion_Beaufort[n] = len(idx_prom_beaufort)
        SeaIce_MYI_import_Beaufort[n] = len(idx_import_beaufort_MYI)
        SeaIce_MYI_export_Beaufort[n] = len(idx_export_beaufort_MYI)
        
    return SeaIce_MYI_aniversary, SeaIce_MYI_aniv_Beaufort, SeaIce_Promotion_Arctic, SeaIce_Promotion_Beaufort, SeaIce_MYI_import_Beaufort, SeaIce_MYI_export_Beaufort

def LITS_Age2(xstart, ystart, startyear, startweek, endyear, endweek, eps = 0.15): 
    
    """
        Lagrangian Ice Tracker System (LITS) 

        The LITS advects a passive tracer from initial xstart,ystart EASE grid 
        coordinate(s), from start date (startyear,startweek) to end date
        (endyear,endweek). In this version of the LITS, we use weekly sea ice 
        drifts from Polar Pathfinder (NSIDC).
        
        This version was modified to track the age of each sea ice parcels. 
        This will also be modified to track the region in which each tracer will be present at each week. 

        [SeaIce_Tracers, SeaIce_Tracers_melted, SeaIce_Melt, SeaIce_Active] = 
            LITS_Age2(xstart, ystart, startyear, startweek, endyear, endweek)

        Where: 

            xstart and ystart are initial EASE Grid coordinates.

            startyear and startweek are the calendar year and week to initialize 
            the track at.
            endyear and endweek define the end date of the track.

            XLag and YLag are the weekly positions of the tracked "parcels".

            ULag and VLag are the weekly velocities of the parcels. 

            ActiveFlag defines when the parcel is still active (vs. having already
            melted).
            
            ActiveGeo defines where the active parcels are located at each week. 
            
            But in this version, we store all of those information in a matrix : SeaIce_Tracers
            
            Dimensions of SeaIce_Tracers : 
                    0 : XLag
                    1 : YLag 
                    2 : ActiveFlag
                    3 : Birth, when the open water became sea ice (jweek)
                    4 : Melt, when the ice floe became open water (jweek)
                    5 : Age, counter of the ice floes age (add 1/52 each week) 
                    6 : AgeRound, ceil the Age dimesions
                    7 : jweek, which jweek we are at
                    8 : week, which week of the year
                    9 : Year, which Year we are


        Additional notes :

            The inputs xstart,ystart are vectors, with one element for each parcel 
            of simulated ice to be tracked forward (end date > start date) 
            or backward (end dates < start date).

            ActiveFlag = 1 is an active parcel.
            ActiveFlag = 0 is a melted parcel.

            Assumes that the velocity vectors from Polar Pathfinder have been 
            pre-processed into a set of interpolant functions (Fu and Fv).

            Calculations are done in the 25km EASE Grid.
            
            ActiveGeo = 1 is in Beaufort Sea
            ActiveGeo = 2 is in Chukchi Sea 
            ActiveGeo = 3 is in East Siberian Sea
            ActiveGeo = 4 is in Laptev Sea
            ActiveGeo = 5 is in the Central Arctic Ocean
            ActiveGeo = 0 is anywhere else 

            
            The initial LITS package was modified to fit the needs of Félix St-Denis. I mainly 
            used the code of Rachel Kim and applied modifications to it. Some things are done 
            trough matlab wrappers. 
    """
    
    SICPath = '/aos/home/fstdenis/LITS/SIC_V4/'
    GridDir = '/aos/home/fstdenis/LITS/Grid/Poly_Masks/'


    ArcOceanMask = scipy.io.loadmat('./ArcticOceanMask4IceEdge.mat')
    ArcOceanMask = ArcOceanMask['ArcticOceanMask4IceEdge'] 
    
    f_1 = open(GridDir+'maskGeo.pckl', 'rb')
    maskGeo = pickle.load(f_1)
    f_1.close()

    dt      = 7*86400  #7days, 86400sec/day # The sea ice velocities are in units of cm/s.
    fac     =  dt * .01/25000 ; # 25000 meters/EASEgrid cell; 100 cm/m 

    #Finding where to end and where to start in Jweek
    jweekstart  = YearWeek2JWeek(startyear, startweek)
    jweekend    = YearWeek2JWeek(endyear, endweek)

    #initialize the matrix that will be used
    SeaIce_Tracers = np.zeros((np.size(xstart), 11)) #Active Flags
    Age_grid_list = []
    SeaIce_Tracers_melted_total = np.zeros((1, 11))  #Melted Ice floes
    SeaIce_MYI_aniversary = np.zeros((jweekend - jweekstart)).astype(int) #Anniversaries/week in Arctic
    SeaIce_MYI_aniv_Beaufort = np.zeros((jweekend - jweekstart)).astype(int) #Anniversaries/week in Beaufort
    SeaIce_Promotion_Arctic = np.zeros((jweekend - jweekstart)).astype(int)
    SeaIce_Promotion_Beaufort = np.zeros((jweekend - jweekstart)).astype(int)
    SeaIce_MYI_Melt_Beaufort = np.zeros((jweekend - jweekstart)).astype(int) #Tracers melted/week in Beaufort
    SeaIce_MYI_Melt_Arctic = np.zeros((jweekend - jweekstart)).astype(int) #Tracers melted/week in Arctic
    SeaIce_MYI_import_Beaufort = np.zeros((jweekend - jweekstart)).astype(int) #Import in Beaufort/week
    SeaIce_MYI_export_Beaufort = np.zeros((jweekend - jweekstart)).astype(int) #Export from Beaufort/week
    SeaIce_Total_Melt = np.zeros((jweekend - jweekstart)).astype(int) #Melt of all sea ice
    OpenWater2FYI_prom_Beaufort = np.zeros((jweekend - jweekstart)).astype(int)
    MYI_extent_Beaufort = np.zeros((jweekend -  jweekstart)).astype(int)
    
    X, Y, i, j, ActiveFlag, Melt, Age, Ageround, ActiveGeo, nWeek, TracerNum \
        = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    
    #-------------------------------------------
    #------ Initialisation
    SeaIce_Tracers[:, X] = xstart 
    SeaIce_Tracers[:, Y] = ystart
    SeaIce_Tracers[:, i] = np.around(xstart)
    SeaIce_Tracers[:, j] = np.around(ystart)
    SeaIce_Tracers[:,ActiveGeo] = maskGeo[xstart, ystart]
    SeaIce_Tracers[:, ActiveFlag] = int(1)
    SeaIce_Tracers[:, nWeek] = 0
    SeaIce_Tracers[:,TracerNum] = np.arange(0, len(xstart)).astype(int)
    SeaIce_Tracers[:, Ageround] = 1
    
    Age_grid_start = np.zeros((361, 361))*np.nan
    Age_grid_start[xstart, ystart] = 1
    
    Age_grid_list.append(Age_grid_start)
    
    IJ_Beaufort_idx = np.where(maskGeo == 1)
    
    IJ_Beaufort = np.vstack((IJ_Beaufort_idx[0], IJ_Beaufort_idx[1])).T
    
    end_melt_season_weeks = list(range(39, 2223+52, 52)) #* for september
    
    x361 = np.arange(0,361,1) #defining 361x361 grid for SIV, SIC interpolation
    y361 = np.arange(0,361,1)
    RGI = spint.RegularGridInterpolator


    if jweekstart < jweekend:
        timearrow = 1 # forward tracking
        julianend = jweekend
    
    else:
        timearrow = -1 # backward tracking
        julianend = jweekend+1
    #-------------------------------------------
        
    #Main loop, iterate trough the weeks until the jweekend
    for jweek in range(jweekstart, julianend, timearrow) :
        
        n = abs(jweekstart - jweek)
        
        if timearrow == 1: # forward tracking
        
            lookupyear, lookupweek  = JWeek2YearWeek(jweek) #converting jweek to Year, Week
            stringyear = str(lookupyear)
            
            lookupnextyear, looknextweek = JWeek2YearWeek(jweek+1)
            stringnextyear = str(lookupnextyear)
        
        
            if lookupweek < 10 : # Add a 0 in front of a week smaller than 10 to follow naming convention.
                stringweek = '0'+str(lookupweek)  
            
                
            else:
                stringweek = str(lookupweek)
                
            if looknextweek < 10 : 
                stringnextweek = '0'+str(looknextweek)
                
            else : 
                stringnextweek = str(looknextweek)
                
                
        #---------------------------------------------------------------------
        # Where we miss some data because of a change of satellite
        # So we just keep the same positions for each time steps but we add +1/52 to the Age
        
        if (stringyear == '2021' and stringweek == '52') : 

            SIC = scipy.io.loadmat(SICPath +'CDR_SIC_'+ stringyear + '_' + stringweek + '.mat') #load .mat file
            SIC_next = 0
        
        else : 
                
                SIC = scipy.io.loadmat(SICPath +'CDR_SIC_'+ stringyear + '_' + stringweek + '.mat') #load .mat file
                SIC_next = scipy.io.loadmat(SICPath +'CDR_SIC_'+ stringnextyear + '_' + stringnextweek + '.mat') #load .mat file
        #---------------------------------------------------------------------

        #last week of the record at the time
        if (stringyear == '2021' and stringweek == '52') : #!
            SIC = SIC['sicEASE'] #storing SIC values from .mat file
            
        else : 
            SIC = SIC['sicEASE'] #storing SIC values from .mat file
            SIC_next = SIC_next['sicEASE']
        
        print('n = ',n, 'stringyear', stringyear, 'stringweek', stringweek)
        
        # #finding the tracers for the week n 
        idx_week = np.where(SeaIce_Tracers[:, nWeek] == n)
        
        #---------------------------------------------------------------------
        # Update the Age 
        if n == 0 :
            # at first step, everything is zero
            print('Age = 0')
            SeaIce_Tracers[idx_week[0], Ageround] = 1 # for september 
            
        else  : 

            SeaIce_Tracers[idx_week[0], ActiveGeo] = maskGeo[SeaIce_Tracers[idx_week[0], i].astype(int), SeaIce_Tracers[idx_week[0], j].astype(int)]
            
            Age_grid_n = np.zeros((361, 361))*np.nan
            Age_grid_n[SeaIce_Tracers[idx_week[0], i].astype(int), SeaIce_Tracers[idx_week[0], j].astype(int)] = SeaIce_Tracers[idx_week[0], Ageround]


            for l in range(len(Pos_Same)) : 

                Age_cell = Age_Same[l]

                Age_grid_n[Pos_Same[l][0].astype(int), Pos_Same[l][1].astype(int)] = Age_cell
 

            Age_grid_list.append(Age_grid_n)
            
            # age distribution in the Beaufort Sea
            
            Age_Beaufort_n = Age_grid_n[IJ_Beaufort_idx]
            
            idx_MYI_beaufort = np.where(Age_Beaufort_n >= 2)
            
            MYI_ext_Beaufort = Age_Beaufort_n[idx_MYI_beaufort]
            
            MYI_extent_Beaufort[n] = len(MYI_ext_Beaufort)
            
        
        if n == 311 : #!
            print('end')
            
        else : 
            
            #-------------------------------------------------------------------------
            # Advection
            SeaIce_Tracers_nextweek = np.zeros((len(idx_week[0]), 11))
            
            XLag_matlab = matlab.double(SeaIce_Tracers[idx_week[0], X].tolist())
            YLag_matlab = matlab.double(SeaIce_Tracers[idx_week[0], Y].tolist())
        

            # #matlab wrapper for velocity interpolation
            uwork, vwork = eng.interp_vel(XLag_matlab, YLag_matlab, stringweek, stringyear, nargout = 2)
           
            uwork = np.asarray(uwork)
            vwork = np.asarray(vwork)
            
            #* for V4
            vwork = np.flipud(vwork)    

            #Making the ice move
            SeaIce_Tracers_nextweek[:, X] = SeaIce_Tracers[idx_week[0], X] + timearrow*vwork*fac
            SeaIce_Tracers_nextweek[:, Y] = SeaIce_Tracers[idx_week[0], Y] + timearrow*uwork*fac
            
            #updating the tracers values 
            SeaIce_Tracers_nextweek[:,i] = np.around(SeaIce_Tracers_nextweek[:, X])
            SeaIce_Tracers_nextweek[:,j] = np.around(SeaIce_Tracers_nextweek[:, Y])
            SeaIce_Tracers_nextweek[:, nWeek] = n+1
            SeaIce_Tracers_nextweek[:, ActiveFlag] = 1 
            SeaIce_Tracers_nextweek[:, Ageround] = SeaIce_Tracers[idx_week[0], Ageround] #* september
            SeaIce_Tracers_nextweek[:, TracerNum] = SeaIce_Tracers[idx_week[0], TracerNum]
            

            #*september
            end_melt_season = end_melt_season_weeks.count(n+1)
            if end_melt_season > 0 : 
                SeaIce_Tracers_nextweek[:, Ageround] += 1
                print('end of melt')
            
            else : 
                # SeaIce_Tracers_nextweek[:, Ageround] = SeaIce_Tracers[idx_week[0], Ageround]
                dummy = 0
                
                
            # Elimination/Finding same IJ
            I_J = np.vstack([SeaIce_Tracers_nextweek[:, i], SeaIce_Tracers_nextweek[:, j]]).T.astype(int)

            unique_IJ, counts_values_position = np.unique(I_J, axis = 0, return_counts = True)
            SameXY = unique_IJ[counts_values_position > 1]
            
            repeated_idx = 0
            idx_delete = []
            Pos_Same = []
            Age_Same = []

            k = 0
            for pos in SameXY:
                repeated_idx = np.argwhere(np.all(I_J == pos, axis=1)).ravel()

                Ages_repeated = SeaIce_Tracers_nextweek[repeated_idx, Ageround]

                AgeMax_idx = np.argmax(Ages_repeated)
                Age_Max = Ages_repeated[AgeMax_idx]

                
                Pos_Same.append(pos)
                Age_Same.append(Age_Max) 
                
                repeated_idx = np.delete(repeated_idx, AgeMax_idx)

                
                len_repeated = len(repeated_idx)
                for l in range(0, len_repeated) : 
                    idx_delete.append(repeated_idx[l].astype(int))  
                    

            # Elimination or not
            SeaIce_Tracers_nextweek = np.delete(SeaIce_Tracers_nextweek, idx_delete, axis = 0)


            #---------------------------------------------------------------------
            # Melt
            rgi_melt = RGI(points=(x361, y361), values = SIC_next, method = 'linear')
            ActiveFlag_out = ((rgi_melt((SeaIce_Tracers_nextweek[:,X], SeaIce_Tracers_nextweek[:,Y]))) > eps)
            
            #Updating the Active Flags in main matrix
            SeaIce_Tracers_nextweek[:, ActiveFlag] = np.asarray(ActiveFlag_out).astype(int)

            # #here we look at the ActiveFlags that became 0, so that the Ice floes have melted
            Ind_melt = np.where(SeaIce_Tracers_nextweek[:, ActiveFlag] == 0) #finding the melted one
            Ind_active = np.where(SeaIce_Tracers_nextweek[:, ActiveFlag] == 1)
            
            # #store the melted in a new matrix with only the melted ice floes
            SeaIce_Tracers_melted = SeaIce_Tracers_nextweek[Ind_melt[0],:]
            SeaIce_Tracers_nextweek_elim_melted = SeaIce_Tracers_nextweek[Ind_active[0],:]
            
            idx_geo = np.where((SeaIce_Tracers_melted[:, X] >=0) & (SeaIce_Tracers_melted[:, X] < 361))

            SeaIce_Tracers_melted[idx_geo[0], ActiveGeo] = maskGeo[SeaIce_Tracers_melted[idx_geo[0],i].astype(int), SeaIce_Tracers_melted[idx_geo[0],j].astype(int)]   
            
            #finding the MYI that melted
            idx_myi_melt = np.where(SeaIce_Tracers_melted[idx_geo[0], Ageround] >= 2.)
            
            Melted_MYI = SeaIce_Tracers_melted[idx_myi_melt[0], :]
            Domain_MYI_Melted = np.where(Melted_MYI[:, ActiveGeo] > 0)
            SeaIce_Domain_MYI_Melted = Melted_MYI[Domain_MYI_Melted[0], :]
            
            idx_Beaufort = np.where(SeaIce_Domain_MYI_Melted[:, ActiveGeo] == 1)

            SeaIce_MYI_Melt_Beaufort[n+1] = len(idx_Beaufort[0])
            SeaIce_MYI_Melt_Arctic[n+1] = len(Domain_MYI_Melted[0])
            SeaIce_Total_Melt[n+1] = len(Ind_melt[0])
            
            SeaIce_Tracers_melted_total = np.append(SeaIce_Tracers_melted_total, SeaIce_Tracers_melted, axis = 0)
            
            #---------------------------------------------------------------------
                 
                 
            #---------------------------------------------------------------------       
            # Formation
            SIC_idc = (np.where(SIC_next*ArcOceanMask > eps))
            SIC_idx, SIC_idy = SIC_idc[0], SIC_idc[1]
            
            #Comparing where there's advected ice and the SIC

            XLag_YLag = np.vstack((SeaIce_Tracers_nextweek_elim_melted[:, i], SeaIce_Tracers_nextweek_elim_melted[:, j])).T
            SIC_x_y = np.vstack((SIC_idx, SIC_idy)).T
    
            Existing_Tracers = (SIC_x_y[:, None] == XLag_YLag).all(-1).any(-1) #finding rows that overlap between existing tracer and SIC map
            
            new = SIC_x_y[~Existing_Tracers] #deleting overlapped xlag,ylag pairs from SIC to identify newly formed ice
            new_XLag = new[:,0]
            new_YLag = new[:,1]
            
            length_new = len(new_XLag)
            num_max = SeaIce_Tracers[:,TracerNum].max()
            
            SeaIce_Tracers_new = np.zeros((np.size(new_XLag), 11))
            SeaIce_Tracers_new[:, ActiveFlag] = int(1)
            SeaIce_Tracers_new[:, X], SeaIce_Tracers_new[:, Y] = new_XLag, new_YLag
            SeaIce_Tracers_new[:, i], SeaIce_Tracers_new[:, j] = new_XLag, new_YLag
            SeaIce_Tracers_new[:, ActiveGeo] = maskGeo[SeaIce_Tracers_new[:, i].astype(int), SeaIce_Tracers_new[:, j].astype(int)]
            SeaIce_Tracers_new[:, nWeek] = n+1
            SeaIce_Tracers_new[:, Ageround] = 1
            SeaIce_Tracers_new[:, TracerNum] = np.arange(num_max+1, num_max + length_new + 1).astype(int)
            
            idx_OW2FYI_Beaufort = np.where(SeaIce_Tracers_new[:, ActiveGeo] == 1)
            OpenWater2FYI_prom_Beaufort[n] = len(idx_OW2FYI_Beaufort[0])
            
            SeaIce_Tracers_nextweek_elim_melted = np.append(SeaIce_Tracers_nextweek_elim_melted, SeaIce_Tracers_new, axis = 0)
            SeaIce_Tracers = np.append(SeaIce_Tracers, SeaIce_Tracers_nextweek_elim_melted, axis = 0 )
            
            #---------------------------------------------------------------------
            

        SeaIce_MYI_aniversary, SeaIce_MYI_aniv_Beaufort, SeaIce_Promotion_Arctic, SeaIce_Promotion_Beaufort, SeaIce_MYI_import_Beaufort, SeaIce_MYI_export_Beaufort = \
            Prom_ExpImp(SeaIce_Tracers, SeaIce_MYI_aniversary, SeaIce_MYI_aniv_Beaufort, SeaIce_Promotion_Arctic, SeaIce_Promotion_Beaufort, 
                        SeaIce_MYI_import_Beaufort, SeaIce_MYI_export_Beaufort, n)
        
    return SeaIce_Tracers, Age_grid_list, SeaIce_Tracers_melted_total, SeaIce_MYI_Melt_Beaufort,  SeaIce_MYI_aniv_Beaufort, SeaIce_MYI_Melt_Arctic, SeaIce_MYI_aniversary, \
            SeaIce_MYI_import_Beaufort, SeaIce_MYI_export_Beaufort, SeaIce_Promotion_Arctic, SeaIce_Promotion_Beaufort, \
                MYI_extent_Beaufort, OpenWater2FYI_prom_Beaufort