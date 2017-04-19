# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 12:23:01 2010

@author: Darren

Updated: Sat Feb. 9 2013
"""


#the data inside each parcel is : [temperature, salinity, nitrate, co2, max_depth, u, v, w_velocity, number of cells here]
    
#trying to make the 'environmental framework' to get the cells into a 3 dimensional framework

#8/21/2013 - changed the number of cells variable to become the PAR data variable


import math
#import copy
import netCDF4
import numpy as np
#import time
import NODC_nitrate
import pickle
from numpy import array as np_array
from numpy import isnan as np_isnan
#from scipy.stats.stats import nanmean
from numpy import nanmean
from scipy.ndimage import map_coordinates

world_latitude = []
world_longitude = []
HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
HYCOM_depth_dictionary = {}
MODIS_PAR_data = [0, 0, 0]  #will be year, day, data
MODIS_PAR_grid_points = []
##filepath_HYCOM = '/home/darren/ibm_model/HYCOM_GOM_data_converted/'
##filepath_MODIS = '/home/darren/ibm_model/MODIS_PAR_data_converted/'
filepath_HYCOM = 'D:/CJunk/HYCOM_GOM_data_converted/'
filepath_HYCOM_future = 'D:/CJunk/HYCOM_GOM_data_future_converted/'
filepath_MODIS = 'D:/CJunk/MODIS_PAR_data_converted/'

for indiv_depth in xrange(201):
    valuetoadd = 1
    while indiv_depth >= HYCOM_depth_values[valuetoadd]:
        valuetoadd += 1
        if valuetoadd > 16:
            valuetoadd = 17
            break
    HYCOM_depth_dictionary[indiv_depth] = valuetoadd - 1



def get_HYCOM_data(year, day, hour, future=False):
    #base_data_location = 'G:/HYCOM_GOM_data_converted/archv.'
    #base_data_location = 'D:/CJunk/HYCOM_GOM_data_converted/archv.'
    base_data_location = filepath_HYCOM + 'archv.'  #this is for campbell13
    #base_data_location = 'J:/HYCOM_GOM_data_converted/archv.'
    #base_data_location = 'C:/Documents and Settings/Lisa/My Documents/Darren_IFCB6/Kbrevis_model/HYCOM_data/archv.'
    #base_data_location = 'H:/temp_HYCOM_data/Converted_HYCOM_data/archv.'
    if future:
        base_data_location = filepath_HYCOM_future + 'archv.'
    #add year
    base_data_location += str(year) + '_'
    #add_day
    if day < 10:
        base_data_location += '0'
    if day < 100:
        base_data_location += '0'
        
    base_data_location += str(day) + '_'
    #add hour
    if hour < 10:
        base_data_location += '0'
    base_data_location += str(hour) + '_3z_adj.nc'
    #print base_data_location
    try:
        dataset = netCDF4.Dataset(base_data_location, mode='r')
    except:
        dataset = None

    return dataset


         

def interpolate_nitrate_values_version2(world, temp_max_depth):
    #the depth values from HYCOM are not evenly spaced and values in-between must be interpolated
    NODC_depth_values = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200]
    world_depth = len(world[0][0])
    val1 = -1
    val2 = 0
    world = world.T
    for zdepth in range(world_depth):
        if zdepth in NODC_depth_values:
            val1 += 1
            val2 += 1
        elif zdepth < temp_max_depth:
            temp_val1 = NODC_depth_values[val1]
            temp_val2 = NODC_depth_values[val2]
            total_difference = temp_val2 - temp_val1
            dist_to_1 = zdepth - temp_val1
            dist_to_2 = temp_val2 - zdepth
            scaled = [1-(dist_to_1 / float(total_difference)), 1 - (dist_to_2 / float(total_difference))]
            if world_depth >= NODC_depth_values[val2]:
                world[2][zdepth] = (world[2][temp_val1] * scaled[0]) + (world[2][temp_val2] * scaled[1])
            else:
                world[2][zdepth] = (world[2][temp_val1] * scaled[0]) + (world[2][temp_val1] * scaled[1])
    world = world.T             
    return world         

def update_nitrate_values_version2(world, temp_max_depth):
    #the depth values from HYCOM are not evenly spaced and values in-between must be interpolated
    NODC_depth_values = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200]
    world_depth = len(world[0][0])
    val1 = -1
    val2 = 0
    world = world.T
    for zdepth in range(world_depth):
        if zdepth in NODC_depth_values:
            val1 += 1
            val2 += 1
        elif zdepth < temp_max_depth:
            temp_val1 = NODC_depth_values[val1]
            temp_val2 = NODC_depth_values[val2]
            total_difference = temp_val2 - temp_val1
            dist_to_1 = zdepth - temp_val1
            dist_to_2 = temp_val2 - zdepth
            scaled = [1-(dist_to_1 / float(total_difference)), 1 - (dist_to_2 / float(total_difference))]
            if world_depth >= NODC_depth_values[val2]:
                world[2][zdepth] = (world[2][temp_val1] * scaled[0]) + (world[2][temp_val2] * scaled[1])
            else:
                world[2][zdepth] = (world[2][temp_val1] * scaled[0]) + (world[2][temp_val1] * scaled[1])
    world = world.T             
    return world 
    
##def nitrate_gradient(location, world):
##    if location[2] < -1:
##        no3_1 = world[location[0]][location[1]][location[2]+1]
##    else:
##        no3_1 = world[location[0]][location[1]][location[2]]
##    if location[2] > max_depth+1:
##        no3_2 = world[location[0]][location[1]][location[2]-1]
##    else:
##        no3_2 = world[location[0]][location[1]][location[2]]
##        
##    return [no3_1, no3_2]

def create_thermocline(depth_of_thermocline, temp_of_thermocline, world):
    for x in range(len(world)):
        for y in range(len(world[0])):
            world[x][y][depth_of_thermocline][0] = float(temp_of_thermocline)

        
def initialize_MODIS_grid(environ_PAR_shape, origin_offset, world_shape):
    #this function will create a list of points mapping the world grid to the MODIS grid points; world grid points are integer and
    #the resulting MODIS grid point locations of those integers are floats
    global MODIS_PAR_grid_points
    print "MODIS_PAR_grid_points initializing..."
    PAR_data_grid_file = filepath_MODIS + 'MODIS_PAR_data_grid.pck'  #this is for running on campbell13
    #PAR_data_grid_file = 'D:/CJunk/MODIS_PAR_data_grid.pck'
    try:
        grid_file = open(PAR_data_grid_file)
        MODIS_PAR_grid_points = pickle.load(grid_file)
        grid_file.close()
         

    except:
        for x in range(450):
            for y in range(450):
                #convert the location to lat/lon
                kbr_x, kbr_y = convert_indexes_to_lat_lon([x, y])
                modis_x = 312 - (((90-kbr_x) / 0.04166666) - 1416)        #original formula: 90 - ((1416 + x) * 0.04166666)
                modis_y = ((kbr_y + 180) / 0.04166666) - 1968   #original formula: -180 + ((1968 + y) * 0.04166666)
                
                if kbr_x != -1 and kbr_y != -1 :
                    MODIS_PAR_grid_points.append([[x, y], [modis_x, modis_y]])
        grid_file = open(PAR_data_grid_file, 'w')
        pickle.dump(MODIS_PAR_grid_points, grid_file)
        grid_file.close()
    print "MODIS_PAR before trimming:", len(MODIS_PAR_grid_points)
    #now trim the data points to the size of the model domain being used for this run
    keepers = []
##    print
##    print "origin_offset", origin_offset
##    print "MODIS_grid", MODIS_PAR_grid_points[0]
##    print "environ shape", environ_PAR_shape
    for data_point in MODIS_PAR_grid_points:
        if 0 < data_point[0][0] - origin_offset[1] < world_shape[1] and 0 < data_point[0][1] - origin_offset[0] < world_shape[0]:
            keepers.append(data_point)
            #adjust the values to the origin offset so the lat/lons match up
            keepers[-1][0][0] -= origin_offset[1]
            keepers[-1][0][1] -= origin_offset[0]
        
    MODIS_PAR_grid_points = keepers

        
    print "Done!", "Length PAR grid:", len(MODIS_PAR_grid_points)
    
def convert_MODIS_grid_to_KbrModel(new_data, environ_PAR_shape):   
    global MODIS_PAR_grid_points
    #this will convert the model grid to the MODIS grid
    #MODIS grid has 312 vertical points (Latitude) and 432 horizontal points (Longitude)
    #model grid can vary but has max 350 vertical and 450 horizontal points
    temp_PAR = np.zeros(environ_PAR_shape)#.T
    for data_point in MODIS_PAR_grid_points:
            temp_val = map_coordinates(new_data, [[data_point[1][0]],[data_point[1][1]]], order=1)
            if temp_val > 0:
                temp_PAR[data_point[0][0]][data_point[0][1]] = temp_val
            else:
                temp_val = map_coordinates(new_data, [[data_point[1][0]],[data_point[1][1]]], order=0)
                if temp_val > 0:
                    temp_PAR[data_point[0][0]][data_point[0][1]] = temp_val
    #temp_PAR = temp_PAR.T        
    temp_PAR[temp_PAR == 0] = np.nan
    return temp_PAR



    

def PAR_data_MODIS(year, day, environ_PAR, origin_offset, initialize=False, world_shape=False, reverse=False):
    #this will use the PAR info from MODIS to reduce sunlight at certain locations across the Gulf of Mexico
    #it will add the PAR data to the environment variable so that it's only calculated every time the environment updates
    #other option was to have each cell find it but with increased cells this results in more function calls than simply updating the environment
    global MODIS_PAR_data, MODIS_PAR_grid_points
    if MODIS_PAR_data[0] != year:
        f = netCDF4.Dataset(filepath_MODIS + 'MODIS_PAR_' + str(year) + '.cdf', mode='r')
        MODIS_PAR_data[2] = f.variables['PAR'][:]
        f.close()
        MODIS_PAR_data[0] = year
    if MODIS_PAR_data[1] != day:
        MODIS_PAR_data[1] = day
    #need to convert from MODIS grid to model grid
    if initialize:
        initialize_MODIS_grid(environ_PAR.shape, origin_offset, world_shape)
    if reverse:
        temp_new_PAR = convert_MODIS_grid_to_KbrModel(MODIS_PAR_data[2][day-1], environ_PAR.shape)
        environ_PAR = nanmean([environ_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR], axis=0)
    else:
        temp_new_PAR = convert_MODIS_grid_to_KbrModel(MODIS_PAR_data[2][day], environ_PAR.shape)
        environ_PAR = nanmean([environ_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR, temp_new_PAR], axis=0)
    print 'Environ_PAR:', environ_PAR.shape
    return environ_PAR

          
def nitrate_from_NODC(julian_day):
    #this should load the nitrate data from the NODC files
    #it will return a list of functions, one for each of 7 depths (0, 10, 20, 30, 50, 75, 100)
    #the function is passed the longitude and latitude and it will return an array with the nitrate concentration at that depth
    #e.g. w(-88, 26) will return array[[0.123456]]
    environ_nitrate = NODC_nitrate.get_monthly_nitrate(julian_day)
    return environ_nitrate



    
def nitrate_concentrations(source, new_or_modify, world, mixed_layer_depth = None, cell_depth = None, cell_loc = None):
    #start placing the nitrate into the system
        #(-0.01*z + ((0.75*(1.0 + math.tanh(-z-13.0))) / 2.0))  <-- this is the original formula from Liu et al 2001 for increasing nitrate as the depth increases
        #i'm going to use this formula for depth increases and the x domain; so one side of the bracket can/will have higher concentrations
        #this makes a concentration; i need to convert this concentration to a stock value so i can reduce the nitrate accordingly
    background_value = 0.05
    scale_factor = 0.2
    if new_or_modify == 'new':
##        
        #world = world.T
        world[2] = background_value
        ###scale according to environmental parameters
        #world[2] *= (36. / world[1]) #scale according to salinity #old way pre Nov 2016
        world[2] += ((-19.13*np.log(world[1])) + 68.952)*0.5  #calculated from the MS06 and MS07 cruises off LA using a logarithimic trendline in Excel
        world[2][world[2] < 0] = background_value
        #world[2] *= 10 * (30. / world[0])**3 #scale according to temperature
        
        
#        ###this uses the mixed_layer_depth to apportion nitrate in the water column
#        temp_nitrate = (HYCOM_depth_values - world[7].T)
#        temp_tanh = np.vectorize(math.tanh)
#        temp_nitrate = temp_tanh(temp_nitrate)
#        temp_nitrate = (scale_factor *(-0.01*(-np.array(HYCOM_depth_values)) + (0.75*(1.01 + temp_nitrate))))  #was 0.1 as the first value
#        world[2] += temp_nitrate.T
#        ####end of MLD 
        
#        world[2][world[4] < 51] += 0.05 #this adds a little extra nitrate for regions closer to shore
#        world[2][world[4] < 41] += 0.1 #this adds a little extra nitrate for regions closer to shore
#        world[2][world[4] < 31] += 0.2 #this adds a little extra nitrate for regions closer to shore
#        world[2][world[4] < 21] += 0.5 #this adds a little extra nitrate for regions closer to shore
#        world[2][world[4] < 11] += 1.0 #this adds a little extra nitrate for regions closer to shore
#        world[2][world[4] < 5]  += 2.0 #this adds a little extra nitrate for regions closer to shore
         
        #this one will have nitrate controlled by temp and salinity with temperature more dominant in deeper areas and salinity dominant in shallow areas
        #gradients = [[0.2, 0.8], [0.4, 0.6], [0.5, 0.5], [0.6, 0.4], [0.7, 0.3], [0.8, 0.2], [0.9, 0.1], [1,0]]
        #for temp_depth, threshold_factor in zip([50, 40, 30, 25, 20, 15, 10, 5], gradients):
        #    world[2][world[4] <= temp_depth] = background_value
        #    world[2][world[4] <= temp_depth] += threshold_factor[0] * ((-19.13*np.log(world[1][world[4] <= temp_depth])) + 68.952)  #calculated from the MS06 and MS07 cruises off LA using a logarithimic trendline in Excel
        #    world[2][world[2] <= 0] = background_value
        #    world[2][world[4] <= temp_depth] += threshold_factor[1] * (0.05 * (30. / world[0][world[4] <= temp_depth])**3) #scale according to temperature
        #world[2][world[4] > 51] += 0.05 * (30. / world[0][world[4] > 51])**3 #scale according to temperature
        
        
    elif new_or_modify == 'gradient_check':
        ##this was the original way to check the gradient in the water column, doesn't work when I modify the N profile differently
        conc_above = background_value + scale_factor *(-0.01*(-(cell_depth-1)) + ((0.75*(1.0 + math.tanh((cell_depth-1) - mixed_layer_depth))) * 0.5)) #was 0.1 for first value
        conc_below = background_value + scale_factor *(-0.01*(-(cell_depth+1)) + ((0.75*(1.0 + math.tanh((cell_depth+1) - mixed_layer_depth))) * 0.5)) #was 0.1 for first value
                
        return conc_above, conc_below

    
    #the conversion to pmol of nitrate (for 1 cubic meter of water): take the concentration
    #concen/1000000, *1000, *(10**12) <--makes it into picomols
    #subtract out how many picomols were used and then to take it back into the concentration:
    #picomols/(10**12), /1000, *1000000
    #the short way to do so is: concen*(10**9) --> to picomols
    #picomols / 10**9 --> to concen
    #nitrate above here
    

def convert_lat_lon(indata):
    #this will take a lat/lon location as input [lat, lon] and return the integer indexes closest to that location
    temp_lat = indata[1]
    temp_lon = indata[0]
    x_index = 0
    y_index = 0
    try:
        while temp_lon > world_longitude[x_index]:
            x_index += 1
            if x_index == len(world_longitude):
                x_index = -1
                break
        while temp_lat > world_latitude[y_index]:
            y_index += 1
            if y_index == len(world_latitude):
                y_index = -1
                break
    except:
        x_index = -1
        y_index = -1
        
    return [x_index, y_index]

def convert_indexes_to_lat_lon(indata):
    #this will take a integer indexes as input [x, y] and return the lat/lon closest to that location
    temp_x = indata[0]
    temp_y = indata[1]
    x_index = int(temp_x)
    y_index = int(temp_y)
    try:
        final_x = world_latitude[x_index]
        final_y = world_longitude[y_index]
        residual_x = (world_latitude[x_index+1] - world_latitude[x_index]) * (temp_x - x_index)
        residual_y = (world_longitude[y_index+1] - world_longitude[y_index]) * (temp_y - y_index)
        final_x += residual_x
        final_y += residual_y
        
    except:
        final_x = -1
        final_y = -1
    
    return [final_x, final_y]

def initialize_environment_version3(size_x_dimension, size_y_dimension, size_z_dimension, year, day, hour, origin_offset, nitrate_source):
    #now I need to make the returned world only the 17 depths of the model
    #intervening points will be interpolated on the fly, this should cut down the memory footprint considerably
    #update the environmental parameters based on the time of day and date using HYCOM model data
    #the data inside each parcel is : [temperature, salinity, nitrate, co2, max_depth, u, v, mixed_layer_depth, daily PAR]
    #8/21/2013 changing the parcel to have PAR data instead of number of cells here
    
    global world_latitude, world_longitude, MODIS_PAR_grid_points
    HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
    new_data = get_HYCOM_data(year, day, hour)
    print "Initializing world..."
    
    temp_temp = new_data.variables['temperature'][0]
    temp_sal = new_data.variables['salinity'][0]
    temp_u = new_data.variables['u'][0]
    temp_v = new_data.variables['v'][0]
    
    #temp_depths = new_data.variables['Depth'][:]
    temp_lats = new_data.variables['Latitude'][:]
    temp_lons = new_data.variables['Longitude'][:]
    temp_mld = new_data.variables['mixed_layer_depth'][0]
    temp_mld = np_array([temp_mld] * 17)

    world_latitude = temp_lats[origin_offset[1]:origin_offset[1]+size_y_dimension]
    world_longitude = temp_lons[origin_offset[0]:origin_offset[0]+size_x_dimension]
    
    zero_float = np.zeros((17, 350, 450), dtype=np.float)
    zero_int1 = np.zeros((17, 350, 450), dtype = np.int)

    zero_int2 = np.zeros((350, 450), dtype = np.float)
    zero_int2[zero_int2 == 0] = np.nan
    zero_int2 = np_array([PAR_data_MODIS(year, day, zero_int2, origin_offset, initialize=True, world_shape=(size_x_dimension, size_y_dimension, size_z_dimension))] * 17)
    ##zero_int2 is to keep track of cell counts, i'm going to add it to the salinity values to get the nans and then convert all other values back to zero
    ##zero_int2 = zero_int2 + temp_sal
    ##zero_int2[zero_int2 < 1000] = 0
    ##zero_int2[zero_int2 > 1000000] = np.nan
    
    #convert missing data
##    temp_temp[temp_temp > 1E7] = np.nan
##    temp_sal[temp_sal > 1E7] = np.nan
##    temp_u[temp_u > 1E7] = np.nan
##    temp_v[temp_v > 1E7] = np.nan
##    temp_mld[temp_mld > 1E7] = np.nan

    
    co2 = np.zeros((17, 350, 450), dtype=np.float)
    co2[co2 == 0] = 1
    #begin creation of world by combining arrays, reordering axes, and adding in missing depth levels
    temp_world = np_array([temp_temp, temp_sal, zero_float, co2, zero_int1, temp_u, temp_v, temp_mld, zero_int2])
    temp_world = temp_world.T
    print temp_world.shape
    temp_world = np.rollaxis(temp_world, 2, 0) #move the depth (z-axis) to the front
    print temp_world.shape
    print "created the basic world framework"
    
    temp_world = np.rollaxis(temp_world, 2, 0)
    temp_world = np.rollaxis(temp_world, 2, 0)
    print temp_world.shape
    #end of world creation

    #begin world trimming to size and location requested
    if size_x_dimension < 450:
        temp_world = temp_world[origin_offset[0]:size_x_dimension+origin_offset[0]]
    if size_y_dimension < 350:
        temp_world = np.rollaxis(temp_world, 1, 0)
        temp_world = temp_world[origin_offset[1]:size_y_dimension+origin_offset[1]]
        temp_world = np.rollaxis(temp_world, 1, 0)
    print temp_world.shape

    #####starting nitrate placement; doing it here so that the nitrate values will be placed in the world before interpolation
    if nitrate_source == 'uniform':
        #add in a nitrate value
        print "Adding nitrate concentrations",
        nitrate_concentrations(1, 'new', temp_world.T)     #placed the nitrate gradient maker into a function so that the gradient can be modified in the middle of a run
        #temp_world = temp_world.T
        print "Done!"
        
    elif nitrate_source == 'NODC':
        print "Loading in nitrate concentrations from NODC data...",
        nitrate_func = nitrate_from_NODC(day) #should return seven functions for depths 0, 10, 20, 30, 50, 75, 100
        print "nitrate functions loaded"
        depths_nitrate = [0,10,20,30,50,75,100,125,150,200]
            
        for x in range(size_x_dimension):
            for y in range(size_y_dimension):
                for z in depths_nitrate:
                    if z < size_z_dimension:
                        if temp_world[x][y][z][0] < 1E6:
                            temp_nitrate_value = 1. * (nitrate_func[depths_nitrate.index(z)](temp_lons[x+origin_offset[0]], temp_lats[y+origin_offset[1]])[0][0])
                            if temp_nitrate_value < 0:
                                temp_nitrate_value = 0.0
                            elif temp_nitrate_value > 50:
                                temp_nitrate_value = 50.0
                            temp_world[x][y][z][2] = temp_nitrate_value
        temp_world = interpolate_nitrate_values_version2(temp_world, size_z_dimension)
        print "Finished placing nitrate values"
    #####
    
    
    world = temp_world
    ######update max depth for each location
    #HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
    for x in range(size_x_dimension):
        for y in range(size_y_dimension):
            maximum_depth_here = 0
            if np_isnan(world[x][y][0][0]):
                #for z in range(size_z_dimension):
                for z in range(17):
                    world[x][y][z][4] = np.nan
                    world[x][y][z][8] = np.nan
            else:
                for z in range(17):
                    if not np_isnan(world[x][y][z][0]):
                        maximum_depth_here = HYCOM_depth_values[z]
                        
                for z in range(17):
                    world[x][y][z][4] = maximum_depth_here
    ########
                    
    print "Finished intializing"
    
    
    new_data.close()
    return world


def update_environment_version3(world, year, day, hour, origin_offset, nitrate_source, future_date=False, climate_change=False):
    #update the environmental parameters based on the time of day and date using HYCOM model data
    #the data inside each parcel is : [temperature, salinity, nitrate, co2, max_depth, u, v, w_velocity, number of cells here]
    if check_future(year, day, future_date):
        new_data = get_HYCOM_data(year, day, hour, True)
    else:
        new_data = get_HYCOM_data(year, day, hour)
        
    if not new_data == None:
        #print "Updating environment..."
        temp_temp = new_data.variables['temperature'][0]
        temp_sal = new_data.variables['salinity'][0]
        temp_u = new_data.variables['u'][0]
        temp_v = new_data.variables['v'][0]
        temp_mld = new_data.variables['mixed_layer_depth'][0]
        #temp_depths = new_data.variables['Depth'][:]
        temp_lats = new_data.variables['Latitude'][:]
        temp_lons = new_data.variables['Longitude'][:]
##        temp2_mld = []
##        for x in range(17):
##            temp2_mld.append(temp_mld)
##        temp_mld = np_array(temp2_mld)
        temp_mld = np_array([temp_mld]*17)
        #convert missing data
##        temp_temp[temp_temp > 1E7] = np.nan
##        temp_sal[temp_sal > 1E7] = np.nan
##        temp_u[temp_u > 1E7] = np.nan
##        temp_v[temp_v > 1E7] = np.nan
##        temp_mld[temp_mld > 1E7] = np.nan

        #extract nitrate, co2, max_depth, and number of cells
##        size_x_dimension = len(world)
##        size_y_dimension = len(world[0])
##        size_z_dimension = len(world[0][0])
        size_x_dimension, size_y_dimension, size_z_dimension, num_variables = world.shape
        
        world = world.T
        world_nitrate = world[2]
        world_co2 = world[3]
        world_depth = world[4]
        
        
        #update the sunlight if necessary
        if day != MODIS_PAR_data[1]:
            if check_future(year, day, future_date): #this will make the PAR the same value across the model domain
                #what kind of skies do you want to simulate for the future (1500=65=sunny; 400=10=cloudy)
                sunlight = 50.  #this value is in different units (Einsteins/ m-2 d-1) while most of my model uses microeinsteins/ m-2 s-1
                temp_PAR = world[8][0]
                temp_PAR[temp_PAR < sunlight] = sunlight  #super sunny skies
                temp_PAR[temp_PAR > sunlight] = sunlight
                world_PAR = np_array([temp_PAR]*17)
            else: #assume it's running the past and PAR data is available
                temp_PAR = PAR_data_MODIS(year, day, world[8][0], origin_offset)
                world_PAR = np_array([temp_PAR]*17)
        else:
            world_PAR = world[8]
        
        #add the empty data to new_data
        temp_data = np_array([temp_temp, temp_sal, temp_u, temp_v, temp_mld])
        #print "Temp_data:", temp_data.shape

        #####starting nitrate placement; doing it here so that the nitrate values will be placed in the world before interpolation
        if nitrate_source == 'NODC':
            #print "Loading in nitrate concentrations from NODC data...",
            if hour == 0: #hard-coded to only update the nitrogen once per day
                nitrate_func = nitrate_from_NODC(day) #should return seven functions for depths 0, 10, 20, 30, 50, 75, 100
                #print "nitrate functions loaded"
                depths_nitrate = [0,10,20,30,50,75,100,125,150,200]
                for x in range(size_x_dimension):
                    for y in range(size_y_dimension):
                        for z in depths_nitrate:
                            if z < size_z_dimension:
                                if world_nitrate[z][y][x] < 1E6:
                                    temp_nitrate_value = 1.0 * (nitrate_func[depths_nitrate.index(z)](temp_lons[x+origin_offset[0]], temp_lats[y+origin_offset[1]])[0][0])
                                    if temp_nitrate_value < 0:
                                        temp_nitrate_value = 0.0
                                    elif temp_nitrate_value > 50:
                                        temp_nitrate_value = 50.0
                                    world_nitrate[z][y][x] = (world_nitrate[z][y][x] + temp_nitrate_value) * 0.5
                world_nitrate = update_nitrate_values_version2(world_nitrate, size_z_dimension)
            #print "Finished placing nitrate values"
        
        if nitrate_source == 'uniform':
            nitrate_concentrations(1, 'new', world)        
            world_nitrate = world[2]
        
        #####
##        if len(temp_data) > size_z_dimension:
##            for x in range(len(temp_data)-1, size_z_dimension-1, -1):
##                temp_data = np.delete(temp_data, x, 0)
        temp_data = np.rollaxis(temp_data, 2, 0) #move the y dimension to the front
        if len(temp_data) > size_y_dimension:
            temp_data = temp_data[origin_offset[1]:size_y_dimension+origin_offset[1]]  #cut to the region of interest
        temp_data = np.rollaxis(temp_data, 3, 0) #move the x dimension to the front
        if len(temp_data) > size_x_dimension:
            temp_data = temp_data[origin_offset[0]:size_x_dimension+origin_offset[0]]  #cut to the region of interest

        temp_data = np.rollaxis(temp_data, 3, 2)

        #print "Temp_data:", temp_data.shape
        temp_temp, temp_sal, temp_u, temp_v, temp_mld = np.split(temp_data.T, 5, axis=0)
        #print temp_temp.shape
        #print world.shape
        #print temp_sal.shape
        #print temp_mld.shape
        world = np_array([temp_temp[0], temp_sal[0], world_nitrate, world_co2, world_depth, temp_u[0], temp_v[0], temp_mld[0], world_PAR])
        #world[0][world[4] < 5] += 5 #trying to increase the temperature of shallow areas to avoid cells getting caught, remove/kill them instead
        world = world.T
        #print "World", world.shape
        #print "Interpolating values"
        #time1 = time.time()
        #world = interpolate_remaining_values_version2(world)
        #print "Time to interpolate:", time.time()-time1
        
        if climate_change:
            world = invoke_climate_change(world)
            
        new_data.close()
        
    return world

def check_future(year, day, future_date):
    #this will check the current date to see if it's past the 'future' point and should be predictive
    if year > future_date[0]:
        return True
    elif year < future_date[0]:
        return False
    else:
        if day >= future_date[1]:
            return True
    return False
    
def invoke_climate_change(world):
    #this function will increase the temperature of the water by a variable amount
    #higher increase at the surface and progessively less of an increase the deeper it goes
    amount_of_change = 2 #degrees C
    
    world = world.T
    for x in range(world.shape[1]):
        world[0, x, :, :] += (amount_of_change/(HYCOM_depth_values[x]+1.))
    world = world.T
    return world

                              

def interpolate_location_values_version3(world, locx, locy, locz):
    #the depth values from HYCOM are not evenly spaced and values in-between must be interpolated
    #the incoming locations must be in float format so that i can calculate the interpolated values
    #the data inside each parcel is : [temperature, salinity, nitrate, co2, max_depth, u, v, mixed_layer_depth, number of cells here]
    #HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
    global HYCOM_depth_values, HYCOM_depth_dictionary
    depth1 = HYCOM_depth_dictionary[int(locz)]
    depth2 = depth1 + 1
    if depth2 > 16:
        depth2 = 16

    #interpolate the depth    
    if depth1 == depth2:
        scaledz = [1, 0]
    else:
        temp_val1 = HYCOM_depth_values[depth1]
        temp_val2 = HYCOM_depth_values[depth2]
        total_difference = float(temp_val2 - temp_val1)
        scaledz = [1-((locz - temp_val1) / total_difference), 1 - ((temp_val2 - locz) / total_difference)]
    
    #interpolate the x range: the int(x), the int(x)+1, at int(y) and int(y)+1 if they exist
    locx_int = int(locx)
    locy_int = int(locy)
    if locx_int < 0 or locy_int < 0:
        locx_int = 1000     #set to high number in order to crash the program, a value less than zero should have been taken care of by now
        locy_int = 1000
    locx_int2 = locx_int + 1
    locy_int2 = locy_int + 1
    if locx_int2 >= len(world):
        locx_int2 = locx_int
    if locy_int2 >= len(world[0]):
        locy_int2 = locy_int
    
    loc000 = [locx_int, locy_int, depth1]
    loc001 = [locx_int2, locy_int, depth1]
    loc010 = [locx_int, locy_int2, depth1]
    loc011 = [locx_int2, locy_int2, depth1]
    loc100 = [locx_int, locy_int, depth2]
    loc101 = [locx_int2, locy_int, depth2]
    loc110 = [locx_int, locy_int2, depth2]
    loc111 = [locx_int2, locy_int2, depth2]

    try:
        if world[loc100[0], loc100[1], loc100[2]][4] == 0 and world[loc111[0], loc111[1], loc111[2]][4] == 0:
            locations = [loc000, loc001, loc010, loc011, loc000, loc001, loc010, loc011]
            #print "using one depth"
        else:
            locations = [loc000, loc001, loc010, loc011, loc100, loc101, loc110, loc111]
    except:
        print loc100, "    :    ", loc111
        if world[loc100[0], loc100[1], loc100[2]][4] == 0 and world[loc111[0], loc111[1], loc111[2]][4] == 0:
            locations = [loc000, loc001, loc010, loc011, loc000, loc001, loc010, loc011]
        
   #check the locations to see if they exist
    #valid_locations = [False,False,False,False,False,False,False,False]
    valid_locations = [not world[locations[indiv_loc][0], locations[indiv_loc][1], locations[indiv_loc][2]][0] < 1E7 for indiv_loc in xrange(8)]
    #for indiv_loc in xrange(8):
    #    #temp_loc = world[locations[indiv_loc][0], locations[indiv_loc][1], locations[indiv_loc][2]]
    #    if not world[locations[indiv_loc][0], locations[indiv_loc][1], locations[indiv_loc][2]][0] < 1E7:
    #        valid_locations[indiv_loc] = True

         
    scaledx1 = [1 - (locx - loc000[0]), 1 - (1 - (locx - loc000[0]))] #000, 001
    scaledy1 = [1 - (locy - loc000[1]), 1 - (1 - (locy - loc000[1]))] #how far is it from y1, y2
    scaledy2 = [1 - (locy - loc100[1]), 1 - (1 - (locy - loc100[1]))] #how far is it from y3, y4

    
    #interpolate the x's
    if not True in valid_locations: #will run this if all data points exist
        x1 = (world[locations[0][0], locations[0][1], locations[0][2]] * scaledx1[0]) + (world[locations[1][0], locations[1][1], locations[1][2]] * scaledx1[1])  #loc000+loc001 at depth 1
        x2 = (world[locations[2][0], locations[2][1], locations[2][2]] * scaledx1[0]) + (world[locations[3][0], locations[3][1], locations[3][2]] * scaledx1[1])  #loc010+loc011 at depth 1
        x3 = (world[locations[4][0], locations[4][1], locations[4][2]] * scaledx1[0]) + (world[locations[5][0], locations[5][1], locations[5][2]] * scaledx1[1])  #loc000+loc001 at depth2
        x4 = (world[locations[6][0], locations[6][1], locations[6][2]] * scaledx1[0]) + (world[locations[7][0], locations[7][1], locations[7][2]] * scaledx1[1])  #loc010+loc011 at depth2

    else: #if one of the data points is missing this section will run
        bad_data_point = np_array([0.,0.,0.,0.,0.,0.,0.,0.,0.])
        if not valid_locations[0]:
            if not valid_locations[1]:
                x1 = (world[locations[0][0], locations[0][1], depth1] * scaledx1[0]) + (world[locations[1][0], locations[1][1], depth1] * scaledx1[1])  #loc000+loc001 at depth 1
            else:
                x1 = world[locations[0][0], locations[0][1], depth1]
        else:
            if not valid_locations[1]:
                x1 = world[locations[1][0], locations[1][1], depth1]
            else:
                x1 = bad_data_point
                scaledy1[0] = 0
                scaledy1[1] = 1.

        if not valid_locations[2]:
            if not valid_locations[3]:
                x2 = (world[locations[2][0], locations[2][1], depth1] * scaledx1[0]) + (world[locations[3][0], locations[3][1], depth1] * scaledx1[1])  #loc010+loc011 at depth 1
            else:
                x2 = world[locations[2][0], locations[2][1], depth1]
        else:
            if not valid_locations[3]:
                x2 = world[locations[3][0], locations[3][1], depth1]
            else:
                x2 = bad_data_point
                scaledy1[1] = 0
                if scaledy1[0] > 0:
                    scaledy1[0] = 1.
                    
        if not valid_locations[4]:
            if not valid_locations[5]:
                x3 = (world[locations[4][0], locations[4][1], depth2] * scaledx1[0]) + (world[locations[5][0], locations[5][1], depth2] * scaledx1[1])  #loc000+loc001 at depth2
            else:
                x3 = world[locations[4][0], locations[4][1], depth2]
        else:
            if not valid_locations[5]:
                x3 = world[locations[5][0], locations[5][1], depth2]
            else:
                x3 = bad_data_point
                scaledy2[0] = 0
                scaledy2[1] = 1.

        if not valid_locations[6]:
            if not valid_locations[7]:
                x4 = (world[locations[6][0], locations[6][1], depth2] * scaledx1[0]) + (world[locations[7][0], locations[7][1], depth2] * scaledx1[1])  #loc010+loc011 at depth2
            else:
                x4 = world[locations[6][0], locations[6][1], depth2]
        else:
            if not valid_locations[7]:
                x4 = world[locations[7][0], locations[7][1], depth2]
            else:
                x4 = bad_data_point
                scaledy2[1] = 0
                if scaledy2[0] > 0:
                    scaledy2[0] = 1.

    
    
        #check to see if it's a bad location (land)
        if x1[0] == x2[0] == x3[0] == x4[0] == 0:
            #x1[0] = x2[1] = x3[2] = x4[3] = np.nan
            interp_location = [np.nan]*9
            #print "herehereherehereherehereherehereherehereherehereherehere"
            return interp_location, np.nan, np.nan

    #print 'loc:', locx, locy, locz
    #print 'x1', x1
    #print 'x2', x2
    #print 'x3', x3
    #print 'x4', x4

    #by now missing values should have been compensated for by using the complete value of the data point that was available
    y1 = (x1 * scaledy1[0]) + (x2 * scaledy1[1])  # value for loc000+loc001
    y2 = (x3 * scaledy2[0]) + (x4 * scaledy2[1])  # value for loc010+loc011
    #print 'y1', y1
    #print 'y2', y2
    #print 'scaled', scaledx1, scaledy1, scaledz

    if 0 < y1[0] < 1E7:
        if 0 < y2[0] < 1E7:
            interp_location = (y1 * scaledz[0]) + (y2 * scaledz[1])  #value for interpolated depth; this is the main location, below here is the gradient for nitrate
        else:
            interp_location = y1  #value for interpolated depth; this is the main location, below here is the gradient for nitrate
    else:
        
        if 0 < y2[0] < 1E7:
            interp_location = y2  #value for interpolated depth; this is the main location, below here is the gradient for nitrate
        else:
            #print "blah blah blah blah blah"
            interp_location = [np.nan]*9

    #get nitrate gradient here, the new way of calculating the nitrate concentration based on mixed layer depth allows me to calculate it easier
    #conc_above, conc_below = nitrate_concentrations(1, 'gradient_check', world, interp_location[7], locz, interp_location)
    #the previous method of calculating the nitrate gradient was based on a fixed formula; the method below should account for differences in nitrate fields
    conc_above = y1[2]
    conc_below = y2[2]
    #print interp_location, conc_above, conc_below
    return interp_location, conc_above, conc_below                    



def create_world_version3(env_width, env_length, depth, year, day, hour, origin_offset, nitrate_source, reverse, culture):
    depth = depth * -1
    if any(culture):
        temp_world = initialize_culture_environment(depth, culture)
    else:
        if reverse == 'forward':
            temp_world = initialize_environment_version3(env_width, env_length, depth, year, day, hour, origin_offset, nitrate_source)
        else:
            temp_world = initialize_environment_version3_reverse(env_width, env_length, depth, year, day, hour, origin_offset, nitrate_source)
    return temp_world, world_latitude, world_longitude

def update_world(world, year, day, hour, origin_offset, nitrate_source, reverse, culture, future_date=False, climate_change=False):
    #this function will be the sole function called from outside to update the environment based on what is being simulated
    #this should allow me to add a culture function without too much trouble
    if any(culture):
        update_culture_environment(culture)
    else:
        if reverse == 'forward':
            world = update_environment_version3(world, year, day, hour, origin_offset, nitrate_source, future_date, climate_change)
        else:
            world = update_environment_version3_reverse(world, year, day, hour, origin_offset, nitrate_source, climate_change)
    return world

def initialize_culture_environment(depth, culture_info):
    max_light, env_salinity, env_temperature, env_nutrients, lights_on_off = culture_info
    temp_depth = [env_temperature, env_salinity, env_nutrients, 0, depth, 0, 0, 0, max_light]
    world = np.array([[[temp_depth]*depth]])
    return world

def update_culture_environment(culture_info):
    pass
    
##def world_creator_version3(width, length, depth, year, day, hour, origin_offset, nitrate_source):
##    xdim = []
##    ydim = []
##    zdim = []
##    size_x_dimension = width
##    size_y_dimension = length
##    size_z_dimension = depth
##    #start the model with HYCOM model data
##    world = initialize_environment_version3(size_x_dimension, size_y_dimension, size_z_dimension, year, day, hour, origin_offset, nitrate_source)
##    return world, world_latitude, world_longitude
##
##def world_creator_version3_reverse(width, length, depth, year, day, hour, origin_offset, nitrate_source):
##    xdim = []
##    ydim = []
##    zdim = []
##    size_x_dimension = width
##    size_y_dimension = length
##    size_z_dimension = depth
##    #start the model with HYCOM model data
##    world = initialize_environment_version3_reverse(size_x_dimension, size_y_dimension, size_z_dimension, year, day, hour, origin_offset, nitrate_source)
##    return world, world_latitude, world_longitude

def initialize_environment_version3_reverse(size_x_dimension, size_y_dimension, size_z_dimension, year, day, hour, origin_offset, nitrate_source):
    #now I need to make the returned world only the 17 depths of the model
    #intervening points will be interpolated on the fly, this should cut down the memory footprint considerably
    #update the environmental parameters based on the time of day and date using HYCOM model data
    #the data inside each parcel is : [temperature, salinity, nitrate, co2, max_depth, u, v, mixed_layer_depth, daily PAR]
    global world_latitude, world_longitude
    HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
    new_data = get_HYCOM_data(year, day, hour)
    print "Initializing world..."
    
    temp_temp = new_data.variables['temperature'][0]
    temp_sal = new_data.variables['salinity'][0]
    temp_u = new_data.variables['u'][0]
    temp_v = new_data.variables['v'][0]
    
    #temp_depths = new_data.variables['Depth'][:]
    temp_lats = new_data.variables['Latitude'][:]
    temp_lons = new_data.variables['Longitude'][:]
    temp_mld = new_data.variables['mixed_layer_depth'][0]

    world_latitude = temp_lats[origin_offset[1]:origin_offset[1]+size_y_dimension]
    world_longitude = temp_lons[origin_offset[0]:origin_offset[0]+size_x_dimension]

    
    #temp2_mld = []
    #for x in range(17):
    #    temp2_mld.append(temp_mld)
    temp_mld = np_array([temp_mld]*17)
    
    zero_float = np.zeros((17, 350, 450), dtype=np.float)
    zero_int1 = np.zeros((17, 350, 450), dtype = np.int)
    zero_int2 = np.zeros((350, 450), dtype = np.float)
    zero_int2[zero_int2 == 0] = np.nan
    zero_int2 = np_array([PAR_data_MODIS(year, day, zero_int2, origin_offset, initialize=True, world_shape=(size_x_dimension, size_y_dimension, size_z_dimension))] * 17)
    ##zero_int2 is to keep track of cell counts, i'm going to add it to the salinity values to get the nans and then convert all other values back to zero
    ##zero_int2 = zero_int2 + temp_sal
    ##zero_int2[zero_int2 < 1000] = 0
    ##zero_int2[zero_int2 > 1000000] = np.nan
    
    #reverse the currents
    temp_u *= -1
    temp_v *= -1
    
    #convert missing data
##    temp_temp[temp_temp > 1E7] = np.nan
##    temp_sal[temp_sal > 1E7] = np.nan
##    temp_u[temp_u > 1E7] = np.nan
##    temp_v[temp_v > 1E7] = np.nan
##    temp_mld[temp_mld > 1E7] = np.nan

    
    co2 = np.zeros((17, 350, 450), dtype=np.float)
    co2[co2 == 0] = 1
    #begin creation of world by combining arrays, reordering axes, and adding in missing depth levels
    temp_world = np_array([temp_temp, temp_sal, zero_float, co2, zero_int1, temp_u, temp_v, temp_mld, zero_int2])
    temp_world = temp_world.T
    print temp_world.shape
    temp_world = np.rollaxis(temp_world, 2, 0) #move the depth (z-axis) to the front
    print temp_world.shape
    print "created the basic world framework"
    
    temp_world = np.rollaxis(temp_world, 2, 0)
    temp_world = np.rollaxis(temp_world, 2, 0)
    print temp_world.shape
    #end of world creation

    #begin world trimming to size and location requested
    if size_x_dimension < 450:
        temp_world = temp_world[origin_offset[0]:size_x_dimension+origin_offset[0]]
    if size_y_dimension < 350:
        temp_world = np.rollaxis(temp_world, 1, 0)
        temp_world = temp_world[origin_offset[1]:size_y_dimension+origin_offset[1]]
        temp_world = np.rollaxis(temp_world, 1, 0)
    print temp_world.shape

##    #####starting nitrate placement; doing it here so that the nitrate values will be placed in the world before interpolation
##    if nitrate_source == 'NODC':
##        print "Loading in nitrate concentrations from NODC data...",
##        nitrate_func = nitrate_from_NODC(day) #should return seven functions for depths 0, 10, 20, 30, 50, 75, 100
##        print "nitrate functions loaded"
##        depths_nitrate = [0,10,20,30,50,75,100,125,150,200]
##        for x in range(size_x_dimension):
##            for y in range(size_y_dimension):
##                for z in depths_nitrate:
##                    if z < temp_max_depth:
##                        if temp_world[x][y][z][0] < 1E6:
##                            temp_nitrate_value = 1. * (nitrate_func[depths_nitrate.index(z)](temp_lons[x+origin_offset[0]], temp_lats[y+origin_offset[1]])[0][0])
##                            if temp_nitrate_value < 0:
##                                temp_nitrate_value = 0.0
##                            elif temp_nitrate_value > 50:
##                                temp_nitrate_value = 50.0
##                            temp_world[x][y][z][2] = temp_nitrate_value
##        temp_world = interpolate_nitrate_values_version2(temp_world, temp_max_depth)
##        print "Finished placing nitrate values"
##    #####
    
    if nitrate_source == 'uniform':
        #add in a nitrate value
        print "Adding nitrate concentrations",
        nitrate_concentrations(1, 'new', temp_world)     #placed the nitrate gradient maker into a function so that the gradient can be modified in the middle of a run
        print "Done!"
    world = temp_world
    ######update max depth for each location
    #HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
    for x in range(size_x_dimension):
        for y in range(size_y_dimension):
            maximum_depth_here = 0
            if np_isnan(world[x][y][0][0]):
                #for z in range(size_z_dimension):
                for z in range(17):
                    world[x][y][z][4] = np.nan
                    world[x][y][z][8] = np.nan
            else:
                for z in range(17):
                    if not np_isnan(world[x][y][z][0]):
                        maximum_depth_here = HYCOM_depth_values[z]
                        
                for z in range(17):
                    world[x][y][z][4] = maximum_depth_here
    ########
                    
    print "Finished intializing"
    
    new_data.close()
    return world


def update_environment_version3_reverse(world, year, day, hour, origin_offset, nitrate_source, climate_change=False):
    #update the environmental parameters based on the time of day and date using HYCOM model data
    #the data inside each parcel is : [temperature, salinity, nitrate, co2, max_depth, u, v, w_velocity, number of cells here]
    new_data = get_HYCOM_data(year, day, hour)
    if not new_data == None:
        #print "Updating environment..."
        temp_temp = new_data.variables['temperature'][0]
        temp_sal = new_data.variables['salinity'][0]
        temp_u = new_data.variables['u'][0]
        temp_v = new_data.variables['v'][0]
        temp_mld = new_data.variables['mixed_layer_depth'][0]
        #temp_depths = new_data.variables['Depth'][:]
        #temp_lats = new_data.variables['Latitude'][:]
        #temp_lons = new_data.variables['Longitude'][:]
        temp_mld = np_array([temp_mld]*17)
        new_data.close()
        
        #reverse the currents
        temp_u *= -1
        temp_v *= -1
        
        #convert missing data
##        temp_temp[temp_temp > 1E7] = np.nan
##        temp_sal[temp_sal > 1E7] = np.nan
##        temp_u[temp_u > 1E7] = np.nan
##        temp_v[temp_v > 1E7] = np.nan
##        temp_mld[temp_mld > 1E7] = np.nan

        #extract nitrate, co2, max_depth, and number of cells
        #size_x_dimension = len(world)
        #size_y_dimension = len(world[0])
        #size_z_dimension = len(world[0][0])
        size_x_dimension, size_y_dimension, size_z_dimension, num_variables = world.shape
        
        world = world.T
        world_nitrate = world[2]
        world_co2 = world[3]
        world_depth = world[4]

        #update the sunlight if necessary
        if day != MODIS_PAR_data[1]:
            temp_PAR = PAR_data_MODIS(year, day, world[8][0], origin_offset, reverse=True)
            world_PAR = np_array([temp_PAR]*17)
        else:
            world_PAR = world[8]
            
        #add the empty data to new_data
        temp_data = np_array([temp_temp, temp_sal, temp_u, temp_v, temp_mld])
        #print "Temp_data:", temp_data.shape
        #insert blank depths
##        empty = np.zeros((450, 350, 5), dtype=np.float)
##        temp_data = temp_data.T
##        temp_data = np.rollaxis(temp_data, 2, 0)
##        temp_max_depth = 200
##        if 125 < size_z_dimension <= 150:
##            temp_max_depth = 150
##        elif 100 < size_z_dimension <= 125:
##            temp_max_depth = 125
##        elif 75 < size_z_dimension <= 100:
##            temp_max_depth = 100
##        elif 50 < size_z_dimension <= 75:
##            temp_max_depth = 75
##        elif 30 < size_z_dimension <= 50:
##            temp_max_depth = 50
##        elif 20 < size_z_dimension <= 30:
##            temp_max_depth = 30
##        elif 10 < size_z_dimension <= 20:
##            temp_max_depth = 20
##        elif 0 < size_z_dimension <= 10:
##            temp_max_depth = 10
##        for x in range(temp_max_depth): #insert the intervening depths; used to be only the z dimension given but with nitrate being added I have to do all 100 at first now
##            if x not in temp_depths:
##                temp_data = np.insert(temp_data, x, empty, 0)

##        #####starting nitrate placement; doing it here so that the nitrate values will be placed in the world before interpolation
        if nitrate_source == 'NODC':
            #print "Loading in nitrate concentrations from NODC data...",
            if hour == 0: #hard-coded to only update the nitrogen once per day
                nitrate_func = nitrate_from_NODC(day) #should return seven functions for depths 0, 10, 20, 30, 50, 75, 100
                #print "nitrate functions loaded"
                depths_nitrate = [0,10,20,30,50,75,100,125,150,200]
                for x in range(size_x_dimension):
                    for y in range(size_y_dimension):
                        for z in depths_nitrate:
                            if z < size_z_dimension:
                                if world_nitrate[z][y][x] < 1E6:
                                    temp_nitrate_value = 1.0 * (nitrate_func[depths_nitrate.index(z)](temp_lons[x+origin_offset[0]], temp_lats[y+origin_offset[1]])[0][0])
                                    if temp_nitrate_value < 0:
                                        temp_nitrate_value = 0.0
                                    elif temp_nitrate_value > 50:
                                        temp_nitrate_value = 50.0
                                    world_nitrate[z][y][x] = (world_nitrate[z][y][x] + temp_nitrate_value) * 0.5
                world_nitrate = update_nitrate_values_version2(world_nitrate, size_z_dimension)
            #print "Finished placing nitrate values"
        
        if nitrate_source == 'uniform':
            nitrate_concentrations(1, 'new', world)        
            world_nitrate = world[2]
        
        #####
##        if len(temp_data) > size_z_dimension:
##            for x in range(len(temp_data)-1, size_z_dimension-1, -1):
##                temp_data = np.delete(temp_data, x, 0)
        temp_data = np.rollaxis(temp_data, 2, 0) #move the y dimension to the front
        if len(temp_data) > size_y_dimension:
            temp_data = temp_data[origin_offset[1]:size_y_dimension+origin_offset[1]]  #cut to the region of interest
        temp_data = np.rollaxis(temp_data, 3, 0) #move the x dimension to the front
        if len(temp_data) > size_x_dimension:
            temp_data = temp_data[origin_offset[0]:size_x_dimension+origin_offset[0]]  #cut to the region of interest

        temp_data = np.rollaxis(temp_data, 3, 2)

        #print "Temp_data:", temp_data.shape
        temp_temp, temp_sal, temp_u, temp_v, temp_mld = np.split(temp_data.T, 5, axis=0)
        #print temp_temp.shape
        #print temp_sal.shape
        #print temp_mld.shape
        
        world = np_array([temp_temp[0], temp_sal[0], world_nitrate, world_co2, world_depth, temp_u[0], temp_v[0], temp_mld[0], world_PAR])
        world = world.T
        #print "World", world.shape
        #print "Interpolating values"
        #time1 = time.time()
        #world = interpolate_remaining_values_version2(world)
        #print "Time to interpolate:", time.time()-time1
    
    
        
    return world
