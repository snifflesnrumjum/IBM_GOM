# this is a working model of a single cell swimming vertically and being exposed to different nitrate and sunlight values
# depending on depth. the cell swims up or down depending on whether it needs more carbon or nitrate. if it has enough resources
# available, when the time is right (midnight in this case) the cell divides. right now it's set to make this cell the "poor daughter"
# cell with the minimum resources. it's presumed that the other cell gets the leftovers.
# i'm going to try to extend this model, in this file, to a bloom patch right now and see how it goes. from there i hope to incorporate a 3-dimensional
# component to check what happens in the field or test tube. ideally i would like to replicate the features seen in a test tube in
# the lab.  9/3/2009

#i would like to call this program as a module in another one and send each cell through this one

#3/8/2010 I am now converting this into a class. A cell will now be a separate class that can be called. This will allow for greater
# flexibility in program design as well as help to reduce the amount of work needed in order to add/remove a component

#v3 - adding the individual cell sunrise/sunset calculator to the model. This way the day length will vary over the year but also, sunrise and sunset for cells in very
#     different positions will be more accurate now. 

#import copy
#import pyximport; pyximport.install()  #for cython use
import random
import math
import time
import numpy
#import cellmodel_environment_HYCOM_v2 as cellmodel_environment_HYCOM  #to get the world update function
import cellmodel_environment_HYCOM_v4_predictive# as cellmodel_environment_HYCOM_v2
import ephem 
from scipy.integrate import quad
get_daylength = ephem.Observer()
#data_outputfile = open('D:/CJunk/Kbr_Model_data/temp/single_cell_data.txt', 'w')
#data_outputfile.close()
class Kbrevis_bloom:
    def define_bloom(self, information):
        self.bloom_id = information['bloom_id']
        self.cells = information['cells']
        self.parentbloom = information['parentbloom']
        self.when_bloom_created = information['when_bloom_created']
        self.location = [0,0,0]
        self.cell_growth = information['cell_growth']
        self.mutation_occurs = information['mutation_occurs']
        self.mutation_rate_cell_properties = information['mutation_rate_cell_properties']
        self.simulation_time_interval = information['simulation_time_interval']
        self.cell_death = information['cell_death']
   
    def avg_bloom_stats(self):
        avg_cn = []
        avg_cnmin = []
        avg_cnmax = []
        avg_n = []
        avg_nmin = []
        avg_nmax = []
        avg_ic = []
        avg_izt_threshold = []
        avg_pl = []
        avg_pma = []
        avg_pmb = []
        avg_rm = []
        avg_s250 = []
        avg_d = []
        avg_hc = []
        avg_growthrate = []
        avg_c_divide = []
        avg_n_divide = []
        avg_salinity_acc = []
        avg_salinity_pref = []
        avg_temperature_acc = []
        avg_temperature_pref = []
        
        
        for cell in self.cells:
            avg_cn.append(cell.cn)
            avg_cnmin.append(cell.cnmin)
            avg_cnmax.append(cell.cnmax)
            avg_n.append(cell.n)
            avg_nmin.append(cell.nmin)
            avg_nmax.append(cell.nmax)
            avg_ic.append(cell.Ic)
            avg_izt_threshold.append(cell.izt_threshold)
            avg_pl.append(cell.pl)
            avg_pma.append(cell.pma)
            avg_pmb.append(cell.pmb)
            avg_rm.append(cell.rm)
            avg_s250.append(cell.s250)
            avg_d.append(cell.d)
            avg_hc.append(cell.hc)
            avg_growthrate.append(cell.growthrate)
            avg_c_divide.append(cell.c_divide_threshold)
            avg_n_divide.append(cell.n_divide_threshold)
            avg_salinity_acc.append(cell.salinity_acclimated)
            avg_temperature_acc.append(cell.temperature_acclimated)
            avg_salinity_pref.append(cell.salinity_preferred)
            avg_temperature_pref.append(cell.temperature_preferred)
        
        avg_cn = numpy.mean(avg_cn)
        avg_cnmin = numpy.mean(avg_cnmin)
        avg_cnmax = numpy.mean(avg_cnmax)
        avg_n = numpy.mean(avg_n)
        avg_nmin = numpy.mean(avg_nmin)
        avg_nmax = numpy.mean(avg_nmax)
        avg_ic = numpy.mean(avg_ic)
        avg_izt_threshold = numpy.mean(avg_izt_threshold)
        avg_pl = numpy.mean(avg_pl)
        avg_pma = numpy.mean(avg_pma)
        avg_pmb = numpy.mean(avg_pmb)
        avg_rm = numpy.mean(avg_rm)
        avg_s250 = numpy.mean(avg_s250)
        avg_d = numpy.mean(avg_d)
        avg_hc = numpy.mean(avg_hc)
        avg_growthrate = numpy.mean(avg_growthrate)
        avg_c_divide = numpy.mean(avg_c_divide)
        avg_n_divide = numpy.mean(avg_n_divide)
        avg_salinity_acc = numpy.mean(avg_salinity_acc)
        avg_salinity_pref = numpy.mean(avg_salinity_pref)
        avg_temperature_acc = numpy.mean(avg_temperature_acc)
        avg_temperature_pref = numpy.mean(avg_temperature_pref)
        
        return [avg_cn, avg_cnmin, avg_cnmax, avg_n, avg_nmin, avg_nmax, avg_ic, avg_izt_threshold, avg_pl, avg_pma, avg_pmb, avg_rm, avg_s250, avg_d, avg_hc, avg_growthrate, avg_c_divide,
                avg_n_divide, avg_salinity_acc, avg_salinity_pref, avg_temperature_acc, avg_temperature_pref]
        
    def what_month_is_it(self, year, julian_day):
        if year % 4 == 0:
            month_ranges = [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        else:
            month_ranges = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
        month_to_return = 0
        while julian_day > month_ranges[month_to_return]:
            month_to_return += 1
        if julian_day > 31:
            day_of_month = julian_day - month_ranges[month_to_return]
        else:
            day_of_month = julian_day
        #adding 1 to the return value for non zero index used in datetime module for the month
        return str(month_to_return + 1) + '/' + str(day_of_month)
            
    def process_one_day(self, bloom_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, world_update_hours, origin_offset, nitrate_source, world_latitude, world_longitude, reverse, culture, simulation_time_interval, carrying_capacity, future_date): #used to have env_temperature and env_salinity instead of world
        HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
        time1 = time.time()
        time_of_day = 3.00 #starts at 3am, adds 0.05 (3 min) during each time step until it cycles back to 3am; CHANGED to accept different time steps (2/17/2014)
        time_of_day2 = int(round(time_of_day*100))
        orig_cell_conc = 0
        day_night = 1 # 0 == day; 1 == night
        current_hour -= 1 #because the first time step below is going to increment the current hour up to 3
        cell_conc_above = [0]*abs(max_depth)
        size_x_dimension = len(world)
        size_y_dimension = len(world[0])
        temp_xloc, temp_yloc, temp_zloc = [-1, -1, -1]   
        for time_step in xrange(int(1440 / (simulation_time_interval / 60.))):
            #time2 = time.time()        
            if time_of_day >= 23.96:
                time_of_day = 0.00
                current_hour = 0
                current_day += 1
                if current_day > 365:
                    current_year += 1
                    current_day = 1
            time_of_day2 = int(round(time_of_day*100))
            if time_of_day2 % 100 == 0:
                if time_of_day > 0.0:
                    current_hour += 1
                if current_hour in world_update_hours:
                    world = cellmodel_environment_HYCOM_v4_predictive.update_world(world, current_year, current_day, current_hour, origin_offset, nitrate_source, reverse, culture, future_date, climate_change=False)
                    
##            if 5.96 < time_of_day < 18.0:    #check to see what time it is; if it's daytime, have the day_night variable indicate this
##                day_night = 0
##            else:
##                day_night = 1
            day_night = None
            
            current_date = str(current_year) + '/' + self.what_month_is_it(current_year, current_day)
            for cell in self.cells:
                
                if 2.96 < time_of_day < 3.04:
                    cell.liveordie(carrying_capacity, False)
                
                if cell.cn <= 0:
                    continue
                
                orig_xloc = cell.x
                orig_yloc = cell.y
                orig_zloc = cell.z
                try:
                    temp_environ_info = cellmodel_environment_HYCOM_v4_predictive.interpolate_location_values_version3(world, orig_xloc, orig_yloc, orig_zloc)
                    environ_info = temp_environ_info[0][:5]
                    if True in numpy.isnan(environ_info):
                        #print "first one"
                        #print cell.cn, orig_xloc, orig_yloc, orig_zloc, '   ', environ_info
                        break
                        cell.cn = 0
                        #print orig_xloc, orig_yloc, orig_zloc, environ_info
                        continue
                    #print "Actual:", temp_environ_info
                    #print "Surface:", cellmodel_environment_HYCOM_v2.interpolate_location_values_version3(world, orig_xloc, orig_yloc, 0)
                    #print
                except:
                    temp_environ_info = cellmodel_environment_HYCOM_v4_predictive.interpolate_location_values_version3(world, orig_xloc, orig_yloc, orig_zloc)
                    #print "hi"
                    cell.cn = 0
                    continue
                environ_info = list(environ_info)
                environ_info.extend([temp_environ_info[1], temp_environ_info[2]])
                #environ_info.append(temp_environ_info[2])
                environ_info.extend(list(temp_environ_info[0][5:]))
                
                cell.time_step_environment([time_of_day, day_night, cell_conc_above, orig_cell_conc, max_depth, max_light, environ_info], world_latitude, world_longitude, current_date, culture)
                temp_xloc = cell.x
                temp_yloc = cell.y
                temp_zloc = cell.z
                
                
                if temp_xloc >= size_x_dimension:
                    #cell.x = len(world)-0.1 
                    cellx_int = int(temp_xloc) - 1 
                elif temp_xloc <= 0:
                    #cell.x = 0.01
                    cellx_int = 0
                else:
                    cellx_int = int(temp_xloc)
                    
                if temp_yloc >= size_y_dimension:
                    #cell.y = len(world[0])-0.1
                    celly_int = int(temp_yloc) - 1
                elif temp_yloc <= 0:
                    #cell.y = 0.01
                    celly_int = 0
                else:
                    celly_int = int(temp_yloc)

                cellz_int = int(cell.z)
                max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
                if max_depth_at_x_y_location < 1E7: #not np_isnan(max_depth_at_x_y_location):
                    if temp_zloc < abs(max_depth) and 0 < temp_zloc < max_depth_at_x_y_location:
                            cellz_int = int(cell.z)       
                    elif temp_zloc >= max_depth_at_x_y_location and max_depth_at_x_y_location != 0:
                        cell.z = max_depth_at_x_y_location - 0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc > max_depth_at_x_y_location and max_depth_at_x_y_location == 0:
                        if cell.z < max_depth_at_x_y_location + 4:
                            cellz_int = int(cell.z)
                        else:
                            cell.z = max_depth_at_x_y_location + 4
                            cellz_int = int(cell.z)
                    elif temp_zloc >= abs(max_depth):
                        cell.z = abs(max_depth)-0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc < 0:
                        cell.z = 0.01
                        cellz_int = 0
                    else:
                        cellz_int = int(cell.z)
                
                zdepth = 0
                while cellz_int >= HYCOM_depth_values[zdepth+1]:
                    zdepth += 1
                cellz_int = zdepth

                orig_zdepth = 0
                while orig_zloc >= HYCOM_depth_values[orig_zdepth+1]:
                    orig_zdepth += 1
                
                if not world[cellx_int, celly_int, cellz_int][0] < 1E7: #checking for NaN here 
                    try:
                        if world[cellx_int, celly_int, orig_zdepth][0] < 1E7:
                            cell.z = orig_zloc
                        elif max_depth_at_x_y_location < 1E7:
                            if world[cellx_int][celly_int][HYCOM_depth_values.index(max_depth_at_x_y_location)][0] < 1E7: #checking for NaN here
                                cell.z = max_depth_at_x_y_location
                                print "Cell could not move:", orig_xloc, orig_yloc, orig_zloc
                                print "New locations:", cellx_int, celly_int, cellz_int
                            cell.same_location += 1
                        elif world[cellx_int][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.y = orig_yloc
                            if cell.z > world[cellx_int][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[cellx_int][int(orig_yloc)][cellz_int][4]
                            cell.same_location = 0
                        elif world[int(orig_xloc)][celly_int][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            if cell.z > world[int(orig_xloc)][celly_int][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][celly_int][cellz_int][4]
                            cell.same_location = 0
                        elif world[int(orig_xloc)][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            if cell.z > world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]
                            cell.same_location = 0
                        else:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            cell.z = orig_zloc
                            cell.same_location += 1
                            
                    except:
                        print cellx_int, celly_int, cellz_int, cell.x, cell.y, cell.z, orig_xloc, orig_yloc, orig_zloc
                        
                if not any(culture):            
                    if cell.x >= size_x_dimension - 1:
                        cell.cn = -9
                    elif cell.x <= 0:
                        cell.cn = -9
                    elif cell.y >= size_y_dimension - 1:
                        cell.cn = -9
                    elif cell.y <= 0:
                        cell.cn = -9                    
                        #world[temp_xloc][temp_yloc][temp_zloc][8]-= 1
                        #world[int(cell.x)][int(cell.y)][int(cell.z)][8]+=1
            if time_of_day2 % 100 == 0:
                print time_of_day,
                print round(temp_xloc, 3), round(temp_yloc, 3), round(temp_zloc, 3) #round(self.cells[0].cn, 3),
                #print round(self.cells[-1].n, 3),  round(self.cells[-1].max_light_PAR_value, 3), round(self.cells[-1].e3, 3), round(self.cells[-1].hc, 3),
                #print round(self.cells[-1].Ic, 3), round(self.cells[-1].Ih, 3), round(self.cells[-1].pm, 3), round(self.cells[-1].pmd, 3), round(self.cells[-1].A, 3), round(self.cells[-1].Q, 3), round(self.cells[-1].cnphoto, 3)
                #data_outputfile = open('D:/CJunk/Kbr_Model_data/culture/single_cell_data.txt', 'a')
                #data_outputfile.write(str([time_of_day, temp_zloc*-1, self.cells[0].cn, self.cells[-1].n, self.cells[-1].I, self.cells[-1].e3, self.cells[-1].hc, self.cells[-1].Ih])[:-1] + '\n')
                #data_outputfile.close()
                
            time_of_day += simulation_time_interval / 3600. #used to be 0.05
        ######this section will run the liveordie code to determine which cells to remove or reproduce
        addcells = []
        remove_index = 0
        global cell_growth
        for cell in self.cells:
            if cell.cn <= 0 and self.cell_death:
                self.cells[self.cells.index(cell)] = -9
            elif cell.cn <= 0 and not self.cell_death:
                cell.cn = cell.cnmin * 1.1 #reanimate the dead cell with a small amount of carbon, but not enough to divide 
            elif cell.divide:
                if reverse == 'forward' and self.cell_growth:
                    cell.divide = False
                    tempcell = self.Kbrcell(cell, 1, cell.species_number)
                    if self.mutation_occurs == True:
                        #mutation(tempcell)  #this puts the new cell through the mutation process and leaves the original one alone
                        tempcell.mutate_loci()
                    addcells.append(tempcell)
                    
        self.cells.sort()        
        while self.cells[remove_index-1] == -9:
            remove_index -= 1
        print "Dead cells:", self.cells.count(-9), remove_index
        if remove_index < 0 and self.cell_death:
            self.cells = self.cells[:remove_index]
        
        self.cells.extend(addcells)
        print "New cells: ", len(addcells)
        #####################################end liveordie section
        time4 = time.time()
        print "One day took: ", time4-time1
        return current_year, current_day, current_hour, world
        
    def process_one_day_multi(self, bloom_slice, bloom_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, world_update_hours, origin_offset, nitrate_source, world_latitude, world_longitude, reverse, culture, simulation_time_interval, carrying_capacity, future_date): #used to have env_temperature and env_salinity instead of world
        HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
        time1 = time.time()
        time_of_day = 3.00 #starts at 3am, adds 0.05 (3 min) during each time step until it cycles back to 3am
        time_of_day2 = int(round(time_of_day*100))
        orig_cell_conc = 0
        cell_conc_above = [0]*abs(max_depth)
        size_x_dimension = len(world)
        size_y_dimension = len(world[0])
        day_night = 1 # 0 == day; 1 == night
        current_hour -= 1 #because the first time step below is going to increment the current hour up to 3
        temp_xloc, temp_yloc, temp_zloc = [-1, -1, -1] 
        for time_step in range(int(1440 / (simulation_time_interval / 60.))):
            if time_of_day >= 23.96:
                time_of_day = 0.00
                current_hour = 0
                current_day += 1
                if current_day > 365:
                    current_year += 1
                    current_day = 1
            time_of_day2 = int(round(time_of_day*100))
            if time_of_day2 % 100 == 0:
                if time_of_day > 0.0:
                    current_hour += 1
                if current_hour in world_update_hours:
                    world = cellmodel_environment_HYCOM_v4_predictive.update_world(world, current_year, current_day, current_hour, origin_offset, nitrate_source, reverse, culture, future_date, climate_change=False)
                    
            day_night = None
            current_date = str(current_year) + '/' + self.what_month_is_it(current_year, current_day)
            #bloom_slice = [self.process_indiv_cell(cell, time_of_day, carrying_capacity) for cell in bloom_slice]
            for cell in bloom_slice:
                if 2.96 < time_of_day < 3.04:
                    cell.liveordie(carrying_capacity, False) 
                if cell.cn <= 0:
                    continue
                orig_xloc = cell.x
                orig_yloc = cell.y
                orig_zloc = cell.z
                try:
                    temp_environ_info = cellmodel_environment_HYCOM_v4_predictive.interpolate_location_values_version3(world, orig_xloc, orig_yloc, orig_zloc)
                    environ_info = temp_environ_info[0][:5]
                    if True in numpy.isnan(environ_info):
                        cell.cn = 0
                        #print orig_xloc, orig_yloc, orig_zloc, environ_info
                        continue
                except:
                    cell.cn = 0
                    continue
                environ_info = list(environ_info)
                environ_info.extend([temp_environ_info[1], temp_environ_info[2]])
                #environ_info.append(temp_environ_info[2])
                environ_info.extend(list(temp_environ_info[0][5:]))

                cell.time_step_environment([time_of_day, day_night, cell_conc_above, orig_cell_conc, max_depth, max_light, environ_info], world_latitude, world_longitude, current_date, culture)
                temp_xloc = cell.x
                temp_yloc = cell.y
                temp_zloc = cell.z
                
                
                if temp_xloc >= size_x_dimension:
                    #cell.x = len(world)-0.1 
                    cellx_int = int(temp_xloc) - 1
                elif temp_xloc <= 0:
                    #cell.x = 0.01
                    cellx_int = 0
                else:
                    cellx_int = int(temp_xloc)
                    
                if temp_yloc >= size_y_dimension:
                    #cell.y = len(world[0])-0.1
                    celly_int = int(temp_yloc) - 1
                elif temp_yloc <= 0:
                    #cell.y = 0.01
                    celly_int = 0
                else:
                    celly_int = int(temp_yloc)

                cellz_int = int(cell.z)
                try:
                    max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
                except:
                    print cellx_int, celly_int
                    max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
                if max_depth_at_x_y_location < 1E7: #not np_isnan(max_depth_at_x_y_location):
                    if temp_zloc < abs(max_depth) and 0 < temp_zloc < max_depth_at_x_y_location:
                            cellz_int = int(cell.z)       
                    elif temp_zloc >= max_depth_at_x_y_location and max_depth_at_x_y_location != 0:
                        cell.z = max_depth_at_x_y_location - 0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc > max_depth_at_x_y_location and max_depth_at_x_y_location == 0:
                        if cell.z < max_depth_at_x_y_location + 4:
                            cellz_int = int(cell.z)
                        else:
                            cell.z = max_depth_at_x_y_location + 4
                            cellz_int = int(cell.z)
                    elif temp_zloc >= abs(max_depth):
                        cell.z = abs(max_depth)-0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc < 0:
                        cell.z = 0.01
                        cellz_int = 0
                    else:
                        cellz_int = int(cell.z)
                
                zdepth = 0
                while cellz_int >= HYCOM_depth_values[zdepth+1]:
                    zdepth += 1
                cellz_int = zdepth

                orig_zdepth = 0
                while orig_zloc >= HYCOM_depth_values[orig_zdepth+1]:
                    orig_zdepth += 1
                
                if not world[cellx_int, celly_int, cellz_int][0] < 1E7: #checking for NaN here 
                    try:
                        if world[cellx_int, celly_int, orig_zdepth][0] < 1E7:
                            cell.z = orig_zloc
                            cell.same_location = 0
                        elif max_depth_at_x_y_location < 1E7:
                            if world[cellx_int][celly_int][HYCOM_depth_values.index(max_depth_at_x_y_location)][0] < 1E7: #checking for NaN here
                                cell.z = max_depth_at_x_y_location
                            cell.same_location = 0
                        elif world[cellx_int][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.y = orig_yloc
                            if cell.z > world[cellx_int][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[cellx_int][int(orig_yloc)][cellz_int][4]
                            cell.same_location = 0
                        elif world[int(orig_xloc)][celly_int][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            if cell.z > world[int(orig_xloc)][celly_int][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][celly_int][cellz_int][4]
                            cell.same_location = 0
                        elif world[int(orig_xloc)][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            if cell.z > world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]
                            cell.same_location = 0
                        else:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            cell.z = orig_zloc
                            cell.same_location += 1
                    except:
                        print cellx_int, celly_int, cellz_int, cell.x, cell.y, cell.z, orig_xloc, orig_yloc, orig_zloc
                        
                            
                if cell.x >= size_x_dimension - 1:
                    cell.cn = -9
                elif cell.x <= 0:
                    cell.cn = -9
                elif cell.y >= size_y_dimension - 1:
                    cell.cn = -9
                elif cell.y <= 0:
                    cell.cn = -9                    
                    #world[temp_xloc][temp_yloc][temp_zloc][8]-= 1
                    #world[int(cell.x)][int(cell.y)][int(cell.z)][8]+=1
            #if time_of_day2 % 200 == 0:
                #print time_of_day,
                #print temp_xloc, temp_yloc, temp_zloc,
                #print cell_conc_above, current_day, current_hour    
            time_of_day += (simulation_time_interval / 3600.)
        ######this section will run the liveordie code to determine which cells to remove or reproduce
        addcells = []
        remove_index = 0
        global cell_growth
        for cell in bloom_slice:
            if cell.cn <= 0 and self.cell_death:
                bloom_slice[bloom_slice.index(cell)] = -9
            elif cell.cn <= 0 and not self.cell_death:
                cell.cn = cell.cnmin * 1.1 #reanimate the dead cell with a small amount of carbon, but not enough to divide 
            elif cell.divide:
                if reverse == 'forward' and self.cell_growth:
                    cell.divide = False
                    tempcell = self.Kbrcell(cell, 1, cell.species_number)
                    if self.mutation_occurs == True:
                        #mutation(tempcell)  #this puts the new cell through the mutation process and leaves the original one alone
                        tempcell.mutate_loci()
                    addcells.append(tempcell)
                    
        bloom_slice.sort()        
        while bloom_slice[remove_index-1] == -9:
            remove_index -= 1
        print "Dead cells:", bloom_slice.count(-9), remove_index
        if remove_index < 0:
            bloom_slice = bloom_slice[:remove_index]
        bloom_slice.extend(addcells)
        print "New cells: ", len(addcells)
        #####################################end liveordie section
        time4 = time.time()
        print "One day took: ", time4-time1
        return bloom_slice
    
    def process_indiv_cell(self, cell, time_of_day, carrying_capacity): #still needs work but might be a faster way to process cells
        if 2.96 < time_of_day < 3.04:
            cell.liveordie(carrying_capacity, False) 
        if cell.cn <= 0:
            pass
        else:
            orig_xloc = cell.x
            orig_yloc = cell.y
            orig_zloc = cell.z
            try:
                temp_environ_info = cellmodel_environment_HYCOM_v4_predictive.interpolate_location_values_version3(world, orig_xloc, orig_yloc, orig_zloc)
                environ_info = temp_environ_info[0][:5]
                if True in numpy.isnan(environ_info):
                    cell.cn = 0
                    #print orig_xloc, orig_yloc, orig_zloc, environ_info
                    #continue
            except:
                cell.cn = 0
                #continue
            environ_info = list(environ_info)
            environ_info.extend([temp_environ_info[1], temp_environ_info[2]])
            #environ_info.append(temp_environ_info[2])
            environ_info.extend(list(temp_environ_info[0][5:]))

            cell.time_step_environment([time_of_day, day_night, cell_conc_above, orig_cell_conc, max_depth, max_light, environ_info], world_latitude, world_longitude, current_date, culture)
            temp_xloc = cell.x
            temp_yloc = cell.y
            temp_zloc = cell.z
            
            
            if temp_xloc >= size_x_dimension:
                #cell.x = len(world)-0.1 
                cellx_int = int(temp_xloc) - 1
            elif temp_xloc <= 0:
                #cell.x = 0.01
                cellx_int = 0
            else:
                cellx_int = int(temp_xloc)
                
            if temp_yloc >= size_y_dimension:
                #cell.y = len(world[0])-0.1
                celly_int = int(temp_yloc) - 1
            elif temp_yloc <= 0:
                #cell.y = 0.01
                celly_int = 0
            else:
                celly_int = int(temp_yloc)

            cellz_int = int(cell.z)
            try:
                max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
            except:
                print cellx_int, celly_int
                max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
            if max_depth_at_x_y_location < 1E7: #not np_isnan(max_depth_at_x_y_location):
                if temp_zloc < abs(max_depth) and 0 < temp_zloc < max_depth_at_x_y_location:
                        cellz_int = int(cell.z)       
                elif temp_zloc >= max_depth_at_x_y_location and max_depth_at_x_y_location != 0:
                    cell.z = max_depth_at_x_y_location - 0.01
                    cellz_int = int(cell.z)
                elif temp_zloc > max_depth_at_x_y_location and max_depth_at_x_y_location == 0:
                    if cell.z < max_depth_at_x_y_location + 4:
                        cellz_int = int(cell.z)
                    else:
                        cell.z = max_depth_at_x_y_location + 4
                        cellz_int = int(cell.z)
                elif temp_zloc >= abs(max_depth):
                    cell.z = abs(max_depth)-0.01
                    cellz_int = int(cell.z)
                elif temp_zloc < 0:
                    cell.z = 0.01
                    cellz_int = 0
                else:
                    cellz_int = int(cell.z)
            
            zdepth = 0
            while cellz_int >= HYCOM_depth_values[zdepth+1]:
                zdepth += 1
            cellz_int = zdepth

            orig_zdepth = 0
            while orig_zloc >= HYCOM_depth_values[orig_zdepth+1]:
                orig_zdepth += 1
            
            if not world[cellx_int, celly_int, cellz_int][0] < 1E7: #checking for NaN here 
                try:
                    if world[cellx_int, celly_int, orig_zdepth][0] < 1E7:
                        cell.z = orig_zloc
                        cell.same_location = 0
                    elif max_depth_at_x_y_location < 1E7:
                        if world[cellx_int][celly_int][HYCOM_depth_values.index(max_depth_at_x_y_location)][0] < 1E7: #checking for NaN here
                            cell.z = max_depth_at_x_y_location
                        cell.same_location = 0
                    elif world[cellx_int][int(orig_yloc)][cellz_int][0] < 1E7:
                        cell.y = orig_yloc
                        if cell.z > world[cellx_int][int(orig_yloc)][cellz_int][4]:
                            cell.z = world[cellx_int][int(orig_yloc)][cellz_int][4]
                        cell.same_location = 0
                    elif world[int(orig_xloc)][celly_int][cellz_int][0] < 1E7:
                        cell.x = orig_xloc
                        if cell.z > world[int(orig_xloc)][celly_int][cellz_int][4]:
                            cell.z = world[int(orig_xloc)][celly_int][cellz_int][4]
                        cell.same_location = 0
                    elif world[int(orig_xloc)][int(orig_yloc)][cellz_int][0] < 1E7:
                        cell.x = orig_xloc
                        cell.y = orig_yloc
                        if cell.z > world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]:
                            cell.z = world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]
                        cell.same_location = 0
                    else:
                        cell.x = orig_xloc
                        cell.y = orig_yloc
                        cell.z = orig_zloc
                        cell.same_location += 1
                except:
                    print cellx_int, celly_int, cellz_int, cell.x, cell.y, cell.z, orig_xloc, orig_yloc, orig_zloc
                    
                        
            if cell.x >= size_x_dimension - 1:
                cell.cn = -9
            elif cell.x <= 0:
                cell.cn = -9
            elif cell.y >= size_y_dimension - 1:
                cell.cn = -9
            elif cell.y <= 0:
                cell.cn = -9
        return cell
    def process_one_day_reverse(self, bloom_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, world_update_hours, origin_offset, nitrate_source, world_latitude, world_longitude, reverse, culture, simulation_time_interval, carrying_capacity): #used to have env_temperature and env_salinity instead of world
        HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
        time1 = time.time()
        time_of_day = 3.00 #starts at 3am, adds 0.05 (3 min) during each time step until it cycles back to 3am
        time_of_day2 = int(round(time_of_day*100))
        orig_cell_conc = 0
        cell_conc_above = [0]*abs(max_depth)
        size_x_dimension = len(world)
        size_y_dimension = len(world[0])
        day_night = 1 # 0 == day; 1 == night
        current_hour -= 1 #because the first time step below is going to increment the current hour up to 3
        temp_xloc, temp_yloc, temp_zloc = [-1, -1, -1] 
        for time_step in range(int(1440 / (simulation_time_interval / 60.))):
                  
            if time_of_day < 0.0:
                time_of_day = 23.95
                current_hour = 23
                current_day -= 1
                if current_day < 1:
                    current_year -= 1
                    current_day = 365
            time_of_day2 = int(round(time_of_day*100))
            if time_of_day2 % 100 == 0:
                if time_of_day < 23.95:
                    current_hour -= 1
                if current_hour in world_update_hours:
                    current_hr_index = world_update_hours.index(current_hour) - 1
                    if current_hr_index < 0:
                        world = cellmodel_environment_HYCOM_v4_predictive.update_world(world, current_year, current_day-1, world_update_hours[current_hr_index], origin_offset, nitrate_source, reverse, culture)
                    else:
                        world = cellmodel_environment_HYCOM_v4_predictive.update_world(world, current_year, current_day, world_update_hours[current_hr_index], origin_offset, nitrate_source, reverse, culture)    
##            if 5.96 < time_of_day < 18.0:    #check to see what time it is; if it's daytime, have the day_night variable indicate this
##                day_night = 0
##            else:
##                day_night = 1
            day_night = None
            
            current_date = str(current_year) + '/' +  self.what_month_is_it(current_year, current_day)
            for cell in self.cells:
                if 2.96 < time_of_day < 3.04:
                    cell.liveordie(carrying_capacity, reverse) 
                if cell.cn <= 0:
                    continue
                orig_xloc = cell.x
                orig_yloc = cell.y
                orig_zloc = cell.z
                try:
                    temp_environ_info = cellmodel_environment_HYCOM_v4_predictive.interpolate_location_values_version3(world, orig_xloc, orig_yloc, orig_zloc)
                    environ_info = temp_environ_info[0][:5]
                    if True in numpy.isnan(environ_info):
                        cell.cn = 0
                        #print orig_xloc, orig_yloc, orig_zloc, environ_info
                        continue
                    #if cell.cellid == 2:
                        #print "Cell location:", orig_xloc, orig_yloc, orig_zloc
                        #print "Actual:", temp_environ_info
                        #print "Surface:", cellmodel_environment_HYCOM_v2.interpolate_location_values_version3(world, orig_xloc, orig_yloc, 0)
                        #print
                except:
                    cell.cn = 0
                    continue
                environ_info = list(environ_info)
                environ_info.extend([temp_environ_info[1], temp_environ_info[2]])
                #environ_info.append(temp_environ_info[2])
                environ_info.extend(list(temp_environ_info[0][5:]))
                #environ_info should have [temperature, salinity, nitrate, co2, max_depth, no3_above, no3_below, u, v, mixed layer depth, number cells] 
                cell.time_step_environment([time_of_day, day_night, cell_conc_above, orig_cell_conc, max_depth, max_light, environ_info], world_latitude, world_longitude, current_date, culture, True)
                temp_xloc = cell.x
                temp_yloc = cell.y
                temp_zloc = cell.z
                
                
                if temp_xloc >= size_x_dimension:
                    #cell.x = len(world)-0.1 
                    cellx_int = int(temp_xloc) - 1
                elif temp_xloc <= 0:
                    #cell.x = 0.01
                    cellx_int = 0
                else:
                    cellx_int = int(temp_xloc)
                    
                if temp_yloc >= size_y_dimension:
                    #cell.y = len(world[0])-0.1
                    celly_int = int(temp_yloc) - 1
                elif temp_yloc <= 0:
                    #cell.y = 0.01
                    celly_int = 0
                else:
                    celly_int = int(temp_yloc)

                cellz_int = int(cell.z)
                max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
                if max_depth_at_x_y_location < 1E7: #not np_isnan(max_depth_at_x_y_location):
                    if temp_zloc < abs(max_depth) and 0 < temp_zloc < max_depth_at_x_y_location:
                            cellz_int = int(cell.z)       
                    elif temp_zloc >= max_depth_at_x_y_location and max_depth_at_x_y_location != 0:
                        cell.z = max_depth_at_x_y_location - 0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc > max_depth_at_x_y_location and max_depth_at_x_y_location == 0:
                        if cell.z < max_depth_at_x_y_location + 4:
                            cellz_int = int(cell.z)
                        else:
                            cell.z = max_depth_at_x_y_location + 4
                            cellz_int = int(cell.z)
                    elif temp_zloc >= abs(max_depth):
                        cell.z = abs(max_depth)-0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc < 0:
                        cell.z = 0.01
                        cellz_int = 0
                    else:
                        cellz_int = int(cell.z)
                
                zdepth = 0
                while cellz_int >= HYCOM_depth_values[zdepth+1]:
                    zdepth += 1
                cellz_int = zdepth

                orig_zdepth = 0
                while orig_zloc >= HYCOM_depth_values[orig_zdepth+1]:
                    orig_zdepth += 1
                
                if not world[cellx_int, celly_int, cellz_int][0] < 1E7: #checking for NaN here 
                    try:
                        if world[cellx_int, celly_int, orig_zdepth][0] < 1E7:
                            cell.z = orig_zloc
                            cell.same_location = 0
                        elif max_depth_at_x_y_location < 1E7:
                            if world[cellx_int][celly_int][HYCOM_depth_values.index(max_depth_at_x_y_location)][0] < 1E7: #checking for NaN here
                                cell.z = max_depth_at_x_y_location
                                cell.same_location = 0
                        elif world[cellx_int][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.y = orig_yloc
                            if cell.z > world[cellx_int][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[cellx_int][int(orig_yloc)][cellz_int][4]
                            cell.same_location = 0
                        elif world[int(orig_xloc)][celly_int][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            if cell.z > world[int(orig_xloc)][celly_int][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][celly_int][cellz_int][4]
                            cell.same_location = 0
                        elif world[int(orig_xloc)][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            if cell.z > world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]
                            cell.same_location = 0
                        else:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            cell.z = orig_zloc
                            cell.same_location += 1
                    except:
                        print cellx_int, celly_int, cellz_int, cell.x, cell.y, cell.z, orig_xloc, orig_yloc, orig_zloc
                        
                            
                if cell.x >= size_x_dimension - 1:
                    cell.cn = -9
                elif cell.x <= 0:
                    cell.cn = -9
                elif cell.y >= size_y_dimension - 1:
                    cell.cn = -9
                elif cell.y <= 0:
                    cell.cn = -9                    
                    #world[temp_xloc][temp_yloc][temp_zloc][8]-= 1
                    #world[int(cell.x)][int(cell.y)][int(cell.z)][8]+=1
                #if cell.x == orig_xloc and cell.y == orig_yloc and cell.z == orig_zloc:
                #    print orig_xloc, orig_yloc, orig_zloc, environ_info
            if time_of_day2 % 200 == 0:
                print time_of_day,
                print temp_xloc, temp_yloc, temp_zloc, self.cells[-1].cn, self.cells[-1].n, self.cells[-1].noon, self.cells[-1].max_light_PAR_value
                #print cell_conc_above, current_day, current_hour    
            time_of_day -= (simulation_time_interval / 3600.)
        ######this section will run the liveordie code to determine which cells to remove or reproduce
        addcells = []
        remove_index = 0
        global cell_growth
        for cell in self.cells:
            if cell.cn <= 0:
                self.cells[self.cells.index(cell)] = -9
            elif cell.divide:
                if reverse == 'forward' and self.cell_growth:
                    cell.divide = False
                    tempcell = self.Kbrcell(cell, 1, cell.species_number)
                    if self.mutation_occurs == True:
                        #mutation(tempcell)  #this puts the new cell through the mutation process and leaves the original one alone
                        tempcell.mutate_loci()
                    addcells.append(tempcell)
                    
        self.cells.sort()        
        while self.cells[remove_index-1] == -9:
            remove_index -= 1
        print "Dead cells:", self.cells.count(-9), remove_index
        if remove_index < 0:
            self.cells = self.cells[:remove_index]
        self.cells.extend(addcells)
        print "New cells: ", len(addcells)
        #####################################end liveordie section
        time4 = time.time()
        print "One day took: ", time4-time1
        return current_year, current_day, current_hour, world
        
    def process_one_day_multi_reverse(self, bloom_slice, bloom_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, world_update_hours, origin_offset, nitrate_source, world_latitude, world_longitude, reverse, culture, simulation_time_interval, carrying_capacity): #used to have env_temperature and env_salinity instead of world
        HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]
        time1 = time.time()
        time_of_day = 3.00 #starts at 3am, adds 0.05 (3 min) during each time step until it cycles back to 3am
        time_of_day2 = int(round(time_of_day*100))
        #orig_cell_conc = len(self.cells)# + (x*0.25)
        orig_cell_conc = 0
        cell_conc_above = [0]*abs(max_depth)
        size_x_dimension = len(world)
        size_y_dimension = len(world[0])
        day_night = 1 # 0 == day; 1 == night
        current_hour -= 1 #because the first time step below is going to increment the current hour up to 3
        temp_xloc, temp_yloc, temp_zloc = [-1, -1, -1] 
        for time_step in range(int(1440 / (simulation_time_interval / 60.))):
            #time2 = time.time()        
            if time_of_day < 0.0:
                time_of_day = 23.95
                current_hour = 23
                current_day -= 1
                if current_day < 1:
                    current_year -= 1
                    current_day = 365
            time_of_day2 = int(round(time_of_day*100))
            if time_of_day2 % 100 == 0:
                if time_of_day > 0.0:
                    current_hour += 1
                if current_hour in world_update_hours:
                    current_hr_index = world_update_hours.index(current_hour) - 1
                    if current_hr_index < 0:
                        world = cellmodel_environment_HYCOM_v4_predictive.update_world(world, current_year, current_day-1, world_update_hours[current_hr_index], origin_offset, nitrate_source, reverse, culture)
                    else:
                        world = cellmodel_environment_HYCOM_v4_predictive.update_world(world, current_year, current_day, world_update_hours[current_hr_index], origin_offset, nitrate_source, reverse, culture)    
##            if 5.96 < time_of_day < 18.0:    #check to see what time it is; if it's daytime, have the day_night variable indicate this
##                day_night = 0
##            else:
##                day_night = 1
            day_night = None
            
            current_date = str(current_year) + '/' + self.what_month_is_it(current_year, current_day)
            for cell in bloom_slice:
                if 2.96 < time_of_day < 3.04:
                    cell.liveordie(carrying_capacity, reverse) 
                if cell.cn <= 0:
                    continue
                orig_xloc = cell.x
                orig_yloc = cell.y
                orig_zloc = cell.z
                try:
                    temp_environ_info = cellmodel_environment_HYCOM_v4_predictive.interpolate_location_values_version3(world, orig_xloc, orig_yloc, orig_zloc)
                    environ_info = temp_environ_info[0][:5]
                    if True in numpy.isnan(environ_info):
                        cell.cn = 0
                        #print orig_xloc, orig_yloc, orig_zloc, environ_info
                        continue
                except:
                    cell.cn = 0
                    continue
                environ_info = list(environ_info)
                environ_info.extend([temp_environ_info[1], temp_environ_info[2]])
                #environ_info.append(temp_environ_info[2])
                environ_info.extend(list(temp_environ_info[0][5:]))

                cell.time_step_environment([time_of_day, day_night, cell_conc_above, orig_cell_conc, max_depth, max_light, environ_info], world_latitude, world_longitude, current_date, culture, True)
                temp_xloc = cell.x
                temp_yloc = cell.y
                temp_zloc = cell.z
                
                
                if temp_xloc >= size_x_dimension:
                    #cell.x = len(world)-0.1 
                    cellx_int = int(temp_xloc) - 1
                elif temp_xloc <= 0:
                    #cell.x = 0.01
                    cellx_int = 0                
                else:
                    cellx_int = int(temp_xloc)
                    
                if temp_yloc >= size_y_dimension:
                    #cell.y = len(world[0])-0.1
                    celly_int = int(temp_yloc) - 1
                elif temp_yloc <= 0:
                    #cell.y = 0.01
                    celly_int = 0
                else:
                    celly_int = int(temp_yloc)

                cellz_int = int(cell.z)
                max_depth_at_x_y_location = abs(world[cellx_int, celly_int, 0][4])
                if max_depth_at_x_y_location < 1E7: #not np_isnan(max_depth_at_x_y_location):
                    if temp_zloc < abs(max_depth) and 0 < temp_zloc < max_depth_at_x_y_location:
                            cellz_int = int(cell.z)       
                    elif temp_zloc >= max_depth_at_x_y_location and max_depth_at_x_y_location != 0:
                        cell.z = max_depth_at_x_y_location - 0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc > max_depth_at_x_y_location and max_depth_at_x_y_location == 0:
                        if cell.z < max_depth_at_x_y_location + 4:
                            cellz_int = int(cell.z)
                        else:
                            cell.z = max_depth_at_x_y_location + 4
                            cellz_int = int(cell.z)
                    elif temp_zloc >= abs(max_depth):
                        cell.z = abs(max_depth)-0.01
                        cellz_int = int(cell.z)
                    elif temp_zloc < 0:
                        cell.z = 0.01
                        cellz_int = 0
                    else:
                        cellz_int = int(cell.z)
                
                zdepth = 0
                while cellz_int >= HYCOM_depth_values[zdepth+1]:
                    zdepth += 1
                cellz_int = zdepth

                orig_zdepth = 0
                while orig_zloc >= HYCOM_depth_values[orig_zdepth+1]:
                    orig_zdepth += 1
                
                if not world[cellx_int, celly_int, cellz_int][0] < 1E7: #checking for NaN here 
                    try:
                        if world[cellx_int, celly_int, orig_zdepth][0] < 1E7:
                            cell.z = orig_zloc
                        elif max_depth_at_x_y_location < 1E7:
                            if world[cellx_int][celly_int][HYCOM_depth_values.index(max_depth_at_x_y_location)][0] < 1E7: #checking for NaN here
                                cell.z = max_depth_at_x_y_location
                        elif world[cellx_int][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.y = orig_yloc
                            if cell.z > world[cellx_int][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[cellx_int][int(orig_yloc)][cellz_int][4]
                        elif world[int(orig_xloc)][celly_int][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            if cell.z > world[int(orig_xloc)][celly_int][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][celly_int][cellz_int][4]
                        elif world[int(orig_xloc)][int(orig_yloc)][cellz_int][0] < 1E7:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            if cell.z > world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]:
                                cell.z = world[int(orig_xloc)][int(orig_yloc)][cellz_int][4]
                        else:
                            cell.x = orig_xloc
                            cell.y = orig_yloc
                            cell.z = orig_zloc
                    except:
                        print cellx_int, celly_int, cellz_int, cell.x, cell.y, cell.z, orig_xloc, orig_yloc, orig_zloc
                        
                            
                if cell.x >= size_x_dimension - 1:
                    cell.cn = -9
                elif cell.x <= 0:
                    cell.cn = -9
                elif cell.y >= size_y_dimension - 1:
                    cell.cn = -9
                elif cell.y <= 0:
                    cell.cn = -9                    
                    #world[temp_xloc][temp_yloc][temp_zloc][8]-= 1
                    #world[int(cell.x)][int(cell.y)][int(cell.z)][8]+=1
            #if time_of_day2 % 200 == 0:
                #print time_of_day,
                #print temp_xloc, temp_yloc, temp_zloc,
                #print cell_conc_above, current_day, current_hour    
            time_of_day -= (simulation_time_interval / 3600.)
        ######this section will run the liveordie code to determine which cells to remove or reproduce
        addcells = []
        remove_index = 0
        global cell_growth
        for cell in bloom_slice:
            if cell.cn <= 0:
                bloom_slice[bloom_slice.index(cell)] = -9
            elif cell.divide:
                if reverse == 'forward' and self.cell_growth:
                    cell.divide = False
                    tempcell = self.Kbrcell(cell, 1, cell.species_number)
                    if self.mutation_occurs == True:
                        #mutation(tempcell)  #this puts the new cell through the mutation process and leaves the original one alone
                        tempcell.mutate_loci()
                    addcells.append(tempcell)
                    
        bloom_slice.sort()        
        while bloom_slice[remove_index-1] == -9:
            remove_index -= 1
        print "Dead cells:", bloom_slice.count(-9), remove_index
        if remove_index < 0:
            bloom_slice = bloom_slice[:remove_index]
        bloom_slice.extend(addcells)
        print "New cells: ", len(addcells)
        #####################################end liveordie section
        time4 = time.time()
        print "One day took: ", time4-time1
        return bloom_slice

    
    
    def Kbrcell(self, cell, newcell, species):
        #this function used to be in the main Kbr_Growthmodel program but I moved it here to improve speed when adding/removing cells
        cellname = -1
        growthrate = cell.growthrate
        cell_origin = cell.origin
        divide = False
        e3 = cell.e3
        hc = cell.hc
        mutation_rate_cell_properties = self.mutation_rate_cell_properties
        other_loci = cell.msat_loci[:] #12 extra loci right now; added 9/26/2010
        other_loci_mut_rates = cell.msat_mut_rate[:]
        if random.random() < mutation_rate_cell_properties:
            #this will set all parameters to the parent cell
            cnmax = cell.cnmax
            cnmin = cell.cnmin
            nmax = cell.nmax
            nmin = cell.nmin
            rm = cell.rm
            pl = cell.pl
            pma = cell.pma
            pmb = cell.pmb
            s250 = cell.s250
            d = cell.d
            ic = cell.Ic
            izt_threshold = cell.izt_threshold
            hc_threshold = cell.hc_threshold
            c_divide_threshold = cell.c_divide_threshold
            n_divide_threshold = cell.n_divide_threshold
            correction_factor = cell.correction_factor
            salinity_preferred= cell.salinity_preferred
            temperature_preferred = cell.temperature_preferred
            
            #below here you can uncomment something to have it vary throughout the simulation
            if random.random() < mutation_rate_cell_properties:
                cnmax = round(random.randint(int((cell.cnmax*100)-50), int((cell.cnmax*100)+50))/100., 3)
            else:
                cnmax = cell.cnmax
            if random.random() < mutation_rate_cell_properties:
                cnmin = round(random.randint(int((cell.cnmin*100)-50), int((cell.cnmin*100)+50))/100., 3)
            else:
                cnmin = cell.cnmin
            if random.random() < mutation_rate_cell_properties:
                nmax = round(random.randint(int((cell.nmax*100)-50), int((cell.nmax*100)+50))/100., 3)
            else:
                nmax = cell.nmax
            if random.random() < mutation_rate_cell_properties:
                nmin = round(random.randint(int((cell.nmin*100)-50), int((cell.nmin*100)+50))/100., 3)
            else:
                nmin = cell.nmin
            #if random.random() < mutation_rate_cell_properties:
            #    rm = round(random.randint(int((cell.rm*1000)-50), int((cell.rm*1000)+50))/1000., 3)
            #else:
            #    rm = cell.rm
            #if random.random() < mutation_rate_cell_properties:
            #    pl = round(random.randint(int((cell.pl*1000)-50), int((cell.pl*1000)+50))/1000., 3)
            #else:
            #    pl = cell.pl
            #if random.random() < mutation_rate_cell_properties:
            #    pma = round(random.randint(int((cell.pma*1000)-50), int((cell.pma*1000)+50))/1000., 3)
            #else:
            #    pma = cell.pma
            #if random.random() < mutation_rate_cell_properties:
            #    pmb = round(random.randint(int((cell.pmb*1000)-50), int((cell.pmb*1000)+50))/1000., 3)
            #else:
            #    pmb = cell.pmb
            #if random.random() < mutation_rate_cell_properties:
            #    s250 = round(random.randint(int((cell.s250*1000)-50), int((cell.s250*1000)+50))/1000., 3)
            #else:
            #    s250 = cell.s250
            #if random.random() < mutation_rate_cell_properties:
            #    d = round(random.randint(int((cell.d*1000)-50), int((cell.d*1000)+50))/1000., 3)
            #else:
            #    d = cell.d
            if random.random() < mutation_rate_cell_properties:
                ic = round(random.randint(int((cell.Ic*100)-50), int((cell.Ic*100)+50))/100., 3)
            else:
                ic = cell.Ic
            if random.random() < mutation_rate_cell_properties:
                izt_threshold = round(random.randint(int((cell.izt_threshold*1000)-50), int((cell.izt_threshold*1000)+50))/1000., 3)
            else:
                izt_threshold = cell.izt_threshold
            if random.random() < mutation_rate_cell_properties:
                hc_threshold = round(random.randint(int((cell.hc_threshold*1000)-50), int((cell.hc_threshold*1000)+50))/1000., 3)
            else:
                hc_threshold = cell.hc_threshold
            if random.random() < mutation_rate_cell_properties:
                c_divide_threshold = round(random.randint(int((cell.c_divide_threshold*1000)-50), int((cell.c_divide_threshold*1000)+50))/1000., 3)
            else:
                c_divide_threshold = cell.c_divide_threshold         
            if random.random() < mutation_rate_cell_properties:
                n_divide_threshold = round(random.randint(int((cell.n_divide_threshold*1000)-50), int((cell.n_divide_threshold*1000)+50))/1000., 3)
            else:
                n_divide_threshold = cell.n_divide_threshold

            #environmental characteristics
            if random.random() < mutation_rate_cell_properties:
                temperature_preferred = round(random.randint(int((cell.temperature_preferred*100)-50), int((cell.temperature_preferred*100)+50))/100., 3)
            else:
                temperature_preferred = cell.temperature_preferred
            if random.random() < mutation_rate_cell_properties:
                salinity_preferred = round(random.randint(int((cell.salinity_preferred*100)-50), int((cell.salinity_preferred*100)+50))/100., 3)
            else:
                salinity_preferred = cell.salinity_preferred
            if random.random() < mutation_rate_cell_properties:
                correction_factor = round(random.randint(int((cell.correction_factor*1000)-50), int((cell.correction_factor*1000)+50))/1000., 3)
            else:
                correction_factor = cell.correction_factor
            
        else:
            cnmax = cell.cnmax
            cnmin = cell.cnmin
            nmax = cell.nmax
            nmin = cell.nmin
            rm = cell.rm
            pl = cell.pl
            pma = cell.pma
            pmb = cell.pmb
            s250 = cell.s250
            d = cell.d
            ic = cell.Ic
            izt_threshold = cell.izt_threshold
            hc_threshold = cell.hc_threshold
            c_divide_threshold = cell.c_divide_threshold
            n_divide_threshold = cell.n_divide_threshold
            salinity_preferred = cell.salinity_preferred
            temperature_preferred = cell.temperature_preferred
            correction_factor = cell.correction_factor
            
        if newcell == 1:
            cn = cnmin
            n = nmin
            z = cell.z
        else:
            cn = cell.cn
            n = cell.n
            z = cell.z
            ic = 6

        cell_data = [cellname, cell.cellid, growthrate, cell.locus1allele, cell.locus2allele, cell.locus3allele, cell.locus4allele, cell.locus5allele, cell.locus6allele,
                     divide, cn, cnmax, cnmin, d, e3, hc, ic, n, nmax, nmin, pl, pma, pmb, rm, s250, cell.x, cell.y, z, c_divide_threshold, hc_threshold, n_divide_threshold, izt_threshold,
                     cell.salinity_acclimated, cell.salinity_stress_level, cell.temperature_acclimated, cell.temperature_stress_level, correction_factor, salinity_preferred, temperature_preferred,
                     other_loci, other_loci_mut_rates, cell_origin, self.simulation_time_interval]

        species_dictionary = {0 : Kbrevis_cell(),
                                  1 : Dinophysis_cell(),
                                  2 : Ptexanum_cell(),
                                  3 : Pminimum_cell(),
                                  4 : Asterionellopsis_cell(),
                                  5 : Thalassionema_cell(),

                                  }
        new_Kbr_cell = species_dictionary[species]
        new_Kbr_cell.define_cell(cell_data)
        return new_Kbr_cell

class Kbrevis_cell:
    def define_cell(self, information):
        #everything here is designed to model a Karenia brevis cell
        #other species can be modeled by modifying the values upon initializing the cell (see different species classes at the bottom)
        self.cellid = information[0]
        self.parentcell = information[1]
        self.growthrate = information[2]
        self.locus1allele = information[3]
        self.locus2allele = information[4]
        self.locus3allele = information[5]
        self.locus4allele = information[6]
        self.locus5allele = information[7]
        self.locus6allele = information[8]
        self.divide = False #information[9]
        self.cn = information[10]
        self.cnmax = information[11]
        self.cnmin = information[12]
        self.d = information[13]
        self.e3 = information[14]
        self.hc = information[15]
        self.Ic = information[16] #light compensation threshold; minimum value of PAR needed for net photosynthesis
        self.n = information[17] #current amount of internal nitrogen
        self.nmax = information[18] #maximum possible internal nitrogen
        self.nmin = information[19] #minimum internal nitrogen; below this the cell dies
        self.pl = information[20] 
        self.pma = information[21] # a constant for determining Pm; Liu et al 2001
        self.pmb = information[22] # a constant for determining Pm; Liu et al 2001
        self.rm = information[23] #dark carbon respiration rate; Liu et al 2001
        self.s250 = information[24] #asymptotic swimming speed acclimated to a light intensity of 250umol quanta m-2 s-1; Liu et al 2001
        self.x = information[25] #cell location x
        self.y = information[26] #cell location y
        self.z = information[27] #cell location z

        
        
        #can add below here more initialization parameters (different C/N thresholds for division, different swim rates, different protein production rates)
        self.c_divide_threshold = information[28] #minimum threshold for carbon in order to allow division
        self.hc_threshold = information[29] #threshold for swimming based on photoinhibition
        self.n_divide_threshold = information[30] #minimum threshold for nitrogen in order to allow division
        self.izt_threshold = information[31] #minimum light threshold to control swimming
        self.no3z_threshold = 1.26 #threshold of nitrate to control swimming
        self.msat_loci = information[39] #microsatellite variables
        self.msat_mut_rate = information[40] #mutation rates for the different microsatellites
        self.origin = information[41] #parent cell that this cell originated from
        self.dt = information[42] #the time step in # seconds (dT == change in time); 180seconds was used in Liu et al 2001; now comes from the simulation_time_interval variable in the main program
        self.dt_minutes = self.dt / 60.
        self.dt_hours = self.dt / 3600. #timestep in hours
        self.kn = 0.42 #half-saturation constant; Liu et al 2001
        self.kq = 8.75 #minimum cellular nitrogen quota for protein synthesis to take place; Liu et al 2001
        self.a = 46.0 #constant in formula for Pm in umol quanta m-2 s-2; Liu et al 2001
        self.b = 17.0 #constant in formula for Pm in umol quanta m-2 s-1; Liu et al 2001
        self.pmc = 3.33 #maximum increment of the diel photosynthesis variation; Liu et al 2001
        self.fi = 0 # ??
        self.t = 1 #time in seconds from Liu et al 2001; not sure this is used at all in this program
        self.Q = 0
        self.A = 0
        self.td = 12 #length of daylight hours from Liu et al 2001; new version of the program makes this variable irrelevant as time is calculated based on date/location

        self.culture = False
        self.cnfull = (self.cnmax-self.cnmin)*0.95
        self.cnhigh = (self.cnmax-self.cnmin)*0.6
        #division_rate calculation
        self.days_since_division = 0.01

        #environmental characteristics/stressors
        self.correction_factor = information[36]    #0.75 for SP3; this will likely be different for every culture
        self.salinity_preferred = information[37] #35. for SP3
        self.salinity_acclimated = information[32]
        self.salinity_stress_level = information[33]
        self.temperature_preferred = information[38] #25. for SP3
        self.temperature_acclimated = information[34]
        self.temperature_stress_level = information[35]
        self.temp_cn = 0
        self.cnprotein = 0
        self.daylength = ephem.Observer() #this is the object used to calculate sunrise/sunset for the cell's location 
        self.sunrise = 6.0 #this will be overwritten by the ephem.Observer
        self.sunset = 18.0 #this will be overwritten by the ephem.Observer
        self.noon = 6.0 #this will be calculated from the sunrise and sunset times
        self.midnight = 6.0 #this will be calculated from the sunset/sunrise times and used in the swimming behavior
        self.max_light_PAR_value = 0 #this will come from the PAR values from the MODIS satellites
        self.max_depth_here = 10.
        #multispecies variables
        self.species = 'Kbrevis'
        self.type = 'dinoflagellate'
        self.species_number = 0
        self.chlorophyll_content = 0.0000000000425
        self.myself()
        self.correction_factor = 1.0
        self.same_location = 0 #this will keep track of how many times a cell remains in the exact same location; if more than 2 days then cell dies because it's likely caught at the coast        
        
    def myself(self):
        pass
    def izt_func_integration(self, xy, orig_cell_conc, weight_of_each_cell, chla, ek0):
        #this is an attempt to calculate the integration equation wiht a pure python function instead of my FOR loop
        cell_conc_above = orig_cell_conc*xy*weight_of_each_cell
        cell_conc_above += weight_of_each_cell
        ektotal = ((ek0 + 0.054*((chla*cell_conc_above)**(0.6666666)) + 0.0088*(chla*cell_conc_above))) * 0.1
        return ektotal
        
    def izt_func(self, time_of_day_temp, z, cell_conc_above, orig_cell_conc, I, salinity=None):  #trying to make the PAR intensity function from Lui et al 2001b
        #ek_range = self.ekrange_func(z, cell_conc_above, orig_cell_conc)  #used to have this as another function only called here so I moved here for performance
        ektotal = 0.0
        #print time_of_day_temp,
        ###########old ek_range_function here
        ek0 = 0.1
        if self.max_depth_here < 200.:  #changed this to use the max depth to determine how close it is to the coast; coastal areas will have more sediment and more attenuation of light
            ek0 += 0.05   #used to be 0.25 here; 7/21/11 changed to 0.1; can range from 0.04 - 0.40; 0.1 is what Lui et al 2001 used. this is the light attenuation due to water (and sediments in the water)
        if self.max_depth_here < 50.:  #changed this to use the max depth to determine how close it is to the coast; coastal areas will have more sediment and more attenuation of light
            ek0 += 0.10   #used to be 0.25 here; 7/21/11 changed to 0.1; can range from 0.04 - 0.40; 0.1 is what Lui et al 2001 used. this is the light attenuation due to water (and sediments in the water)
        if self.max_depth_here < 10.:
            ek0 += 0.15 #higher attenuation value for coastal areas; uses this when the max depth is less than 10m
        if orig_cell_conc <= 0:
            orig_cell_conc = 1
        chla = self.chlorophyll_content #0.0000000000425     #units: pg chla cell-1
        weight_of_each_cell = 1000
        
        ####trying to make the ek0 value salinity dependent 
        if salinity:
            ek0 = 0.1 + ((36. / salinity) / 2.)**5
        
#        for xy in xrange(0, int(abs(z)*10), 1):
#            cell_conc_above = orig_cell_conc*xy*weight_of_each_cell
#            cell_conc_above += weight_of_each_cell
#            ektotal += ((ek0 + 0.054*((chla*cell_conc_above)**(0.6666666)) + 0.0088*(chla*cell_conc_above))) * 0.1
        #trying the python integration function here:
        ektotal, ekerror = quad(self.izt_func_integration, 0, int(abs(z)*10), args=(orig_cell_conc, weight_of_each_cell, chla, ek0))
        ###########
        if 0 <= time_of_day_temp < self.noon:
            time_of_day_temp += 24
        ek_range = ektotal
        
        #if time_of_day_temp - (self.sunrise + self.noon) < 0:
        #    time_of_day_temp = self.sunrise + self.noon
        if time_of_day_temp - (self.sunrise) == 0:
            time_of_day_temp = self.sunrise + 0.01
        
        if self.culture:
            izt0 = I
            izt = izt0 * math.exp(-ek_range)
        else:
            izt0 = I*math.sin(((time_of_day_temp-(self.sunrise))*math.pi)/self.td)
            izt = abs(izt0 * math.exp(-ek_range))
        if izt > I:
            izt = I
        elif izt < 0:
            izt = 0
        #print 'izt:', izt, -ek_range
        return izt


        
    def dcnphoto_func(self, I, Ic, Ik, pmd, pl, hc, e3, time_of_day, pma, pmb, pmc, a, b, t, fi, td, rm, Ih, pm):
        dcnphoto = 0
        #Q = self.dQ_func(pmd, hc, pl, I, Ih)
        Q = pmd + hc*(pl-pmd)  #this is the line from the dQ_func; moved here for performance improvement, eliminated a function call
        #A = self.dA_func(e3, Ic, Q, rm)
        A = ((Q*math.tanh(e3))/(rm*e3)) * (Ic*Ic) + Ic - 1 #this is the line from the dA_func; moved here for performance improvement
        if A < 1.0:
            A = 1.0
        if I < Ic:
            dcnphoto_low = rm*(-1 + (I/A) + (1-(Ic/A))*math.exp(((A+1)/Ic)*(I-Ic)))  #units: pmol C;the change in net photosynthesis carbon; only used during daylight when PAR intensity is >= 0 but less than Ic
            dcnphoto = dcnphoto_low
            #print "low", dcnphoto, I, Ic, time_of_day, pmd, Q
        elif I >= Ic:
            dcnphoto_high = Q * math.tanh(e3) * (1.0 - math.exp((-I/Ik) * math.tanh(I-Ic)))   #units: pmol C; the change in net photosynthesis carbon; only used during daylight when PAR intensity (I) is >= the PAR compensation threshold (Ic)    
            dcnphoto = dcnphoto_high 
            #print "high", dcnphoto
            self.Q = Q
            self.A = A
        return dcnphoto

    def pm_func(self, pma, pmb, e3, a, b):
        pm = pma + pmb * math.tanh((e3-a)/b)
        return pm

    def pmd_func(self, pm, pmc, time_of_day, fi, td):
        if self.culture:
            t = 6
        else:
            t = ((time_of_day-self.sunrise)) #used to have -6
        pmd = pm + pmc * (math.sin(((t+fi) * math.pi) / td)**3)
        return pmd


    def ik_func(self, e3):
        ik = e3
        return ik

    def ih_func(self, e3):
        ih = (2 * e3) + 5
        return ih

    def ve3_func(self, e3, I):    #the governing equation for the 3day PAR exposure --> e3
        t3 = (3 * 24 * 60 * 60.) / self.dt
        ve3 = e3 + ((I-e3) / t3)
        if ve3 < 0:
            ve3 = 0.0
        return ve3

    def hi_func(self, I, Ih):
        I = I*1.0
        if I < Ih:
            hi = 0
        elif I >= Ih:
            hi = 1.0-math.exp(-(I-Ih)/Ih)
            self.hi = hi
        return hi

    def vn_func(self, kn, no3_conc, vmax):   #the governing equation for the internal cellular nitrogen
        #no3z = float(self.no3_func(z, max_depth, cell_conc_above)) #old way; now i use the world
        no3z = no3_conc
        vn = ((vmax * (no3z / (kn + no3z)))  / 60.) * self.dt_minutes
    
        return vn
    def vmax_func(self, n):
        vmax = (5.46 * math.exp(-0.186 * n)) #5.46 for Kbrevis
        return vmax
    
    def vhc_func(self, I, Ih, hc):    #the governing equation for the cumulative photoinhibition
        hi = self.hi_func(I, Ih)
        gamma = 60. / self.dt_minutes        #units: h; time scale of induction and recovery of the photoinhibition
        vhc = (hi - hc) / gamma
        return vhc

    def hc_func(self, I, Ih, hc):
        hc_new = hc + self.vhc_func(I, Ih, hc)
        return hc_new
    
    def saccli_func(self, s250, I):     #the light-acclimated swimming speed
        #used to have d here equal to 0.26; in an effort speed things up I have eliminated the variable and hard-coded it
        s250_2 = s250 * 277.77777777#/ 60. / 60. * 1000000
        alpha = 0.55
        #saccli = s250 * (1 + (d * (((math.tanh((alpha*I)/(d*s250_2))) - math.tanh((alpha*250)/(d*s250_2))))))
        saccli = s250 * (1 + (0.26 * (((math.tanh((alpha*I)/(0.26*s250_2))) - math.tanh((alpha*250)/(0.26*s250_2))))))
        return saccli

    def swimming_velocity(self, s250, I, cn, z, n, day_night, time_of_day, izt, hc, cnhigh, kn, cnmax, cnmin, nmax, nmin, max_depth, cell_conc_above, no3_conc, no3_conc_above, no3_conc_below, saccli):    #these are the different swimming velocity functions
        #these functions help determine how far to move a particular cell during a given time slot
        #saccli = self.saccli_func(s250, izt)
        swim_dir, swim_speed = self.swimming_orientation(cn, n, day_night, time_of_day, z, izt, hc, kn, cnmax, cnmin, nmax, nmin, max_depth, cell_conc_above, no3_conc, no3_conc_above, no3_conc_below)
        #no3z = self.no3_func(z, max_depth, cell_conc_above) #old way; now i use the world
        no3z = no3_conc
        vz = 0
        if swim_speed <> 0:
            if swim_speed == 1:
                vz = saccli * (((cnmax-cn)/(cnmax-cnmin))*((cnmax-cn)/(cnmax-cnmin))) #swim speed 1
            elif swim_speed == 2:
                vz = -saccli #swim speed 2
            elif swim_speed == 3:
                vz = saccli * (0.2 + (0.8 * (1.0 - (no3z/(kn+no3z))))) * (((cnmax-cn)/(cnmax-cnhigh))**0.5) #swim speed 3
            elif swim_speed == 4:
                vz = -saccli * (0.2 + (0.8 * (1.0 - (no3z/(kn+no3z))))) * (((cn-cnhigh)/(cnmax-cnhigh))**0.5)#swim speed 4
            elif swim_speed == 5 and swim_dir == 1:
                vz = saccli * ((0.2 + 0.8 * (1- (no3z/(kn+no3z)))) * (0.515 *( 1+ math.tanh(12-(13.33 * ((n-nmin)/(nmax-nmin)))))) - 0.069) #swim speed 5_1
            elif swim_speed == 5 and swim_dir == -1:
                vz = -saccli * ((0.2 + 0.8 * (1- (no3z/(kn+no3z)))) * (0.515 *( 1+ math.tanh(12-(13.33 * ((n-nmin)/(nmax-nmin)))))) - 0.069)# swim speed 5_2
            elif swim_speed == 6:
                vz = 0.25 * saccli * ((cnmax-cn)/(cnmax-cnmin)) #swim speed 6
            
        else:    
        #elif swim_speed == 0:
            vz = 0 #swim speed 0
            
        return vz                                 


    def swimming_orientation(self, cn, nn, day_night, time_of_day, z, izt, hc, kn, cnmax, cnmin, nmax, nmin, max_depth, cell_conc_above, no3_conc, no3_conc_above, no3_conc_below):
        c = (cn-cnmin) * (1./(cnmax-cnmin))         #this scales the value to the range (max/min) where 90 = 1.0 and 36 = 0.0
        swim_orientation = 0
        n = (nn-nmin) * (1./(nmax-nmin))
        #no3_1 = self.no3_func(z, max_depth, cell_conc_above) #<-- old way; now i'm using the world
        no3_1 = no3_conc
        no3 = round(no3_1, 2)
        swim_speed = 0
        sal_stress = self.salinity_stress_level
        temp_stress = self.temperature_stress_level
        izt_thresh = self.izt_threshold
        no3z1 = no3_conc_above
        no3z2 = no3_conc_below
        no3z = ((no3z1-no3_1)+(no3_1-no3z2)) * 0.5 #((no3z1-no3_1)+(no3_1-no3z2))/2.0 #it's computationally faster to multiply by 0.5 than to divide by 2.
        if day_night == 0:                  #going to try incorporating the temp and salinity stress here to see if i can get them to respond to the environment (NOT cross thermoclines) 7/21/11
            
            if hc < self.hc_threshold: #used to be 0.8 here
                if c < 0.6 and sal_stress < 0.4 and temp_stress < 0.4:
                    swim_orientation = 1
                    swim_speed = 1
                elif 0.6 <= c < 0.95 and izt < izt_thresh and sal_stress < 0.2 and temp_stress < 0.2:
                    swim_orientation = 1
                    swim_speed = 1
                elif 0.6 <= c < 0.8 and n >= 0.9 and izt >= izt_thresh and sal_stress < 0.2 and temp_stress < 0.2:
                    swim_orientation = 0
                
                elif 0.6 <= c < 0.8 and n < 0.9 and izt >= izt_thresh and sal_stress < 0.2 and temp_stress < 0.2:
                    if no3 == self.no3z_threshold:
                        swim_orientation = 0
                    elif no3 > self.no3z_threshold:
                        swim_orientation = -1
                        swim_speed = 4
                    elif no3 < self.no3z_threshold and no3z >= 0.01:
                        swim_orientation = 1
                        swim_speed = 3
                    else:
                        swim_orientation = -1
                        swim_speed = 4
                elif 0.8 <= c < 0.95 and izt >= 17.5 and sal_stress < 0.2 and temp_stress < 0.2:
                    if no3 == self.no3z_threshold:
                        swim_orientation = 0
                    elif no3 > self.no3z_threshold:
                        swim_orientation = -1
                        swim_speed = 4
                    elif no3 < self.no3z_threshold and no3z >= 0.01:
                        swim_orientation = 1
                        swim_speed = 3
                    else:
                        swim_orientation = -1
                        swim_speed = 4
                elif c >= 0.95 and sal_stress < 0.1 and temp_stress < 0.1:
                    if no3 == self.no3z_threshold:
                        swim_orientation = 0
                    elif no3 > self.no3z_threshold:
                        swim_orientation = -1
                        swim_speed = 4
                    elif no3 < self.no3z_threshold and no3z >= 0.01:
                        swim_orientation = 1
                        swim_speed = 3
                    else:
                        swim_orientation = -1
                        swim_speed = 4
            else:
            #elif hc >= self.hc_threshold: #used to be 0.8 here
                swim_orientation = -1
                swim_speed = 2
        else:
            
            if self.sunrise - self.midnight <= time_of_day < self.sunrise - 2:
                if no3 >= self.no3z_threshold and sal_stress < 0.2 and temp_stress < 0.2:
                    swim_orientation = 0
                elif no3 < self.no3z_threshold and no3z >= 0.01 and sal_stress < 0.2 and temp_stress < 0.2:
                    swim_orientation = 1
                    swim_speed = 5
                else: 
                    swim_orientation = -1
                    swim_speed = 5
            elif self.sunrise - 2 <= time_of_day < self.sunrise:
                if n >= 0.9 or c < 0.1 and sal_stress < 0.2 and temp_stress < 0.2:
                    swim_orientation = 1
                    swim_speed = 6
                
                else:
                    swim_orientation = 0
            else: #if temp_sunset <= time_of_day <= temp_sunset + 6:
                if no3z >= 0.01 and sal_stress < 0.2 and temp_stress < 0.2:
                    swim_orientation = 1
                    swim_speed = 5
                else:# no3z < 0.01:
                    swim_orientation = -1
                    swim_speed = 5
        return [swim_orientation, swim_speed]

    def salinity_stress(self, env_salinity):
        sal_accl = self.salinity_acclimated
        sal_pref = self.salinity_preferred
        sal_accl = sal_accl + ((env_salinity-sal_accl) / (2 * 1440 / self.dt_minutes))#0.0010416666) #/ 960.) #used to have 480 here for both salinity and temperature 480=one day in the 3min timesteps I use, 960 = 2days
        sal_stress_level = (abs(env_salinity - sal_accl)/sal_accl) - (1-(abs(sal_accl-sal_pref)/sal_pref)) / (2 * 1440 / self.dt_minutes)#* 0.0010416666 #/960.
        
        if 0 <= sal_stress_level <= 1:
            self.salinity_stress_level = sal_stress_level
        elif sal_stress_level > 1.0:
            self.salinity_stress_level = 1.0
        else:
            self.salinity_stress_level = 0.0
        self.salinity_acclimated = sal_accl
        
    
    def temperature_stress(self, env_temperature):
        temp_accl = self.temperature_acclimated
        temp_pref = self.temperature_preferred
        temp_accl = temp_accl + ((env_temperature - temp_accl) / (2 * 1440 / self.dt_minutes))#* 0.0010416666) #/ 960.) multiplying by a decimal is a lot faster than division so that's why the funny number is here
        temp_stress_level = ((abs(env_temperature - temp_accl)/temp_accl) - (1-(abs(temp_accl-temp_pref)/temp_pref)) / (2 * 1440 / self.dt_minutes))#* 0.0010416666) #/960.
        
        if 0 <= temp_stress_level <= 1:
            self.temperature_stress_level = temp_stress_level
        elif temp_stress_level > 1:
            self.temperature_stress_level = 1.0
        else:
            self.temperature_stress_level = 0.0
        self.temperature_acclimated = temp_accl
        
    def swim_x_direction(self, u, world_latitude, world_longitude):
        try:
            longdist_at_current_lat = ((math.pi * 6378137 * math.cos(math.radians(world_latitude[int(self.y)]))) /
                                        (180 * (1-(0.00669437 * (math.sin(math.radians(world_latitude[int(self.y)]))**2)))**.5))
            longdist_at_current_lat *= 0.04
        except:
            print self.x, len(world_latitude)
            longdist_at_current_lat = ((math.pi * 6378137 * math.cos(math.radians(world_latitude[int(self.y)]))) /
                                        (180 * (1-(0.00669437 * (math.sin(math.radians(world_latitude[int(self.y)]))**2)))**.5))
        self.x += ((u + numpy.random.normal(0, 0.01))*180)/longdist_at_current_lat #incoming value is in m/s and my time steps are 3 min (180 sec.); resolution of model data is ~3.5-4km, changes depending on how far north it is
        self.x += (u*self.dt)/longdist_at_current_lat #incoming value is in m/s and my time steps are 3 min (180 sec.); resolution of model data is ~3.5-4km, changes depending on how far north it is
                             
    def swim_y_direction(self, v):
        self.y += ((v + numpy.random.normal(0, 0.01))*180)/4430.  #incoming value is in m/s and my time steps are 3 min; about 4.43km in 1/25 degree
        self.y += (v*self.dt)/4430.  #incoming value is in m/s and my time steps are 3 min; about 4.43km in 1/25 degree
        
    def calculate_sun_values(self):
        self.sunrise = self.daylength.next_rising(ephem.Sun()).datetime()
        self.sunset = self.daylength.next_setting(ephem.Sun()).datetime()
        self.sunrise = self.sunrise.hour + ((self.sunrise.minute / self.dt_minutes) * self.dt_hours)#this converts the datetime object in the 24.00 hour type I'm using to track time used to be: self.sunrise.hour + ((self.sunrise.minute / 3.) * 0.05)
        self.sunset = self.sunset.hour + ((self.sunset.minute / self.dt_minutes) * self.dt_hours) #this converts the datetime object in the 24.00 hour type I'm using to track time, used to be: self.sunset.hour + ((self.sunrise.minute / 3.) * 0.05)
        if self.sunset < self.sunrise:
            self.noon = (24 - (self.sunrise - self.sunset)) / 2. #this will give the hours between sunrise and noon; to get the actual time of noon, add the sunrise time to this value
        else:
            self.noon = (self.sunset - self.sunrise) / 2. #this will give the hours between sunrise and noon; to get the actual time of noon, add the sunrise time to this value
        self.midnight = (24 - (self.noon * 2)) * 0.5 #this will provide an estimate of midnight; it's the number of hours after sunset until the midnight
        self.td = self.noon * 2.        

    def convert_and_calc_PAR(self, max_light_PAR):
        #this will take the PAR value from MODIS (given in E / m**-2 / d**-1) and convert it to uE / m**-2 / s**-1
        result = (max_light_PAR * 1000000) / ((self.noon * 2) * 3600.)
        return result

    def mutate_loci(self):
        for locus in range(len(self.msat_loci)):
            if random.random() < self.msat_mut_rate[locus]:
                up_or_down = random.random()
                if up_or_down < 0.5:
                    if self.msat_loci[locus] > 1:
                        self.msat_loci[locus] -= 1
                else:
                    if self.msat_loci[locus] < 49 :
                        self.msat_loci[locus] += 1

    
    def time_step_environment(self, variables, world_latitude, world_longitude, current_date, culture, reverse_run=False):
        time_of_day = variables[0]
        day_night = variables[1]        
        cell_conc_above = variables[2]
        orig_cell_conc = variables[3]
        max_depth = variables[4]
        max_light = variables[5]
        #may want to add a temperature/salinity variable here at some point to take into account differences under different conditions
        env_temperature = variables[6][0] #used to be just 6
        env_salinity = variables[6][1] #used to be just a 7
        no3_conc = float(variables[6][2])
        co2_conc = float(variables[6][3])
        self.max_depth_here = variables[6][4]
        no3_conc_above = float(variables[6][5])
        no3_conc_below = float(variables[6][6])
        u_vel = variables[6][7]
        v_vel = variables[6][8]
        mixed_layer_depth = variables[6][9]
        max_light_PAR = variables[6][10]
        
        
        self.z *= -1
        if time_of_day < .03:
            self.days_since_division += 1
        if not any(culture):
            if 3.0 < time_of_day <= 3.05:
                self.daylength.lat = str(world_latitude[int(self.y)])
                self.daylength.lon = str(world_longitude[int(self.x)])
                self.daylength.date = current_date + ' 03:03'
                self.calculate_sun_values()
            max_light_PAR = self.convert_and_calc_PAR(max_light_PAR)  #find out the value for max light based on cloud cover and PAR from merged MODIS Aqua and Terra data; gaps will be filled by using the previous days values
            self.max_light_PAR_value = max_light_PAR #for bookkeeping purposes
        else:
            #self.culture = True
            max_light_PAR = culture[0]
            self.max_light_PAR_value = max_light_PAR
            self.sunrise = culture[4][0]
            self.sunset = culture[4][1]
            if self.sunset < self.sunrise:
                self.noon = (24 - (self.sunrise - self.sunset)) * 0.5 #this will give the hours between sunrise and noon; to get the actual time of noon, add the sunrise time to this value
            else:
                self.noon = (self.sunset - self.sunrise) * 0.5 #this will give the hours between sunrise and noon; to get the actual time of noon, add the sunrise time to this value
            self.midnight = (24 - (self.noon * 2)) * 0.5 #this will provide an estimate of midnight; it's the number of hours after sunset until the midnight
            self.td = self.noon * 2.  
        #indiv cell values above here, everything below here is either calculated from the above data or constant
        
        #I = 0    
##        if day_night == 1:      #change the I value according to whether it's day or night
##            I = 0
##        elif day_night == 0:    #full sun is reportedly 1500; we grow ours at 70 and the high light experiment was 140 (this is microeinsteins)
##            I = max_light
        if self.sunrise <= time_of_day < self.sunset or (self.sunrise <= time_of_day and self.sunset < self.sunrise) or (time_of_day < self.sunset and self.sunset < self.sunrise): #this should turn the sun on when the current time is between sunrise and sunset
            if numpy.isnan(max_light_PAR):
                I = 0#max_light
                day_night = 0
            else:
                I = max_light_PAR
                day_night = 0
        else:
            I = 0
            day_night = 1
        #print time_of_day, I, self.sunrise, self.sunset, '         ',
        if I == 0:
            izt = 0
        else:            
            izt = self.izt_func(time_of_day, self.z, cell_conc_above, orig_cell_conc, I, env_salinity)  #use this for normal field conditions
        self.I = izt
        self.e3 = self.ve3_func(self.e3, izt)
        
##        if self.e3 < 0:
##            self.e3 = 0
        self.Ih = self.ih_func(self.e3)
        self.Ik = self.e3
        #self.hi = self.hi_func(izt, self.Ih)
        self.hc = self.hc_func(izt, self.Ih, self.hc)
        self.vmax = self.vmax_func(self.n)
##        try:
##            self.vmax = (5.46 * math.exp(-0.186 * self.n))
##        except:
##            print self.n
##            self.vmax = (5.46 * math.exp(-0.186 * self.n))
        #pm = pma + pmb * math.tanh((e3-a)/b)
        self.pm = self.pm_func(self.pma, self.pmb, self.e3, self.a, self.b)
        #pmd = pm + pmc * math.sin(((t+fi) * math.pi) / td)  #units: pmol C h-1; the dark adapted production rate
        self.pmd = self.pmd_func(self.pm, self.pmc, time_of_day, self.fi, self.td)   #use this for field conditions; sun rising and setting
        
        #environmental stressors
        self.salinity_stress(env_salinity)
        self.temperature_stress(env_temperature)
        
        #internal carbon/nitrogen variables
        if day_night == 0:
            self.cnphoto = self.dcnphoto_func(izt, self.Ic, self.Ik, self.pmd, self.pl, self.hc, self.e3, time_of_day, self.pma, self.pmb, self.pmc, self.a, self.b, self.t, self.fi, self.td, self.rm, self.Ih, self.pm)
        else: #if day_night == 1:
            self.cnphoto = 0
        temp_n = self.vn_func(self.kn, no3_conc, self.vmax)  #returns pmol/3min; no3_conc is in micromolar units
        #this is the original stress function
        temp_n *= ((1 - (0.02*self.salinity_stress_level)) * (1 - (0.02*self.temperature_stress_level))*
                             (math.exp(-((env_salinity-self.salinity_preferred)**2)*0.006666666))*
                             (math.exp(-((env_temperature-self.temperature_preferred)**2)*0.006666666))*
                             self.correction_factor)
        
        #print "temp_n:", temp_n, self.vmax, no3_conc
        if not reverse_run:
            if (no3_conc*(1E6)) >= temp_n:
                self.n += temp_n #self.vn_func(self.z, self.kn, self.n, max_depth, cell_conc_above)
                #no3_conc = ((no3_conc*(1E6))-temp_n)*(1E-6)
            else:
                self.n += (no3_conc * (1E6))
                #no3_conc = 0.0
        else:
            if (no3_conc*(1E6)) >= temp_n:
                self.n -= temp_n #self.vn_func(self.z, self.kn, self.n, max_depth, cell_conc_above)
                #no3_conc = ((no3_conc*(1E6))-temp_n)*(1E-6)
            else:
                self.n -= (no3_conc * (1E6))
                #no3_conc = 0.0
            if self.n < 5:
            #    self.cn = 10
                self.n = 20
        #self.cnprotein = self.dcnprotein_func(self.kq, self.n, self.vmax)  #units returned are pmols*N/cell*h
        dcnprotein = 0.814 * self.vmax * (1 - (self.kq / (0.869*self.n)))
        if dcnprotein < 0:
            dcnprotein = 0
        self.cnprotein = dcnprotein
        #self.cn += self.vcn_func(day_night, self.cnphoto, self.cnprotein, self.rm) *
        #            (1 - self.salinity_stress_level * (abs(self.salinity_acclimated-self.salinity_preferred)/self.salinity_preferred)) *
        #            (1 - self.temperature_stress_level * (abs(self.temperature_acclimated-self.temperature_preferred)/self.temperature_preferred))

        
        if day_night == 0:
            vcn = ((self.cnphoto - dcnprotein) / 60.) * self.dt_minutes #/ 20.  #light
        else:
            vcn = ((-self.rm - dcnprotein) / 60.) * self.dt_minutes #/ 20.       #dark
        
        #self.temp_cn = self.vcn_func(day_night, self.cnphoto, self.cnprotein, self.rm)

        if vcn > 0:
            #self.temp_cn = vcn * ((1 - (0.02*self.salinity_stress_level)) * (1 - (0.02*self.temperature_stress_level))*
            self.temp_cn = vcn * (
                             (math.exp(-((env_salinity-self.salinity_preferred)**2)*0.02))*   #was 0.006666666 here instead of 0.02
                             (math.exp(-((env_temperature-self.temperature_preferred)**2)*0.02))*
                             self.correction_factor)
        else:
            self.temp_cn = vcn
        #self.temp_cn = vcn   
                                
        if not reverse_run:
            self.cn += self.temp_cn
        else:
            if self.temp_cn > 0:
                self.cn -= self.temp_cn
            else:
                self.cn += self.temp_cn
                
        if self.cn > self.cnmax:
            self.cn = self.cnmax
        if self.n > self.nmax:
            self.n = self.nmax
        
        #if np_isnan(self.salinity_stress_level) or np_isnan(self.temperature_stress_level) or np_isnan(self.salinity_preferred) or np_isnan(self.temperature_preferred):
        #    self.cn = 0
        #swimming orientation/speed variables
        
        if self.type == 'dinoflagellate':
            saccli = self.saccli_func(self.s250, izt)
            self.vz = self.swimming_velocity(1, I, self.cn, self.z, self.n, day_night, time_of_day, izt, self.hc, self.cnhigh, self.kn, self.cnmax, self.cnmin, self.nmax, self.nmin, max_depth, cell_conc_above, no3_conc, no3_conc_above, no3_conc_below, saccli)
            #if self.vz < 0:  #this line is used for the 'swim down only' simulation runs
            #    self.vz = 0
            if not reverse_run:
                self.z += (self.vz / 60.) * self.dt_minutes #((self.vz / 60.) *3)
            else:
                self.z += (self.vz / 60.) * self.dt_minutes #((self.vz / 60.) *3)
        elif self.type == 'diatom':
            if not reverse_run:
                if self.z >= -mixed_layer_depth:
                    self.z += numpy.random.normal(0.0001 ,0.001) # this line is to simulate random "swimming" up or down and at varying speeds
                else:
                    self.z -= 0.001 #arbitrarily chose 0.001 value for sinking rate
            else:
                if self.z >= -mixed_layer_depth:
                    self.z -= numpy.random.normal(0.0001 ,0.001) # this line is to simulate random "swimming" up or down and at varying speeds
                else:
                    self.z += 0.001 #arbitrarily chose 0.001 value for sinking rate
        #self.z += w_vel   #this adds the model vertical velocity to the cell's movements; Rob suggested removing it due to it being unreliable
                               
        if self.z >= 0: 
            self.z = -0.01
        elif self.z <= max_depth:
            self.z = max_depth + 0.01
        if not any(culture):
            self.swim_x_direction(u_vel, world_latitude, world_longitude)
            #self.x += (u_vel*180)*0.00025 #/4000. #incoming value is in m/s and my time steps are 3 min (180 sec.); resolution of model data is ~3.5-4km 
            self.swim_y_direction(v_vel)
            #self.y += (v_vel*180)*0.00025 #/4000.  #incoming value is in m/s and my time steps are 3 min; resolution of model data is ~3.5-4km    
        #resources_used = copy.deepcopy([env_temperature, env_salinity, no3_conc, co2_conc, number_of_cells_here])
        resources_used = [env_temperature*1.0, env_salinity*1.0, no3_conc*1.0, co2_conc*1.0, 1]
        
        self.z *= -1
        
        #print self.cn, self.n, vcn, self.cnphoto, dcnprotein, self.e3, izt, time_of_day
        return resources_used

    def liveordie(self, carrying_capacity, reverse_run):
        
        cellliveordie = random.random()
        if self.same_location >= 100:
            self.cn = -1        
        elif cellliveordie < carrying_capacity:
            self.cn = -1
        elif self.cn < self.cnmin *.9:
            self.cn = -1
        elif self.cn > (self.cnmax*self.c_divide_threshold) and self.n > (self.nmax*self.n_divide_threshold):# and cell.divide == 1:# and replicatecell < 1.9999:
            self.divide = True
            self.cn -= self.cnmin
            self.n -= self.nmin
            self.growthrate = 1./self.days_since_division
            self.days_since_division = 0.01
        
        if reverse_run:
            if self.cn < self.cnmin:
                self.cn = self.cnmax
            if self.n < self.nmin:
                self.n = self.nmax
        
        return
##################################################### below here are old designations and program stuff, use for reference
#variable declaration begins below here
#cell_id = cellid        #a unique value to define this cell for later divisions and keeping track of it
#parent_cell = parent    #the cell id of the cell that divided and created this cell
#genotype = [0,0,0,0,0,0]    #the genotype of this cell; it comes from the parent cell after having been passed through the mutation option/function
#growth_rate = 0.4       # maximum growth rate for this cell; not sure if i'm going to give each cell and indep rate or not
#dt = 180.0              #units: s; the time step for each calculation = 3 min
#e3_initial = 2                #units: umol quanta m-2 s-1; an initializing value; variable between [0,522] the sun-shade photoacclimation (3day PAR exposure)
#e3 = 1#e3_initial + ve3_func(e3_initial, I)    #units: umol quanta m-2 s-1; variable between [0,522] the sun-shade photoacclimation (3day PAR exposure)
#ek0 = 0.1               #units: m-1; the PAR attenuation coefficient due to water alone
#t3 = 3                  #units: days; the time scale of the sun-shade acclimated parameters
#Ih = ih_func(e3)
#Ik = e3
#a = 46.0                #units: umol quanta m-2 s-1
#b = 17.0                #units: umol quanta m-2 s-1
#t = 1                   #units: s; time
#fi = 0
#td = 12                 #units: h; length of the daylight hours
#s250 = 1
#d = 0.26
#defining the CARBON characteristics for this cell; also contains some PHOTOSYNTHESIS stuff
#cn = 90.0  #units: pmol C; the internal/cellular carbon for this cell; the value will be set to the lowest value indicating this cell was the least daughter
#cnmax = 90.0 #units: pmol C; the maximum internal cellular carbon
#cnmin = 36.0 #units: pmol C; the minimum internal cellular carbon
#cnfull = 87.3 #units: pmol C; the internal carbon threshold for swimming orientation
#cnhigh = 68.4 #units: pmol C; the internal carbon threshold for swimming orientation
#rm = 3* (0.333/60.)      #units: pmol C; the dark carbon respiration rate; i may need to add another "light carbon respiration rate" for daytime carbon fixation
#cnphoto = 4.0      #units: pmol C; net photosynthetic rate
#cnprotein = 3.0    #units: pmol C; carbon in the cellular protein
#hc = 0.31
#pl = 0.25       #units: pmol C h-1; the light-adapted production rate
#pma = 0.67      #units: pmol C h-1; a constant for determining pm
#pmb = 0.25      #units: pmol C h-1; a constant for determining pm
#pmc = 3.33      #units: pmol C h-1;  max increment of the diel photosynthesis variation
#pm = pma + pmb * math.tanh((e3-a)/b)    #units: pmol C h-1; the sun-shade acclimated maximum photosynthetic rate
#pmd = pm + pmc * math.sin(((t+fi) * math.pi) / td)  #units: pmol C h-1; the dark adapted production rate
#Q = pmd + hc_func(I, Ih, hc)*(pl-pmd)
#A = ((Q*math.tanh(e3))/(rm*e3)) * (Ic**2) + Ic - 1
#imax = 1500     #units: umol quanta m-2 s-1
#izt = izt_func(day_night, time_of_day, z)    
# all CARBON related things should be above here if possible


# defining the NITROGEN characteristics for this cell
#n = 7.0         #units: pmol N; the internal cellular nitrogen
#nmax = 23.3     #units: pmol N; the maximum internal nitrogen
#nmin = 6.32     #units: pmol N; the minimum internal nitrogen
#kn = 0.42       #units: uM NO3-N; half-saturation constant
#kq = 8.75       #units: pmol N; the minimum cellular nitrogen quota for protein synthesis to take place
#no3th = 1.26    #units: uM NO3-N; the ambient nitrate concentration threshold for swimming orientation control
    
    
    #all NITROGEN related things should be above here if possible


class Dinophysis_cell(Kbrevis_cell):
    print "UhOh!"
    def myself(self):
        self.species = 'Dinophysis'
        self.type = 'dinoflagellate'
        self.species_number = 1
        self.pl = 0.798 #(.58 from Harding 1988)
        self.rm = 0.5
        self.Ic = 5.0
        self.cnmax = 92.8
        self.cnmin = 27.8
        self.nmax = 20.17
        self.nmin = 6.05
        self.s250 = 0.6  #m/s
        self.kn = 3.0            
        self.kq = 7.57          
        self.no3z_threshold = 2.0 #Kbrevis has a lower Kn (Ks) value and is more sensitive to NO3; Pmin is less sensitive so the concentration to control swimming should be higher right?
                            
        
class Ptexanum_cell(Kbrevis_cell):
    #should run the define cell function first and then the following lines will change the values of the important parameters
    def myself(self):
        self.species = 'Ptexanum'
        self.type = 'dinoflagellate'
        self.species_number = 2
        self.pl = 0.798  #0.798
        self.rm = 0.42
        self.Ic = 5.0
        self.cnmax = 92.8
        self.cnmin = 27.8 #(~30% of the cnmax; this was just an estimate made by me)
        self.nmax = 20.17 #(based on the 4.6 C:N ratio given by Falkowski et al 1985)
        self.nmin = 6.05 #(~30% of the nmax; just an estimate made by me)
        self.s250 = 0.50 #m/s; maximum swim speed at light intensity of 250um
        self.kn = 0.42
        self.kq = 7.57
        self.no3z_threshold = 1.0
        self.chlorophyll_content = 0.0000000000425
    
class Pminimum_cell(Kbrevis_cell):

    def myself(self):
        self.species = 'Pminimum'
        self.type = 'dinoflagellate'
        self.species_number = 3
        self.pl = 0.083 #(.58 from Harding 1988)
        self.rm = 0.1
        self.Ic = 3.5
        self.cnmax = 12.5
        self.cnmin = 3.75
        self.nmax = 6.52
        self.nmin = 1.95
        self.s250 = 0.38  #m/s
        self.kn = 5.0            
        self.kq = 2.44          
        self.no3z_threshold = 3.0 #Kbrevis has a lower Kn (Ks) value and is more sensitive to NO3; Pmin is less sensitive so the concentration to control swimming should be higher right?
                
                
      
class Asterionellopsis_cell(Kbrevis_cell):
    def myself(self):
        self.species = 'Asterionellopsis'
        self.type = 'diatom'
        self.species_number = 4
        self.pl = .7   #see Chan 1980 - values were ~0.48 for Cylindrotheca fusiformis and ~0.39 for Thalassiosira eccentrica (these values are pmol C/(pg chlA * h)-->this assumes one pg chlA
        self.rm = 0.05  #according to Chan 1980 the 'dark uptake rate' was ~1-1.5% of the max photosyn rate 
        self.Ic = 2.5
        self.pma = 2.5 #arbitrary choice to make Q approach a max value of ~10
        self.pmb = 2.5 #arbitrary choice to make Q approach a max value of ~10
        self.pmc = 5.0 #arbitrary choice to make Q approach a max value of ~10
        
        self.cnmax = 73.0 #an approx 7.3:1 C:N ratio per Sarthou et al 2005; C.fusiformis had 32.5 pgC/cell and T. eccentrica had 1750pgC/cell from Chan 1980
        self.cnmin = 15.0
        self.nmax = 10.0 # 
        self.nmin = 2.0
        self.s250 = 0.0  #m/s
        self.chlorophyll_content = 49.7E-12 #pg chlA in a cell per Chan 1980 for Thalassiosira eccentrica; it's only 0.98 for Cylindrotheca fusiformis; and 42.5 for Kbrevis (Liu et al 2001)
        self.kn = 1.6   #1.6 per Sarthou et al 2005     
        self.kq = 2.44          
        self.no3z_threshold = 1.0 #Kbrevis has a lower Kn (Ks) value and is more sensitive to NO3; Pmin is less sensitive so the concentration to control swimming should be higher right?
        self.correction_factor = 1.0
    def vmax_func(self, n):
        vmax = (7.0 * math.exp(-0.186 * n)) #5.46 for Kbrevis
        return vmax
         
class Thalassionema_cell(Kbrevis_cell):
    def myself(self):
        self.species = 'Thalassionema'
        self.type = 'diatom'
        self.species_number = 5
        self.pl = .7   #see Chan 1980 - values were ~0.48 for Cylindrotheca fusiformis and ~0.39 for Thalassiosira eccentrica (these values are pmol C/(pg chlA * h)-->this assumes one pg chlA
        self.rm = 0.05  #according to Chan 1980 the 'dark uptake rate' was ~1-1.5% of the max photosyn rate 
        self.Ic = 2.5
        self.pma = 2.5 #arbitrary choice to make Q approach a max value of ~10
        self.pmb = 2.5 #arbitrary choice to make Q approach a max value of ~10
        self.pmc = 5.0 #arbitrary choice to make Q approach a max value of ~10
        
        self.cnmax = 73.0 #an approx 7.3:1 C:N ratio per Sarthou et al 2005; C.fusiformis had 32.5 pgC/cell and T. eccentrica had 1750pgC/cell from Chan 1980
        self.cnmin = 15.0
        self.nmax = 10.0 # 
        self.nmin = 2.0
        self.s250 = 0.0  #m/s
        self.chlorophyll_content = 49.7E-12 #pg chlA in a cell per Chan 1980 for Thalassiosira eccentrica; it's only 0.98 for Cylindrotheca fusiformis; and 42.5 for Kbrevis (Liu et al 2001)
        self.kn = 1.6   #1.6 per Sarthou et al 2005     
        self.kq = 2.44          
        self.no3z_threshold = 1.0 #Kbrevis has a lower Kn (Ks) value and is more sensitive to NO3; Pmin is less sensitive so the concentration to control swimming should be higher right?
        self.correction_factor = 1.0
    def vmax_func(self, n):
        vmax = (7.0 * math.exp(-0.186 * n)) #5.46 for Kbrevis
        return vmax    
    
