import random
import pylab
from pylab import contourf, contour
import cellmodel_fieldbloom_envir_HYCOM_v4_multispp as cellmodel_bloom_cell_class
import cellmodel_environment_HYCOM_v4_predictive as cellmodel_environment
import time
import copy
import numpy
import pp
import IBM_sample_information_NOAA

#6/15/09 - added a function that will make the first split equal so i can compare two different populations (1 vs 2)
#        - also made it 20 blooms total and added a size of each bloom output for each day and a popgen output for every 5th day

#6/23/09 - changing the different source setting (first - so that there is one), and now once the first two blooms split they will undergo
#        - a changeable number of mutation cycles in order to cause the allele freq distrib to be different in both blooms
#8/3/09 - added the genepop function to have it make genepop files from the model data as well; i know create popgen, arlequin, genepop
#       - add a bottleneck function in which a large portion of a population is removed? how would this be biologically relevant?

#2/11/2013 - converting to HYCOM model for use with HYCOM model data for the environmental data

#6/3/09 - adding a mating/crossing over function to make it more 'realistic'

#added changes to make the world use real data from the HYCOM Gulf of Mexico dataset Feb. 2013

#going to change the cell so that it starts with mostly full carbon stores in order to allow it to "burn-in" without dying immediately, also eliminated the live_or_die step for the first ten days

#4/29/2014 - changing the program to incorporate multiple species at once, each species having a different chlorophyll content and swimming speed

#09/09/2016 - adding the functionality to use forecast currents and forecast where cells will be going

mating_occurs = False
freq_of_mating = .001   #this is the frequency of mating events
mutation_occurs = False
mutation_rates = [0.001, 0.001, 0.0001, 0.0001, 0.00001, 0.00001]
migration = False
migration_rate = 0.001
skew_freqs = False                    
skewlevel = 4
max_pop_size = 100000000

mutation_rate_cell_properties = 0.05
predictive = False 
future_date = [2016, 252, 3]  #this will be used for predicting where cells will go based on forecast currents from HYCOM [year, julian day, hour]

start_hour = 3     #this is the hour that the world data will load first
start_day = 182    #this is the day that the world data will load first
start_year = 2011   #this is the year that the world data will load first

current_hour = 3    #this is the hour the simulation will start on; should probably always be a 3
current_day = 182   #this is the day the simulation will start on
current_year = 2011 #this is the year the simulation will start on

environment_update_interval = [0,3,6,9,12,15,18,21] #make this a list with the hours from each day that I want to use to update the environment
nitrate_source = 'uniform' #can be uniform or NODC
run_forward_or_backward = 'forward'  #can be forward or backward
cell_death = False
cell_growth = True
track_cell_locations = True
cell_loc_frequency = 1
add_real_sample_data = False
max_num_cells_per_sample = 500
use_multiple_species = False
species_to_model = ['Karenia']  #Kbrevis, Ptexanum, Pminimum, Dinophysis
init_azim = -60
init_elev = 30


length_of_run = 200
simulation_time_interval = 180 #number of seconds in a time step for each day, probably best to keep it in multiples of 60; originally set to 180
#give the x,y,z start location , how many cells at this location, and the concentration at this location;
#multiple starts can be given, just give multiple lists
start_locations_and_concentrations = []
#for x in range(150, 250, 5):
#    for y in range(275, 310, 5):
#        for z in range(10, 60, 5):
#            start_locations_and_concentrations.append([x, y, z, 1, 2])


environ_width = 300  #this plus offset[0] cannot be greater than 450
environ_length = 350    #this plus offset[1] cannot be greater than 350
environ_origin_offset = (0,0) # (0,200) puts it starting at Brownsville
max_depth = -100 #has to be either 1, 6, 11, 16, 21, 26, 31, 41, 51, 61, 71, 81, 91, 101, 126, 151, 201  based on depths given by HYCOM model

HYCOM_depth_values = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200]

identical_cells = False
###########################culture variables
simulate_a_culture = False #making this true will create a world with one grid point, 0m deep
max_light = 160  #this is the light power being given in the simulation 1500=sunlight 70=our incubator
env_salinity = 20.0 #what salinity should the culture environment be
env_temperature = 15.0 #what temperature should the culture environment be
env_nutrients = 1.00 #what nutrient concentration (Nitrate) should be in the culture medium
num_cells_to_start = 1 #how many culture cells to start with
lights_on_off = (6.00, 18.00)  #what hour to turn the lights on/off (ie sunrise/sunset)

if simulate_a_culture:
    culture = [max_light, env_salinity, env_temperature, env_nutrients, lights_on_off]
    for x in range(num_cells_to_start):
        start_locations_and_concentrations.append([0, 0, 1, 2])
else:
    culture = [None]

###########################end culture variables

##########################
#this is code for the parallel python library to use
use_threads = True    #do you want to use multiple processors?
num_processors_to_use = 4  #how many processors do you want to use?
if use_threads:
    ppservers = ()
    job_server = pp.Server(num_processors_to_use,  ppservers=ppservers)#,  secret = 'pop') #this should limit the processing to only 4 cores/threads
##########################

######temp warnings
import warnings
warnings.filterwarnings(action='ignore', message= '.*less.*', category=RuntimeWarning)
warnings.filterwarnings(action='ignore', message= '.*greater.*', category=RuntimeWarning)
######

#trying a new way where the locations are given as lat/lon with concentration and the function in the world library converts it to integers
extra_bloom_dates = IBM_sample_information_NOAA.extra_bloom_dates
cell_id_number = 1

######variables for plotting cells and environmental info
cell_concs_time = [0]
temp_time = [0]
sal_time = [0]



print "creating world...."

world, world_latitude, world_longitude = cellmodel_environment.create_world_version3(environ_width, environ_length, max_depth, start_year, start_day, start_hour,
                                                                                     environ_origin_offset, nitrate_source, run_forward_or_backward, culture)

print "World lengths: ", len(world), len(world[0]), len(world[0][0])
print 'done!'


write_log_files = True  #if false then no files will be written
if write_log_files == True:
    bloommodelpath = "D:/CJunk/Kbr_Model_data/temp/"          #windows path
    #bloommodelpath = 'C:/Documents and Settings/Lisa/My Documents/Darren_IFCB6/Kbrevis_model/outfiles/' #campbell8 path    
    #bloommodelpath = '/host/CJunk/Kbr_Model_data/temp/'    #linux path
    #bloommodelpath = '/home/darren/ibm_model/temp/'  #campbel13 path
    
    bloomsizepath = bloommodelpath + "bloom_sizes.txt"
    bloom_size_file = open(bloomsizepath, 'w')
    #bloom_cell_loc_str = bloommodelpath + 'cellmodel.csv'
    #bloom_cell_location = open(bloom_cell_loc_str, 'w')
    cell_information_str = bloommodelpath + "cell_info.txt"
    cell_info_file = open(cell_information_str, 'w')
    cell_location_file = open(bloommodelpath + 'cell_locations.txt', 'w')
    cell_location_file.write('Cell locations\n')
    cell_location_file.close() #clears the file at the beginning of every run; I append to this file later in the run
    model_run_info = open(bloommodelpath + 'run_info.txt', 'w')
    
def Kbrcell(cell, newcell, species):
    global cell_id_number, mutation_rate_cell_properties
    cellname = cell_id_number
    growthrate = cell.growthrate
    cell_origin = cell.origin
    divide = 1
    e3 = cell.e3
    hc = cell.hc
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
        #z = float(random.randint(max_depth+1, -1))
        z = cell.z
        ic = 6#float(random.randint(5,17))
        #c_divide_threshold = random.randint(0,100)/100.
        #n_divide_threshold = random.randint(0,100)/100.
    cell_data = [cellname, cell.cellid, growthrate, cell.locus1allele, cell.locus2allele, cell.locus3allele, cell.locus4allele, cell.locus5allele, cell.locus6allele,
                 divide, cn, cnmax, cnmin, d, e3, hc, ic, n, nmax, nmin, pl, pma, pmb, rm, s250, cell.x, cell.y, z, c_divide_threshold, hc_threshold, n_divide_threshold, izt_threshold,
                 cell.salinity_acclimated, cell.salinity_stress_level, cell.temperature_acclimated, cell.temperature_stress_level, correction_factor, salinity_preferred, temperature_preferred,
                 other_loci, other_loci_mut_rates, cell_origin, simulation_time_interval]
    #cell= [ 0     ,     1     ,      2      ,       3,           4      ,      5      ,       6      ,     7      ,   8  ,
    #        9,   10,   11  , 12,  13, 14, 15, 16, 17,  18,19,   20 , 21 , 22, 23  ,24, 25,26,27,        28,               29,           30,              31,
    #       32, 33, 34, 35, 36, 37, 38,
    #       39, 40, 41, 42] 
    species_dictionary = {0 : cellmodel_bloom_cell_class.Kbrevis_cell(),
                              1 : cellmodel_bloom_cell_class.Dinophysis_cell(),
                              2 : cellmodel_bloom_cell_class.Ptexanum_cell(),
                              3 : cellmodel_bloom_cell_class.Pminimum_cell(),
                              4 : cellmodel_bloom_cell_class.Asterionellopsis_cell(),
                              5 : cellmodel_bloom_cell_class.Thalassionema_cell(),

                              }
    new_Kbr_cell = species_dictionary[species]
    #new_Kbr_cell = cellmodel_bloom_cell_class.Kbrevis_cell()
    
    new_Kbr_cell.define_cell(cell_data)
    if cellname % 100 == 0:
        cell_info_file.write(str(cell_data))
        cell_info_file.write("\n")
    cell_id_number += 1
    return new_Kbr_cell
def addonesubtractone():
    randomvariable = random.random()
    addorsubtract = 0
    if 1 > randomvariable > .5:
        addorsubtract = 1 #this means to add one repeat
        #if it stays 0 then it means to subtract one repeat
    return addorsubtract

def indiv_locus_mutation(locusmutation_rate, originalcell):
    returncell = originalcell   
    addorsubtract = addonesubtractone()
    if addorsubtract == 0:
        if originalcell > 0:
            returncell = originalcell - 1
    elif addorsubtract == 1:
        if originalcell < 50:
            returncell = originalcell + 1
    elif addorsubtract == 2:
        if originalcell < 48:
            returncell = originalcell + 3
    elif addorsubtract == 3:
        if originalcell > 4:
            returncell = originalcell - 3
    return returncell

def mutation(cell):
    temp_mutation_count = 0
    mutatelocus1 = random.random() / mutation_rates[0]
    if mutatelocus1 <= 1:
        cell.locus1allele = indiv_locus_mutation(mutatelocus1, cell.locus1allele)
        temp_mutation_count += 1
    mutatelocus2 = random.random() / mutation_rates[1]
    if mutatelocus2 <= 1:
        cell.locus2allele = indiv_locus_mutation(mutatelocus2, cell.locus2allele)
        temp_mutation_count += 1
    mutatelocus3 = random.random() / mutation_rates[2]
    if mutatelocus3 <= 1:
        cell.locus3allele = indiv_locus_mutation(mutatelocus3, cell.locus3allele)
        temp_mutation_count += 1
    mutatelocus4 = random.random() / mutation_rates[3]
    if mutatelocus4 <= 1:
        cell.locus4allele = indiv_locus_mutation(mutatelocus4, cell.locus4allele)
        temp_mutation_count += 1
    mutatelocus5 = random.random() / mutation_rates[4]
    if mutatelocus5 <= 1:
        cell.locus5allele = indiv_locus_mutation(mutatelocus5, cell.locus5allele)
        temp_mutation_count += 1
    mutatelocus6 = random.random() / mutation_rates[5]
    if mutatelocus6 <= 1:
        cell.locus6allele = indiv_locus_mutation(mutatelocus6, cell.locus6allele)
        temp_mutation_count += 1
    

def newgrowthrate(bloom):
    #cell.growthrate = round(random.randint(15, 40)/100., 5)
    for cell in bloom.cells:
        if cell.days_since_division > 0.05:
            cell.growthrate = 1./cell.days_since_division
        else:
            cell.growthrate = 0.000


def mutation_firsttime(cell):
    mutatelocus1 = random.random() / (0.01)
    if mutatelocus1 <= 1:
        cell.locus1allele = indiv_locus_mutation(mutatelocus1, cell.locus1allele) 
    mutatelocus2 = random.random() / (0.01)
    if mutatelocus2 <= 1:
        cell.locus2allele = indiv_locus_mutation(mutatelocus2, cell.locus2allele)
    mutatelocus3 = random.random() / (0.001)
    if mutatelocus3 <= 1:
        cell.locus3allele = indiv_locus_mutation(mutatelocus3, cell.locus3allele)
    mutatelocus4 = random.random() / (0.001)
    if mutatelocus4 <= 1:
        cell.locus4allele = indiv_locus_mutation(mutatelocus4, cell.locus4allele)
    mutatelocus5 = random.random() / (0.001)
    if mutatelocus5 <= 1:
        cell.locus5allele = indiv_locus_mutation(mutatelocus5, cell.locus5allele)
    mutatelocus6 = random.random() / (0.025)
    if mutatelocus6 <= 1:
        cell.locus6allele = indiv_locus_mutation(mutatelocus6, cell.locus6allele)
    cell.mutate_loci()
#Kbrcell(parentcell, growthrateinit, locus1, locus2, locus3, locus4, locus5, locus6, cn, cnmax, cnmin, d, e3, hc, ic, n, nmax, nmin, pl,
#            pma, pmb, rm, s250, x, y, z, newcell, izt_threshold, hc_threshold, c_divide_threshold, n_divide_threshold):
        
def createfirstbloom(start_locs_concs):
    global world
    bloomtuple = []
    other_loci = [25,25,25,25,25,25,25,25,25,25,25,25] #12 extra loci right now; added 9/26/2010
    other_loci_mut_rates = [0.001,0.001,0.001,0.0001,0.0001,0.0001,0.0005,0.0005,0.0005,0.002,0.002,0.002]
    species_dictionary = {0 : cellmodel_bloom_cell_class.Kbrevis_cell(),
                          1 : cellmodel_bloom_cell_class.Dinophysis_cell(),
                          2 : cellmodel_bloom_cell_class.Ptexanum_cell(),
                          3 : cellmodel_bloom_cell_class.Pminimum_cell(),
                          4 : cellmodel_bloom_cell_class.Asterionellopsis_cell(),
                          5 : cellmodel_bloom_cell_class.Thalassionema_cell(),
                          }
    for bloom_start_info in start_locs_concs:
        world_starting_location, temp_nitrate_above, temp_nitrate_below = cellmodel_environment.interpolate_location_values_version3(world,bloom_start_info[0],bloom_start_info[1],bloom_start_info[2])
        if world_starting_location[0] < 1E7:
            #print bloom_start_info, world_starting_location
            if not identical_cells:
                ######################
                #
                #this one makes each cell unique to begin with
                for cell in range(bloom_start_info[3]):
                    master_cell = species_dictionary[bloom_start_info[4]]
                    cell = [0, 0, .3, 30, 30, 30, 30, 30, 30, 0, float(random.randint(30, 50)), 90, 36, 0.26, 500, .1, 20., float(random.randint(10, 23)), 23.3,
                            6.32, 0.25, 0.67, 0.25, 0.333, 1., float(bloom_start_info[0]), float(bloom_start_info[1]), float(bloom_start_info[2]), 0.9, 0.8, 0.95, 17.5,
                            world_starting_location[1], 0.0, world_starting_location[0], 0.0, random.random(), world_starting_location[1],
                            world_starting_location[0], other_loci, other_loci_mut_rates, (bloom_start_info[0], bloom_start_info[1], bloom_start_info[2]), simulation_time_interval]
             
                    master_cell.define_cell(cell)
                    bloomtuple.append(Kbrcell(master_cell, 0, bloom_start_info[4]))
                    bloomtuple[-1].days_since_division = 1
                #######################
            elif not any(culture):
                ######################identical cells but not a culture
                #this one makes all cells identical to each other at the start
                master_cell = species_dictionary[bloom_start_info[4]]
                cell = [0, 0, .3, 30, 30, 30, 30, 30, 30, 0, 66, 90, 36, 0.26, 1, .1, 6., 16.32, 23.3,
                        6.32, 0.25, 0.67, 0.25, 0.333, 1., float(bloom_start_info[0]), float(bloom_start_info[1]), float(bloom_start_info[2]), 0.9, 0.8, 0.75, 17.5,
                        env_salinity, 0.0, env_temperature, 0.0, random.random(), random.randint(25, 40), random.randint(15, 35),
                        other_loci, other_loci_mut_rates, (bloom_start_info[0], bloom_start_info[1], bloom_start_info[2]), simulation_time_interval]
                master_cell.define_cell(cell)
                for cell in range(bloom_start_info[3]):  #use this when you want all cells identical
                    bloomtuple.append(Kbrcell(master_cell, 0, bloom_start_info[4]))
                    
                    bloomtuple[-1].days_since_division = 1
                #####################
            else:
                ######################identical cells AND it's a culture
                #this one makes all cells identical to each other at the start
                master_cell = species_dictionary[bloom_start_info[4]]
                cell = [0, 0, .3, 30, 30, 30, 30, 30, 30, 0, 66, 90, 36, 0.26, 1, .1, 6., 16.32, 23.3,
                        6.32, 0.25, 0.67, 0.25, 0.333, 1., float(bloom_start_info[0]), float(bloom_start_info[1]), float(bloom_start_info[2]), 0.9, 0.8, 0.75, 17.5,
                        env_salinity, 0.0, env_temperature, 0.0, 0.250, 35, 25,
                        other_loci, other_loci_mut_rates, (bloom_start_info[0], bloom_start_info[1], bloom_start_info[2]), simulation_time_interval]
                master_cell.define_cell(cell)
                for cell in range(bloom_start_info[3]):  #use this when you want all cells identical
                    bloomtuple.append(Kbrcell(master_cell, 0, bloom_start_info[4]))
                    
                    bloomtuple[-1].days_since_division = 1
                #####################

                    
    print "Randomizing first bloom's allele distribution"
    firstbloomx = 0
    while firstbloomx < 50:
        for cell in bloomtuple:
            mutation_firsttime(cell)
            
        firstbloomx += 1
        if firstbloomx % 500 == 0:
            print firstbloomx, 
    print "Done!"
    return bloomtuple


def add_cells_to_bloom(inbloom, start_locs_concs): 
    global world, HYCOM_depth_values, species_to_model
    species_list = {0:'Kbrevis',
                    1: 'Dinophysis',
                    2: 'Ptexanum',
                    3: 'Pminimum',
                    4: 'Asterionellopsis',
                    5: 'Thalassionema',
                }
    bloomtuple = inbloom.cells
    temp_start_locs = []
    for bloom_start in start_locs_concs:
        #print bloom_start
        if species_list[bloom_start[4]] in species_to_model:
            temp_indexes = cellmodel_environment.convert_lat_lon([bloom_start[0]+0.02, bloom_start[1] - 0.02])
            print temp_indexes, bloom_start
            if -1 in temp_indexes or 0 in temp_indexes:
                continue #-1 means the location is outside of the current map. this line just tells the program to skip it
            zdepth = 0
            while bloom_start[2] > HYCOM_depth_values[zdepth+1]:
                zdepth += 1
            for x in range(temp_indexes[0], temp_indexes[0]+1):
                for y in range(temp_indexes[1], temp_indexes[1]+1):
                    world_starting_location, temp_nitrate_above, temp_nitrate_below = cellmodel_environment.interpolate_location_values_version3(world,x,y,bloom_start[2])
                    if world_starting_location[0] < 1E7:
                        if bloom_start[2] > world_starting_location[4]:
                            temp_start_locs.append([x,y,world_starting_location[4], bloom_start[3], bloom_start[4]])
                        else:
                            temp_start_locs.append([x,y,bloom_start[2], bloom_start[3], bloom_start[4]])
                    else:
                        temp_x, temp_y = adjust_world_lat_lon_location([x,y])
                        world_starting_location, temp_nitrate_above, temp_nitrate_below = cellmodel_environment.interpolate_location_values_version3(world,temp_x,temp_y,bloom_start[2])
                        if world_starting_location[0] < 1E7:
                            temp_start_locs.append([temp_x,temp_y,bloom_start[2], bloom_start[3], bloom_start[4]])
        else:
            print "Wrong species:", bloom_start
    print "Extra cells:", len(temp_start_locs), len(bloomtuple)
    print temp_start_locs
    print
    other_loci = [25,25,25,25,25,25,25,25,25,25,25,25] #12 extra loci right now; added 9/26/2010
    other_loci_mut_rates = [0.001,0.001,0.001,0.0001,0.0001,0.0001,0.0005,0.0005,0.0005,0.002,0.002,0.002]
    species_dictionary = {0 : cellmodel_bloom_cell_class.Kbrevis_cell(),
                              1 : cellmodel_bloom_cell_class.Dinophysis_cell(),
                              2 : cellmodel_bloom_cell_class.Ptexanum_cell(),
                              3 : cellmodel_bloom_cell_class.Pminimum_cell(),
                              4 : cellmodel_bloom_cell_class.Asterionellopsis_cell(),
                              5 : cellmodel_bloom_cell_class.Thalassionema_cell(),
                              }
    for bloom_start_info in temp_start_locs:
        print bloom_start_info
        world_starting_location, temp_nitrate_above, temp_nitrate_below = cellmodel_environment.interpolate_location_values_version3(world,bloom_start_info[0],bloom_start_info[1],bloom_start_info[2])
        if world_starting_location[0] < 1E7:
            if identical_cells == False:
                ######################
                #
                #this one makes each cell unique to begin with
                for cell in range(min(bloom_start_info[3], max_num_cells_per_sample)):
                    master_cell = species_dictionary[bloom_start_info[4]]
                    cell = [0, 0, .3, 30, 30, 30, 30, 30, 30, 0, 90, random.randint(60,90), 36, 0.26, 1, .1, 6., float(random.randint(10,20)), 23.3,
                            6.32, 0.25, 0.67, 0.25, 0.333, 1., float(bloom_start_info[0]), float(bloom_start_info[1]), float(bloom_start_info[2]), 0.9, 0.8, 0.95, 17.5,
                            env_salinity, 0.0, env_temperature, 0.0, random.random(), random.randint(25, 40),
                            random.randint(15, 35), other_loci, other_loci_mut_rates, (bloom_start_info[0], bloom_start_info[1], bloom_start_info[2]), simulation_time_interval]
             
                    master_cell.define_cell(cell)
                    bloomtuple.append(Kbrcell(master_cell, 0, bloom_start_info[4]))
                    bloomtuple[-1].days_since_division = 1
                #######################
            else:
                ######################
                #this one makes all cells identical to each other at the start
                master_cell = species_dictionary[bloom_start_info[4]]
                cell = [0, 0, .3, 30, 30, 30, 30, 30, 30, 0, 66, 90, 36, 0.26, 1, .1, 6., 16.32, 23.3,
                        6.32, 0.25, 0.67, 0.25, 0.333, 1., float(bloom_start_info[0]), float(bloom_start_info[1]), float(bloom_start_info[2]), 0.9, 0.8, 0.75, 17.5,
                        env_salinity, 0.0, env_temperature, 0.0, random.random(), random.randint(25, 40), random.randint(15, 35),
                        other_loci, other_loci_mut_rates, (bloom_start_info[0], bloom_start_info[1], bloom_start_info[2]), simulation_time_interval]
                master_cell.define_cell(cell)
                for cell in range(min(bloom_start_info[3], max_num_cells_per_sample)):  #use this when you want all cells identical
                    bloomtuple.append(Kbrcell(master_cell, 0, bloom_start_info[4]))
                    
                    bloomtuple[-1].days_since_division = 1
                #####################

    print "length after:", len(bloomtuple)
    return bloomtuple


def adjust_world_lat_lon_location(inloc):
    #try to find a spot near the original location in which place the cells
    #try one tick to the east/west depending on which side of the gulf it is
    if inloc[0] < -90:  #this set is for the locations near Texas, it will check to the east first, then south, then north, but not west
        if not numpy.isnan(world[inloc[0]+1, inloc[1], 0][0]):
            return [inloc[0]+1, inloc[1]]
        elif not numpy.isnan(world[inloc[0], inloc[1]-1, 0][0]):
            return [inloc[0], inloc[1]-1]
        elif not numpy.isnan(world[inloc[0], inloc[1]+1, 0][0]):
            return [inloc[0], inloc[1]+1]
        else:
            return inloc
    else:   #this set is for the Florida locations, it will check west first, then south, then north, but not east
        if not numpy.isnan(world[inloc[0]-1, inloc[1], 0][0]):
            return [inloc[0]-1, inloc[1]]
        elif not numpy.isnan(world[inloc[0], inloc[1]-1, 0][0]):
            return [inloc[0], inloc[1]-1]
        elif not numpy.isnan(world[inloc[0], inloc[1]+1, 0][0]):
            return [inloc[0], inloc[1]+1]
        else:
            return inloc
    

        
def mating_in_bloom(sourcebloom):  # i need to rewrite this function to incorporate the fact that I'm using a cell class instead of a tuple
    temp_exchanged_alleles_count = 0
    cells_to_mate = random.sample(sourcebloom.cells, 2)
    for locus in range(2,8):
        if random.random() < .0: #temporarily changed to 0 while needing to be rewritten
            tempcell1 = cells_to_mate[0][locus]
            tempcell2 = cells_to_mate[1][locus]
            cells_to_mate[0][locus] = tempcell2
            cells_to_mate[1][locus] = tempcell1
            temp_exchanged_alleles_count += 1
    if temp_exchanged_alleles_count > 3:
        newgrowthrate(cells_to_mate[0])
        newgrowthrate(cells_to_mate[1])


def write_bloom_sizes(field_blooms):
    string_to_write = "Day " + str(x) +":\n"
    bloom_size_file.write(string_to_write)
    xtemp = 1
    for bloom in field_blooms:
        string = str(xtemp) + ": " + str(len(bloom.cells))+ "\n"
        bloom_size_file.write(string)
        xtemp += 1


def save_avg_data():
    f = open(bloommodelpath + "junk_data.pck", 'w')
    import pickle
    d = [growth_trend, cnmax_trend, nmax_trend, c_divide_trend, n_divide_trend, ic_avg_trend,
	     cnmin_trend, nmin_trend, rm_trend, pl_trend, pma_trend, pmb_trend, s250_trend, d_trend, izt_threshold_trend,
	     hc_threshold_trend, cn_trend, n_trend, salinity_trend, temperature_trend, depth_trend]
    pickle.dump(d, f)
    f.close()
    f = open(bloommodelpath + "junk_data.csv", 'w')
    f.write('Growth rate, Cnmax, Nmax, C divide, N divide, Ic, Cnmin, Nmin, Rm, Pl, Pma, Pmb, s250, D, Izt, Hc, Cn, N, Salinity, Temperature, Depth, \n')
    for datapoint in range(len(d[0])):
        for datatype in d:
            f.write(str(datatype[datapoint]))
            f.write(',')
        f.write('\n')
    f.close()
    
def avg_growth_rate(bloom):
    #y = field_blooms.index(bloom)
    #print concen_array[y]
    locus1_counts = []
    gwthrates = []
    cnmaxrates = []
    nmaxrates = []
    c_dividerates = []
    n_dividerates = []
    ic_avg = []
    cnmin_temp = []
    nmin_temp = []
    rm_temp = []
    pl_temp = []
    pma_temp = []
    pmb_temp = []
    s250_temp = []
    d_temp = []
    izt_threshold_temp = []
    hc_threshold_temp = []
    cn_temp = []
    n_temp = []
    sal_stress = []
    temp_stress = []
    sal_prefer = []
    temp_prefer = []
    depth_temp = []
    for cell in bloom.cells:
        
        if cell.locus1allele not in locus1_counts:
            locus1_counts.append(cell.locus1allele)
        gwthrates.append(cell.growthrate)
        cnmaxrates.append(cell.cnmax)
        nmaxrates.append(cell.nmax)
        c_dividerates.append(cell.c_divide_threshold)
        n_dividerates.append(cell.n_divide_threshold)
        ic_avg.append(cell.Ic)
        cnmin_temp.append(cell.cnmin)
        nmin_temp.append(cell.nmin)
        rm_temp.append(cell.rm)
        pl_temp.append(cell.pl)
        pma_temp.append(cell.pma)
        pmb_temp.append(cell.pmb)
        s250_temp.append(cell.s250)
        d_temp.append(cell.d)
        izt_threshold_temp.append(cell.izt_threshold)
        hc_threshold_temp.append(cell.hc_threshold)
        cn_temp.append(cell.cn)
        n_temp.append(cell.n)
        sal_prefer.append(cell.salinity_preferred)
        temp_prefer.append(cell.temperature_preferred)
        sal_stress.append(cell.salinity_stress_level)
        temp_stress.append(cell.temperature_stress_level)
        depth_temp.append(-cell.z)
    if generation % 1 ==0:
        
        locus1_trend.append(locus1_counts)
        #print locus1_trend
        growth_trend.append(numpy.mean(gwthrates))
        cnmax_trend.append(numpy.mean(cnmaxrates))
        nmax_trend.append(numpy.mean(nmaxrates))
        c_divide_trend.append(numpy.mean(c_dividerates))
        n_divide_trend.append(numpy.mean(n_dividerates))
        ic_avg_trend.append(numpy.mean(ic_avg))
        cnmin_trend.append(numpy.mean(cnmin_temp))
        nmin_trend.append(numpy.mean(nmin_temp))
        rm_trend.append(numpy.mean(rm_temp))
        pl_trend.append(numpy.mean(pl_temp))
        pma_trend.append(numpy.mean(pma_temp))
        pmb_trend.append(numpy.mean(pmb_temp))
        s250_trend.append(numpy.mean(s250_temp))
        d_trend.append(numpy.mean(d_temp))
        izt_threshold_trend.append(numpy.mean(izt_threshold_temp))
        hc_threshold_trend.append(numpy.mean(hc_threshold_temp))
        cn_trend.append(numpy.mean(cn_temp))
        n_trend.append(numpy.mean(n_temp))
        salinity_trend.append(numpy.mean(sal_stress))
        temperature_trend.append(numpy.mean(temp_stress))
        depth_trend.append(numpy.mean(depth_temp))
    
    print "Avg growth rate:        ", round(numpy.mean(gwthrates), 4)
    print "Avg cnmax rate:         ", round(numpy.mean(cnmaxrates), 4)
    print "Avg cn:                 ", round(numpy.mean(cn_temp), 4)
    print "Avg n:                  ", round(numpy.mean(n_temp), 4)
    print "Avg sal_prefer:         ", round(numpy.mean(sal_prefer), 4)
    print "Avg temp_prefer:        ", round(numpy.mean(temp_prefer), 4)
    print "Avg ic_rate:            ", round(numpy.mean(ic_avg),4)
    print "Avg stress (sal):       ", round(numpy.mean(sal_stress), 4)
    print "Avg stress (temp):      ", round(numpy.mean(temp_stress), 4)
    print "Avg depth :             ", round(numpy.mean(depth_temp), 4)
    print round(bloom.cells[0].salinity_acclimated, 3), round(bloom.cells[0].temperature_acclimated, 3), round(bloom.cells[0].temp_cn, 5), round(bloom.cells[0].cn, 3), round(bloom.cells[0].n, 3),
    print round(bloom.cells[0].cnprotein, 3)

    model_run_info.write("Avg growth rate:     " + str(round(numpy.mean(gwthrates), 4)) + '\n')
    model_run_info.write("Avg cnmax rate:      " + str(round(numpy.mean(cnmaxrates), 4)) + '\n')
    model_run_info.write("Avg cn:              " + str(round(numpy.mean(cn_temp), 4)) + '\n')
    model_run_info.write("Avg n:               " + str(round(numpy.mean(n_temp), 4)) + '\n')
    model_run_info.write("Avg sal_prefer:      " + str(round(numpy.mean(sal_prefer), 4)) + '\n')
    model_run_info.write("Avg temp_prefer:     " + str(round(numpy.mean(temp_prefer), 4)) + '\n')
    model_run_info.write("Avg ic_rate:         " + str(round(numpy.mean(ic_avg),4)) + '\n')
    model_run_info.write("Avg stress (sal):    " + str(round(numpy.mean(sal_stress), 4)) + '\n')
    model_run_info.write("Avg stress (temp):   " + str(round(numpy.mean(temp_stress), 4)) + '\n')
    model_run_info.write("Avg depth :          " + str(round(numpy.mean(depth_temp), 4)) + '\n')


def process_bloom_slice(bloom_slice, environ_data):
    for cell in bloom_slice:
        cell.time_step_environment(environ_data)
    return bloom_slice

def skew_allele_freqs(bloom):
    global skewlevel
    amt_skew = skewlevel
    
    for cell in bloom.cells:
        if cell.locus1allele < 49-amt_skew:
            cell.locus1allele += amt_skew
        if cell.locus2allele < 49-amt_skew:
            cell.locus2allele += amt_skew
        if cell.locus3allele < 49-amt_skew:
            cell.locus3allele += amt_skew
        if cell.locus4allele < 49-amt_skew:
            cell.locus4allele += amt_skew
        #if cell.locus5allele < 49-amt_skew:
        #    cell.locus5allele += amt_skew
        #if cell.locus6allele < 49-amt_skew:
        #    cell.locus6allele += amt_skew
        
def modify_concen_array(concen):
    #num_blooms = len(concen_array)
    new_concen_array = []
    for bloom in range(1):
        tempbloom = []
        for x in range(len(concen[0][0])):
            tempbloom.append([])
        for y in range(len(concen[0])):
            for x in range(len(concen[0][0])):
                tempbloom[x].append(concen[0][y][x])
        new_concen_array.append(tempbloom)
        
    return new_concen_array


def percent_concen_array(concen):
    #this function scales the concentrations so that later dates with more cells don't swamp out the earlier dates
    #num_blooms = len(concen_array)
    new_concen_array = []
    for bloom in range(1):
        tempbloom = []
        for x in range(len(concen[0][0])):
            tempbloom.append([])
        for y in range(len(concen[0])):
            total = sum(concen[0][y])
            for x in range(len(concen[0][0])):
                temp_pct = int((concen[0][y][x]/float(total))* 100)
                tempbloom[x].append(temp_pct)
        new_concen_array.append(tempbloom)
        
    return new_concen_array



def start_first_cells_in_world(firstbloom):
    global world
    for cell in firstbloom:
        world[int(cell.x)][int(cell.y)][int(cell.z)][8] += 1

def make_figure():
    import pylab
    from pylab import plot as plt
    
    a = plt(growth_trend)
    pylab.savefig('c:/cjunk/kbr_model_data/growth_trend.png', dpi = 100)
    del a
    b = plt(ic_avg_trend)
    pylab.savefig('c:/cjunk/kbr_model_data/ic_avg_trend.png', dpi = 100)
    del b
    c = plt(salinity_trend)
    pylab.savefig('c:/cjunk/kbr_model_data/salinity_trend.png', dpi = 100)
    del c
    d = plt(temperature_trend)
    pylab.savefig('c:/cjunk/kbr_model_data/temperature_trend.png', dpi = 100)
    del d
    e = modify_concen_array(concen_array)
    f = pylab.contourf(e[0])
    pylab.savefig('c:/cjunk/kbr_model_data/cell_concentrations.png', dpi=100)
    del f


    
def make_figure2(cell_concs, temper, salin, nitrate, ssh, scatter_locs, PAR_img, mixed_layer_depth, blank_mask):
    
    global current_year, current_day, current_hour, output_figure, max_depth, bloommodelpath, HYCOM_depth_values
    namestr = str(current_year)
    if current_day < 100:
        namestr += '0'
    if current_day < 10:
        namestr += '0'
    namestr += str(current_day)
    if current_hour < 10:
        namestr += '0'
    namestr += str(current_hour)
    #for depthrange in range(0, 17, 1):
    #    output_figure = contourf(cell_concs[0+depthrange:1+depthrange].sum(0), range(6))
    #    pylab.colorbar()
    #    pylab.title('Cell concentration ' + namestr + ' ' + str(HYCOM_depth_values[depthrange]) + 'm', size=20)
    #    pylab.savefig(bloommodelpath + 'cell_concentration' + namestr + '_' + str(HYCOM_depth_values[depthrange]) + 'm.png', dpi = 250)
    #    pylab.clf()
     
    output_figure = contourf(cell_concs.sum(0), range(11))
    pylab.title(''.join(['Cell concentration ', namestr, ' ', 'All depths']), size=20)        
    pylab.colorbar()   
    pylab.savefig(''.join([bloommodelpath, 'cell_concentration', namestr, '_all.png']), dpi = 150)
    pylab.clf()

    output_figure = contourf(blank_mask, range(6), alpha=.1)
    pylab.title(''.join(['Cell scatter ', namestr, ' ', 'All depths']), size=20)
    pylab.xlim(0, environ_width+1)
    pylab.ylim(0, environ_length+1)
    temp_data = split_locations_by_species(scatter_locs)
    for species_plot in temp_data:
        if len(temp_data[species_plot][0]) > 0:
            pylab.scatter(temp_data[species_plot][0], temp_data[species_plot][1], c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
    #pylab.scatter(scatter_locs[0], scatter_locs[1], c=scatter_locs[2], s=2, alpha=0.9, edgecolor = '', cmap='autumn', vmin=0.0, vmax=75.0)
    #pylab.colorbar()
    pylab.savefig(''.join([bloommodelpath, 'cell_scatter', namestr, '_all.png']), dpi = 150)
    pylab.clf()
    
#    #trying 3d scatter plot here
#    #print 'Making 3d plot',
#    output_figure = contourf(blank_mask, range(6), alpha=.1)
#    output_figure = pylab.figure()
#    ax = output_figure.add_subplot(111, projection='3d')
#    pylab.title(''.join(['Cell scatter ', namestr, ' ', 'All depths']), size=20)
#    pylab.xlim(0, environ_width+1)
#    pylab.ylim(0, environ_length+1)
#    temp_data = split_locations_by_species(scatter_locs)
#    for species_plot in temp_data:
#        if len(temp_data[species_plot][0]) > 0:
#            depth_locs = numpy.array(temp_data[species_plot][2])*-1
#            #pylab.scatter(temp_data[species_plot][0], temp_data[species_plot][1], c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
#            ax.scatter(temp_data[species_plot][0], temp_data[species_plot][1], depth_locs, c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
#    ax.set_zlim(max_depth, 0)
#    #ax.azim = -46 #change the view/azimuth
#    #ax.elev = -163 #change the view/elevation
#    
#    #pylab.colorbar()
#    pylab.savefig(''.join([bloommodelpath, 'cell_3d_scatter', namestr, '_all.png']), dpi = 150)
#    pylab.clf()
#    #print 'Done!'
#    ###end 3d scatter
    
    
    
    output_figure = contourf(salin[0], numpy.arange(10, 41, 0.25))
    pylab.colorbar(ticks=range(10, 51, 5))
    pylab.title(''.join(['Salinity ', namestr, ' 0m']), size=20)
    pylab.savefig(''.join([bloommodelpath, 'salinity', namestr, '.png']), dpi = 150)
    pylab.clf()
    
    output_figure = contourf(temper[0], numpy.arange(5, 35.25, 0.25))
    pylab.colorbar(ticks=range(5, 41, 5))
    pylab.title(''.join(['Temperature ', namestr, ' 0m']), size=20)
    pylab.savefig(''.join([bloommodelpath, 'temperature', namestr, '.png']), dpi = 150)
    pylab.clf()

    output_figure = contourf(mixed_layer_depth[0], numpy.arange(0, 100, 1))
    pylab.colorbar(ticks=range(0, 101, 5))
    pylab.title(''.join(['Mixed Layer Depth ', namestr, ' ']), size=20)
    pylab.savefig(''.join([bloommodelpath, 'mld', namestr, '.png']), dpi = 150)
    pylab.clf()
    
    output_figure = contourf(PAR_img[0], range(70))
    pylab.colorbar(ticks=range(0, 71, 5))
    pylab.title(''.join(['PAR ', namestr, ' 0m']), size=20)
    pylab.savefig(''.join([bloommodelpath, 'PAR', namestr, '.png']), dpi = 150)
    pylab.clf()
    
    if not ssh == -1:
        output_figure = contourf(cell_concs.sum(0), range(6))
        pylab.twinx()
        pylab.twiny()
        output_figure = contour(ssh, numpy.arange(-45, 75, 5), cmap = pylab.get_cmap('PRGn'))
        #pylab.colorbar(ticks=range(-45, 75, 10))
        pylab.title('SSH+Cell conc ' + namestr + ' 0m\n', size=20)
        pylab.savefig(bloommodelpath + 'SSH_cellconc' + namestr +'.png', dpi = 150)
        pylab.clf()

        output_figure = contourf(ssh, numpy.arange(-45, 75, .1), cmap = pylab.get_cmap('PRGn'))
        pylab.colorbar(ticks=range(-45, 76, 10))
        pylab.title('SSH ' + namestr + '0m', size=20)
        pylab.savefig(bloommodelpath + 'SSH ' + namestr +'.png', dpi = 150)
        pylab.clf()

    for depthrange in range(5):
        output_figure = contourf(nitrate[depthrange], numpy.arange(0.0, 10.05, 0.05))
        pylab.colorbar()
        pylab.title(''.join(['Nitrate concen ', namestr, ' ', str(HYCOM_depth_values[depthrange]), 'm']), size=20)
        pylab.savefig(''.join([bloommodelpath, 'nitrate', namestr, '_', str(HYCOM_depth_values[depthrange]), 'm.png']), dpi = 200)
        pylab.clf()
    pylab.close('all')
 
def make_figure3(cell_concs, temper, salin, nitrate, ssh, scatter_locs, PAR_img, mixed_layer_depth, blank_mask):
    #this will make a single figure with all the plots on it
    global current_year, current_day, current_hour, output_figure, max_depth, bloommodelpath, HYCOM_depth_values, init_azim, init_elev
    namestr = str(current_year)
    if current_day < 100:
        namestr += '0'
    if current_day < 10:
        namestr += '0'
    namestr += str(current_day)
    if current_hour < 10:
        namestr += '0'
    namestr += str(current_hour)
    overall_figure = pylab.figure(dpi=600) 
    #pylab.title(''.join(['Model day ', namestr]))
#    output_figure = overall_figure.add_subplot(234)
#    pylab.contourf(cell_concs.sum(0), range(51))
#    pylab.title(''.join(['Concentration']), size=10)        
#    #pylab.colorbar(orientation='horizontal')
#    pylab.xticks(size=5)
#    pylab.yticks(size=5)
#    
    output_figure = overall_figure.add_subplot(221)
    pylab.contourf(blank_mask, range(6), alpha=.1)
    pylab.title(''.join(['Scatter ', namestr]), size=10)
    pylab.xlim(0, environ_width+1)
    pylab.ylim(0, environ_length+1)
    temp_data = split_locations_by_species(scatter_locs)
    for species_plot in temp_data:
        if len(temp_data[species_plot][0]) > 0:
            pylab.scatter(temp_data[species_plot][0], temp_data[species_plot][1], c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
    pylab.xticks(size=5)
    pylab.yticks(size=5)
    
    ax = overall_figure.add_subplot(224, projection='3d')
    depth_0 = get_world_mask(0)
    pylab.contourf(depth_0[0], depth_0[1], depth_0[2], colors='b', alpha=.1)
    depth_25 = get_world_mask(25)
    pylab.contourf(depth_25[0], depth_25[1], depth_25[2], colors='b', alpha=.1)
    depth_50 = get_world_mask(50)
    pylab.contourf(depth_50[0], depth_50[1], depth_50[2], colors='b', alpha=.1)
    if current_day % 1 == 0:
        init_azim += 1.
        init_elev += 0.15
    ax.azim = init_azim #change the view/azimuth
    ax.elev = init_elev #change the view/elevation
    pylab.title(''.join(['Scatter 3D']), size=10)
    pylab.xlim(0, environ_width+1)
    pylab.ylim(0, environ_length+1)
    temp_data = split_locations_by_species(scatter_locs)
    for species_plot in temp_data:
        depth_locs = numpy.array(temp_data[species_plot][2])*-1
        #print 'depth_locs', depth_locs, temp_data['Ptexanum'][2], temp_data[species_plot][2]
        if len(temp_data[species_plot][0]) > 0:
            #pylab.scatter(temp_data[species_plot][0], temp_data[species_plot][1], c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
            ax.scatter(temp_data[species_plot][0], temp_data[species_plot][1], depth_locs, c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
    ax.set_zlim(max_depth, 0)
    ax.set_zticklabels(range(max_depth, 1, 10), size=5)
    pylab.xticks(size=5)
    pylab.yticks(size=5)
    
    ###adding in a scatter of the xloc by depth and yloc by depth
        
    output_figure = overall_figure.add_subplot(223)
    pylab.xlim(0, environ_width+1)
    pylab.ylim(max_depth, 0)
    #temp_data = split_locations_by_species(scatter_locs)
    for species_plot in temp_data:
        if len(temp_data[species_plot][0]) > 0:
            pylab.scatter(temp_data[species_plot][0], depth_locs, c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
    pylab.xticks(size=5)
    pylab.yticks(size=5)
    
    output_figure = overall_figure.add_subplot(222)
    pylab.xlim(0, (max_depth*-1) + 1)
    pylab.ylim(0, environ_length+1)
    #temp_data = split_locations_by_species(scatter_locs)
    for species_plot in temp_data:
        depth_locs = numpy.array(temp_data[species_plot][2])*-1
        depth_locs_pos = depth_locs * -1
        if len(temp_data[species_plot][0]) > 0:
            pylab.scatter(depth_locs_pos, temp_data[species_plot][1], c=temp_data[species_plot][3], s=20, alpha=0.4, edgecolor = '')
    #pylab.xticks(range(0, (max_depth*-1)+1, 10), [str(xtic) for xtic in range(0, -51, -10)], size=5)
    pylab.xticks(size=5)
    pylab.yticks(size=5)
    
#    output_figure = overall_figure.add_subplot(235)
#    pylab.contourf(salin[0], numpy.arange(20, 41, 0.25))
#    #pylab.colorbar(ticks=range(20, 41, 5), orientation='horizontal')
#    pylab.title(''.join(['Salinity']), size=10)
#    pylab.xticks(size=5)
#    pylab.yticks(size=5)
#    
#    output_figure = overall_figure.add_subplot(236)
#    pylab.contourf(temper[0], numpy.arange(5, 40.25, 0.25))
#    #pylab.colorbar(ticks=range(5, 41, 5), orientation='horizontal')
#    pylab.title(''.join(['Temperature']), size=10)
#    pylab.xticks(size=5)
#    pylab.yticks(size=5)
    
#    output_figure = overall_figure.add_subplot(259)
#    pylab.contourf(mixed_layer_depth[0], numpy.arange(0, 100, 1))
#    #pylab.colorbar(ticks=range(0, 101, 5), orientation='horizontal')
#    pylab.title(''.join(['Mixed Layer Depth']), size=10)
#    pylab.xticks(size=5)
#    pylab.yticks(size=5)
#    
#    output_figure = overall_figure.add_subplot(2,5,10)
#    pylab.contourf(PAR_img[0], range(70))
#    #pylab.colorbar(ticks=range(0, 71, 5), orientation='horizontal')
#    pylab.title(''.join(['PAR']), size=10)
#    pylab.xticks(size=5)
#    pylab.yticks(size=5)
    
    
    pylab.savefig(''.join([bloommodelpath, 'Master', namestr, '.png']), dpi = 300)
    pylab.clf()
    pylab.close('all')
    
def get_world_mask(in_depth):
    ###this will provide a mask for plotting a contour map with the land masked at certain depths
    temp_world = world.T
    blank_depth = numpy.zeros((len(temp_world[0][0]), len(temp_world[0][0][0])))
    blank_depth[blank_depth == 0] = in_depth * -1
    blank_mask = numpy.ma.masked_array(blank_depth, numpy.isnan(temp_world[0][HYCOM_depth_values.index(in_depth)]))
    X, Y = numpy.meshgrid(range(len(temp_world[0][0][0])), range(len(temp_world[0][0])))
    return X, Y, blank_mask
    
    
def split_locations_by_species(in_loc_data):
    #print 'in_loc_data', in_loc_data
    #this function will split the location list into separate lists based on species
    loclists = { 'Kbrevis':    [[],[],[], 'r'],
                 'Dinophysis': [[],[],[], 'y'],
                 'Ptexanum':   [[],[],[], 'k'],
                 'Pminimum':   [[],[],[], 'c'],
                 'Asterionellopsis': [[],[],[], 'm'],
                 'Thalassionema': [[],[],[], 'g'],
                 }
    remove_these_spp = []
    for spp in range(len(in_loc_data[3])):
        loclists[in_loc_data[3][spp]][0].append(in_loc_data[0][spp])
        loclists[in_loc_data[3][spp]][1].append(in_loc_data[1][spp])
        loclists[in_loc_data[3][spp]][2].append(in_loc_data[2][spp])
    for spp in loclists.keys():
        if len(loclists[spp][0]) == 0:
            remove_these_spp.append(spp)
    for spp in remove_these_spp:
        loclists.pop(spp)
    #print 'returned', loclists
    return loclists

def poll_cell_locations(blooms, world):
    #in order to use multiprocessing I'm going to have to move my cell tracking to an external function (here)
    #and only update cell locations once per day
    global max_depth, HYCOM_depth_values, cell_id_number
    size_x_dimension = len(world)
    size_y_dimension = len(world[0])
    size_z_dimension = len(world[0][0])
    
    zeros = numpy.zeros((size_z_dimension, size_y_dimension, size_x_dimension), dtype = numpy.int)
    xlocs = []
    ylocs = []
    zlocs = []
    origins = []
    cellids = []
    species_list = []
    for bloom in blooms:
        for cell in bloom.cells:
            if cell.cellid == -1:
                cell.cellid = cell_id_number
                cell_id_number += 1
            temp_xloc = round(cell.x)
            temp_yloc = round(cell.y)
            temp_zloc = round(cell.z)
            
            if temp_xloc >= size_x_dimension:
                temp_xloc = int(cell.x) - 1
            if temp_yloc >= size_y_dimension:
                temp_yloc = int(cell.y) - 1
            if temp_zloc >= size_z_dimension:
                temp_zloc = int(cell.z) - 1
            zdepth = 0
            while temp_zloc > HYCOM_depth_values[zdepth+1]:
                zdepth += 1
            zeros[int(zdepth), int(temp_yloc), int(temp_xloc)] += 1
            xlocs.append(temp_xloc)
            ylocs.append(temp_yloc)
            zlocs.append(temp_zloc)
            origins.append(cell.origin)
            cellids.append(cell.cellid)
            species_list.append(cell.species)
    if track_cell_locations:
        if generation % cell_loc_frequency == 0:
            track_cell_locs(cellids, origins, xlocs, ylocs, zlocs, species_list)
        
    #world = world.T
    #old_locations = world[8][:]
    #world[8] = zeros
    #world = world.T
    old_locations = zeros
    return old_locations, [xlocs, ylocs, zlocs, species_list]

def world_map_maker(world, cell_locations, scatter_locs):
    global temp_time, sal_time, current_day, current_year
    temp_world = world.T
    blank_mask = numpy.ma.masked_array(numpy.zeros((len(temp_world[0][0]), len(temp_world[0][0][0]))), numpy.isnan(temp_world[0][0]))
    conc_mask = numpy.ma.masked_array(cell_locations, numpy.isnan(temp_world[0]))
    temp_mask = numpy.ma.masked_array(temp_world[0], numpy.isnan(temp_world[1]))
    salin_mask = numpy.ma.masked_array(temp_world[1], numpy.isnan(temp_world[0]))
    nitrate_mask = numpy.ma.masked_array(temp_world[2], numpy.isnan(temp_world[0]))
    PAR_mask = numpy.array(temp_world[8]) #numpy.ma.masked_array(temp_world[8], numpy.isnan(temp_world[1]))
    mixed_layer_depth = numpy.ma.masked_array(temp_world[7], numpy.isnan(temp_world[0]))
    #cell_concs_time[0] = conc_mask.sum(0)
    temp_time[0] = temp_mask[0]
    sal_time[0] = salin_mask[0]
    if bloommodelpath == "C:/CJunk/Kbr_Model_data/temp/":
        day_str = ''
        if current_day < 100:
            day_str += '0'
        if current_day < 10:
            day_str += '0'
        day_str += str(current_day)
        #ssh = netCDF4.Dataset('C:/CJunk/Satellite_imagery/FTP Downloads/SSH/ssha_gom_wrt_model_'+str(current_year)+day_str+'.nc').variables['ssh']
        ssh = -1
    else:
        ssh = -1
    make_figure2(conc_mask, temp_time, sal_time, nitrate_mask, ssh, scatter_locs, PAR_mask, mixed_layer_depth, blank_mask)
    make_figure3(conc_mask, temp_time, sal_time, nitrate_mask, ssh, scatter_locs, PAR_mask, mixed_layer_depth, blank_mask)


def track_cell_locs(cellids, origins, xlocs, ylocs, zlocs, species_list):
    cell_location_file = open(bloommodelpath + 'cell_locations.txt', 'r') #read the file in to get the cell ids and lines
    cell_dictionary = {}
    header_line = cell_location_file.next()
    for line in cell_location_file:
        if len(line) > 1:
            cellname, celllocs = line[:-1].split(',', 1)
            cell_dictionary[int(cellname)] = celllocs

    cell_location_file.close()
    cell_location_file = open(bloommodelpath + 'cell_locations.txt', 'w') #rewrite the file with new data added
    cell_location_file.write(header_line[:-1] + ',' + str(generation) + '\n')
    for cell_index, indiv_cell in enumerate(cellids):
        if indiv_cell in cell_dictionary: #this will add the  latest location data to the cells already present in the file
            cell_location_file.write(''.join([str(indiv_cell), ',', cell_dictionary.pop(indiv_cell), ',', str([xlocs[cell_index], ylocs[cell_index], zlocs[cell_index]]), ',', species_list[cell_index], '\n']))

        else: #this will add new cells to the file and skip the data points that they missed
            cell_location_file.write(''.join([str(indiv_cell), ',']))
            cell_location_file.write(''.join([',' for numgensskipped in range(generation / cell_loc_frequency)])) #this line replaces the two below it; should be faster 
            #for numgensskipped in range(generation / cell_loc_frequency):
            #    cell_location_file.write(',')
            cell_location_file.write(''.join([str([xlocs[cell_index], ylocs[cell_index], zlocs[cell_index]]), '\n']))
            
    #for indiv_cell in cell_dictionary: #this should write all the dead cells at the end of the file
    #    cell_location_file.write(str(indiv_cell) + ',' + cell_dictionary[indiv_cell] + '\n')
    dead_cell_file = open(bloommodelpath + 'cell_locations_dead_cells.txt', 'a')  #append the dead cells to this file to avoid having to rewrite them every time
    [dead_cell_file.write(''.join([str(indiv_cell), ',', cell_dictionary[indiv_cell], '\n'])) for indiv_cell in cell_dictionary]  #this line replaces the two above it
    cell_location_file.close()
    dead_cell_file.close()

def log_start_parameters(run_info):
    run_info.write('Mating_occurs: '+str(mating_occurs)+'\n')
    run_info.write('Freq of mating: '+str(freq_of_mating) + '\n')
    run_info.write('Mutation occurs: '+ str(mutation_occurs) + '\n')
    run_info.write('Msat mutation rates: ' +str(mutation_rates) + '\n')
    run_info.write('Migration: '+str(migration) + '\n')
    run_info.write('Migration rate: ' +str(migration_rate) + '\n')
    run_info.write('Skew freqs: ' +str(skew_freqs) + '\n')
    run_info.write('Skew level: '+str(skewlevel) + '\n')
    run_info.write('Max pop size: '+str(max_pop_size) + '\n')
    run_info.write('Mutation rate cell properties: ' + str(mutation_rate_cell_properties) + '\n')
    run_info.write('Start hour: ' + str(start_hour) + '\n')
    run_info.write('Start day: ' + str(start_day) + '\n')
    run_info.write('Start year: ' +str(start_year) + '\n')
    run_info.write('Environment update interval: ' +str(environment_update_interval) + '\n')
    run_info.write('Nitrate source: ' +str(nitrate_source) + '\n')
    run_info.write('Run forward or backward: ' +str(run_forward_or_backward) + '\n')
    run_info.write('Cell death: ' +str(cell_death) + '\n')
    run_info.write('Cell growth: ' +str(cell_growth) + '\n')
    run_info.write('Track cell locs: ' + str(track_cell_locations) + '\n')
    run_info.write('Cell loc freq: ' +str(cell_loc_frequency) + '\n')
    run_info.write('Add real sample data: ' +str(add_real_sample_data) + '\n')
    run_info.write('Use multiple species: ' +str(use_multiple_species) + '\n')
    run_info.write('Species to model: ' +str(species_to_model) + '\n') 
    run_info.write('Length of run: '+str(length_of_run) + '\n')
    run_info.write('Simulation time interval: '+str(simulation_time_interval) + '\n')
    run_info.write('Environ width: ' +str(environ_width) + '\n')
    run_info.write('Environ length: ' +str(environ_length) + '\n')
    run_info.write('Environ origin offset: ' +str(environ_origin_offset) + '\n')
    run_info.write('Max depth: ' +str(max_depth) + '\n')
    run_info.write('Identical cells: ' +str(identical_cells) + '\n')
    run_info.write('Simulate a culture: ' + str(simulate_a_culture) + '\n')
    run_info.write('Max light: ' +str(max_light) + '\n')
    run_info.write('Env_salinity culture: ' +str(env_salinity) + '\n')
    run_info.write('Env temp culture: ' + str(env_temperature) + '\n')
    run_info.write('Env nutrients culture: ' + str(env_nutrients) + '\n')
    run_info.write('Number of cells to start: ' +str(num_cells_to_start) + '\n')
    run_info.write('Lights on/off:' + str(lights_on_off) + '\n')
    run_info.write('Future date:' + str(future_date) + '\n')
    return
##############################################################################################################################
##############################################################################################################################
#####Below here starts the program
# functions are defined by this point or in another library
numberofblooms = 0
field_blooms=[0]
allblooms = ["originalbloom"]



field_blooms[0] = cellmodel_bloom_cell_class.Kbrevis_bloom()
field_blooms[0].define_bloom({'bloom_id':0, 'cells':[], 'parentbloom':0, 'when_bloom_created':0, 'cell_growth':cell_growth, 'mutation_occurs':mutation_occurs,
                              'mutation_rate_cell_properties':mutation_rate_cell_properties, 'simulation_time_interval':simulation_time_interval,
                              'cell_death':cell_death})
field_blooms[0].cells = createfirstbloom(start_locations_and_concentrations)
numberofblooms = len(field_blooms)
#start_first_cells_in_world(field_blooms[0].cells)
print "Length of initial bloom:", len(field_blooms[0].cells)
print "bloom1 created"
runs = 60
bottlenecked_yet = 0
x = 1
y = 0
generation = 0
cell_conc_above = []
concen_array = []
locus1_trend = []
growth_trend = []
cnmax_trend = []
nmax_trend = []
c_divide_trend = []
n_divide_trend = []
ic_avg_trend = []
cnmin_trend = []
nmin_trend = []
rm_trend = []
pl_trend = []
pma_trend = []
pmb_trend = []
s250_trend = []
d_trend = []
izt_threshold_trend = []
hc_threshold_trend = []
cn_trend = []
n_trend = []
salinity_trend = []
temperature_trend = []
depth_trend = []
blm_number = 0
cell_x_vals = []
cell_y_vals = []
cell_z_vals = []
locs = [cell_x_vals, cell_y_vals, cell_z_vals]
##for blm in range(2):
##    concen_array.append([])
##for xtemp in range(int(max_depth), 1):  #old way
##    cell_conc_above.append(0)
    
output_figure = pylab.figure()
log_start_parameters(model_run_info)
max_pop_size = float(max_pop_size)
def main():
    global generation, x, mutate_original_two, freq_new_bloom, max_blooms, bottlenecked_yet, skew_freqs, field_blooms
    global migration_rate, migration_paths, gen_to_produce_new_blooms, sample_frequency, temp_bloom_creator_count
    global current_year, current_day, current_hour, world, environment_update_interval, environ_origin_offset
    global extra_bloom_dates, extra_bloom_locations
    while generation < length_of_run:
        cell_locs_plot, scatter_plot_locs = poll_cell_locations(field_blooms, world)
        if not simulate_a_culture:
            #print 'scatter_plot_locs', scatter_plot_locs[:10]
            world_map_maker(world, cell_locs_plot, scatter_plot_locs)
        
        y = 0
        bloomcount = 1
        blm_number = 0
        print current_year, current_day, current_hour
        model_run_info.write(str(current_year) +','+str(current_day) +', '+str(current_hour) + '\n')
        for bloomx in range(len(field_blooms)):
            bloom = field_blooms[bloomx]
            if (current_year, current_day, current_hour) in extra_bloom_dates and add_real_sample_data:
                print "Adding extra cells!", len(bloom.cells)
                bloom.cells = add_cells_to_bloom(bloom, extra_bloom_dates[(current_year, current_day, current_hour)])
            avg_growth_rate(bloom)
##            if generation > -1:
##                liveordie(bloom)
            carrying_capacity = len(bloom.cells) / max_pop_size
            if carrying_capacity > 0.9:
                carrying_capacity = 0.9
            print "Bloom", bloomcount, ":"
            bloomcount += 1
            if len(bloom.cells) < 100 or use_threads == False:
                if run_forward_or_backward == 'forward':
                    current_year, current_day, current_hour, world = bloom.process_one_day(blm_number, concen_array, max_depth, max_light, world, current_year, current_day,
                                                                                           current_hour, environment_update_interval, environ_origin_offset, nitrate_source,
                                                                                           world_latitude, world_longitude, run_forward_or_backward, culture, simulation_time_interval,
                                                                                           carrying_capacity, future_date,)
                else:
                    current_year, current_day, current_hour, world = bloom.process_one_day_reverse(blm_number, concen_array, max_depth, max_light, world, current_year, current_day,
                                                                                                   current_hour, environment_update_interval, environ_origin_offset, nitrate_source,
                                                                                                   world_latitude, world_longitude, run_forward_or_backward, culture, simulation_time_interval,
                                                                                                   carrying_capacity)
                current_hour = 3
            else:
                cells_to_be_processed = []
                time1 = time.time()
                if run_forward_or_backward == 'forward':
                    
                    for processor in range(num_processors_to_use):
                        if processor + 1 != num_processors_to_use:  #this one slices into the appropriate sizes
                            cells_to_be_processed.append(job_server.submit(bloom.process_one_day_multi,
                                                    (bloom.cells[(processor*len(bloom.cells)/num_processors_to_use):((processor+1)*len(bloom.cells)/num_processors_to_use)],
                                                     blm_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, environment_update_interval,
                                                     environ_origin_offset, nitrate_source, world_latitude, world_longitude, run_forward_or_backward, culture, simulation_time_interval,
                                                     carrying_capacity, future_date,),(Kbrcell, ),
                                                    ("cellmodel_fieldbloom_envir_HYCOM_v4_multispp","random","numpy","ephem",), globals=globals()))
                        else:                   #this one cleans up the end slice
                            cells_to_be_processed.append(job_server.submit(bloom.process_one_day_multi,
                                                    (bloom.cells[(processor*len(bloom.cells)/num_processors_to_use):],
                                                     blm_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, environment_update_interval,
                                                     environ_origin_offset, nitrate_source, world_latitude, world_longitude, run_forward_or_backward, culture, simulation_time_interval,
                                                     carrying_capacity, future_date,),(Kbrcell, ), 
                                                    ("cellmodel_fieldbloom_envir_HYCOM_v4_multispp","random","numpy","ephem",), globals=globals()))
                    current_hour = 3
                    current_day += 1
                    if current_day > 365:
                        current_year += 1
                        current_day = 1
                else:
                    for processor in range(num_processors_to_use):
                        if processor + 1 != num_processors_to_use:  #this one slices into the appropriate sizes
                            cells_to_be_processed.append(job_server.submit(bloom.process_one_day_multi_reverse,
                                                    (bloom.cells[(processor*len(bloom.cells)/num_processors_to_use):((processor+1)*len(bloom.cells)/num_processors_to_use)],
                                                     blm_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, environment_update_interval, environ_origin_offset, nitrate_source, world_latitude, world_longitude, run_forward_or_backward, culture, simulation_time_interval, carrying_capacity),(Kbrcell, ),
                                                    ("cellmodel_fieldbloom_envir_HYCOM_v4_multispp","random","numpy","ephem",), globals=globals()))
                        else:                   #this one cleans up the end slice
                            cells_to_be_processed.append(job_server.submit(bloom.process_one_day_multi_reverse,
                                                    (bloom.cells[(processor*len(bloom.cells)/num_processors_to_use):],
                                                     blm_number, concen_array, max_depth, max_light, world, current_year, current_day, current_hour, environment_update_interval, environ_origin_offset, nitrate_source, world_latitude, world_longitude, run_forward_or_backward, culture, simulation_time_interval, carrying_capacity),(Kbrcell, ),
                                                    ("cellmodel_fieldbloom_envir_HYCOM_v4_multispp","random","numpy","ephem",), globals=globals()))
                    current_hour = 3
                    current_day -= 1
                    if current_day < 1:
                        current_year -= 1
                        current_day = 365
                
                
                bloom_slices = []
                for cell in cells_to_be_processed:
                    bloom_slices.append(cell())
                temp_bloom = []
                
                for x1 in range(len(bloom_slices)):
                    if len(bloom_slices[x1]) > 0:
                        for y1 in range(len(bloom_slices[x1])):
                           temp_bloom.append(bloom_slices[x1][y1])
                field_blooms[y].cells = copy.deepcopy(temp_bloom)
                world = cellmodel_environment.update_world(world, current_year, current_day, 0, environ_origin_offset, nitrate_source, 'forward', culture, future_date, climate_change=False)
                print "One day took:", time.time()-time1
            #print 'Size of bloom ', y+1, ': ', len(temp_bloom)    
            #field_blooms[y].cells = copy.deepcopy(temp_bloom)
            
            blm_number += 1
            if x % 3 == 0:
                newgrowthrate(bloom)
            if mating_occurs == True:
                to_mate_or_not = int(len(bloom.cells)*freq_of_mating)
                if to_mate_or_not > random.random():
                    for reproduce in range(to_mate_or_not):
                        mating_in_bloom(bloom)
              
                        
            y += 1
        
            
        write_bloom_sizes(field_blooms)
        
        if x % 1 == 0:
            print "x:",x, "Number of blooms: ", len(field_blooms), "Bloom size:", len(field_blooms[0].cells)
        
        x += 1
        generation += 1    
    save_avg_data()    
    cell_location_file.close()

main()
#import cProfile
#cProfile.run('main()')
