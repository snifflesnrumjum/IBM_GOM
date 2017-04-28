#one off program to convert the NODC Nitrate csv files from the Gulf of Mexico into usable files for my Karenia model

#need to interpolate horizontally and vertically
#depths are : 0, 10, 20, 30, 50, 75, 100, 125, etc
#i plan to use scipy.interpolate.interp2d to return a function that takes x, y coordinates and returns a function

import numpy
import scipy.interpolate

def get_monthly_nitrate(julian_day):
    
    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    current_month = what_month_is_it(julian_day)
    
    indata = []
    infile = open('D:/CJunk/HYCOM_GOM_data_converted/NODC_nitrate/n' + months[current_month] + 'an01.csv')
    infile.next()
    infile.next()
    for line in infile:
        indata.append(line[:-1].split(','))
        for inval in range(len(indata[-1])):
            if indata[-1][inval] != '':
                indata[-1][inval] = float(indata[-1][inval])
                if inval == 0 or inval == 1:
                    indata[-1][inval] = str(indata[-1][inval])
                    
            else:
                indata[-1][inval] = 0
    lat = {'14.5':0, '15.5':1, '16.5':2, '17.5':3, '18.5':4, '19.5':5, '20.5':6, '21.5':7, '22.5':8, '23.5':9, '24.5':10, '25.5':11, '26.5':12,
           '27.5':13, '28.5':14, '29.5':15, '30.5':16, '31.5':17, '32.5':18, '33.5':19}
    lon = {'-98.5':0, '-97.5':1, '-96.5':2, '-95.5':3, '-94.5':4, '-93.5':5, '-92.5':6, '-91.5':7, '-90.5':8,  '-89.5':9, '-88.5':10, '-87.5':11, '-86.5':12,
           '-85.5':13, '-84.5':14, '-83.5':15, '-82.5':16, '-81.5':17, '-80.5':18, '-79.5':19}
    nitrate_data = numpy.zeros((len(lat), len(lon), 10), dtype=float)
    for indiv_lat in range(len(indata)):
        for depth in range(2, 12):
            if indata[indiv_lat][0] in lat and indata[indiv_lat][1] in lon:
                nitrate_data[lat[indata[indiv_lat][0]], lon[indata[indiv_lat][1]], depth-2] = indata[indiv_lat][depth]
    nitrate_interpolations = []
    nitrate_data = nitrate_data.T
    xlons = lon.keys()
    ylats = lat.keys()
    for s in range(len(lon.keys())):
        xlons[s] = float(xlons[s])
        ylats[s] = float(ylats[s])
    xlons.sort()
    ylats.sort()
    xlons = numpy.array(xlons)
    ylats = numpy.array(ylats)
    
    for depth in range(10):
        w = scipy.interpolate.RectBivariateSpline(xlons, ylats, nitrate_data[depth], s=0)
        nitrate_interpolations.append(w)
    return nitrate_interpolations
                               


def what_month_is_it(julian_day):
    year = 2009
    if year % 4 == 0:
        month_ranges = [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
    else:
        month_ranges = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    month_to_return = 0
    while julian_day > month_ranges[month_to_return]:
        month_to_return += 1

    return month_to_return



              
                

