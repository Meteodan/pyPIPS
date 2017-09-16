# pyPIPS_merge
#
# This script is based on Sean Waugh's MATLAB script "POD_dataMerger.m"
# Original description follows:
#
# This program is written to combine and process the retrieved SD card data
# from the mobile PIPS. It will merge the two data file types together to
# create a single, line matched record of both information streams. It also
# derives the ambient wind direction, dewpoint, and RH (using the method
# outlined by Richardson et al. (1998).

# Richardson, S. J., S. E. Fredrickson, F. V. Brock, and J. A. Brotzge, 1998: 
# Combination temperature and relative humidity probes: avoiding large air 
# temperature errors and associated relative humidity errors. Preprints, 
# 10th Symposium on Meteorological Observations and Instrumentation, 
# Phoenix, AZ, USA, American Meteorological Society, 278-283.

# Original MATLAB script written by Sean Waugh: 1/25/2016
# Contact: sean.waugh@noaa.gov or 405-312-7585
#
# Python version written by Dan Dawson: 09/13/2017
# Contact: dandawson@purdue.edu

# Import required modules

import numpy as N
import os
import glob
import sys
import re
import csv
from datetime import datetime,timedelta

# Location of input files. You should only need to change these for your particular
# dataset.

PIPS_data_dir = '/Users/ddawson/Dropbox/PIPS_data/Irma2017/PIPS2B/converted/'
output_filename = 'PIPS_merged.txt'

# Function definitions

def process_onesec_record(row):
    """Processes a single row from a 1-s data file"""
    for field,value in row.iteritems():
        value = value.strip()
        # First, strip out all quotes from each field
        value = value.replace('"','')
        # Convert strings to numeric values where appropriate
        if(field != 'GPSDate' and field != 'GPSTime'):
            try:
                value = N.float(value)
                # This part isn't really needed, but is added for consistency with output
                # from the original Matlab script
                if(N.isnan(value)):
                    value = 'NaN'
                if(N.int(value) == value):
                    value = N.int(value)
            except:
                pass
        row[field] = value
    # Parse GPS latitude and longitude
    tokens = row['GPSLat'].strip().split()
    try:
        row['GPSLat'] = '{:7.4f}'.format(N.float(tokens[0])/100.)
        row['GPSLatHem'] = tokens[1]
    except:
        row['GPSLat'] = 'NaN'
        row['GPSLatHem'] = ''
    tokens = row['GPSLon'].strip().split()
    try:
        row['GPSLon'] = '{:7.4f}'.format(N.float(tokens[0])/100.)
        row['GPSLonHem'] = tokens[1]
    except:
        row['GPSLon'] = 'NaN'
        row['GPSLonHem'] = ''
    
    return row

def process_tensec_record(row):
    """Processes a single row from a 10-s data file"""
    for field,value in row.iteritems():
        # First, strip out all quotes from each field
        value = value.replace('"','')
        # Convert strings to numeric values where appropriate
        try:
            value = N.float(value)
            # This part isn't really needed, but is added for consistency with output
            # from the original Matlab script
            if(N.isnan(value)):
                value = 'NAN'
            if(N.int(value) == value):
                value = N.int(value)
        except:
            pass
        row[field] = value
    
    return row
    
def parseTimeStamp(timestring):
    """Parses a logger timestamp string and returns a datetime object"""
    date = timestring[0] # .strip('-')
    time = timestring[1] # .strip(':')
    year = N.int(date[:4])
    month = N.int(date[5:7])
    day = N.int(date[8:])
    hour = N.int(time[:2])
    min = N.int(time[3:5])
    sec = N.int(time[6:])
    return datetime(year,month,day,hour,min,sec)
   

# The following are taken from https://stackoverflow.com/questions/
# 5967500/how-to-correctly-sort-a-string-with-a-number-inside?noredirect=1&lq=1

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

# From https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression?noredirect=1&lq=1

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result
    
fieldnames_onesec = ['TIMESTAMP', 'RECORD', 'BattV', 'PTemp_C', 'WindDir', 'WS_ms', 
                     'WSDiag', 'FastTemp', 'SlowTemp', 'RH', 'Pressure(1)', 'FluxDirection', 
                     'GPSTime', 'GPSStatus', 'GPSLat', 'GPSLon', 'GPSSpd', 'GPSDir', 
                     'GPSDate', 'GPSMagVar', 'GPSAlt']
fieldnames_tensec = ['TIMESTAMP', 'RECORD', 'ParsivelStr']

fieldnames_output = ['TIMESTAMP','RECORD','BattV','PTemp_C','WindDir','WS_ms','WSDiag',
                     'FastTemp','SlowTemp','RH','Pressure(1)','FluxDirection','GPSTime',
                     'GPSStatus','GPSLat','GPSLatHem','GPSLon','GPSLonHem','GPSSpd','GPSDir',
                     'GPSDate','GPSMagVar','GPSAlt','WindDirAbs','Dewpoint','RHDer','ParsivelStr']

# Get lists of the 1-s and 10-s data files
filelist_onesec = glob.glob(PIPS_data_dir+'*_OneHz*')
filelist_tensec = glob.glob(PIPS_data_dir+'*_TenHz*')

# Sort the lists in natural order
filelist_onesec.sort(key=natural_keys)
filelist_tensec.sort(key=natural_keys)

# Check sizes of the lists. There should be the same number of files in both lists.
# If there isn't, something is wrong and the user needs to figure that out before running 
# this script.
# Actually could probably handle this by finding the union of the two lists,
# but I'm too lazy to mess with that right now.

nfiles = len(filelist_onesec)
if(nfiles != len(filelist_tensec)):
    sys.exit("Number of files for 1-s and 10-s data do not match! Fix this and run the script again!")

# Start reading in data. Start with the 1-s files

dict_onesec = {}
for i,file in enumerate(filelist_onesec):
    print "Reading 1-s data file: ",os.path.basename(file)   
    with open(file) as f:
        f.next() # Read and discard first header line
        fieldnames = f.next().strip().replace('"','').split(',') # The field names are contained in the second header line
        if (fieldnames != fieldnames_onesec):
            sys.exit("Something's wrong with this file, aborting!")
        #print fieldnames
        f.next() # Read and discard third header line
        f.next() # Read and discard fourth header line
        
        # The remaining lines contain all the data. Read and parse them into a dictionary
        # using the field names as keys
        # See https://stackoverflow.com/questions/14091387/creating-a-dictionary-from-a-csv-file
        reader = csv.DictReader(f,fieldnames=fieldnames)
        for row in reader:
            #print row
            # Process the dictionary of values in each row
            row = process_onesec_record(row)
            for column, value in row.iteritems():
                dict_onesec.setdefault(column, []).append(value)
                    
# Now read the 10-s files

dict_tensec = {}
for i,file in enumerate(filelist_tensec):
    print "Reading 10-s data file: ",os.path.basename(file)
    with open(file) as f:
        f.next() # Read and discard first header line
        fieldnames = f.next().strip().replace('"','').split(',') # The field names are contained in the second header line
        if (fieldnames != fieldnames_tensec):
            sys.exit("Something's wrong with this file, aborting!")
        #print fieldnames
        f.next() # Read and discard third header line
        f.next() # Read and discard fourth header line

        # The remaining lines contain all the data. Read and parse them into a dictionary
        # using the field names as keys
        # See https://stackoverflow.com/questions/14091387/creating-a-dictionary-from-a-csv-file
        reader = csv.DictReader(f,fieldnames=fieldnames)
        for row in reader:
            #print row
            # Process the dictionary of values in each row
            row = process_tensec_record(row)
            for column, value in row.iteritems():
                dict_tensec.setdefault(column, []).append(value)
        #print dict_tensec
        
# Now for the fun part! Merge all the records into a single 1-s file, matching the Parsivel
# records with the corresponding 1-s record as we go

numrecords = len(dict_onesec['TIMESTAMP']) # These are the number of lines (records) we
                                           # have to work with.

PIPS_outputfile = os.path.join(PIPS_data_dir,output_filename)

with open(PIPS_outputfile, 'w') as f:
    writer = csv.DictWriter(f,fieldnames=fieldnames_output)
    writer.writeheader()
    j = 0 # the j-index is for the Parsivel (10-s) records
    #firstRecord = True
    for i in xrange(numrecords): # Loop through the 1-s records
        # Parse the timestamp for the 1-s records and create a datetime object out of it
        timestring_onesec = dict_onesec['TIMESTAMP'][i].strip().split()
        datetime_onesec = parseTimeStamp(timestring_onesec)
        
        # Fill in known 1-s values into output row dictionary
        outputrow = {key: value[i] for key,value in dict_onesec.iteritems()}
        # Derive absolute wind direction, dewpoint, and RH
        # Note, the output of the string "NaN" instead of just letting it output
        # the float version is not really needed, but is done just for consistency with 
        # the original Matlab script (the float version outputs as "nan").
        try:
            winddirabs = N.mod(dict_onesec['FluxDirection'][i] + \
                               dict_onesec['WindDir'][i],360.)
            if(N.isnan(winddirabs)):
                outputrow['WindDirAbs'] = 'NaN'
            else: 
                outputrow['WindDirAbs'] = '{:4.1f}'.format(winddirabs).strip().rstrip('0').rstrip('.')
        except:
            outputrow['WindDirAbs'] = 'NaN'
            
        RH = dict_onesec['RH'][i]
        SlowT = dict_onesec['SlowTemp'][i]
        FastT = dict_onesec['FastTemp'][i]
        
        try:
            dewpoint = 243.04*(N.log(RH/100.)+((17.625*SlowT)/(243.04+SlowT)))/ \
                              (17.625-N.log(RH/100.)-((17.625*SlowT)/(243.04+SlowT)))
            if(N.isnan(dewpoint)):
                outputrow['Dewpoint'] = 'NaN'
            else:
                outputrow['Dewpoint'] = '{:7.4f}'.format(dewpoint).strip().rstrip('0').rstrip('.')
        except:
            outputrow['Dewpoint'] = 'NaN'
        
        try:
            RHDer = 100.*(N.exp((17.625*dewpoint)/(243.04+dewpoint))/ \
                          N.exp((17.625*FastT)/(243.04+FastT)))
            if(N.isnan(RHDer)):
                outputrow['RHDer'] = 'NaN'
            else:
                outputrow['RHDer'] =  '{:7.4f}'.format(RHDer).strip().rstrip('0').rstrip('.')
        except:
            outputrow['RHDer'] = 'NaN'
        
        # Next, we need to fill in the Parsivel data for the records that need it.
        nextParsivelTimeStamp = dict_tensec['TIMESTAMP'][j]
        # Parse the timestamp for the 10-s records and create a datetime object out of it
        timestring_tensec = nextParsivelTimeStamp.strip().split()
        datetime_tensec = parseTimeStamp(timestring_tensec)
        
        # To match the Parsivel timestamps with the appropriate 1-s timestamps, we
        # exploit the fact that the time is always increasing in the dataset (or it should
        # be; if not, something is wrong! This version of the script doesn't catch such 
        # errors, but it probably should be added at some point).
        # We basically loop one at a time through the 1-s records, and try to match the
        # timestamp up with the next Parsivel timestamp. If we find a match, great! Add it 
        # to the end of the 1-s record, and then go on to the next Parsivel timestamp. 
        # This way, we only sweep through both the 1-s and 10-s records once, instead of checking
        # the entire list of 1-s times for every 10-s record, which saves a lot of 
        # time (i.e search is of ~O(N) instead of O(N*n) where, N is the number of 1-s records
        # and n is the number of 10-s records).
        
        # This shouldn't really happen if the logger is working correctly, but in case the
        # First n Parsivel times are *earlier* than the first 1-s time, we have to increment
        # the j index until we get our first match.
        while(datetime_onesec > datetime_tensec):
            j += 1
            nextParsivelTimeStamp = dict_tensec['TIMESTAMP'][j]
            timestring_tensec = nextParsivelTimeStamp.strip().split()
            datetime_tensec = parseTimeStamp(timestring_tensec)
        
        if(datetime_onesec == datetime_tensec): # Found a match! Add Parsivel string to end
                                                # of 1-s record, and then increment the j counter.
            outputrow['ParsivelStr'] = dict_tensec['ParsivelStr'][j]
            j += 1
        else:   # No match, put a NaN at the end of the 1-s record.
            outputrow['ParsivelStr'] = 'NaN'
                
        #Write the output record
        print "Writing output record # ",i+1," of ",numrecords
        writer.writerow(outputrow)