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
# Latest update 08/30/2020
# Contact: dandawson@purdue.edu

# Import required modules

import os
import glob
from pyPIPS.pips_io import correct_PIPS
import sys
import re
import csv
from datetime import datetime
import numpy as np

# Location of input files. You should only need to change these for your particular
# dataset.

_PIPS_data_dir = '/Users/ddawson/Dropbox/PIPS_data/SPOTTR2019/051919_test/converted/PIPS1A/'
_output_filename = 'PIPS_merged.txt'

fieldnames_onesec = ['TIMESTAMP', 'RECORD', 'BattV', 'PTemp_C', 'WindDir', 'WS_ms',
                     'WSDiag', 'FastTemp', 'SlowTemp', 'RH', 'Pressure', 'FluxDirection',
                     'GPSTime', 'GPSStatus', 'GPSLat', 'GPSLon', 'GPSSpd', 'GPSDir',
                     'GPSDate', 'GPSMagVar', 'GPSAlt']
# In newer versions of the file, the "RECORD" field has been removed.
fieldnames_onesecv2 = fieldnames_onesec[:]
fieldnames_onesecv2.remove('RECORD')
fieldnames_onesec_TriPIPS = fieldnames_onesec[:]
fieldnames_onesec_TriPIPS.remove('FastTemp')
fieldnames_onesec_TriPIPSv2 = fieldnames_onesecv2[:]
fieldnames_onesec_TriPIPSv2.remove('FastTemp')

fieldnames_tensec = ['TIMESTAMP', 'RECORD', 'ParsivelStr']
fieldnames_tensecv2 = fieldnames_tensec[:]
fieldnames_tensecv2.remove('RECORD')

fieldnames_output = ['TIMESTAMP', 'BattV', 'PTemp_C', 'WindDir', 'WS_ms', 'WSDiag',
                     'FastTemp', 'SlowTemp', 'RH', 'Pressure', 'FluxDirection', 'GPSTime',
                     'GPSStatus', 'GPSLat', 'GPSLatHem', 'GPSLon', 'GPSLonHem', 'GPSSpd', 'GPSDir',
                     'GPSDate', 'GPSMagVar', 'GPSAlt', 'WindDirAbs', 'Dewpoint', 'RHDer',
                     'ParsivelStr']

fieldnames_output_TriPIPS = fieldnames_output[:]
fieldnames_output_TriPIPS.remove('FastTemp')

# Function definitions


def process_onesec_record(row):
    """Processes a single row from a 1-s data file"""
    # First, remove the "RECORD" item if it exists
    if 'RECORD' in row:
        row.pop('RECORD')
    for field, value in row.items():
        value = value.strip()
        # First, strip out all quotes from each field
        value = value.replace('"', '')
        # Convert strings to numeric values where appropriate
        if field not in ('GPSDate', 'GPSTime'):
            try:
                value = np.float(value)
                # This part isn't really needed, but is added for consistency with output
                # from the original Matlab script
                if np.isnan(value):
                    value = 'NaN'
                if np.int(value) == value:
                    value = np.int(value)
            except Exception:
                pass
        row[field] = value
    # Parse GPS latitude and longitude
    try:
        tokens = row['GPSLat'].strip().split()
        try:
            row['GPSLat'] = '{:7.4f}'.format(np.float(tokens[0]) / 100.)
            row['GPSLatHem'] = tokens[1]
        except Exception:
            row['GPSLat'] = 'NaN'
            row['GPSLatHem'] = ''
        tokens = row['GPSLon'].strip().split()
        try:
            row['GPSLon'] = '{:7.4f}'.format(np.float(tokens[0]) / 100.)
            row['GPSLonHem'] = tokens[1]
        except Exception:
            row['GPSLon'] = 'NaN'
            row['GPSLonHem'] = ''
    except Exception:
        row['GPSLat'] = 'NaN'
        row['GPSLatHem'] = ''
        row['GPSLon'] = 'NaN'
        row['GPSLonHem'] = ''
    return row


def process_tensec_record(row):
    """Processes a single row from a 10-s data file"""
    for field, value in row.items():
        # First, strip out all quotes from each field
        value = value.replace('"', '')
        value = value.replace('\n', '')
        # Convert strings to numeric values where appropriate
        try:
            value = np.float(value)
            # This part isn't really needed, but is added for consistency with output
            # from the original Matlab script
            if np.isnan(value):
                value = 'NAN'
            if np.int(value) == value:
                value = np.int(value)
        except Exception:
            pass
        row[field] = value

    return row


def parseTimeStamp(timestring):
    """Parses a logger timestamp string and returns a datetime object"""
    date = timestring[0]  # .strip('-')
    time = timestring[1]  # .strip(':')
    year = np.int(date[:4])
    month = np.int(date[5:7])
    day = np.int(date[8:])
    hour = np.int(time[:2])
    minute = np.int(time[3:5])
    sec = np.int(time[6:8])
    return datetime(year, month, day, hour, minute, sec)


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
    return [atoi(c) for c in re.split(r'(\d+)', text)]


# From https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
def sortby(X, Y):
    """Sorts list X by values in list Y"""
    return [x for _, x in sorted(zip(Y, X), key=lambda pair: pair[0])]

# From
# https://stackoverflow.com/questions/38987/
# how-to-merge-two-dictionaries-in-a-single-expression?noredirect=1&lq=1


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def readData(data_dir):
    """Reads the data from a list of 1-s and 10-s files and stores them in two dictionaries"""
    # Get lists of the 1-s and 10-s data files
    filelist_onesec = glob.glob(os.path.join(data_dir, '*_OneHz*'))
    filelist_tensec = glob.glob(os.path.join(data_dir, '*_TenHz*'))

    # Sort the lists in natural order
    filelist_onesec.sort(key=natural_keys)
    filelist_tensec.sort(key=natural_keys)

    # Check sizes of the lists. There should be the same number of files in both lists.
    # If there isn't, something is wrong and the user needs to figure that out before running
    # this script.
    # Actually could probably handle this by finding the union of the two lists,
    # but I'm too lazy to mess with that right now.

#     nfiles = len(filelist_onesec)
#     if(nfiles != len(filelist_tensec)):
#         sys.exit("Number of files for 1-s and 10-s data do not match! Fix this and run the script
#                  again!")

    # Start reading in data. Start with the 1-s files
    TriPIPS = False
    dict_onesec = {}
    for filename in filelist_onesec:
        print("Reading 1-s data file: ", os.path.basename(filename))
        with open(filename) as f:
            try:
                next(f)  # Read and discard first header line
            except Exception:
                print('File has no data. Skipping!')
                continue
            # The field names are contained in the second header line
            fieldnames = next(f).strip().replace('"', '').split(',')
            # Some corrupt files have split the first header line into two. Try to detect that
            # here and discard the extra line
            if 'TIMESTAMP' not in fieldnames[0]:
                fieldnames = next(f).strip().replace('"', '').split(',')
            # Some versions of the file have 'Pressure(1)' instead of 'Pressure'. Just change it
            # to 'Pressure'
            fieldnames = [field if field != 'Pressure(1)' else 'Pressure'
                          for field in fieldnames]
            if fieldnames == fieldnames_onesec or fieldnames == fieldnames_onesecv2:
                # print("We are dealing with an original PIPS data file!")
                TriPIPS = False
            elif (fieldnames == fieldnames_onesec_TriPIPS
                  or fieldnames == fieldnames_onesec_TriPIPSv2):
                # print("We are dealing with a TriPIPS data file (no FastTemp)!")
                TriPIPS = True
            else:
                print(fieldnames, '\n',
                      fieldnames_onesec)
                # sys.exit("Something's wrong with this file, aborting!")
                print("Something's wrong with this file, skipping!")
                continue
            # print fieldnames
            next(f)  # Read and discard third header line
            next(f)  # Read and discard fourth header line

            # The remaining lines contain all the data. Read and parse them into a dictionary
            # using the field names as keys
            # See
            # https://stackoverflow.com/questions/14091387/creating-a-dictionary-from-a-csv-file
            reader = csv.DictReader(f, fieldnames=fieldnames)
            for row in reader:
                # Process the dictionary of values in each row
                row = process_onesec_record(row)
                for column, value in row.items():
                    dict_onesec.setdefault(column, []).append(value)

    # Now read the 10-s files

    dict_tensec = {}
    for filename in filelist_tensec:
        print("Reading 10-s data file: ", os.path.basename(filename))
        with open(filename) as f:
            try:
                next(f)  # Read and discard first header line
            except Exception:
                print('File has no data. Skipping!')
                continue
            # The field names are contained in the second header line
            fieldnames = next(f).strip().replace('"', '').split(',')
            # Some corrupt files have split the first header line into two. Try to detect that
            # here and discard the extra line
            if 'TIMESTAMP' not in fieldnames[0]:
                fieldnames = next(f).strip().replace('"', '').split(',')
            if fieldnames not in (fieldnames_tensec, fieldnames_tensecv2):
                # print(fieldnames, fieldnames_tensec, fieldnames_tensecv2)
                # sys.exit("Something's wrong with this file, aborting!")
                print("Something's wrong with this file, skipping!")
                continue
            # print fieldnames
            next(f)  # Read and discard third header line
            next(f)  # Read and discard fourth header line

            # The remaining lines contain all the data. Read and parse them into a dictionary
            # using the field names as keys
            # See
            # https://stackoverflow.com/questions/14091387/creating-a-dictionary-from-a-csv-file
            reader = csv.DictReader(f, fieldnames=fieldnames)
            for row in reader:
                # Process the dictionary of values in each row
                row = process_tensec_record(row)
                for column, value in row.items():
                    dict_tensec.setdefault(column, []).append(value)
            # print dict_tensec
    return dict_onesec, dict_tensec, TriPIPS


def mergeData(data_dir, out_filename, verbose=False):

    # First, read the data from the files
    dict_onesec, dict_tensec, TriPIPS = readData(data_dir)

    if TriPIPS:
        output_fields = fieldnames_output_TriPIPS
    else:
        output_fields = fieldnames_output

    # Now for the fun part! Merge all the records into a single 1-s file, matching the Parsivel
    # records with the corresponding 1-s record as we go
    try:
        numrecords = len(dict_onesec['TIMESTAMP'])  # These are the number of lines (records) we
    except:
        numrecords = 0
    # have to work with.
    try:
        numparsivelrecords = len(dict_tensec['TIMESTAMP'])  # Number of Parsivel records
    except:
        numparsivelrecords = 0

    PIPS_outputfile = os.path.join(data_dir, out_filename)

    with open(PIPS_outputfile, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=output_fields)
        writer.writeheader()
        j = 0  # the j-index is for the Parsivel (10-s) records
        # firstRecord = True
        datetime_onesec_list = []
        for i in range(numrecords):  # Loop through the 1-s records
            # Parse the timestamp for the 1-s records and create a datetime object out of it
            timestring_onesec = dict_onesec['TIMESTAMP'][i].strip().split()
            # print('i', 'timestring_onesec', i, timestring_onesec)
            datetime_onesec = parseTimeStamp(timestring_onesec)
            datetime_onesec_list.append(datetime_onesec)

        # Sort the records in increasing order of time. (Some files are corrupted and have some
        # of their records out of order)
        indices = list(range(numrecords))
        indices_sorted = sortby(indices, datetime_onesec_list)
        if not indices == indices_sorted:
            print("Times out of order! Need to sort!")
            datetime_onesec_list.sort()
            for field, records in dict_onesec.copy().items():
                dict_onesec[field] = [records[i] for i in indices_sorted]

        # Now sort the 10-s data
        datetime_tensec_list = []
        for i in range(numparsivelrecords):
            timestring_tensec = dict_tensec['TIMESTAMP'][i].strip().split()
            datetime_tensec = parseTimeStamp(timestring_tensec)
            datetime_tensec_list.append(datetime_tensec)

        indices = list(range(numparsivelrecords))
        indices_sorted = sortby(indices, datetime_tensec_list)
        if not indices == indices_sorted:
            print("Times out of order! Need to sort!")
            datetime_tensec_list.sort()
            for field, records in dict_tensec.copy().items():
                dict_tensec[field] = [records[i] for i in indices_sorted]

        for i in range(numrecords):
            datetime_onesec = datetime_onesec_list[i]
            if i > 0:
                delta_t = (datetime_onesec_list[i] - datetime_onesec_list[i - 1]).total_seconds()
                if delta_t > 1.:
                    print("Time gap detected between {} and {}".format(
                        datetime_onesec_list[i - 1].strftime('%Y-%m-%d %H:%M:%S'),
                        datetime_onesec_list[i].strftime('%Y-%m-%d %H:%M:%S')))
            # Fill in known 1-s values into output row dictionary
            outputrow = {key: value[i] for key, value in dict_onesec.items()}
            # Derive absolute wind direction, dewpoint, and RH
            # Note, the output of the string "NaN" instead of just letting it output
            # the float version is not really needed, but is done just for consistency with
            # the original Matlab script (the float version outputs as "nan").
            try:
                winddirabs = np.mod(dict_onesec['FluxDirection'][i] +
                                    dict_onesec['WindDir'][i], 360.)
                if np.isnan(winddirabs):
                    outputrow['WindDirAbs'] = 'NaN'
                else:
                    outputrow['WindDirAbs'] = '{:4.1f}'.format(
                        winddirabs).strip().rstrip('0').rstrip('.')
            except Exception:
                outputrow['WindDirAbs'] = 'NaN'

            RH = dict_onesec['RH'][i]
            SlowT = dict_onesec['SlowTemp'][i]
            if not TriPIPS:
                FastT = dict_onesec['FastTemp'][i]

            try:
                dewpoint = 243.04 * (np.log(RH / 100.) + ((17.625 * SlowT) / (243.04 + SlowT))) / \
                    (17.625 - np.log(RH / 100.) - ((17.625 * SlowT) / (243.04 + SlowT)))
                if np.isnan(dewpoint):
                    outputrow['Dewpoint'] = 'NaN'
                else:
                    outputrow['Dewpoint'] = '{:7.4f}'.format(
                        dewpoint).strip().rstrip('0').rstrip('.')
            except Exception:
                outputrow['Dewpoint'] = 'NaN'

            try:
                RHDer = 100. * (np.exp((17.625 * dewpoint) / (243.04 + dewpoint)) /
                                np.exp((17.625 * FastT) / (243.04 + FastT)))
                if np.isnan(RHDer):
                    outputrow['RHDer'] = 'NaN'
                else:
                    outputrow['RHDer'] = '{:7.4f}'.format(RHDer).strip().rstrip('0').rstrip('.')
            except Exception:
                outputrow['RHDer'] = 'NaN'

            if verbose:
                print("Processing one-sec timestamp {}".format(
                    datetime_onesec.strftime('%Y-%m-%d %H:%M:%S')))

            if j < numparsivelrecords:  # Only do the following if we still have Parsivel data
                datetime_tensec = datetime_tensec_list[j]

                # Next, we need to fill in the Parsivel data for the records that need it.
                # nextParsivelTimeStamp = dict_tensec['TIMESTAMP'][j]
                # Parse the timestamp for the 10-s records and create a datetime object out of it
                # timestring_tensec = nextParsivelTimeStamp.strip().split()
                # datetime_tensec = parseTimeStamp(timestring_tensec)

                # To match the Parsivel timestamps with the appropriate 1-s timestamps, we
                # exploit the fact that the time is always increasing in the dataset (or it should
                # be; if not, something is wrong! This version of the script doesn't catch such
                # errors, but it probably should be added at some point).
                # We basically loop one at a time through the 1-s records, and try to match the
                # timestamp up with the next Parsivel timestamp. If we find a match, great! Add it
                # to the end of the 1-s record, and then go on to the next Parsivel timestamp.
                # This way, we only sweep through both the 1-s and 10-s records once, instead of
                # checking the entire list of 1-s times for every 10-s record, which saves a lot of
                # time (i.e search is of ~O(N) instead of O(N*n) where, N is the number of 1-s
                # records and n is the number of 10-s records).

                # This shouldn't really happen if the logger is working correctly, but in case the
                # First n Parsivel times are *earlier* than the first 1-s time, we have to increment
                # the j index until we get our first match.
                while datetime_onesec > datetime_tensec:
                    if verbose:
                        print("Initial syncing of one-sec and ten-sec times:")
                        print("onesec = {}, tensec = {}".format(
                            datetime_onesec.strftime('%Y-%m-%d %H:%M:%S'),
                            datetime_tensec.strftime('%Y-%m-%d %H:%M:%S')))
                        print(j)
                    j += 1
                    try:
                        datetime_tensec = datetime_tensec_list[j]
                    except:
                        break
#                     nextParsivelTimeStamp = dict_tensec['TIMESTAMP'][j]
#                     timestring_tensec = nextParsivelTimeStamp.strip().split()
#                     datetime_tensec = parseTimeStamp(timestring_tensec)

                # Found a match! Add Parsivel string to end
                # of 1-s record, and then increment the j
                # counter.
                if datetime_onesec == datetime_tensec:
                    if verbose:
                        print("Found match for one-sec and ten-sec data!")
                        print("onesec = {}, tensec = {}".format(
                            datetime_onesec.strftime('%Y-%m-%d %H:%M:%S'),
                            datetime_tensec.strftime('%Y-%m-%d %H:%M:%S')))
                    outputrow['ParsivelStr'] = dict_tensec['ParsivelStr'][j]
                    j += 1
                else:   # No match, put a NaN at the end of the 1-s record.
                    outputrow['ParsivelStr'] = 'NaN'
            else:
                outputrow['ParsivelStr'] = 'NaN'

            # Write the output record
            print("Writing output record # ", i + 1, " of ", numrecords)
            writer.writerow(outputrow)


if __name__ == "__main__":

    if len(sys.argv) > 1:
        PIPS_data_dir = sys.argv[1]
        output_filename = sys.argv[2]
        if len(sys.argv) > 3:
            serialnum = sys.argv[3]
            print(serialnum)
    else:
        PIPS_data_dir = _PIPS_data_dir
        output_filename = _output_filename
        serialnum = None

    mergeData(PIPS_data_dir, output_filename)
    if serialnum:
        print("Correcting offset Parsivel strings!")
        output_path = os.path.join(PIPS_data_dir, output_filename)
        corrected_filename = 'corrected_' + output_filename
        corrected_output_path = os.path.join(PIPS_data_dir, corrected_filename)
        correct_PIPS(serialnum, output_path, corrected_output_path)
