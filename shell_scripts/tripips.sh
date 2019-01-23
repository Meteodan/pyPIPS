#!/bin/bash

# creates a datestamped (year,month,day,hour,minute of the script run time) folder, cds into it
# mkdir ~/Desktop/tridata/$(date +'oneHz%y%m%d%H%M') && cd $_

# grabs specified number of historic records from the query table, creates <filename>.html
# wget -nH --cut-dirs=100 "http://10.163.29.26/?command=TableDisplay&table=One_Hz&records=100" -O tablesnap.html && cd ..

# moves newly-created, datestamped folder to a permanent location
# mv one* ~/Desktop/exampledepotmountpoint/

# grabs specified number of historic records from the query table, creates <filename>.html
wget -nH --cut-dirs=100 "http://10.163.29.26/?command=TableDisplay&table=One_Hz&records=100" -O ~/sshfs_mounts/depot/data/Projects/TriPIPS/webdata/onesec/$(date +'oneHz%y%m%d%H%M.html')
