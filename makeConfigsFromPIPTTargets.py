"""

Script to generate config files from a PIPT format targets file

"""

import os, sys
import astropy.table as atpy

CONFIG="""# Parameters file for use with makeSALTSlitMaskFiles.py

# Cluster/field name/centre coords
# Output from the script will be written to MaskFiles_<name>/
name: $NAME
RADeg: $RADEG
decDeg: $DECDEG

# The format of the galaxy catalog file - either "FetchDR8", "FelipeS82", "FelipeDR8", "MattFITS", or "zCluster"
catalogFormat: "zCluster"

# Path to the galaxy catalog file - OR the database, if catalogFormat is 'zCluster'
galaxyCatalogFile: "DECaLSDR9"

# Below specify arbitrary cuts on catalog for top priority targets to put slits on
targetCuts:
    - 'r > 17'
    - 'r < 22'
    - 'r-z < 1'

# Cuts below are used to specify secondary targets, these are used when no more of
# the galaxies that pass targetCuts can be placed in a mask
fillerCuts:
    - 'r > 17'
    - 'r < 23'

# Slit dimensions
slitLengthArcsec: 9.0
slitWidthArcsec: 1.5

# Minimum gap to leave between slits
safetyArcsec: 1.0

# List of object IDs in the galaxy catalog to include in every mask (e.g., the BCG)
alwaysIncludeIDs: []

# Slit length for the above - if this was a BCG, might want a longer slit, for example
alwaysSlitLengthArcsec: 10.0

# Reference stars, for slit mask alignment - these should be selected from brightStars.reg
# To make brightStars.reg without making slit mask files, run the script with refStarIDs=[]
refStarIDs:
    - 5923
    - 4972
    - 5389
    - 4225

# Fluffy stuff
CMDCol: "r-z"

# Number of masks to make - objects will not overlap between masks unless they are in alwaysIncludeIDs
numMasks: 6

"""

outDir="configs"
tab=atpy.Table().read("pipt-format-targets.csv")

for row in tab:
    conf=CONFIG.replace("$NAME", row['target name'].replace(" ", "_"))
    conf=conf.replace("$RADEG", str(row['right ascension']))
    conf=conf.replace("$DECDEG", str(row['declination']))
    with open(outDir+os.path.sep+"%s.yml" % (row['target name'].replace(" ", "_")), "w") as outFile:
        outFile.write(conf)
