#!/usr/bin/env python3

"""

For SMERGERS - based on old code last used in 2017

"""

import glob
import astropy.io.fits as pyfits
import os
import sys
import numpy as np
import datetime
import time
import pylab as plt
import operator
from astLib import *
from zCluster import *
import astropy.table as atpy
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import urllib3
import math
import matplotlib.patches as patches
import string
from PIL import Image
import mechanize
from zCluster import retrievers
sys.path.append("rss-proptools")
from finder_chart import finderchart
import slitmask
import zipfile
import IPython
import yaml

#-------------------------------------------------------------------------------------------------------------
def fetchSDSSDR8Image(name, RADeg, decDeg, sizeArcmin = 12.0, JPEGFolder = "SDSSDR8Images", refetch = False):
    """Fetches the SDSS .jpg for the given image size using the casjobs webservice, stores them under
    jpegFolder. makeSDSSPlots loads these jpegs in, and use matplotlib to make them into plots with
    coord axes etc.
    
    """

    if os.path.exists(JPEGFolder) == False:
        os.makedirs(JPEGFolder)
    
    print(">>> Fetching SDSS DR8 image ...")
                 
    outFileName=JPEGFolder+os.path.sep+name.replace(" ", "_")+".jpg"
    
    # new - as reducing size of image, no point in getting them at ~6 Mb a pop
    SDSSWidth=1200.0
    SDSSScale=(sizeArcmin*60.0)/SDSSWidth # 0.396127
    
    if os.path.exists(outFileName) == False or refetch == True:        
        try:
            urlString="http://skyservice.pha.jhu.edu/DR8/ImgCutout/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)
            urlString=urlString+"&scale="+str(SDSSScale)+"&width="+str(int(SDSSWidth))+"&height="+str(int(SDSSWidth))
            urllib.request.urlretrieve(urlString, filename = outFileName)
        except:
            print("eh?")
            IPython.embed()
            sys.exit()
            print("... WARNING: couldn't get SDSS DR8 image ...")
            outFileName=None

#-------------------------------------------------------------------------------------------------------------
def fetchLegacySurveyImage(name, RADeg, decDeg, sizeArcmin = 12, sizePix = 1200, JPEGFolder = "RGBImages",
                           refetch = False, layer = 'ls-dr9', bands = 'grz'):
    """Fetches .jpg cut-out from legacysurvey.org sky viewer. Based on the code in sourcery.

    Valid layers include e.g. decals-dr7, des-dr1 etc..

    """

    outFileName=JPEGFolder+os.path.sep+name.replace(" ", "_")+".jpg"
    os.makedirs(JPEGFolder, exist_ok = True)

    http=urllib3.PoolManager()

    decalsWidth=sizePix
    decalsPixScale=(sizeArcmin*60.0)/float(decalsWidth)
    if os.path.exists(outFileName) == False or refetch == True:
        #http://legacysurvey.org/viewer/jpeg-cutout?ra=52.102810&dec=-21.670020&size=2048&layer=des-dr1&pixscale=0.3809&bands=grz
        urlString="http://legacysurvey.org/viewer/jpeg-cutout?ra=%.6f&dec=%.6f&size=%d&layer=%s&pixscale=%.4f&bands=%s" % (RADeg, decDeg, decalsWidth, layer, decalsPixScale, bands)
        print("... getting RGB image: %s" % (urlString))
        resp=http.request('GET', urlString)
        with open(outFileName, 'wb') as f:
            f.write(resp.data)
            f.close()

#-------------------------------------------------------------------------------------------------------------
def makeSlitsRGBPlot(name, RADeg, decDeg, slitMask, outFileName, JPEGFolder = "RGBImages", sizeArcmin = 12.0, remakePlots = True,
                      noAxes = False, plotClusterPos = True, 
                      figSize = (10, 7)):
    """Makes astPlot plots out of SDSS JPEG image, overplotting slitMask.
    
    If noAxes == True, do not plot coord axes and instead make the plot fill the plot area (for nemomontage).

    """
    
    #print ">>> Making SDSS DR7, DR8 or S82 plots ..."
    plotInteractiveStatus=plt.isinteractive()
    if plotInteractiveStatus == True:
        plt.matplotlib.interactive(False)
    
    sizeDeg=sizeArcmin/60.0
            
    # Load data
    inJPGPath=JPEGFolder+os.path.sep+name.replace(" ", "_")+".jpg"
    if os.path.exists(inJPGPath) == False:
        print("... couldn't find .jpg image - aborting ...")
        IPython.embed()
        sys.exit()
    
    im=Image.open(inJPGPath)
    data=np.array(im)
    try:
        data=np.flipud(data)
        if noAxes == False:             # If we're not plotting coords, don't flip or it will be backwards
            data=np.fliplr(data)
    except:
        "... something odd about image (1d?) - aborting ..."
        IPython.embed()
        sys.exit()
    
    R=data[:, :, 0]
    G=data[:, :, 1]
    B=data[:, :, 2]
    cutLevels=[[R.min(), R.max()], [G.min(), G.max()], [B.min(), B.max()]]
    
    # Make a WCS
    xSizeDeg, ySizeDeg=sizeArcmin/60.0, sizeArcmin/60.0
    xSizePix=R.shape[1]
    ySizePix=R.shape[0]
    xRefPix=xSizePix/2.0
    yRefPix=ySizePix/2.0
    xOutPixScale=xSizeDeg/xSizePix
    yOutPixScale=ySizeDeg/ySizePix
    newHead=pyfits.Header()
    newHead['NAXIS']=2
    newHead['NAXIS1']=xSizePix
    newHead['NAXIS2']=ySizePix
    newHead['CTYPE1']='RA---TAN'
    newHead['CTYPE2']='DEC--TAN'
    newHead['CRVAL1']=RADeg
    newHead['CRVAL2']=decDeg
    newHead['CRPIX1']=xRefPix+1
    newHead['CRPIX2']=yRefPix+1
    newHead['CDELT1']=xOutPixScale
    newHead['CDELT2']=xOutPixScale    # Makes more sense to use same pix scale
    newHead['CUNIT1']='DEG'
    newHead['CUNIT2']='DEG'
    wcs=astWCS.WCS(newHead, mode='pyfits')
    astImages.saveFITS(outFileName.replace(".png", ".fits"), R, wcs)
    
    # Make plot
    fig=plt.figure(figsize = figSize)
    if noAxes == True:
        axes=[0, 0, 1, 1]
        axesLabels=None
    else:
        axes=[0.1,0.1,0.8,0.8]
        axesLabels="sexagesimal"
    p=astPlots.ImagePlot([R, G, B], wcs, cutLevels = cutLevels, title = name, axes = axes, 
                            axesLabels = axesLabels)

    RAs=[]
    decs=[]
    for key in list(slitMask.keys()):
        obj=slitMask[key]
        if obj['type'] == 'object':
            RAs.append(obj['RADeg'])
            decs.append(obj['decDeg'])
    p.addPlotObjects(RAs, decs, 'slitObjects', size = sizeDeg/40.0*3600.0, symbol = 'box', color = "red")

    plt.savefig(outFileName)
    plt.close()
    
    plt.matplotlib.interactive(plotInteractiveStatus)
    
#-------------------------------------------------------------------------------------------------------------
def catalog2DS9(catalog, outFileName, constraintsList = [], addInfo = [], idKeyToUse = 'name', 
                RAKeyToUse = 'RADeg', decKeyToUse = 'decDeg', color = "cyan", includeRSSFoV = False,
                centreRADeg = None, centreDecDeg = None):
    """Converts a catalog containing object dictionaries into a ds9 region file. Objects will be labelled
    in the .reg file according to the idKeyToUse.
    
    If color == 'key', then use 'color' key in object dictionary to set color.
    
    constraintsList works the same way as selectFromCatalog function
    
    """

    # Cut catalog according to constraints
    cutCatalog=selectFromCatalog(catalog, constraintsList) 
    
    outFile=open(outFileName, "w")
    timeStamp=datetime.datetime.today().date().isoformat()
    comment="# DS9 region file\n"
    outFile.write(comment)
    outFile.write('global dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    for obj in cutCatalog:
        if len(addInfo) > 0:
            infoString=""
            for d in addInfo:
                if infoString != "":
                    infoString=infoString+" "
                if obj[d['key']] != None:
                    infoString=infoString+d['fmt'] % (obj[d['key']])
                else:
                    infoString=infoString+"%s" % (str(obj[d['key']]))
            infoString=" ["+infoString+"]"
        else:
            infoString=""
        if color == 'key':
            colorString=obj['color']
        else:
            colorString=color
        outFile.write("fk5;point(%.6f,%.6f) # point=boxcircle color={%s} text={%s%s}\n" \
                    % (obj[RAKeyToUse], obj[decKeyToUse], colorString, obj[idKeyToUse], infoString))
    if includeRSSFoV == True and centreRADeg is not None and centreDecDeg is not None:
        outFile.write('circle(%.6f,%.6f,%.6f") # color=%s\n' % (centreRADeg, centreDecDeg, 4.0*60.0, "green"))
    outFile.close()
        
#-------------------------------------------------------------------------------------------------------------
def selectFromCatalog(catalog, constraintsList):
    """Given a catalog (list of dictionaries representing objects), return a list of objects matching the
    given constraintsList. Each item in constraintsList is a string in the form:
    
    "key < value", "key > value", etc.
    
    Note that the spaces between key, operator (e.g. '<') and value are essential!
    
    """
    
    passedConstraint=catalog
    for constraintString in constraintsList:
        lastCatalog=passedConstraint
        passedConstraint=[]
        for obj in lastCatalog:         
            key, op, value=constraintString.split()
            if eval("obj['%s'] %s %s" % (key, op, value)) == True:
                passedConstraint.append(obj)
    
    return passedConstraint
    
#-------------------------------------------------------------------------------------------------------------
def catalogToSlitDataFile(catalog, slitWidth, slitLength, fileName):
    """Converts given object list into a file that can be imported into RSMT as slit data - doesn't include
    reference stars.
    
    """
    
    header="""# This file contains slit data for the RSS Slitmask Tool.
            #
            # Its column contain the following data:
            #
            # Column 1: right ascension of the slit center (in degrees or as hh:mm:ss.s)
            # Column 2: declination of the slit center (in degrees or as [+-]dd mm ss.s)
            # Column 3: slit width (in arcseconds)
            # Column 4: slit height (in arcseconds)
            # Column 5: slit tilt (in degrees)
            
            """

    outFile=open(fileName, "w")
    for obj in catalog:
        outFile.write("%.6f %.6f %.2f %.2f 0.0\n" % (obj['RADeg'], obj['decDeg'], slitWidth, slitLength))
    outFile.close()

#-------------------------------------------------------------------------------------------------------------
def makeTargetSlitDataFile(fileName, targetCatalog, fillerCatalog, slitWidth, slitLength, centreRADeg, centreDecDeg, safetyArcsec,
                           refStarSlitLength = 5.0, refStarSlitWidth = 5.0, refStarsList = [], alwaysSlitLength = 10.0, alwaysIncludeList =[],
                           previousMasksList = [], DSSImageFileName = None, imageCuts = [200, 400]):
    """This makes a slit data file for import with a sensible set of objects chosen to fill the mask. We're
    assuming the refStarSlitLength based on what it looks like in the tool, though surely that is far too 
    large?
    
    Note the BIG assumption made here is that we're designing masks with decDeg in Y direction - i.e. if the 
    PA is 90 deg or something then this will not work!
    
    If DSSImageFileName is not None, make a finder chart plot with astPlot
    
    Returns the slits info so we can exclude objects from future masks
    
    """

    print(">>> Making slit data file %s" % (fileName))
    
    # This will hold the objects that go into the mask
    slitsDict={}
    
    # First, add slits for ref stars - we can now include these in output
    refCount=0
    for obj in refStarsList:
        refCount=refCount+1
        slitsDict['ref%d' % (refCount)]={'slitLength': refStarSlitLength, 'slitWidth': refStarSlitWidth,
                                         'RADeg': obj['RADeg'], 'decDeg': obj['decDeg'],
                                         'decMin': obj['decDeg']-refStarSlitLength/3600.0/2.0,
                                         'decMax': obj['decDeg']+refStarSlitLength/3600.0/2.0,
                                         'type': 'refStar',
                                         'refStarID': obj['id'],
                                         'mag_r': obj['r']}
                                
    # Now put objects into the mask, starting from nearest centre, we add a flag key for whether object is
    # already included in mask
    catalogs=[targetCatalog, fillerCatalog]
    for catalog in catalogs:
        for obj in catalog:
            obj['rDeg']=astCoords.calcAngSepDeg(obj['RADeg'], obj['decDeg'], centreRADeg, centreDecDeg)
            obj['inMask']=False
            obj['overlapping']=False    # flag, for would be overlapping with object already added to mask
        catalog=sorted(catalog, key=operator.itemgetter('rDeg'))
        
        # If we had supplied a list of objects which are already in mask, this is where we would apply it...
        for maskDict in previousMasksList:
            for key in list(maskDict.keys()):
                slit=maskDict[key]
                for obj in catalog:
                    if obj['RADeg'] == slit['RADeg'] and obj['decDeg'] == slit['decDeg']:
                        obj['inMask']=True
    
    # Ditto for adding in objects to include in every mask
    # Be careful here, we're not checking these for overlaps
    alwaysCount=1
    for selObj in alwaysIncludeList:
        try:
            slitsDict['always%d' % (alwaysCount)]={'slitLength': alwaysSlitLength, 'slitWidth': slitWidth,
                                'RADeg': selObj['RADeg'], 'decDeg': selObj['decDeg'],
                                'decMin': selObj['decDeg']-alwaysSlitLength/3600.0/2.0,
                                'decMax': selObj['decDeg']+alwaysSlitLength/3600.0/2.0,
                                'type': 'object',
                                'catalogID': selObj['id'],
                                'mag_r': selObj['r']}
        except:
            print("Eh?")
            ipshell()
            sys.exit()
        alwaysCount=alwaysCount+1
        
    # Add galaxy targets to mask according to priority
    fillMaskWithTargets(targetCatalog, slitsDict, centreRADeg, centreDecDeg, slitWidth, slitLength, safetyArcsec)
    fillMaskWithTargets(fillerCatalog, slitsDict, centreRADeg, centreDecDeg, slitWidth, slitLength, safetyArcsec)
    
    # Write to file - old:
    #header="""# This file contains slit data for the RSS Slitmask Tool.
##
## Its columns contain the following data.
##
## For a slit:
## Column 1: the string SLI
## Column 2: right ascension of the slit center (in degrees or as hh:mm:ss.s)
## Column 3: declination of the slit center (in degrees or as [+-]dd mm ss.s)
## Column 4: slit width (in arcseconds)
## Column 5: slit height (in arcseconds)
## Column 6: slit tilt (in degrees)
##
## For a reference star:
## Column 1: the string REF
## Column 2: right ascension of the reference star center (in degrees or as hh:mm:ss.s)
## Column 3: declination of the reference star center (in degrees or as [+-]dd mm ss.s )
            #"""           
    #outFile=file(fileName, "w")
    #for key in slitsDict:
        #slit=slitsDict[key]
        #if slit['type'] == 'object':
            #outFile.write("SLI %.6f %.6f %.2f %.2f 0.0\n" % (slit['RADeg'], slit['decDeg'], 
                                                         #slit['slitWidth'], slit['slitLength']))
        #elif slit['type'] == 'refStar':
            #outFile.write("REF %.6f %.6f\n" % (slit['RADeg'], slit['decDeg']))
    #outFile.close()
    
    #---
    #Object Name -- (format=%16s, char. string)  No double quotes allowed. UTF-8 encoded.  Maximum of 16 characters.
    #RA -- (real) hh:mm:ss.sss or an ordinary float (eg., 237.8156)
    #Dec -- (real) sdd:mm:ss.ss or an ordinary float (eg., -45.9345233)
    #Equinox -- (real) eg., 2000; currently there is no provision for precession, but this could be added on short notice if needed
    #Magnitude -- (real) eg., 21.99
    #PassBand -- (char) eg., V; currently truncated to one character
    #Priority_Code -- (int) as follows:
        #1  A value one means the object is pre-selected for the mask. 
        #1-0  Weighting for the priority of an object   
        #0 goes along for the ride (do not use for program targets!); others ignored.
        #-1 indicates a reference star
    # New:
    #header="#Object_Name    RA  Dec Equinox Magnitude   PassBand    Priority_Code\n"
    #outFile=open(fileName, "w")
    #outFile.write(header)
    #for key in slitsDict:
        #slit=slitsDict[key]
        #line="%.6f %.6f 2000 %.3f r" % (slit['RADeg'], slit['decDeg'], slit['mag_r'])
        #if slit['type'] == 'object':
            #line="tID"+str(slit['catalogID'])+" "+line+" 1\n"
        #elif slit['type'] == 'refStar':
            #line="sID"+str(slit['refStarID'])+" "+line+" -1\n"
        #else:
            #raise Exception("didn't understand slit['type'] %s" % (slit['type']))
        #outFile.write(line)
    #outFile.close()
    # Long format
    header="#name targ_ra targ_dec equinox mag band priority width length tilt\n"
    with open(fileName, "w") as outFile:
        outFile.write(header)
        for key in slitsDict:
            slit=slitsDict[key]
            line="%.6f %.6f 2000 %.3f r" % (slit['RADeg'], slit['decDeg'], slit['mag_r'])
            if slit['type'] == 'object':
                line="tID"+str(slit['catalogID'])+" "+line+" 1"
            elif slit['type'] == 'refStar':
                line="sID"+str(slit['refStarID'])+" "+line+" -1"
            else:
                raise Exception("didn't understand slit['type'] %s" % (slit['type']))
            line=line+" %.3f %.3f 0\n" % (slit['slitWidth'], slit['slitLength'])
            outFile.write(line)
    #---

    # Write as DS9 .reg file here, complete with circle for RSS fov etc..
    outFile=open(fileName.replace(".txt", ".reg"), "w")
    outFile.write("# Region file format: DS9 version 3.0\n")
    outFile.write("global font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n")
    for key in slitsDict:
        slit=slitsDict[key]
        if slit['type'] == 'refStar':
            color='white'
        elif slit['type'] == 'object':
            color='cyan'
        outFile.write("fk5;box("+str(slit['RADeg'])+", "+str(slit['decDeg'])+", "+str(slit['slitWidth'])+"\", "+str(slit['slitLength'])+"\", "+str(0.0)+") # text={"+str(key)+"} color="+str(color)+"\n")
    outFile.write('circle(%.6f,%.6f,%.6f") # color=%s\n' % (centreRADeg, centreDecDeg, 4.0*60.0, "green"))
    outFile.close()
    
    
    # Make a finding chart - note we don't draw on slit dimensions here, as that means upgrading ImagePlot
    # We only draw reference stars for now
    if DSSImageFileName != None:
        img=pyfits.open(DSSImageFileName)
        wcs=astWCS.WCS(DSSImageFileName)
        
        clip=astImages.clipImageSectionWCS(img[0].data, wcs, centreRADeg, centreDecDeg, 12/60.0)
        fig=plt.figure(figsize=(12,12))
        p=astPlots.ImagePlot(clip['data'], clip['wcs'], cutLevels = imageCuts, 
                             title = os.path.split(fileName)[-1].replace(".txt", "")+" (PA=180)", colorMapName = "gray_r")
        
        p.addPlotObjects([centreRADeg], [centreDecDeg], 'fov', symbol='circle', size=8.0*60, width=1.0, 
                         color='blue')
                         
        refRAs=[]
        refDecs=[]
        for key in list(slitsDict.keys()):
            slit=slitsDict[key]
            if slit['type'] == 'refStar':
                refRAs.append(slit['RADeg'])
                refDecs.append(slit['decDeg'])
        p.addPlotObjects(refRAs, refDecs, 'refStars', symbol='box', size=10.0, width=1.0, 
                         color='red')
     
        p.draw()

        for key in list(slitsDict.keys()):
            slit=slitsDict[key]
            if slit['type'] == 'object':
                x, y=wcs.wcs2pix(slit['RADeg'], slit['decDeg'])
                widthPix=(slit['slitWidth']/3600.0)/wcs.getPixelSizeDeg()
                lengthPix=(slit['slitLength']/3600.0)/wcs.getPixelSizeDeg()
                c=patches.Rectangle((x-widthPix/2, y-lengthPix/2), widthPix, lengthPix, 
                                 fill=False, edgecolor='green', linewidth=1.0)
                p.axes.add_patch(c)
                      
        
        p.save(fileName.replace(".txt", ".png"))
        plt.close()
        
    # Return this so that we can exclude these slits from the masks we make next
    return slitsDict

#-------------------------------------------------------------------------------------------------------------
def fillMaskWithTargets(catalog, slitsDict, centreRADeg, centreDecDeg, slitWidth, slitLength, safetyArcsec):
    """This adds slits to the mask, prioritising objects closer to the centre and avoiding overlaps.
    Can run this multiple times with different catalogs to prioritise.
    
    """

    # Work out current object count
    objectCount=0
    for key in list(slitsDict.keys()):
        if "obj" in key:
            objID=int(key.split("obj")[-1])
            if objID > objectCount:
                objectCount=objID
    objectCount=objectCount+1
    
    # The way this works is we try to add slits until the mask area is full
    maxRDeg=3.9/60.0 # 4' radius from centre is max fov
    maxXDeg=2.6/60.0 # max distance in X direction from centre (assumes PA = 180.)
    cutCatalog=selectFromCatalog(catalog, ["rDeg < %.6f" % (maxRDeg), 'inMask == False', 'overlapping == False'])
    maskFull=False  # by full, here we may mean of this particular type of target
    while maskFull == False:
        # Mark as overlapping all objects too close to objects already in mask
        for obj in cutCatalog:
            for key in list(slitsDict.keys()):
                rDegDecMin=astCoords.calcAngSepDeg(centreRADeg, obj['decDeg'], 
                                             centreRADeg, slitsDict[key]['decMin'])
                rDegDecMax=astCoords.calcAngSepDeg(centreRADeg, obj['decDeg'], 
                                             centreRADeg, slitsDict[key]['decMax'])
                if min([rDegDecMin, rDegDecMax]) < slitsDict[key]['slitLength']/3600.0/2.0+safetyArcsec/3600.0:
                    obj['overlapping']=True
        # Now add to mask the nearest object to central position that's not overlapping
        cutCatalog=selectFromCatalog(cutCatalog, ["rDeg < %.6f" % (maxRDeg), 'inMask == False', 'overlapping == False'])
        XMin, XMax, blah, blah=astCoords.calcRADecSearchBox(centreRADeg, centreDecDeg, maxXDeg)
        cutCatalog=selectFromCatalog(cutCatalog, ["RADeg > %.6f" % (XMin)])
        cutCatalog=selectFromCatalog(cutCatalog, ["RADeg < %.6f" % (XMax)])
        cutCatalog=sorted(cutCatalog, key=operator.itemgetter('rDeg'))
        if len(cutCatalog) > 0:
            print("... adding object %d ..." % (objectCount))
            selObj=cutCatalog[0]
            selObj['inMask']=True
            slitsDict['obj%d' % (objectCount)]={'slitLength': slitLength, 'slitWidth': slitWidth,
                                        'RADeg': selObj['RADeg'], 'decDeg': selObj['decDeg'],
                                        'decMin': selObj['decDeg']-slitLength/3600.0/2.0,
                                        'decMax': selObj['decDeg']+slitLength/3600.0/2.0,
                                        'type': 'object',
                                        'catalogID': selObj['id'],
                                        'mag_r': selObj['r']}
            objectCount=objectCount+1
        else:
            maskFull=True
    
#-------------------------------------------------------------------------------------------------------------
def makeRefStarsList(starsPhotList, namesList, label, outDir):
    """Given a list of object names, make a list of ref stars for feeding into the mask target selection 
    routine.
    
    label is for appending to .reg file name
    
    """
    
    refStars=[]
    print(">>> Reference star properties:")
    for obj in starsPhotList:
        if obj['id'] in namesList:
            refStars.append(obj)
            print("... %s (r=%.2f) ..." % (obj['id'], obj['r']))
    catalog2DS9(refStars, outDir+os.path.sep+"refStars_%s.reg" % (label), idKeyToUse = 'id', color = 'white')

    return refStars

#-------------------------------------------------------------------------------------------------------------
def SDSSRetriever(RADeg, decDeg, halfBoxSizeDeg = 9.0/60.0, DR = 7, objType = 'star'):
    """Retrieves SDSS main photometry at the given position.
        
    """
    
    CACHE_DIR="SDSSCache"
    if os.path.exists(CACHE_DIR) == False:
        os.makedirs(CACHE_DIR)
    
    if DR == 7:
        url='http://cas.sdss.org/astrodr7/en/tools/search/x_sql.asp'
        outFileName=CACHE_DIR+os.path.sep+"SDSSDR7_%.4f_%.4f_%.4f_%s.csv" % (RADeg, decDeg, halfBoxSizeDeg, objType)
    elif DR == 8:
        url='http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp'
        outFileName=CACHE_DIR+os.path.sep+"SDSSDR8_%.4f_%.4f_%.4f_%s.csv" % (RADeg, decDeg, halfBoxSizeDeg, objType)
    elif DR == 9:
        url='http://skyserver.sdss3.org/dr9/en/tools/search/x_sql.asp'
        outFileName=CACHE_DIR+os.path.sep+"SDSSDR9_%.4f_%.4f_%.4f_%s.csv" % (RADeg, decDeg, halfBoxSizeDeg, objType)
        
    print("... getting SDSS DR%d photometry (file: %s) ..." % (DR, outFileName))

    # First, check if we previously downloaded a catalog and it went wrong - if it did, delete the file so
    # we can try fetching it again
    #if os.path.exists(outFileName) == True:
        #inFile=file(outFileName, "r")
        #lines=inFile.readlines()
        #inFile.close()
        #if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
            #os.remove(outFileName)
        
    if os.path.exists(outFileName) == False:
    
        # Clean galaxy photometry query - note flags for r-band only, may want to change
        # Assuming for limits here on equator, so don't bother with cos(dec)
        #
        # We may want to add something that does multiple queries if we want a bigger area from the standard
        # SDSS query interface
        RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, halfBoxSizeDeg)
        if objType == 'star':
            sql="""SELECT ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z,flags_r,run,extinction_g,extinction_r,extinction_i 
                FROM Star 
                WHERE 
                ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f 
                AND ((flags_r & 0x10000000) != 0) 
                -- detected in BINNED1 
                AND ((flags_r & 0x8100000c00a4) = 0) 
                -- not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP, 
                -- SATURATED, or BAD_COUNTS_ERROR 
                AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2)) 
                -- not DEBLEND_NOPEAK or small PSF error 
                -- (substitute psfmagerr in other band as appropriate) 
                AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0) 
                -- not INTERP_CENTER or not COSMIC_RAY
                """ % (RAMin, RAMax, decMin, decMax)
        elif objType == 'galaxy':
            #sql="""SELECT ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z,flags_r,run,extinction_g,extinction_r,extinction_i
                #FROM Galaxy 
                #WHERE 
                #ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f 
                #AND ((flags_r & 0x10000000) != 0) 
                #-- detected in BINNED1 
                #AND ((flags_r & 0x8100000c00a0) = 0) 
                #-- not NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP, SATURATED, 
                #-- or BAD_COUNTS_ERROR. 
                #-- if you want to accept objects with interpolation problems for PSF mags, 
                #-- change this to: AND ((flags_r & 0x800a0) = 0) 
                #AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2)) 
                #-- not DEBLEND_NOPEAK or small PSF error 
                #-- (substitute psfmagerr in other band as appropriate) 
                #AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0) 
                #-- not INTERP_CENTER or not COSMIC_RAY - omit this AND clause if you want to 
                #-- accept objects with interpolation problems for PSF mags.
                #AND Err_g < 0.2 AND Err_r < 0.2 AND Err_i < 0.2 
                #""" % (RAMin, RAMax, decMin, decMax)
            # Seem to be missing obvious galaxies with the above
            sql="""SELECT ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z,flags_r,run,extinction_g,extinction_r,extinction_i
                FROM Galaxy 
                WHERE 
                ra BETWEEN %.6f and %.6f AND dec BETWEEN %.6f and %.6f 
                """ % (RAMin, RAMax, decMin, decMax)
        else:
            raise Exception("didn't understand object type")
        
        # Filter SQL so that it'll work
        fsql = ''
        for line in sql.split('\n'):
            fsql += line.split('--')[0] + ' ' + os.linesep;
    
        params=urllib.parse.urlencode({'cmd': fsql, 'format': "csv"})
        try:
            response=urllib.request.urlopen(url+'?%s' % (params))
        except:
            print("Connection reset by peer?")
            ipshell()
        lines=response.read()
        lines=lines.split("\n")

        outFile=open(outFileName, "w")
        for line in lines:
            outFile.write(line+"\n")
        outFile.close()
    
    else:
        
        inFile=open(outFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    
    # Parse .csv into catalog
    if lines[0].find("No objects have been found") != -1 or len(lines) > 1 and lines[1][:5] == "ERROR":
        catalog=None
    else:
        catalog=[]
        for line in lines[1:]: # first line always heading
            if len(line) > 3:
                photDict={}
                bits=line.replace("\n", "").split(",")
                try:
                    photDict['RADeg']=float(bits[0])
                except:
                    if lines[1][:46] == '"ERROR: Maximum 60 queries allowed per minute.':
                        print("... exceeded server queries per minute limit - waiting ...")
                        time.sleep(70)
                        os.remove(outFileName)
                        return "retry"
                    else:
                        print("what?")
                        IPython.embed()
                        sys.exit()
                photDict['decDeg']=float(bits[1])
                photDict['u']=float(bits[2])
                photDict['g']=float(bits[3])
                photDict['r']=float(bits[4])
                photDict['i']=float(bits[5])
                photDict['z']=float(bits[6])
                photDict['uErr']=float(bits[7])
                photDict['gErr']=float(bits[8])
                photDict['rErr']=float(bits[9])
                photDict['iErr']=float(bits[10])
                photDict['zErr']=float(bits[11])
                photDict['run']=int(bits[13])
                photDict['extinction_g']=float(bits[14])
                photDict['extinction_r']=float(bits[15])
                photDict['extinction_i']=float(bits[16])                
                catalog.append(photDict)

    # Doing this because for some reason otherwise ID numbers are not consistent between runs,
    # i.e. lines in SDSS catalog file are not necessarily returned in same order if query rerun
    catalog=sorted(catalog, key=operator.itemgetter('decDeg'))
    catalog=sorted(catalog, key=operator.itemgetter('RADeg'))
    idCount=1
    for obj in catalog:
        obj['id']=idCount
        idCount=idCount+1

    catalog2DS9(catalog, outFileName.replace(".csv", ".reg"), idKeyToUse='id', 
                            addInfo = [{'key': 'r', 'fmt': '%.3f'}, {'key': 'rErr', 'fmt': '%.3f'}], color='red')
                                    
    return catalog

#-------------------------------------------------------------------------------------------------------------
def fetchDSSImage(RADeg, decDeg, outFileName, sizeArcmin = 10.0, refetch = False, verbose = True):
    """Fetches DSS image from Skyview
    
    """
    
    sizeDeg=sizeArcmin/60.0
    
    count=0
    sleepiness=10
    retries=1   # don't flog a dead horse
    survey="DSS"
    if os.path.exists(outFileName) == False:# or checkImageDownloadSuccess(outFileName) == False:
        if verbose == True: print(">>> Fetching Skyview image ...")
        success=False
        attempts=0
        #while success == False and attempts < retries:
        urlString="http://skyview.gsfc.nasa.gov/cgi-bin/images?Position="+str(RADeg)+","+str(decDeg)
        urlString=urlString+"&Size="+str(sizeDeg)+"&Pixels=1000&Projection=Tan&Survey="
        urlString=urlString+survey+"&Coordinates=J2000&Return=fits"
        print(urlString)
        urllib.request.urlretrieve(urlString, filename = outFileName)
        #success=checkImageDownloadSuccess(outFileName)
        #if success == False:
            #print("... failed to fetch image - trying again in %d sec ..." % (sleepiness))
            #time.sleep(sleepiness)  # be kind to Skyview
            #attempts=attempts+1

#-------------------------------------------------------------------------------------------------------------
def fetchDSSImageFromESO(RADeg, decDeg, outFileName, sizeArcmin = 10.0, refetch = False, verbose = True):
    """Stupid US government...
    
    """
    
    IPython.embed()
    sys.exit()
        
    # Except this is disallowed by robots.txt
    #br=mechanize.Browser()
    #br.open("http://archive.eso.org/dss/dss")
    #br.select_form(nr=0)
    #br.form.set_value(astCoords.decimal2hms(RADeg, " "), name = "ra")
    #br.form.set_value(astCoords.decimal2dms(decDeg, " "), name = "dec")
    #br.form.set_value(str(sizeArcmin), name = "x")
    #br.form.set_value(str(sizeArcmin), name = "y")
    #response=br.submit()
    
   
            
    
    #http://archive.eso.org/dss/dss
    
#-------------------------------------------------------------------------------------------------------------
def checkImageDownloadSuccess(fileName):
    """Checks if the image fileName was successfully downloaded or not. If it wasn't, we'll have an html 
    file masquerading as a .jpg (this may need tweaking to work with e.g. SDSS, works for Skyview).
    
    """
    
    if os.path.exists(fileName) == True:
        inFile=open(fileName, "r")
        line=inFile.readline()
        inFile.close()
        if line[:10] == 'SIMPLE  = ':
            success=True
        else:
            success=False
    else:
        success=False
    
    return success  

#------------------------------------------------------------------------------------------------------------
def makeRSMTFile(targetDataFileName, imageFileName, targetName, RADeg, decDeg, proposer, proposalCode, maskName, outDir):
    """Write slit mask files in the format SALT uses (XML file and finder chart in ZIP archive).

    The output file name is outDir/maskName.rsmt

    """
    sm=slitmask.SlitMask()

    sm.creator='rss-mask-design'
    sm.proposer=proposer
    sm.proposal_code=proposalCode
    sm.target_name=targetName
    sm.mask_name=maskName

    sm.add_center_ra(RADeg)
    sm.add_center_dec(decDeg)

    if sm.center_dec > -35:
        sm.add_position_angle(180)
    else:
        sm.add_position_angle(0)

    sm.slitlets.readascii(targetDataFileName, form = "long")
    sm.set_MaskPosition()
    sm.outFoV()

    print("\nChecking mask validity in a moment... if your mask is invalid due to not choosing reference stars, use:\n")
    print("    ds9 %s -regions %s\n" % (imageFileName, os.path.dirname(imageFileName)+os.path.sep+"brightStars.reg"))
    print("to inspect the image and choose suitable reference stars (add these to the config file under refStarIDs)\n")

    sm.validate()

    with open("Slitmask.xml", "w") as outFile:
        outFile.write(sm.writexml())

    finderchart("Slitmask.xml", image = imageFileName, outfile = "Slitmask.png")
    plt.close()

    zf=zipfile.ZipFile(RSMTDir+os.path.sep+maskName+".rsmt", mode = "w")
    zf.write("Slitmask.png")
    zf.write("Slitmask.xml")
    zf.close()

    os.remove("Slitmask.xml")
    os.remove("Slitmask.png")

#-------------------------------------------------------------------------------------------------------------
# Main
if len(sys.argv) < 2:
    print("Run: % makeSALTSlitMaskFiles.py <.yml file>")
else:
    
    parFileName=sys.argv[1]
    
    with open(parFileName, "r") as stream:
        clusterDict=yaml.safe_load(stream)

    # Put everthing under here
    outDir="MaskFiles_"+clusterDict['name'].replace(" ", "_")
    os.makedirs(outDir, exist_ok = True)

    RSMTDir="RSMTFiles"
    os.makedirs(RSMTDir, exist_ok = True)

    # Cluster position
    cRADeg=clusterDict['RADeg']
    cDecDeg=clusterDict['decDeg']

    RArcmin=5.0 # just go as wide as we like

    # Load in galaxies catalog, which has bpz photo-z info (z_p and odds)
    catalogFormat=clusterDict['catalogFormat']
    if catalogFormat == 'MattFITS':
        tab=atpy.Table().read(clusterDict['galaxyCatalogFile'])    
        IPython.embed()
        sys.exit()
    elif catalogFormat == 'Mathilde':
        tab=atpy.Table().read(clusterDict['galaxyCatalogFile'])    
        # No ID column, so add one
        tab.add_column(atpy.Column(np.arange(1, len(tab)+1), "id"))
        tab.rename_column("DECDeg", "decDeg")
        # Fudge: make mag column names lower case, so we don't have to go through changing 'i' -> 'I'
        # Actually, easiest if we just copy columns, so then don't have to think in the .par file.
        # We should just add a 'magColumn' parameter to the .par file at some point...
        copyColumns={'RC': 'r', 'IC': 'i'}
        for key in list(copyColumns.keys()):
            tab.add_column(copyColumns[key], tab[key]) 
        tab.add_column(atpy.Column(tab['r'], 'pMag_r')) # yes, another fudge
        starTab=tab[np.where(np.equal(tab['TYPE'], 1))]
        tab=tab[np.where(np.equal(tab['TYPE'], 2))]
    elif catalogFormat == 'FelipeS82':
        tab=atpy.Table().read(clusterDict['galaxyCatalogFile'], type = 'ascii')
        inFile=open(clusterDict['galaxyCatalogFile'])
        columnNames=inFile.readline().lstrip("# ").split()
        inFile.close()
        for i in range(len(columnNames)):
            try:
                tab.rename_column("col%d" % (i+1), columnNames[i])
            except:
                tab.rename_column("col%d" % (i+1), columnNames[i]+"_2")
        tab.rename_column("ra", "RADeg")
        tab.rename_column("dec", "decDeg")
        tab.rename_column("ID", "id")        
    elif catalogFormat == 'FelipeDR8':
        tab=atpy.Table().read(clusterDict['galaxyCatalogFile'], type = 'ascii')
        inFile=open(clusterDict['galaxyCatalogFile'])
        lines=inFile.readlines()
        inFile.close()
        columnNames=[]
        for line in lines:
            if line[0] == "#":
                columnNames.append(line.split(":")[-1].lstrip().rstrip("\n"))
        for i in range(len(columnNames)):
            try:
                tab.rename_column("col%d" % (i+1), columnNames[i])
            except:
                tab.rename_column("col%d" % (i+1), columnNames[i]+"_2")
        tab.rename_column("ra", "RADeg")
        tab.rename_column("dec", "decDeg")
        tab.add_column(atpy.Column(np.arange(len(tab), dtype = int), "id"))
    elif catalogFormat == 'FetchDR8':
        myCatalog=SDSSRetriever(cRADeg, cDecDeg, halfBoxSizeDeg = 5.0/60.0, DR = 8, objType = 'galaxy')
        keys=list(myCatalog[0].keys())
        tab=atpy.Table()
        for k in keys:
            arr=[]
            for obj in myCatalog:
                arr.append(obj[k])
            tab.add_column(atpy.Column(arr, k))
        # NOTE: we don't deredden, but we could
        # Colours useful for J2058
        tab.add_column(atpy.Column(tab['r']-tab['i'], 'r-i'))
        tab.add_column(atpy.Column(tab['g']-tab['r'], 'g-r')) 
        tab.add_column(atpy.Column(tab['g']-tab['i'], 'g-i')) 
        # For compatibility with Felipe catalogs
        tab.add_column(atpy.Column(tab['r'], 'pMag_r'))
        # Debugging
        outFileName=outDir+os.path.sep+"fetchedSDSSDR8GalaxyCatalog.fits"
        if os.path.exists(outFileName) == True:
            os.remove(outFileName)
        tab.table_name="SDSSGalaxies"
        tab.write(outFileName)
    elif catalogFormat == 'zCluster':
        # Use galaxyCatalogFile to select retriever
        database=clusterDict['galaxyCatalogFile']
        retriever, retrieverOptions, passbandSet=retrievers.getRetriever(database)
        if retrieverOptions is None:
            retrieverOptions={}
        retrieverOptions['altCacheDir']="retrieverCache"
        os.makedirs("retrieverCache", exist_ok = True)
        myCatalog=retriever(cRADeg, cDecDeg, halfBoxSizeDeg = 10.0/60.0, optionsDict = retrieverOptions)
        keys=myCatalog[0].keys()
        tab=atpy.Table()
        for k in keys:
            arr=[]
            for objDict in myCatalog:
                arr.append(objDict[k])
            tab[k]=arr
        # Hacky because we want colours for DECaLS
        try:
            tab['r-z']=tab['r']-tab['z']
        except:
            pass
        # Could add photo-z run here if we wanted
    elif catalogFormat == "zClusterOutput":
        tab=atpy.Table().read(clusterDict['galaxyCatalogFile'])
    else:
        raise Exception("didn't understand catalogFormat (%s)" % (catalogFormat))
    # Clunkyness below isn't really needed if we're going zCluster only?
    RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(cRADeg, cDecDeg, RArcmin/60.0)
    RAMask=np.logical_and(np.greater(tab['RADeg'], RAMin), np.less(tab['RADeg'], RAMax))
    decMask=np.logical_and(np.greater(tab['decDeg'], decMin), np.less(tab['decDeg'], decMax))
    posMask=np.logical_and(RAMask, decMask)
    myCatalogTable=tab[posMask]
    myCatalog=[]
    for row in myCatalogTable:
        objDict={}
        for key in list(tab.keys()):
            objDict[key]=row[key]
        myCatalog.append(objDict)
            
    # Get (usually SDSS) star catalog - use this to pick out alignment stars
    if catalogFormat == 'zCluster':# and database == 'CFHTLenS':
        # Assume we already have the retriever set-up above
        retrieverOptions['getStars']=True
        starCatalog=retriever(cRADeg, cDecDeg, halfBoxSizeDeg = 10.0/60.0, optionsDict = retrieverOptions)
        keys=starCatalog[0].keys()
        starTab=atpy.Table()
        for k in keys:
            arr=[]
            for objDict in starCatalog:
                arr.append(objDict[k])
            starTab[k]=arr
        starTab=tab[np.where(tab['r'] < 19)]
        starCatalog=[]
        refStarIDs=[]
        for row in starTab:
            objDict={}
            for key in list(tab.keys()):
                objDict[key]=row[key]
            starCatalog.append(objDict)
            refStarIDs.append(objDict['id'])
        brightStars=selectFromCatalog(starCatalog, ["r < 18.5"])
        catalog2DS9(brightStars, outDir+os.path.sep+"brightStars.reg", addInfo=[{'key': 'r', 'fmt': '%.3f'}], idKeyToUse = 'id',
                    includeRSSFoV = True, centreRADeg = cRADeg, centreDecDeg = cDecDeg)
    elif catalogFormat == 'Mathilde':
        starCatalog=[]
        refStarIDs=[]
        for row in starTab:
            objDict={}
            for key in list(tab.keys()):
                objDict[key]=row[key]
            starCatalog.append(objDict)
            refStarIDs.append(objDict['id'])
    else:
        starCatalog=SDSSRetriever(cRADeg, cDecDeg, halfBoxSizeDeg = 5.0/60.0, DR = 8, objType = 'star')
        # These are just for eyeballing which stars to include, since we'll need ids of these in my star catalogs
        brightStars=selectFromCatalog(starCatalog, ["r < 19"])
        catalog2DS9(brightStars, outDir+os.path.sep+"brightStars.reg", addInfo=[{'key': 'r', 'fmt': '%.3f'}], idKeyToUse = 'id',
                    includeRSSFoV = True, centreRADeg = cRADeg, centreDecDeg = cDecDeg)
        
    # This is just nice
    fetchLegacySurveyImage(clusterDict['name'], clusterDict['RADeg'], clusterDict['decDeg'])

    # Apply cuts on galaxies
    targetCatalog=selectFromCatalog(myCatalog, clusterDict['targetCuts'])

    # This catalog is for when we have areas of the mask with no primary targets in there
    fillerCatalog=selectFromCatalog(myCatalog, clusterDict['fillerCuts'])

    # Write out .reg files
    # Figue out what additional catalog keys to include first
    addInfo_targets=[]
    addedKeys=[]
    for cut in clusterDict['targetCuts']:
        key=cut.split()[0]
        if key not in addedKeys:
            addInfo_targets.append({'key': key, 'fmt': '%.2f'})
            addedKeys.append(key)
    addInfo_filler=[]
    addedKeys=[]
    for cut in clusterDict['fillerCuts']:
        key=cut.split()[0]
        if key not in addedKeys:
            addInfo_filler.append({'key': key, 'fmt': '%.2f'})       
            addedKeys.append(key)
    catalog2DS9(targetCatalog, outDir+os.path.sep+"targetGalaxies.reg", idKeyToUse = 'id', color = 'red',
                addInfo = addInfo_targets)
    catalog2DS9(fillerCatalog, outDir+os.path.sep+"fillerGalaxies.reg", idKeyToUse = 'id', color = 'green',
                addInfo = addInfo_filler)            
    catalog2DS9(starCatalog, outDir+os.path.sep+"starCatalog.reg", idKeyToUse='id', color = 'yellow',
                addInfo=[{'key': 'r', 'fmt': '%.2f'}])

    # List of galaxies to include in every mask - can be in filler catalog (if e.g. BCG missed by phot-z selection)
    alwaysIncludeList=[]
    for obj in myCatalog:
        if obj['id'] in clusterDict['alwaysIncludeIDs']:
            alwaysIncludeList.append(obj)
    for obj in fillerCatalog:
        if obj['id'] in clusterDict['alwaysIncludeIDs']:
            inListAlready=False
            for listObj in alwaysIncludeList:
                if obj['id'] == listObj['id']:
                    inListAlready=True
            if inListAlready == False:
                alwaysIncludeList.append(obj)            
    alwaysSlitLength=clusterDict['alwaysSlitLengthArcsec'] # Option to use bigger slits for e.g. BCGs

    # Now make the files for feeding into saltrsmt
    previousMasksList=[]

    for i in range(1, clusterDict['numMasks']+1):
        
        # Because may need this for setting in saltrsmt
        print(">>> Field centre: %s, %s" % (astCoords.decimal2hms(cRADeg, ":"), astCoords.decimal2dms(cDecDeg, ":")))
        
        DSSImageFileName=outDir+os.path.sep+"DSS_%s.fits" % (clusterDict['name'].replace(" ", "_"))
        fetchDSSImage(cRADeg, cDecDeg, DSSImageFileName, sizeArcmin = 10.0)
        #fetchDSSImageFromESO(cRADeg, cDecDeg, DSSImageFileName, sizeArcmin = 10.0)

        if catalogFormat != "Mathilde":
            try:
                refStarIDs=clusterDict['refStarIDs']
            except:
                refStarIDs=[]
            
        #if refStarIDs == []:
            #print("Need to select refStarIDs from brightStars.reg")
            #sys.exit()
        refStars=makeRefStarsList(starCatalog, refStarIDs, "Mask%d" % (i), outDir)
        
        targetDataFileName=outDir+os.path.sep+"%s_slitData_Mask%d.txt" % (clusterDict['name'].replace(" ", "_"), i)
        maskSlits=makeTargetSlitDataFile(targetDataFileName,
                                        targetCatalog, fillerCatalog, clusterDict['slitWidthArcsec'], clusterDict['slitLengthArcsec'], 
                                        cRADeg, cDecDeg, clusterDict['safetyArcsec'],
                                        refStarSlitLength = 5.0, refStarSlitWidth = 5.0, refStarsList = refStars, previousMasksList = previousMasksList,
                                        alwaysIncludeList = alwaysIncludeList, alwaysSlitLength = alwaysSlitLength,
                                        DSSImageFileName = DSSImageFileName,
                                        imageCuts = [3600, 8800])

        maskName="%sM%d" % (clusterDict['name'].replace(" ", "_"), i)

        makeRSMTFile(targetDataFileName, DSSImageFileName, clusterDict['name'], clusterDict['RADeg'], clusterDict['decDeg'],
                     clusterDict['proposer'], clusterDict['proposalCode'], maskName, RSMTDir)

        previousMasksList.append(maskSlits)

    # Count slits
    slitsCount=0
    for slitMask in previousMasksList:
        for objKey in slitMask:
            if 'catalogID' in list(slitMask[objKey].keys()) and slitMask[objKey]['catalogID'] > 0: # avoid adding arc
                slitsCount=slitsCount+1
    print("Total target galaxy slits placed (ignoring duplicates) = %d" % (slitsCount))

    # Overplot slits on SDSS image
    maskCount=0
    for slitMask in previousMasksList:
        maskCount=maskCount+1
        makeSlitsRGBPlot(clusterDict['name'], clusterDict['RADeg'], clusterDict['decDeg'], slitMask,
                          outDir+os.path.sep+"RGB_SlitsLoc_Mask%d.png" % (maskCount))
        
    # Plot where targeted slits are on the CMD
    parDict=clusterDict
    col=myCatalogTable[parDict['CMDCol']]
    #riErr=np.sqrt(myCatalogTable['errMagAper2_r']**2+myCatalogTable['errMagAper2_i']**2)
    #errMask=np.less(riErr, 0.2)
    r=myCatalogTable['r']
    colSlits=[]
    rSlits=[]
    slitsCount=0
    for slitMask in previousMasksList:
        for objKey in slitMask:
            if 'catalogID' in list(slitMask[objKey].keys()) and slitMask[objKey]['catalogID'] > 0: # avoid adding arc
                objMask=np.equal(myCatalogTable['id'], slitMask[objKey]['catalogID'])
                colSlits.append(myCatalogTable[parDict['CMDCol']][objMask][0])
                rSlits.append(myCatalogTable['r'][objMask][0])
                slitsCount=slitsCount+1
                
    plt.plot(r, col, 'b.')
    plt.plot(rSlits, colSlits, 'ro')
    plt.xlim(15, 25)
    plt.ylim(0, 2)
    plt.xlabel("r (AB)")
    plt.ylabel(parDict['CMDCol'])
    plt.savefig(outDir+os.path.sep+"CMD.png")
    plt.close()


