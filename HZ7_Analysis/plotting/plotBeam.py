#############################################################################################
#                                                                                           #
# Function used to draw the beam using the header of the cube and knowing the plotting axes #
#                                                                                           #
#############################################################################################
import numpy as np 
from matplotlib.patches import Ellipse
from matplotlib.offsetbox import (AnchoredOffsetbox, AuxTransformBox, DrawingArea, TextArea, VPacker)

def drawBeam(cubeHeader,ax):      
    majorBeam, minorBeam, positionalAngle = getBeamInfo(cubeHeader)
    majorBeamPixel, minorBeamPixel = convertArcsecsToPix(cubeHeader,majorBeam*3600),convertArcsecsToPix(cubeHeader,minorBeam*3600)
    aeb = AnchoredEllipseBeam(ax.transData,width=majorBeamPixel,height=minorBeamPixel,angle=positionalAngle)
    ax.add_artist(aeb)

class AnchoredEllipseBeam(AnchoredOffsetbox): #class that I need for plotting ellipse
    def __init__(self, transform, width, height, angle, loc = 'lower left',
                 pad=0.5, borderpad=0.1, prop=None, frameon=False):
        """
        Draw an ellipse the size in data coordinate of the give axes.
        pad, borderpad in fraction of the legend font size (or prop)
        """
        self._box = AuxTransformBox(transform)
        self.ellipse = Ellipse((0, 0), width, height, angle,fill=False,color='k',lw=2)
        self._box.add_artist(self.ellipse)
        super().__init__(loc, pad=pad, borderpad=borderpad,
                         child=self._box, prop=prop, frameon=frameon)

def getBeamInfo(cubeHeader):
    beamMajor = cubeHeader['BMAJ']
    beamMinor = cubeHeader['BMIN']
    beamPositionAngle = cubeHeader['BPA']
    return beamMajor,beamMinor, beamPositionAngle

def convertArcsecsToPix(cubeHeader,arcsecValue):
    degPerPix = cubeHeader['CDELT2']
    arcsecsPerPix = degPerPix * 3600
    return arcsecValue/arcsecsPerPix

#must be in arcseconds
def drawBeamManually(majorBeam,minorBeam,positionalAngle,cdeltValue,ax):
    degPerPix = cdeltValue
    arcsecsPerPix = degPerPix * 3600
    majorBeamPixel, minorBeamPixel = majorBeam/arcsecsPerPix, minorBeam/arcsecsPerPix
    aeb = AnchoredEllipseBeam(ax.transData,width=majorBeamPixel,height=minorBeamPixel,angle=positionalAngle)
    ax.add_artist(aeb)
