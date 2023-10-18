#
# Last modified on Tue Sep  4 16:32:29 PDT 2001 by lindy
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/math/ncoords.py,v 1.2.12.1 2016/02/11 19:04:10 annao Exp $
#

"""ncoords.py - Numeric coordinates

This class is intented to be the base class of a number
of classes which transform and generally operate on lists
of homogeneous coordinates.
"""

import numpy

class Ncoords:
    def __init__(self, refCoords, tolist=1):
        """refCoords is an nx3 list of n points
        
        resultCoords is set up and maintained as homogeneous coords
        if tolist then return the result coords as a python list
        """
        try:
            self.refCoords = numpy.array(numpy.concatenate(
                (refCoords, numpy.ones( (len(refCoords), 1), 'f')), 1))
        except TypeError:
            raise ValueError, "invalid input array"

        self.resultCoords = self.refCoords
        self.tolist = tolist


    def reset(self):
        self.resultCoords = self.refCoords


    def getResultCoords(self):
        """Return the list of result coordinates

        if tolist is set, return an nx3 Python ListType.
        if tolist is not set, return an nx4 numpy array.
        """
        if self.tolist:
            return self.resultCoords[:,:3].tolist()
        else:
            return self.resultCoords


