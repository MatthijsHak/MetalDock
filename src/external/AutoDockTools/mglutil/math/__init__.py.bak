import numpy
import types

def crossProduct (A, B, normal=True):
    """     Return cross product of two vectors A and B
normal: return normalized vector
"""
    res=[ A[1]*B[2] - A[2]*B[1],
          A[2]*B[0] - A[0]*B[2],
          A[0]*B[1] - A[1]*B[0] ]
    if normal:
        return norm(res)
    else:
        return res

def norm (A):
    """     Return normalized vector A.
"""
    if type(A) == types.ListType:
        A=numpy.array(A,'f')
        res= A/numpy.sqrt(numpy.dot(A,A))
        return res.tolist()    
    elif type(A)==numpy.ndarray:    
        return A/numpy.sqrt(numpy.dot(A,A))    
    else:
        print "Need a list or numpy array"
        return None

def getCenter(coords):
    """ get center of all the coords """
    coords=numpy.array(coords, 'f')    
    return (numpy.sum(coords, 0)/len(coords)).tolist()
