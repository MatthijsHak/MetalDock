# 
# The contents of this file are subject to the Mozilla Publics
# License Version 1.1 (the "License"); you may not use this file
# except in compliance with the License. You may obtain a copy of
# the License at http://www.mozilla.org/MPL/
# 
# Software distributed under the License is distributed on an "AS
# IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
# implied. See the License for the specific language governing
# rights and limitations under the License.
# 
# The Original Code is "Java-Python Extension: testplus (JPE-testplus)".
# 
# The Initial Developer of the Original Code is Frederic Bruno Giacometti.
# Portions created by Frederic Bruno Giacometti are
# Copyright (C) 2002 Frederic Bruno Giacometti. All Rights Reserved.
# 
# Contributor: Frederic Giacometti, frederic.giacometti@arakne.com
#

#
#  Michel Sanner   June 28 2002
#
# This version of testplus was modified to run with python interpreter
# version 1.5.2 and older.
# WARNING: the testcollect fucntion will only sort the functions found
#          with version 2.1 and above. For earlier version the order is
#          the order of the hasing in global
#

import re, sys, os, operator
import string, types
from functools import reduce

def chdir(directory=None):
    if directory is None:
        if len(sys.argv):
            directory = os.path.split(sys.argv[0])[0]

    if len(directory):
        os.chdir(directory)


def logfile():
    global __logfile
    try:
        return __logfile
    except NameError:
        # add the regression directory to the path so that data files are
        # found in the local directory
        cwd = os.getcwd()
        direc = os.path.split(sys.argv[0])[0]
        if os.name == 'nt': #sys.platform == 'win32':
            direc = string.replace(direc, '/', '\\')
            
        if cwd[-len(direc):]!=direc:
            d = os.path.join( os.getcwd(), os.path.split(sys.argv[0])[0] )
            __logfile = open( sys.argv[ 0] + '.log', 'w')
        else:
            d = cwd
            __logfile = open( os.path.split(sys.argv[0])[1] + '.log', 'w')
        sys.path.append( d )

        return __logfile

## connect, disconnect and fun are 3 tuple (fun, args, kw) OR just a function

class TestHarness:
    def __init__( self,
                  name,
                  funs = [],
                  dependents = [],
                  connect = (lambda : None, (), {} ),
                  disconnect = (lambda : None, (), {} ) ):
        self.name = name
        self.funs = funs
        self.dependents = dependents
        self.connect = connect
        self.disconnect = disconnect
        self.count = 0


    def __call__( self, fun):
        try:
            if type(fun)==tuple:
                func = fun[0]
                args = fun[1]
                kw = fun[2]
            else:
                func = fun
                args = ()
                kw = {}
            func(*args, **kw)
        except:
            return fun, sys.exc_info()
        else:
            return None


    def __getattr__( self, attr):
        if attr == 'failures':
            result = [_f for _f in self.dependents if _f]
            if not result:
                fail = self( self.connect )
                if fail:
                    result = [fail]
                else:
                    for fun in self.funs:
                        self.count = self.count + 1
                        testname = 'TEST%4i %s.%s ' % (self.count,
                                                       self.name,
                                                       fun.__name__)
                        if len(testname) < 70:
                            testname = testname + ' '*(65-len(testname))
                        print(testname, end=' ')
                        prevstdout = sys.stdout
                        #prevstderr = sys.stderr
                        sys.stdout = logfile() # sys.stderr = logfile()
                        try:
                            print(testname)
                            sys.stdout.flush()
                            res = self( fun)
                            result.append( res )
                            sys.stdout.flush()
                        finally:
                            sys.stdout = prevstdout
                            #sys.stderr = prevstderr
                            if res is None:
                                print('PASSED')
                            else:
                                print('FAILED')
                    result.append( self( self.disconnect))
                    result = [_f for _f in result if _f]
                    
        else:
            raise AttributeError( attr)
        setattr( self, attr, result)
        return result

        
    def __len__( self):
        l1 = len( [x for x in self.failures if not isinstance( x, TestHarness)] )
        l2 = reduce(operator.add,
                    list(map( len, [x for x in self.failures if isinstance( x, TestHarness)] )),
                    0)
        val = l1+l2
        return val

        
    def __str__( self):
        from os import path
        klass = self.__class__
                              
        return string.join( ['\n',
                          'LOGFILE is <%s>' % path.abspath( logfile().name),
                          '\nTEST HARNESS %s: %s error%s'
                          ' out of %i tests:'
                             % (self.name, len( self) or 'SUCCESS - no',
                                1 < len( self) and 's' or '',
                                self.totalcount())]
                           + list(map(lambda x, self=self: re.sub( '\n', '\n    ',
                                   isinstance( x, TestHarness)
                                   and str( x) or self.error2str(*x)),
                                 self.failures))
                           + [''])


    def error2str( self, fun, xxx_todo_changeme):
        (exctype, excvalue, tb) = xxx_todo_changeme
        import traceback
        return '\n%s:\n  <%s>\n' % (exctype, excvalue) + \
               string.join([''] + traceback.format_tb( tb.tb_next))


    def totalcount( self):
        return reduce( operator.add,
                       [x.totalcount() for x in self.dependents],
                       self.count)


def getfile(x):
    return x.__code__.co_filename


def getsourcelines(x):
    return [x.__code__.co_firstlineno]

    
def testcollect( globs,
                 matchfun = lambda x: re.match( 'test', x)):

    result = []
    for x in list(globs.items()):
        if callable( x[ 1]) and matchfun( x[ 0]):
            result.append(x[1])
            result.sort( lambda x, y: cmp( (getfile( x), getsourcelines( x)[ -1]),
                                           (getfile( y), getsourcelines( y)[ -1])))
    return result


def prun( args = sys.argv[ 1:]):
    # win323 loses the stderr for subprocesses ... ?!
    sys.stderr = sys.stdout
    globs = {}
    try:
        for cmd in args:
            exec(cmd, globs)
    except SystemExit:
        raise
    except:
        import traceback
        traceback.print_stack()
        traceback.print_exc()
        sys.exit( 1)

    
