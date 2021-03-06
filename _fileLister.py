# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

#
#
#    ARCS Software 1.0
#
#    COPYRIGHT AND PERMISSION NOTICE
#    Copyright (c) 2006 California Institute of Technology.
#    All rights reserved.
#
#    Permission is hereby granted, free of charge, to any person obtaining a 
#    copy of this software and associated documentation files (the 
#    "Software"), to deal in the Software without restriction, including 
#    without limitation the rights to use, copy, modify, merge, publish, 
#    distribute, and/or sell copies of the Software, and to permit persons 
#    to whom the Software is furnished to do so, provided that the above 
#    copyright notice(s) and this permission notice appear in all copies of 
#    the Software and that both the above copyright notice(s) and this 
#    permission notice appear in supporting documentation.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
#    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT 
#    OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
#    HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL 
#    INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING 
#    FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, 
#    NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION 
#    WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#    Except as contained in this notice, the name of a copyright holder 
#    shall not be used in advertising or otherwise to promote the sale, use 
#    or other dealings in this Software without prior written authorization 
#    of the copyright holder.
#
#    All source code included in this distribution is covered by this notice,
##   unless specifically stated otherwise within each file. See each file within
#    each release for specific copyright holders.
#
#    ARCS is the name of an instrument under construction with U.S. DOE 
#    funding at Oak Ridge National Laboratory.
#

'''
utility functions to list files
'''

def recursiveFileList( path, exts ):
    res = []
    if "" in exts: return _recursiveListByExt( path, "" )
    for ext in exts:
        #print "find files with ext %s in %s ... " % (ext, path)
        l = _recursiveListByExt( path, ext )
        #print "found %s" % ( l, )
        res += l
        pass
    return res


def _listByExt( path, ext ):
    'return a list files all with the same given extension'
    from glob import glob
    import os
    return glob( os.path.join( path, '*%s'%ext ) )


def _recursiveListByExt( path, ext ):
    from os import listdir
    from os.path import join, isdir
    res = _listByExt( path, ext )
    for entry in listdir( path ):
        if ".svn" in entry or "CVS" in entry: continue
        fullpath = join( path, entry )
        if isdir(fullpath):
            #print _recursiveListByExt( fullpath, ext )
            res += _recursiveListByExt( fullpath, ext )
            pass
        pass
    return res

    

#backward compatibility
_recursiveListByExts = recursiveFileList
